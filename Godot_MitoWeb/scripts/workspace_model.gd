class_name WorkspaceModel
extends RefCounted
## Pure data layer for a Mito Web workspace: uploaded sources with their
## import plans, the imported VariantStore, and the saved file format.
## It has no knowledge of browsers or the scene tree, so it is testable
## headlessly.
##
## File format (version 2): the ASCII magic "MITOWEB\n", a u32 format
## version, a u32 uncompressed length, then zstd-compressed var_to_bytes()
## of the workspace dictionary. Version 1 (JSON, sources only) is migrated
## on load.

const FORMAT := "mito-web-workspace"
const FORMAT_VERSION := 2
const MAGIC := "MITOWEB\n"
const FILE_EXTENSION := "mitoweb"

const STATUS_PENDING := "pending"
const STATUS_IMPORTED := "imported"
const STATUS_ERROR := "error"

var data: Dictionary
var store: VariantStore


func _init(initial: Dictionary = {}, initial_store: VariantStore = null) -> void:
	data = initial if not initial.is_empty() else new_workspace_data()
	store = initial_store if initial_store else VariantStore.new()


static func new_workspace_data(name: String = "Untitled project") -> Dictionary:
	var now := _now()
	return {
		"format": FORMAT,
		"format_version": FORMAT_VERSION,
		"name": name,
		"created": now,
		"modified": now,
		"next_id": 1,
		"sources": [],
		"graph": GraphModel.new_graph(),
	}


func get_name() -> String:
	return data.get("name", "")


func set_name(value: String) -> void:
	data["name"] = value
	touch()


func sources() -> Array:
	return data["sources"]


## Workflows saved from graphs in this workspace (see Workflows).
func saved_workflows() -> Array:
	if not data.has("workflows"):
		data["workflows"] = []
	return data["workflows"]


## Live Ensembl results fetched in this workspace: {allele key: {"ensembl": {...}}}.
func live_annotations() -> Dictionary:
	if not data.has("live_annotations"):
		data["live_annotations"] = {}
	return data["live_annotations"]


## Lineage rules ([clade, group] pairs, see Lineages); empty until set.
func lineage_rules() -> Array:
	if not data.has("lineage_rules"):
		data["lineage_rules"] = []
	return data["lineage_rules"]


## The analysis graph (see GraphModel). Created on first use for older files.
func graph() -> Dictionary:
	if not data.has("graph"):
		data["graph"] = GraphModel.new_graph()
	return data["graph"]


func source_by_id(source_id: int) -> Dictionary:
	for source in sources():
		if source["id"] == source_id:
			return source
	return {}


func find_source_by_sha256(sha256: String) -> Dictionary:
	for source in sources():
		if source["sha256"] == sha256:
			return source
	return {}


## Adds an uploaded file and inspects it. Returns [source, is_new]; the same
## content uploaded again returns the existing source.
func add_source(file_name: String, bytes: PackedByteArray) -> Array:
	var sha := sha256_hex(bytes)
	var existing := find_source_by_sha256(sha)
	if not existing.is_empty():
		return [existing, false]
	# The same text in another form (gzipped, CRLF) is also a duplicate.
	var unzipped := GzipReader.decompress(bytes)
	if unzipped["ok"]:
		var content := VcfImporter.content_fingerprint(unzipped["bytes"])
		for other in sources():
			if content_sha256(other) == content:
				return [other, false]
	var source := {
		"id": int(data["next_id"]),
		"name": file_name,
		"size": bytes.size(),
		"sha256": sha,
		"added": _now(),
		"raw": bytes,
		"status": STATUS_PENDING,
		"error": "",
		"info": {},
		"plan": {"build": "", "samples": []},
		"result": {},
	}
	_inspect(source)
	data["next_id"] = int(data["next_id"]) + 1
	sources().append(source)
	touch()
	return [source, true]


## Fingerprint of a source's decompressed text (computed for older files).
func content_sha256(source: Dictionary) -> String:
	var info: Dictionary = source.get("info", {})
	if not info.has("content_sha256"):
		var unzipped := GzipReader.decompress(source["raw"])
		info["content_sha256"] = VcfImporter.content_fingerprint(unzipped["bytes"]) if unzipped["ok"] else source["sha256"]
		source["info"] = info
	return info["content_sha256"]


func _record_keys(source: Dictionary) -> PackedStringArray:
	var info: Dictionary = source.get("info", {})
	if not info.has("record_keys"):
		var unzipped := GzipReader.decompress(source["raw"])
		var reader := VcfReader.parse(unzipped["bytes"].get_string_from_utf8()) if unzipped["ok"] else null
		info["record_keys"] = VcfImporter.record_keys(reader) if reader and reader.ok else PackedStringArray()
		source["info"] = info
	return info["record_keys"]


## Files that repeat data already in the project. Imported files count as
## the originals; otherwise the earlier file does. Returns
## {source_id: {"of": id, "of_name": String, "reason": String}}.
##   same content: identical once decompressed (e.g. a .vcf and its .vcf.gz)
##   same sample: a single-sample file whose sample column matches another
##     file's, e.g. a strict and a filtered call set of one sample; generic
##     column names (e.g. "16") also need half of the calls to be shared
##   same samples: a multi-sample file sharing most sample columns
func duplicate_report() -> Dictionary:
	var order: Array = sources().filter(func(s): return s["status"] == STATUS_IMPORTED)
	order.append_array(sources().filter(func(s): return s["status"] == STATUS_PENDING))
	var by_content := {}
	var by_column := {}
	var out := {}
	for source in order:
		var content := content_sha256(source)
		if by_content.has(content):
			var other: Dictionary = by_content[content]
			out[source["id"]] = {"of": other["id"], "of_name": other["name"], "reason": "same content as %s" % other["name"]}
			continue
		var columns: Array = Array(source.get("info", {}).get("sample_columns", [])).map(func(c): return _sample_key(c)).filter(func(c): return c != "")
		var duplicate := {}
		if columns.size() == 1 and by_column.has(columns[0]):
			var other: Dictionary = by_column[columns[0]]
			var shared := _containment(_record_keys(source), _record_keys(other))
			if _specific(columns[0]) or shared >= 0.5:
				duplicate = {"of": other["id"], "of_name": other["name"],
					"reason": "same sample (%s) as %s; %d%% of calls shared" % [columns[0], other["name"], roundi(shared * 100)]}
		elif columns.size() > 1:
			var counts := {}
			for c in columns:
				if by_column.has(c):
					var other_id: int = by_column[c]["id"]
					counts[other_id] = counts.get(other_id, 0) + 1
			for other_id in counts:
				if counts[other_id] * 2 >= columns.size():
					var other := source_by_id(other_id)
					duplicate = {"of": other_id, "of_name": other["name"],
						"reason": "%d of %d samples also in %s" % [counts[other_id], columns.size(), other["name"]]}
					break
		if duplicate:
			out[source["id"]] = duplicate
			continue
		by_content[content] = source
		for c in columns:
			if not by_column.has(c):
				by_column[c] = source
	return out


## Removes every file the duplicate report flags (and its imported samples).
func remove_duplicates() -> int:
	var report := duplicate_report()
	for source_id in report:
		remove_source(source_id)
	return report.size()


static func _sample_key(column: String) -> String:
	var key := column.strip_edges().to_lower()
	for suffix in [".bam", ".cram", ".sam"]:
		key = key.trim_suffix(suffix)
	return key


## A column name specific enough to identify a sample by itself.
static func _specific(key: String) -> bool:
	return key.length() >= 4 and not key.is_valid_int() and not key in ["sample", "sample1", "tumor", "normal", "unknown", "none"]


static func _containment(a: PackedStringArray, b: PackedStringArray) -> float:
	if a.is_empty() or b.is_empty():
		return 0.0
	var small := a if a.size() <= b.size() else b
	var large := b if small == a else a
	var shared := 0
	for key in small:
		if large.bsearch(key) < large.size() and large[large.bsearch(key)] == key:
			shared += 1
	return float(shared) / small.size()


func remove_source(source_id: int) -> bool:
	var list := sources()
	for i in list.size():
		if list[i]["id"] == source_id:
			store.remove_source(source_id)
			list.remove_at(i)
			touch()
			return true
	return false


## Removes a source's imported samples and returns it to pending, keeping
## its plan so it can be edited and imported again.
func undo_import(source_id: int) -> void:
	var source := source_by_id(source_id)
	if source.is_empty() or source["status"] != STATUS_IMPORTED:
		return
	store.remove_source(source_id)
	source["status"] = STATUS_PENDING
	source["result"] = {}
	touch()


func mark_imported(source_id: int, result: Dictionary) -> void:
	var source := source_by_id(source_id)
	source["status"] = STATUS_IMPORTED
	source["error"] = ""
	source["result"] = result
	touch()


func touch() -> void:
	data["modified"] = _now()


func to_bytes() -> PackedByteArray:
	var payload := data.duplicate()
	payload["store"] = store.to_dict()
	var body := var_to_bytes(payload)
	var out := MAGIC.to_ascii_buffer()
	var header := PackedByteArray()
	header.resize(8)
	header.encode_u32(0, FORMAT_VERSION)
	header.encode_u32(4, body.size())
	out.append_array(header)
	out.append_array(body.compress(FileAccess.COMPRESSION_ZSTD))
	return out


## Parses and validates a saved workspace. Returns {"ok": true, "model":
## WorkspaceModel, "migrated": bool} or {"ok": false, "error": String}.
static func from_bytes(bytes: PackedByteArray) -> Dictionary:
	var magic := MAGIC.to_ascii_buffer()
	if bytes.size() >= magic.size() + 8 and bytes.slice(0, magic.size()) == magic:
		return _from_binary(bytes, magic.size())
	var text := bytes.get_string_from_utf8().strip_edges()
	if text.begins_with("{"):
		return _from_legacy_json(text)
	return _failure("The file is not a Mito Web workspace.")


static func _from_binary(bytes: PackedByteArray, offset: int) -> Dictionary:
	var version := bytes.decode_u32(offset)
	if version > FORMAT_VERSION:
		return _failure(
			"This workspace was saved by a newer version of Mito Web (format %d; this version reads up to %d)."
			% [version, FORMAT_VERSION])
	if version < 2:
		return _failure("The workspace has an invalid format version.")
	var size := bytes.decode_u32(offset + 4)
	var body := bytes.slice(offset + 8).decompress(size, FileAccess.COMPRESSION_ZSTD)
	if body.size() != size:
		return _failure("The workspace file is damaged (it could not be decompressed).")
	var parsed = bytes_to_var(body)
	if typeof(parsed) != TYPE_DICTIONARY or parsed.get("format") != FORMAT:
		return _failure("The workspace file is damaged or is not a Mito Web workspace.")
	if typeof(parsed.get("sources")) != TYPE_ARRAY or typeof(parsed.get("store")) != TYPE_DICTIONARY:
		return _failure("The workspace is missing its sources or variant store.")
	for source in parsed["sources"]:
		var problem := _validate_source(source)
		if problem:
			return _failure(problem)
	if parsed.has("graph"):
		var graph_problem := GraphModel.validate(parsed["graph"])
		if graph_problem:
			return _failure(graph_problem)
	for w in parsed.get("workflows", []):
		if typeof(w) != TYPE_DICTIONARY or GraphModel.validate(w.get("graph")) != "":
			return _failure("A saved workflow in this file is malformed.")
	var restored := VariantStore.from_dict(parsed["store"])
	if not restored["ok"]:
		return _failure(restored["error"])
	parsed.erase("store")
	parsed["format_version"] = FORMAT_VERSION
	return {"ok": true, "model": new(parsed, restored["store"]), "migrated": false}


## Format 1 stored uploaded files as text and nothing else. Each file becomes
## a pending source again.
static func _from_legacy_json(text: String) -> Dictionary:
	var json := JSON.new()
	if json.parse(text) != OK:
		return _failure("The file is not a Mito Web workspace (invalid JSON at line %d)." % json.get_error_line())
	var parsed = json.data
	if typeof(parsed) != TYPE_DICTIONARY or parsed.get("format") != FORMAT:
		return _failure("The file is not a Mito Web workspace.")
	var version = parsed.get("format_version")
	if typeof(version) not in [TYPE_INT, TYPE_FLOAT] or int(version) != 1:
		return _failure("The workspace has an invalid format version.")
	if typeof(parsed.get("sources")) != TYPE_ARRAY:
		return _failure("The workspace is missing its source list.")
	var model := new(new_workspace_data(str(parsed.get("name", "Untitled project"))))
	model.data["created"] = str(parsed.get("created", model.data["created"]))
	for old in parsed["sources"]:
		if typeof(old) != TYPE_DICTIONARY or typeof(old.get("text")) != TYPE_STRING or typeof(old.get("name")) != TYPE_STRING:
			return _failure("The workspace contains a malformed source entry.")
		model.add_source(old["name"], old["text"].to_utf8_buffer())
	return {"ok": true, "model": model, "migrated": true}


static func _validate_source(source) -> String:
	if typeof(source) != TYPE_DICTIONARY:
		return "The workspace contains a malformed source entry."
	var expected := {
		"id": TYPE_INT, "name": TYPE_STRING, "sha256": TYPE_STRING, "raw": TYPE_PACKED_BYTE_ARRAY,
		"status": TYPE_STRING, "plan": TYPE_DICTIONARY, "info": TYPE_DICTIONARY,
	}
	for key in expected:
		if typeof(source.get(key)) != expected[key]:
			return "A source entry has a missing or malformed \"%s\"." % key
	return ""


func _inspect(source: Dictionary) -> void:
	var inspected := VcfImporter.inspect(source["name"], source["raw"])
	if inspected["ok"]:
		source["info"] = inspected["info"]
		source["plan"] = inspected["plan"]
	else:
		source["status"] = STATUS_ERROR
		source["error"] = inspected["error"]


static func sha256_hex(bytes: PackedByteArray) -> String:
	var ctx := HashingContext.new()
	ctx.start(HashingContext.HASH_SHA256)
	ctx.update(bytes)
	return ctx.finish().hex_encode()


static func _failure(message: String) -> Dictionary:
	return {"ok": false, "error": message}


static func _now() -> String:
	return Time.get_datetime_string_from_system(true) + "Z"
