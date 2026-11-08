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
	}


func get_name() -> String:
	return data.get("name", "")


func set_name(value: String) -> void:
	data["name"] = value
	touch()


func sources() -> Array:
	return data["sources"]


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
