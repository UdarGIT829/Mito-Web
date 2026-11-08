class_name WorkspaceModel
extends RefCounted
## Pure data layer for a Mito Web workspace: creation, validation, and
## (de)serialization. It has no knowledge of files, browsers, or the scene
## tree, so it can be tested headlessly.

const FORMAT := "mito-web-workspace"
const FORMAT_VERSION := 1

var data: Dictionary


func _init(initial: Dictionary = {}) -> void:
	data = initial if not initial.is_empty() else new_workspace_data()


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


func find_source_by_sha256(sha256: String) -> Dictionary:
	for source in sources():
		if source["sha256"] == sha256:
			return source
	return {}


## Adds an uploaded file. Returns the existing source instead when the same
## content was already added, so re-uploading a file is harmless.
func add_source(file_name: String, bytes: PackedByteArray) -> Dictionary:
	var sha := sha256_hex(bytes)
	var existing := find_source_by_sha256(sha)
	if not existing.is_empty():
		return existing
	var source := {
		"id": int(data["next_id"]),
		"name": file_name,
		"size": bytes.size(),
		"sha256": sha,
		"added": _now(),
		# Phase 1 keeps the raw text; the import phase replaces this with
		# parsed, columnar sample data.
		"text": bytes.get_string_from_utf8(),
	}
	data["next_id"] = int(data["next_id"]) + 1
	sources().append(source)
	touch()
	return source


func remove_source(source_id: int) -> bool:
	var list := sources()
	for i in list.size():
		if int(list[i]["id"]) == source_id:
			list.remove_at(i)
			touch()
			return true
	return false


func touch() -> void:
	data["modified"] = _now()


func to_json() -> String:
	return JSON.stringify(data, "", false)


## Parses and validates serialized workspace text. Returns
## {"ok": true, "model": WorkspaceModel} or {"ok": false, "error": String}.
static func from_json(text: String) -> Dictionary:
	var json := JSON.new()
	if json.parse(text) != OK:
		return _failure("The file is not a Mito Web workspace (invalid JSON at line %d)." % json.get_error_line())
	var parsed = json.data
	if typeof(parsed) != TYPE_DICTIONARY:
		return _failure("The file is not a Mito Web workspace (not a JSON object).")
	if parsed.get("format") != FORMAT:
		return _failure("The file is not a Mito Web workspace.")
	var version = parsed.get("format_version")
	if typeof(version) not in [TYPE_INT, TYPE_FLOAT] or int(version) < 1:
		return _failure("The workspace has no valid format version.")
	if int(version) > FORMAT_VERSION:
		return _failure(
			"This workspace was saved by a newer version of Mito Web (format %d; this version reads up to %d)."
			% [int(version), FORMAT_VERSION]
		)
	if typeof(parsed.get("sources")) != TYPE_ARRAY:
		return _failure("The workspace is missing its source list.")
	for source in parsed["sources"]:
		if typeof(source) != TYPE_DICTIONARY:
			return _failure("The workspace contains a malformed source entry.")
		for key in ["id", "name", "size", "sha256", "text"]:
			if not source.has(key):
				return _failure("A source entry is missing \"%s\"." % key)
	# JSON numbers come back as floats; restore integer fields.
	parsed["format_version"] = int(version)
	parsed["next_id"] = int(parsed.get("next_id", 1))
	for source in parsed["sources"]:
		source["id"] = int(source["id"])
		source["size"] = int(source["size"])
		parsed["next_id"] = maxi(parsed["next_id"], source["id"] + 1)
	return {"ok": true, "model": new(parsed)}


static func sha256_hex(bytes: PackedByteArray) -> String:
	var ctx := HashingContext.new()
	ctx.start(HashingContext.HASH_SHA256)
	ctx.update(bytes)
	return ctx.finish().hex_encode()


static func _failure(message: String) -> Dictionary:
	return {"ok": false, "error": message}


static func _now() -> String:
	return Time.get_datetime_string_from_system(true) + "Z"
