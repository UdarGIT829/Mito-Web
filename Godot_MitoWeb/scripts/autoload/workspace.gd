extends Node
## Autoload "Workspace": owns the open WorkspaceModel and persists it.
##
## user:// is backed by the browser's IndexedDB in web exports, so the
## autosave survives reloads. Browser storage can still be cleared or evicted,
## which is why the exported workspace file is the durable save.

signal changed
signal error(message: String)

const AUTOSAVE_PATH := "user://workspace.json"
const AUTOSAVE_DELAY_SECONDS := 1.0

var model: WorkspaceModel
var _autosave_timer: Timer


func _ready() -> void:
	_autosave_timer = Timer.new()
	_autosave_timer.one_shot = true
	_autosave_timer.wait_time = AUTOSAVE_DELAY_SECONDS
	_autosave_timer.timeout.connect(save_autosave)
	add_child(_autosave_timer)
	if not _load_autosave():
		model = WorkspaceModel.new()


func is_storage_persistent() -> bool:
	return OS.is_userfs_persistent()


func new_project() -> void:
	model = WorkspaceModel.new()
	_mark_changed()


func rename(value: String) -> void:
	model.set_name(value)
	_mark_changed()


func add_source(file_name: String, bytes: PackedByteArray) -> Dictionary:
	var before := model.sources().size()
	var source := model.add_source(file_name, bytes)
	if model.sources().size() != before:
		_mark_changed()
	return source


func remove_source(source_id: int) -> void:
	if model.remove_source(source_id):
		_mark_changed()


func export_bytes() -> PackedByteArray:
	return model.to_json().to_utf8_buffer()


func export_file_name() -> String:
	var safe := model.get_name().strip_edges().validate_filename().replace(" ", "_")
	if safe.is_empty():
		safe = "workspace"
	return "%s.mitoweb.json" % safe


## Replaces the open workspace with an imported one. Returns an error message,
## or "" on success.
func import_bytes(bytes: PackedByteArray) -> String:
	var result := WorkspaceModel.from_json(bytes.get_string_from_utf8())
	if not result["ok"]:
		return result["error"]
	model = result["model"]
	_mark_changed()
	return ""


func save_autosave() -> void:
	var file := FileAccess.open(AUTOSAVE_PATH, FileAccess.WRITE)
	if file == null:
		error.emit("Autosave failed: %s" % error_string(FileAccess.get_open_error()))
		return
	file.store_string(model.to_json())
	file.close()


func _load_autosave() -> bool:
	if not FileAccess.file_exists(AUTOSAVE_PATH):
		return false
	var result := WorkspaceModel.from_json(FileAccess.get_file_as_string(AUTOSAVE_PATH))
	if not result["ok"]:
		push_warning("Ignoring unreadable autosave: %s" % result["error"])
		return false
	model = result["model"]
	return true


func _mark_changed() -> void:
	changed.emit()
	_autosave_timer.start()
