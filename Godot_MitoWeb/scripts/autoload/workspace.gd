extends Node
## Autoload "Workspace": owns the open WorkspaceModel, runs imports, and
## persists it.
##
## user:// is backed by the browser's IndexedDB in web exports, so the
## autosave survives reloads. Browser storage can still be cleared or evicted,
## which is why the exported workspace file is the durable save.

signal changed
signal error(message: String)
signal import_progress(fraction: float, message: String)
signal import_finished(summary: String)

const AUTOSAVE_PATH := "user://workspace.mitoweb"
const LEGACY_AUTOSAVE_PATH := "user://workspace.json"
const AUTOSAVE_DELAY_SECONDS := 1.0

var model: WorkspaceModel
var importing := false
## Increases whenever imported samples or calls change, so cached analysis
## results know to recompute. Not saved.
var data_revision := 0
var _evaluator := GraphEvaluator.new()
var _last_results := {}
var _importer: VcfImporter
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
	data_revision += 1
	_mark_changed()


func rename(value: String) -> void:
	model.set_name(value)
	_mark_changed()


## Adds uploaded files. Returns [added, duplicates, unreadable].
func add_files(files: Array) -> Array:
	var counts := [0, 0, 0]
	for file in files:
		var added := model.add_source(file["name"], file["bytes"])
		if not added[1]:
			counts[1] += 1
		elif added[0]["status"] == WorkspaceModel.STATUS_ERROR:
			counts[2] += 1
		else:
			counts[0] += 1
	if counts[0] + counts[2] > 0:
		_mark_changed()
	return counts


func remove_sources(source_ids: Array) -> void:
	for source_id in source_ids:
		model.remove_source(source_id)
	data_revision += 1
	_mark_changed()


func undo_imports(source_ids: Array) -> void:
	for source_id in source_ids:
		model.undo_import(source_id)
	data_revision += 1
	_mark_changed()


## Evaluated analysis graph: {node id: {frame, error, summary}}. Cached per
## node, so calling this often is cheap.
func graph_results() -> Dictionary:
	if not importing:
		_last_results = _evaluator.evaluate(model.graph(), model.store, data_revision)
	return _last_results


## Call after editing the analysis graph in place.
func graph_edited() -> void:
	model.touch()
	_mark_changed()


## Call after editing a source's plan in place.
func plan_edited() -> void:
	model.touch()
	_mark_changed()


## Imports the given pending sources one after another.
func import_sources(source_ids: Array) -> void:
	if importing:
		return
	var queue: Array = []
	var missing_build: Array = []
	for source_id in source_ids:
		var source := model.source_by_id(source_id)
		if source.get("status") != WorkspaceModel.STATUS_PENDING:
			continue
		if source["plan"].get("build", "").is_empty():
			missing_build.append(source["name"])
		else:
			queue.append(source)
	if queue.is_empty():
		import_finished.emit(
			"Choose a reference build for: %s." % ", ".join(missing_build) if missing_build
			else "Nothing to import: select pending files.")
		return
	importing = true
	changed.emit()
	var imported := 0
	var calls := 0
	var failures: Array = []
	var cancelled := false
	for k in queue.size():
		var source: Dictionary = queue[k]
		_importer = VcfImporter.new()
		_importer.progress.connect(func(fraction: float, message: String) -> void:
			import_progress.emit((k + fraction) / queue.size(), message))
		var result: Dictionary = await _importer.import_source(source, model.store)
		if result.get("cancelled"):
			cancelled = true
			break
		if result["ok"]:
			model.mark_imported(source["id"], result)
			imported += 1
			calls += result["calls"]
		else:
			failures.append("%s: %s" % [source["name"], result["error"]])
	_importer = null
	importing = false
	data_revision += 1
	_mark_changed()
	var parts: Array = ["Imported %d file(s), %d calls." % [imported, calls]]
	if cancelled:
		parts.append("Cancelled; the file in progress was not imported.")
	if missing_build:
		parts.append("Skipped (no reference build chosen): %s." % ", ".join(missing_build))
	if failures:
		parts.append("Failed: %s" % "; ".join(failures))
	import_finished.emit(" ".join(parts))


func cancel_import() -> void:
	if _importer:
		_importer.cancelled = true


func export_bytes() -> PackedByteArray:
	return model.to_bytes()


func export_file_name() -> String:
	var safe := model.get_name().strip_edges().validate_filename().replace(" ", "_")
	if safe.is_empty():
		safe = "workspace"
	return "%s.%s" % [safe, WorkspaceModel.FILE_EXTENSION]


## Replaces the open workspace with an imported one. Returns an error message,
## or "" on success.
func import_bytes(bytes: PackedByteArray) -> String:
	if importing:
		return "Wait for the current import to finish."
	var result := WorkspaceModel.from_bytes(bytes)
	if not result["ok"]:
		return result["error"]
	model = result["model"]
	data_revision += 1
	_mark_changed()
	return ""


func save_autosave() -> void:
	var file := FileAccess.open(AUTOSAVE_PATH, FileAccess.WRITE)
	if file == null:
		error.emit("Autosave failed: %s" % error_string(FileAccess.get_open_error()))
		return
	file.store_buffer(model.to_bytes())
	file.close()


func _load_autosave() -> bool:
	for path in [AUTOSAVE_PATH, LEGACY_AUTOSAVE_PATH]:
		if not FileAccess.file_exists(path):
			continue
		var result := WorkspaceModel.from_bytes(FileAccess.get_file_as_bytes(path))
		if result["ok"]:
			model = result["model"]
			return true
		push_warning("Ignoring unreadable autosave %s: %s" % [path, result["error"]])
	return false


func _mark_changed() -> void:
	changed.emit()
	if not importing:
		_autosave_timer.start()
