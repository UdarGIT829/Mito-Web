extends Control
## Phase 1 shell: project header, uploaded-source list, workspace
## export/import, and a status bar showing browser storage.

const SOURCE_ACCEPT := ".vcf,.csv,.tsv,.txt"
const WORKSPACE_ACCEPT := ".json"

var _name_edit: LineEdit
var _source_tree: Tree
var _status: Label
var _storage: Label
var _confirm: ConfirmationDialog


func _ready() -> void:
	_build_ui()
	Workspace.changed.connect(_refresh)
	Workspace.error.connect(_set_status)
	WebFiles.files_picked.connect(_on_files_picked)
	WebFiles.storage_estimated.connect(_on_storage_estimated)
	_refresh()
	WebFiles.estimate_storage()
	if not Workspace.is_storage_persistent():
		_set_status("Browser storage is unavailable (private window?). Export your workspace to keep your work.")


func _build_ui() -> void:
	set_anchors_preset(Control.PRESET_FULL_RECT)
	var root := VBoxContainer.new()
	root.set_anchors_preset(Control.PRESET_FULL_RECT)
	root.add_theme_constant_override("separation", 8)
	var margin := MarginContainer.new()
	margin.set_anchors_preset(Control.PRESET_FULL_RECT)
	for side in ["left", "right", "top", "bottom"]:
		margin.add_theme_constant_override("margin_" + side, 12)
	add_child(margin)
	margin.add_child(root)

	var header := HBoxContainer.new()
	header.add_theme_constant_override("separation", 8)
	root.add_child(header)
	var title := Label.new()
	title.text = "Mito Web"
	title.add_theme_font_size_override("font_size", 22)
	header.add_child(title)
	_name_edit = LineEdit.new()
	_name_edit.placeholder_text = "Project name"
	_name_edit.custom_minimum_size.x = 240
	_name_edit.text_submitted.connect(func(_t: String) -> void: _commit_name())
	_name_edit.focus_exited.connect(_commit_name)
	header.add_child(_name_edit)
	var spacer := Control.new()
	spacer.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	header.add_child(spacer)
	_add_button(header, "Upload files…", func() -> void:
		WebFiles.pick_files("sources", SOURCE_ACCEPT, true))
	_add_button(header, "Export workspace", _export_workspace)
	_add_button(header, "Import workspace…", func() -> void:
		_confirm_then("Importing replaces the open project. Export it first if you want to keep it.",
			func() -> void: WebFiles.pick_files("workspace", WORKSPACE_ACCEPT, false)))
	_add_button(header, "New project", func() -> void:
		_confirm_then("Start a new, empty project? The open project will be discarded.", Workspace.new_project))

	_source_tree = Tree.new()
	_source_tree.size_flags_vertical = Control.SIZE_EXPAND_FILL
	_source_tree.columns = 4
	_source_tree.column_titles_visible = true
	_source_tree.hide_root = true
	_source_tree.select_mode = Tree.SELECT_SINGLE
	var widths := [0, 90, 190, 90]
	for i in 4:
		_source_tree.set_column_title(i, ["Source file", "Size", "Added (UTC)", ""][i])
		_source_tree.set_column_expand(i, i == 0)
		_source_tree.set_column_custom_minimum_width(i, widths[i])
	_source_tree.cell_selected.connect(_on_cell_selected)
	root.add_child(_source_tree)

	var footer := HBoxContainer.new()
	root.add_child(footer)
	_status = Label.new()
	_status.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_status.text_overrun_behavior = TextServer.OVERRUN_TRIM_ELLIPSIS
	footer.add_child(_status)
	_storage = Label.new()
	footer.add_child(_storage)

	_confirm = ConfirmationDialog.new()
	add_child(_confirm)


func _add_button(parent: Control, text: String, action: Callable) -> void:
	var button := Button.new()
	button.text = text
	button.pressed.connect(action)
	parent.add_child(button)


func _refresh() -> void:
	if not _name_edit.has_focus():
		_name_edit.text = Workspace.model.get_name()
	_source_tree.clear()
	var tree_root := _source_tree.create_item()
	for source in Workspace.model.sources():
		var item := _source_tree.create_item(tree_root)
		item.set_text(0, source["name"])
		item.set_tooltip_text(0, "SHA-256 %s" % source["sha256"])
		item.set_text(1, String.humanize_size(int(source["size"])))
		item.set_text(2, str(source.get("added", "")).replace("T", " ").trim_suffix("Z"))
		item.set_metadata(0, source["id"])
		item.set_text(3, "Remove")
		item.set_custom_color(3, Color(1.0, 0.55, 0.5))
		item.set_tooltip_text(3, "Remove this file from the project")
	if Workspace.model.sources().is_empty():
		var hint := _source_tree.create_item(tree_root)
		hint.set_text(0, "No files yet. Use “Upload files…” to add VCF or table files; they stay in this browser.")
		for i in 4:
			hint.set_selectable(i, false)


func _commit_name() -> void:
	var value := _name_edit.text.strip_edges()
	if value != Workspace.model.get_name() and not value.is_empty():
		Workspace.rename(value)


func _on_cell_selected() -> void:
	var item := _source_tree.get_selected()
	if item == null or _source_tree.get_selected_column() != 3:
		return
	var source_id = item.get_metadata(0)
	if source_id != null:
		# Rebuilding the Tree inside its own selection signal is ignored, so
		# apply the removal on the next frame.
		_remove_source.call_deferred(int(source_id), item.get_text(0))


func _remove_source(source_id: int, file_name: String) -> void:
	Workspace.remove_source(source_id)
	_set_status("Removed %s." % file_name)


func _on_files_picked(purpose: String, files: Array) -> void:
	if files.is_empty():
		return
	match purpose:
		"sources":
			var added := 0
			for file in files:
				var before := Workspace.model.sources().size()
				Workspace.add_source(file["name"], file["bytes"])
				added += Workspace.model.sources().size() - before
			var skipped := files.size() - added
			_set_status("Added %d file(s)%s." % [added, " (%d already present)" % skipped if skipped else ""])
			WebFiles.estimate_storage(true)
		"workspace":
			var message := Workspace.import_bytes(files[0]["bytes"])
			_set_status(message if message else "Imported %s." % files[0]["name"])
			WebFiles.estimate_storage()


func _export_workspace() -> void:
	_commit_name()
	WebFiles.save_file(Workspace.export_file_name(), Workspace.export_bytes())
	_set_status("Exported %s." % Workspace.export_file_name())


func _on_storage_estimated(usage: int, quota: int, persisted: bool) -> void:
	if usage < 0:
		_storage.text = ""
		return
	_storage.text = "Browser storage: %s of %s%s" % [
		String.humanize_size(usage), String.humanize_size(quota), "" if persisted else " (may be evicted)"]


func _confirm_then(message: String, action: Callable) -> void:
	_confirm.dialog_text = message
	for connection in _confirm.confirmed.get_connections():
		_confirm.confirmed.disconnect(connection["callable"])
	_confirm.confirmed.connect(action, CONNECT_ONE_SHOT)
	_confirm.popup_centered()


func _set_status(message: String) -> void:
	_status.text = message
