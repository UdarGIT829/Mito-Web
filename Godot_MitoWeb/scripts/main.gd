extends Control
## Main screen: project header, Files and Samples tabs, and a status bar with
## import progress and browser storage use.

const WORKSPACE_ACCEPT := ".mitoweb,.json"

var _name_edit: LineEdit
var _tabs: TabContainer
var _status: Label
var _storage: Label
var _progress: ProgressBar
var _cancel: Button
var _confirm: ConfirmationDialog
var _header_buttons: Array = []


func _ready() -> void:
	_build_ui()
	Workspace.changed.connect(_refresh)
	Workspace.error.connect(_set_status)
	Workspace.import_progress.connect(_on_import_progress)
	Workspace.import_finished.connect(_on_import_finished)
	WebFiles.files_picked.connect(_on_files_picked)
	WebFiles.storage_estimated.connect(_on_storage_estimated)
	_refresh()
	WebFiles.estimate_storage()
	if not Workspace.is_storage_persistent():
		_set_status("Browser storage is unavailable (private window?). Export your workspace to keep your work.")


func _build_ui() -> void:
	var margin := MarginContainer.new()
	margin.set_anchors_preset(Control.PRESET_FULL_RECT)
	for side in ["left", "right", "top", "bottom"]:
		margin.add_theme_constant_override("margin_" + side, 12)
	add_child(margin)
	var root := VBoxContainer.new()
	root.add_theme_constant_override("separation", 8)
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
	_add_header_button(header, "Export workspace", _export_workspace)
	_add_header_button(header, "Import workspace…", func() -> void:
		_confirm_then("Importing replaces the open project. Export it first if you want to keep it.",
			func() -> void: WebFiles.pick_files("workspace", WORKSPACE_ACCEPT, false)))
	_add_header_button(header, "New project", func() -> void:
		_confirm_then("Start a new, empty project? The open project will be discarded.", Workspace.new_project))

	_tabs = TabContainer.new()
	_tabs.size_flags_vertical = Control.SIZE_EXPAND_FILL
	root.add_child(_tabs)
	var files := FilesPanel.new()
	files.name = "Files"
	files.status_message.connect(_set_status)
	_tabs.add_child(files)
	var samples := SamplesPanel.new()
	samples.name = "Samples"
	_tabs.add_child(samples)
	var analysis := GraphPanel.new()
	analysis.name = "Analysis"
	analysis.status_message.connect(_set_status)
	_tabs.add_child(analysis)
	var tables := TablePanel.new()
	tables.name = "Tables"
	tables.status_message.connect(_set_status)
	_tabs.add_child(tables)
	analysis.open_table.connect(func(node_id: String) -> void:
		_tabs.current_tab = tables.get_index()
		tables.show_node(node_id))
	var stats := StatsPanel.new()
	stats.name = "Statistics"
	stats.status_message.connect(_set_status)
	_tabs.add_child(stats)
	analysis.open_stats.connect(func(node_id: String) -> void:
		_tabs.current_tab = stats.get_index()
		stats.show_node(node_id))

	var footer := HBoxContainer.new()
	footer.add_theme_constant_override("separation", 8)
	root.add_child(footer)
	_progress = ProgressBar.new()
	_progress.custom_minimum_size = Vector2(200, 0)
	_progress.max_value = 1.0
	_progress.step = 0.001
	_progress.visible = false
	footer.add_child(_progress)
	_cancel = Button.new()
	_cancel.text = "Cancel"
	_cancel.visible = false
	_cancel.pressed.connect(Workspace.cancel_import)
	footer.add_child(_cancel)
	_status = Label.new()
	_status.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_status.text_overrun_behavior = TextServer.OVERRUN_TRIM_ELLIPSIS
	footer.add_child(_status)
	_storage = Label.new()
	footer.add_child(_storage)

	_confirm = ConfirmationDialog.new()
	add_child(_confirm)


func _add_header_button(parent: Control, text: String, action: Callable) -> void:
	var button := Button.new()
	button.text = text
	button.pressed.connect(action)
	parent.add_child(button)
	_header_buttons.append(button)


func _refresh() -> void:
	if not _name_edit.has_focus():
		_name_edit.text = Workspace.model.get_name()
	var busy := Workspace.importing
	_progress.visible = busy
	_cancel.visible = busy
	for button in _header_buttons:
		button.disabled = busy


func _commit_name() -> void:
	var value := _name_edit.text.strip_edges()
	if value != Workspace.model.get_name() and not value.is_empty():
		Workspace.rename(value)


func _on_files_picked(purpose: String, files: Array) -> void:
	if files.is_empty():
		return
	match purpose:
		"sources":
			_set_status(Workspace.add_files(files))
			_tabs.current_tab = 0
			WebFiles.estimate_storage(true)
		"workspace":
			var message := Workspace.import_bytes(files[0]["bytes"])
			_set_status(message if message else "Opened %s." % files[0]["name"])
			WebFiles.estimate_storage()


func _on_import_progress(fraction: float, message: String) -> void:
	_progress.value = fraction
	_set_status(message)


func _on_import_finished(summary: String) -> void:
	_set_status(summary)
	WebFiles.estimate_storage()


func _export_workspace() -> void:
	_commit_name()
	WebFiles.save_file(Workspace.export_file_name(), Workspace.export_bytes(), "application/octet-stream")
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
