class_name WorkflowDialog
extends AcceptDialog
## Lists built-in and saved workflows, shows what each does and needs, asks
## for its settings, and inserts it into the analysis graph. Also saves the
## current graph as a workflow and imports/exports workflow files.

signal insert_requested(graph: Dictionary, replace: bool, title: String)
signal status_message(text: String)

var _list: ItemList
var _title: Label
var _description: Label
var _needs: Label
var _form: GridContainer
var _insert: Button
var _replace: Button
var _export: Button
var _delete: Button
var _save_dialog: ConfirmationDialog
var _save_title: LineEdit
var _save_description: TextEdit
var _entries: Array = []  # workflow dictionaries, parallel to list items (null for headers)
var _inputs := {}  # param id -> Control


func _init() -> void:
	title = "Workflows"
	ok_button_text = "Close"
	min_size = Vector2i(860, 520)
	var root := VBoxContainer.new()
	add_child(root)
	var split := HSplitContainer.new()
	split.size_flags_vertical = Control.SIZE_EXPAND_FILL
	root.add_child(split)
	_list = ItemList.new()
	_list.custom_minimum_size = Vector2(260, 380)
	_list.item_selected.connect(func(_i) -> void: _show_selected())
	split.add_child(_list)
	var detail := VBoxContainer.new()
	detail.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	detail.add_theme_constant_override("separation", 8)
	split.add_child(detail)
	_title = Label.new()
	_title.add_theme_font_size_override("font_size", 18)
	detail.add_child(_title)
	_description = _wrapped()
	detail.add_child(_description)
	_needs = _wrapped()
	_needs.add_theme_color_override("font_color", Color(1, 0.8, 0.4))
	detail.add_child(_needs)
	_form = GridContainer.new()
	_form.columns = 2
	detail.add_child(_form)
	var buttons := HBoxContainer.new()
	detail.add_child(buttons)
	_insert = _button("Insert into graph", func() -> void: _emit_insert(false))
	buttons.add_child(_insert)
	_replace = _button("Replace graph", func() -> void: _emit_insert(true))
	_replace.tooltip_text = "Clears the current graph first."
	buttons.add_child(_replace)
	_export = _button("Export file", _export_selected)
	buttons.add_child(_export)
	_delete = _button("Delete", _delete_selected)
	buttons.add_child(_delete)
	var bottom := HBoxContainer.new()
	root.add_child(bottom)
	bottom.add_child(_button("Save current graph as workflow…", func() -> void:
		_save_title.text = ""
		_save_description.text = ""
		_save_dialog.popup_centered()
		_save_title.grab_focus.call_deferred()))  # typing must go to the name, not press Save
	bottom.add_child(_button("Import workflow file…", func() -> void:
		WebFiles.pick_files("workflow", ".json,.mitoflow", false)))
	_build_save_dialog()
	WebFiles.files_picked.connect(_on_files_picked)


func open() -> void:
	_fill_list()
	popup_centered()


func _fill_list() -> void:
	var previous := _selected()
	_list.clear()
	_entries.clear()
	_header("Built-in")
	for w in Workflows.builtin():
		_list.add_item("  " + w["title"])
		_entries.append(w)
	var saved: Array = Workspace.model.saved_workflows()
	_header("Saved in this project" if saved else "Saved in this project (none yet)")
	for w in saved:
		_list.add_item("  " + w["title"])
		_entries.append(w)
	var select := 1
	for i in _entries.size():
		if _entries[i] != null and previous and _entries[i]["id"] == previous["id"]:
			select = i
	_list.select(select)
	_show_selected()


func _header(text: String) -> void:
	_list.add_item(text, null, false)
	_list.set_item_custom_fg_color(_list.item_count - 1, Color(0.75, 0.85, 1.0))
	_entries.append(null)


func _selected() -> Dictionary:
	var items := _list.get_selected_items()
	if items.is_empty() or _entries[items[0]] == null:
		return {}
	return _entries[items[0]]


func _show_selected() -> void:
	for child in _form.get_children():
		child.queue_free()
	_inputs.clear()
	var w := _selected()
	_title.text = w.get("title", "")
	_description.text = w.get("description", "")
	_needs.text = ("Needs: " + w["needs"]) if w.get("needs") else ""
	var builtin := w.has("build")
	_export.visible = not w.is_empty() and not builtin
	_delete.visible = not w.is_empty() and not builtin
	_insert.disabled = w.is_empty()
	_replace.disabled = w.is_empty()
	for p in w.get("params", []):
		var label := Label.new()
		label.text = p["label"]
		_form.add_child(label)
		var control := _param_control(p)
		_form.add_child(control)
		_inputs[p["id"]] = control


func _param_control(p: Dictionary) -> Control:
	match p["type"]:
		"number":
			var spin := SpinBox.new()
			spin.min_value = p.get("min", 0)
			spin.max_value = p.get("max", 100)
			spin.step = 0.5
			spin.value = p["default"]
			return spin
		"haplogroup":
			var option := OptionButton.new()
			var found := Workflows.haplogroups_in(Workspace.model.store)
			if found.is_empty():
				option.add_item("%s (no haplogroup tags yet)" % p["default"])
				option.set_item_metadata(0, p["default"])
			for pair in found:
				option.add_item("%s (%d samples)" % pair)
				option.set_item_metadata(option.item_count - 1, pair[0])
				if pair[0] == p["default"]:
					option.select(option.item_count - 1)
			return option
		_:  # tag
			var option := OptionButton.new()
			var tags: Array = SampleQuery.tag_counts(Workspace.model.store).map(func(t): return t[0])
			tags.erase(SampleQuery.UNTAGGED)
			if not p["default"] in tags:
				tags.push_front(p["default"])
			for tag in tags:
				option.add_item(tag)
				option.set_item_metadata(option.item_count - 1, tag)
				if tag == p["default"]:
					option.select(option.item_count - 1)
			return option


func _values() -> Dictionary:
	var out := {}
	for id in _inputs:
		var control: Control = _inputs[id]
		if control is SpinBox:
			out[id] = control.value
		elif control is OptionButton:
			out[id] = control.get_item_metadata(control.selected)
	return out


func _emit_insert(replace: bool) -> void:
	var w := _selected()
	if w.is_empty():
		return
	insert_requested.emit(Workflows.build(w, _values()), replace, w["title"])
	hide()


func _export_selected() -> void:
	var w := _selected()
	if w.is_empty() or w.has("build"):
		return
	var name := str(w["title"]).validate_filename().replace(" ", "_") + ".mitoflow.json"
	WebFiles.save_file(name, Workflows.to_file_text(w).to_utf8_buffer(), "application/json")
	status_message.emit("Exported workflow %s." % name)


func _delete_selected() -> void:
	var w := _selected()
	if w.is_empty() or w.has("build"):
		return
	Workspace.model.saved_workflows().erase(w)
	Workspace.graph_edited()
	_fill_list()


func _build_save_dialog() -> void:
	_save_dialog = ConfirmationDialog.new()
	_save_dialog.title = "Save graph as workflow"
	_save_dialog.ok_button_text = "Save"
	var box := VBoxContainer.new()
	_save_dialog.add_child(box)
	box.add_child(_wrapped_text("Saves the current graph (nodes, settings, connections) so it can be inserted again or exported to share."))
	_save_title = LineEdit.new()
	_save_title.placeholder_text = "Name"
	box.add_child(_save_title)
	_save_dialog.register_text_enter(_save_title)  # Enter saves
	_save_description = TextEdit.new()
	_save_description.placeholder_text = "What it does and what it needs (optional)"
	_save_description.custom_minimum_size = Vector2(460, 90)
	box.add_child(_save_description)
	_save_dialog.confirmed.connect(func() -> void:
		var graph := Workspace.model.graph()
		if graph["nodes"].is_empty():
			status_message.emit("The graph is empty; nothing to save.")
			return
		var name := _save_title.text.strip_edges()
		var w := Workflows.from_graph(graph, name if name else "My workflow", _save_description.text.strip_edges())
		Workspace.model.saved_workflows().append(w)
		Workspace.graph_edited()
		status_message.emit("Saved workflow \"%s\"." % w["title"])
		_fill_list())
	add_child(_save_dialog)


func _on_files_picked(purpose: String, files: Array) -> void:
	if purpose != "workflow" or files.is_empty():
		return
	var parsed := Workflows.from_file_text(files[0]["bytes"].get_string_from_utf8())
	if not parsed["ok"]:
		status_message.emit(parsed["error"])
		return
	Workspace.model.saved_workflows().append(parsed["workflow"])
	Workspace.graph_edited()
	status_message.emit("Imported workflow \"%s\"." % parsed["workflow"]["title"])
	if visible:
		_fill_list()


static func _wrapped() -> Label:
	var l := Label.new()
	l.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	l.custom_minimum_size.x = 520  # wrapping labels need a width
	return l


static func _wrapped_text(text: String) -> Label:
	var l := _wrapped()
	l.text = text
	l.custom_minimum_size.x = 460
	return l


static func _button(text: String, action: Callable) -> Button:
	var b := Button.new()
	b.text = text
	b.pressed.connect(action)
	return b
