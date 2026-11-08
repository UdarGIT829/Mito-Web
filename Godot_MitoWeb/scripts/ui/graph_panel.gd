class_name GraphPanel
extends VBoxContainer
## The analysis graph editor. Nodes are GraphNodes built from the GraphModel
## data in the workspace; every edit updates that data, and the graph is
## re-evaluated shortly after (unchanged nodes reuse cached results).

signal status_message(text: String)

const FRAME_COLOR := Color(0.45, 0.8, 1.0)
const EXCLUDE_COLOR := Color(1.0, 0.45, 0.4)
const ERROR_COLOR := Color(1, 0.5, 0.45)
const OK_COLOR := Color(0.7, 0.95, 0.75)
const ADD_ORDER := ["data_provider", "tag_filter", "comparison", "logical_filter", "enrich", "table"]

var _graph_edit: GraphEdit
var _add_menu: MenuButton
var _context_menu: PopupMenu
var _context_position := Vector2.ZERO
var _filter_dialog: FilterDialog
var _table_dialog: TableDialog
var _confirm: ConfirmationDialog
var _evaluator := GraphEvaluator.new()
var _results := {}
var _graph_ref: Dictionary = {}
var _views := {}  # node id -> {"node": GraphNode, "status": Label, ...}
var _tag_pickers: Array = []
var _eval_timer: Timer
var _tags_revision := -1


func _ready() -> void:
	var bar := HBoxContainer.new()
	add_child(bar)
	_add_menu = MenuButton.new()
	_add_menu.text = "Add node"
	_add_menu.flat = false
	for i in ADD_ORDER.size():
		_add_menu.get_popup().add_item(GraphModel.TYPES[ADD_ORDER[i]]["title"], i)
	_add_menu.get_popup().id_pressed.connect(func(id: int) -> void:
		_add_node(ADD_ORDER[id], _visible_center()))
	bar.add_child(_add_menu)
	var example := Button.new()
	example.text = "Insert example"
	example.tooltip_text = "Treated vs untreated: exclusive mutations and those in both, using the tags Treated and Untreated."
	example.pressed.connect(_insert_example)
	bar.add_child(example)
	var clear := Button.new()
	clear.text = "Clear graph"
	clear.pressed.connect(func() -> void:
		_confirm.dialog_text = "Remove every node from the graph?"
		_confirm.popup_centered())
	bar.add_child(clear)
	var hint := Label.new()
	hint.text = "Drag from an output (right) to an input (left). Right-click the canvas to add nodes; select a node and press Delete to remove it."
	hint.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	hint.text_overrun_behavior = TextServer.OVERRUN_TRIM_ELLIPSIS
	hint.add_theme_color_override("font_color", Color(0.7, 0.7, 0.7))
	bar.add_child(hint)

	_graph_edit = GraphEdit.new()
	_graph_edit.size_flags_vertical = Control.SIZE_EXPAND_FILL
	_graph_edit.right_disconnects = true
	_graph_edit.minimap_enabled = false
	_graph_edit.connection_request.connect(_on_connection_request)
	_graph_edit.disconnection_request.connect(_on_disconnection_request)
	_graph_edit.delete_nodes_request.connect(_on_delete_nodes)
	_graph_edit.end_node_move.connect(_on_end_node_move)
	_graph_edit.popup_request.connect(_on_popup_request)
	add_child(_graph_edit)

	_context_menu = PopupMenu.new()
	for i in ADD_ORDER.size():
		_context_menu.add_item("Add " + GraphModel.TYPES[ADD_ORDER[i]]["title"], i)
	_context_menu.id_pressed.connect(func(id: int) -> void: _add_node(ADD_ORDER[id], _context_position))
	add_child(_context_menu)
	_filter_dialog = FilterDialog.new()
	add_child(_filter_dialog)
	_table_dialog = TableDialog.new()
	add_child(_table_dialog)
	_confirm = ConfirmationDialog.new()
	_confirm.confirmed.connect(func() -> void:
		var g := Workspace.model.graph()
		g["nodes"].clear()
		g["edges"].clear()
		_commit(true))
	add_child(_confirm)

	_eval_timer = Timer.new()
	_eval_timer.one_shot = true
	_eval_timer.wait_time = 0.15
	_eval_timer.timeout.connect(_evaluate)
	add_child(_eval_timer)
	Workspace.changed.connect(_on_workspace_changed)
	_rebuild()


func _on_workspace_changed() -> void:
	if not is_same(Workspace.model.graph(), _graph_ref):
		_rebuild()
		return
	if Workspace.data_revision != _tags_revision:
		_refresh_tag_lists()
	_eval_timer.start()


func _rebuild() -> void:
	_graph_ref = Workspace.model.graph()
	_graph_edit.clear_connections()
	for view in _views.values():
		_graph_edit.remove_child(view["node"])
		view["node"].queue_free()
	_views.clear()
	_tag_pickers.clear()
	for n in _graph_ref["nodes"]:
		_create_view(n)
	_sync_connections()
	_refresh_tag_lists()
	_evaluate()


func _commit(rebuild: bool = false) -> void:
	Workspace.graph_edited()
	if rebuild:
		_rebuild()
	else:
		_eval_timer.start()


# --- Editing -----------------------------------------------------------------

func _add_node(type: String, position: Vector2) -> void:
	var n := GraphModel.add_node(Workspace.model.graph(), type, position)
	_create_view(n)
	_refresh_tag_lists()
	_commit()


func _insert_example() -> void:
	var g := Workspace.model.graph()
	var example := GraphModel.example_graph()
	var offset := Vector2(0, 0)
	for n in g["nodes"]:
		offset.y = maxf(offset.y, n["position"].y + 320)
	var renamed := {}
	for n in example["nodes"]:
		var copy := GraphModel.add_node(g, n["type"], n["position"] + offset, n["config"])
		copy["title"] = n["title"]
		renamed[n["id"]] = copy["id"]
	for e in example["edges"]:
		GraphModel.connect_nodes(g, renamed[e["from"]], e["from_port"], renamed[e["to"]], e["to_port"])
	_commit(true)
	status_message.emit("Inserted the example. It uses the tags Treated and Untreated; add them to your samples on the Files tab.")


func _on_connection_request(from_node: StringName, from_port: int, to_node: StringName, to_port: int) -> void:
	var problem := GraphModel.connect_nodes(Workspace.model.graph(), from_node, from_port, to_node, to_port)
	if problem:
		status_message.emit(problem)
		return
	_sync_connections()
	_commit()


func _on_disconnection_request(from_node: StringName, from_port: int, to_node: StringName, to_port: int) -> void:
	GraphModel.disconnect_nodes(Workspace.model.graph(), from_node, from_port, to_node, to_port)
	_sync_connections()
	_commit()


func _on_delete_nodes(names: Array) -> void:
	for node_name in names:
		_remove_node(str(node_name))


func _remove_node(node_id: String) -> void:
	GraphModel.remove_node(Workspace.model.graph(), node_id)
	var view: Dictionary = _views.get(node_id, {})
	if view:
		_tag_pickers.erase(view.get("tags"))
		_graph_edit.remove_child(view["node"])
		view["node"].queue_free()
		_views.erase(node_id)
	_sync_connections()
	_commit()


func _on_end_node_move() -> void:
	for node_id in _views:
		GraphModel.node(Workspace.model.graph(), node_id)["position"] = _views[node_id]["node"].position_offset
	Workspace.graph_edited()


## Canvas position at the middle of the visible area, nudged so repeated
## additions do not stack exactly.
func _visible_center() -> Vector2:
	var center := (_graph_edit.scroll_offset + _graph_edit.size / 2) / _graph_edit.zoom - Vector2(125, 80)
	var nudge := Vector2.ZERO
	while Workspace.model.graph()["nodes"].any(func(n): return n["position"].distance_to(center + nudge) < 20):
		nudge += Vector2(30, 30)
	return center + nudge


func _on_popup_request(at_position: Vector2) -> void:
	_context_position = (at_position + _graph_edit.scroll_offset) / _graph_edit.zoom
	_context_menu.position = Vector2i(get_screen_position() + at_position + Vector2(0, 30))
	_context_menu.popup()


func _sync_connections() -> void:
	_graph_edit.clear_connections()
	for e in Workspace.model.graph()["edges"]:
		_graph_edit.connect_node(e["from"], e["from_port"], e["to"], e["to_port"])


# --- Node views --------------------------------------------------------------

func _create_view(n: Dictionary) -> void:
	var spec: Dictionary = GraphModel.TYPES[n["type"]]
	var config: Dictionary = n["config"]
	var gn := GraphNode.new()
	gn.name = n["id"]
	gn.title = n["title"]
	gn.position_offset = n["position"]
	gn.tooltip_text = spec["help"]
	gn.custom_minimum_size.x = 250
	_apply_opaque_style(gn)
	var close := Button.new()
	close.text = "×"
	close.flat = true
	close.tooltip_text = "Remove this node"
	close.pressed.connect(func() -> void: _remove_node(n["id"]))
	gn.get_titlebar_hbox().add_child(close)
	var view := {"node": gn}
	var rows: Array = []  # Controls in slot order

	match n["type"]:
		"data_provider":
			var tags := CheckListButton.new("All samples")
			tags.set_selected(config["criteria"].get("tags_any", config["criteria"].get("tags_all", [])))
			var match_mode := _option(["with any of these tags", "with all of these tags"], 1 if config["criteria"].has("tags_all") else 0)
			var apply_criteria := func() -> void:
				var chosen := tags.get_selected()
				config["criteria"] = {} if chosen.is_empty() else ({"tags_all": chosen} if match_mode.selected == 1 else {"tags_any": chosen})
				_commit()
			tags.selection_changed.connect(func(_s) -> void: apply_criteria.call())
			match_mode.item_selected.connect(func(_i) -> void: apply_criteria.call())
			var pass_only := CheckBox.new()
			pass_only.text = "Only calls with FILTER = PASS"
			pass_only.button_pressed = config.get("pass_only", false)
			pass_only.toggled.connect(func(on: bool) -> void:
				config["pass_only"] = on
				_commit())
			rows = [_labeled("Samples", tags), match_mode, pass_only]
			view["tags"] = tags
			_tag_pickers.append(tags)
		"tag_filter":
			var tags := CheckListButton.new("Choose tags")
			tags.set_selected(config.get("tags", []))
			tags.selection_changed.connect(func(selected: Array) -> void:
				config["tags"] = selected
				_commit())
			var match_mode := _option(["a carrier has any of them", "every tag is on some carrier"], 1 if config.get("match") == "all" else 0)
			match_mode.item_selected.connect(func(i: int) -> void:
				config["match"] = "all" if i == 1 else "any"
				_commit())
			var restrict := CheckBox.new()
			restrict.text = "Keep only these samples"
			restrict.tooltip_text = "Also narrow the samples (and the presence denominator) to samples with these tags."
			restrict.button_pressed = config.get("restrict_samples", false)
			restrict.toggled.connect(func(on: bool) -> void:
				config["restrict_samples"] = on
				_commit())
			rows = [_labeled("Tags", tags), match_mode, restrict]
			view["tags"] = tags
			_tag_pickers.append(tags)
		"comparison":
			var include := Label.new()
			include.text = "Include"
			var exclude := Label.new()
			exclude.text = "Exclude"
			exclude.add_theme_color_override("font_color", EXCLUDE_COLOR)
			var combine := _option(["Include: in every input", "Include: in any input"], 1 if config.get("combine") == "union" else 0)
			combine.item_selected.connect(func(i: int) -> void:
				config["combine"] = "union" if i == 1 else "intersect"
				_commit())
			rows = [include, exclude, combine]
		"logical_filter":
			var summary := Label.new()
			summary.text = FilterDialog.describe(config.get("filters", {}))
			summary.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
			summary.custom_minimum_size.x = 230
			var edit := Button.new()
			edit.text = "Edit filters…"
			edit.pressed.connect(func() -> void:
				for c in _filter_dialog.saved.get_connections():
					_filter_dialog.saved.disconnect(c["callable"])
				_filter_dialog.saved.connect(func(new_config: Dictionary) -> void:
					config["filters"] = new_config
					summary.text = FilterDialog.describe(new_config)
					_commit(), CONNECT_ONE_SHOT)
				_filter_dialog.edit(config.get("filters", {})))
			var mode := _option(["Remove failing mutations", "Keep all, add a pass/fail column"], 1 if config.get("mode") == "flag" else 0)
			var flag_name := LineEdit.new()
			flag_name.placeholder_text = "column name"
			flag_name.text = config.get("flag_name", "passes")
			flag_name.visible = config.get("mode") == "flag"
			mode.item_selected.connect(func(i: int) -> void:
				config["mode"] = "flag" if i == 1 else "drop"
				flag_name.visible = i == 1
				_commit())
			flag_name.text_changed.connect(func(text: String) -> void:
				config["flag_name"] = text
				_commit())
			rows = [summary, edit, mode, flag_name]
		"enrich":
			var note := Label.new()
			note.text = "Annotation lookups arrive in a later\nphase. Mutations pass through unchanged."
			rows = [note]
		"table":
			var table_name := LineEdit.new()
			table_name.placeholder_text = "Table name"
			table_name.text = config.get("name", "")
			table_name.text_changed.connect(func(text: String) -> void:
				config["name"] = text
				Workspace.graph_edited())
			var open := Button.new()
			open.text = "View table"
			open.pressed.connect(func() -> void:
				var result: Dictionary = _results.get(n["id"], {})
				if result.get("frame") == null:
					status_message.emit(result.get("error", "Nothing to show yet."))
					return
				_table_dialog.show_frame(result["frame"], config.get("name", "")))
			rows = [table_name, open]

	var status := Label.new()
	status.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	status.custom_minimum_size.x = 230
	rows.append(status)
	view["status"] = status
	for row in rows:
		gn.add_child(row)
	var inputs: Array = spec["inputs"]
	for i in inputs.size():
		gn.set_slot_enabled_left(i, true)
		gn.set_slot_type_left(i, 0)
		gn.set_slot_color_left(i, EXCLUDE_COLOR if inputs[i]["name"] == "Exclude" else FRAME_COLOR)
	if spec["outputs"] > 0:
		gn.set_slot_enabled_right(0, true)
		gn.set_slot_type_right(0, 0)
		gn.set_slot_color_right(0, FRAME_COLOR)
	_graph_edit.add_child(gn)
	_views[n["id"]] = view


func _refresh_tag_lists() -> void:
	_tags_revision = Workspace.data_revision
	var tags: Array = SampleQuery.tag_counts(Workspace.model.store).map(func(pair): return pair[0])
	tags.erase(SampleQuery.UNTAGGED)
	for picker in _tag_pickers:
		var items := tags.duplicate()
		for selected in picker.get_selected():
			if not selected in items:
				items.append(selected)  # keep tags not present in the data visible
		picker.set_items(items)


func _evaluate() -> void:
	if Workspace.importing:
		return
	var started := Time.get_ticks_msec()
	_results = _evaluator.evaluate(Workspace.model.graph(), Workspace.model.store, Workspace.data_revision)
	for node_id in _views:
		var result: Dictionary = _results.get(node_id, {})
		var status: Label = _views[node_id]["status"]
		if result.get("error"):
			status.text = result["error"]
			status.add_theme_color_override("font_color", ERROR_COLOR)
		else:
			status.text = result.get("summary", "")
			status.add_theme_color_override("font_color", OK_COLOR)
	var elapsed := Time.get_ticks_msec() - started
	if elapsed > 500:
		status_message.emit("Graph updated in %.1f s." % (elapsed / 1000.0))


static func _option(items: Array, selected: int) -> OptionButton:
	var option := OptionButton.new()
	for item in items:
		option.add_item(item)
	option.select(selected)
	return option


static func _labeled(text: String, control: Control) -> HBoxContainer:
	var box := HBoxContainer.new()
	var label := Label.new()
	label.text = text
	box.add_child(label)
	control.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	box.add_child(control)
	return box


## The default node panel is translucent, which makes overlapping nodes hard
## to read; use an opaque panel with the same shape.
static func _apply_opaque_style(gn: GraphNode) -> void:
	for style_name in ["panel", "panel_selected"]:
		var base := gn.get_theme_stylebox(style_name)
		if base is StyleBoxFlat:
			var opaque: StyleBoxFlat = base.duplicate()
			opaque.bg_color.a = 1.0
			gn.add_theme_stylebox_override(style_name, opaque)
