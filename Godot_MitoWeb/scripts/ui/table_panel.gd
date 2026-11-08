class_name TablePanel
extends VBoxContainer
## The Tables tab: shows the mutations arriving at a Table node of the
## analysis graph. Click a column title to sort or filter that column;
## double-click a row to see the calls behind it.

signal status_message(text: String)

var _picker: OptionButton
var _per_sample: CheckBox
var _info: Label
var _filters_label: Label
var _filter_bar: HBoxContainer
var _clear_filters: Button
var _tree: Tree
var _page_label: Label
var _prev: Button
var _next: Button
var _page_size: OptionButton
var _column_popup: PopupPanel
var _details: AcceptDialog
var _details_tree: Tree

var _node_ids: Array = []
var _node_id := ""
var _model: TableViewModel
var _model_frame: MutationFrame
var _dirty := true
var _lookup: Button
var _ensembl: EnsemblClient
var _model_annotations_revision := -1


func _ready() -> void:
	var bar := HBoxContainer.new()
	add_child(bar)
	bar.add_child(_label("Table"))
	_picker = OptionButton.new()
	_picker.custom_minimum_size.x = 220
	_picker.item_selected.connect(func(i: int) -> void: show_node(_node_ids[i]))
	bar.add_child(_picker)
	_per_sample = CheckBox.new()
	_per_sample.text = "One AF column per sample"
	_per_sample.toggled.connect(func(on: bool) -> void:
		if _model:
			_model.set_per_sample_columns(on)
			_render())
	bar.add_child(_per_sample)
	_info = _label("")
	_info.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_info.text_overrun_behavior = TextServer.OVERRUN_TRIM_ELLIPSIS
	bar.add_child(_info)
	_lookup = _button("Look up missing consequences", _lookup_missing)
	_lookup.tooltip_text = "Asks the Ensembl VEP service about mutations in this table that the bundled gnomAD data does not cover. Results are saved in the workspace."
	bar.add_child(_lookup)
	bar.add_child(_button("Copy all (TSV)", func() -> void:
		if _model:
			DisplayServer.clipboard_set(_model.export_text("\t"))
			status_message.emit("Copied %d rows to the clipboard." % _model.view.size())))
	bar.add_child(_button("Download CSV", func() -> void:
		if _model:
			var file_name := (_current_name() if _current_name() else "table").validate_filename().replace(" ", "_") + ".csv"
			WebFiles.save_file(file_name, _model.export_text(",").to_utf8_buffer(), "text/csv")
			status_message.emit("Downloaded %s (%d rows)." % [file_name, _model.view.size()])))

	_filter_bar = HBoxContainer.new()
	var filter_bar := _filter_bar
	add_child(filter_bar)
	_filters_label = _label("")
	_filters_label.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_filters_label.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	_filters_label.add_theme_color_override("font_color", Color(0.75, 0.85, 1.0))
	filter_bar.add_child(_filters_label)
	_clear_filters = _button("Clear column filters", func() -> void:
		_model.clear_filters()
		_render())
	filter_bar.add_child(_clear_filters)

	_tree = Tree.new()
	_tree.size_flags_vertical = Control.SIZE_EXPAND_FILL
	_tree.hide_root = true
	_tree.column_titles_visible = true
	_tree.select_mode = Tree.SELECT_ROW
	_tree.column_title_clicked.connect(_on_title_clicked)
	_tree.item_activated.connect(_on_row_activated)
	add_child(_tree)

	var pager := HBoxContainer.new()
	add_child(pager)
	_prev = _button("< Previous", func() -> void:
		_model.page -= 1
		_render())
	pager.add_child(_prev)
	_page_label = _label("")
	pager.add_child(_page_label)
	_next = _button("Next >", func() -> void:
		_model.page += 1
		_render())
	pager.add_child(_next)
	pager.add_child(_label("   Rows per page"))
	_page_size = OptionButton.new()
	for size in TableViewModel.PAGE_SIZES:
		_page_size.add_item(str(size))
	_page_size.select(1)
	_page_size.item_selected.connect(func(i: int) -> void:
		if _model:
			_model.page_size = TableViewModel.PAGE_SIZES[i]
			_model.page = 0
			_render())
	pager.add_child(_page_size)
	var hint := _label("   Click a column title to sort or filter. Double-click a row for its calls.")
	hint.add_theme_color_override("font_color", Color(0.7, 0.7, 0.7))
	pager.add_child(hint)

	_column_popup = PopupPanel.new()
	var popup_style := StyleBoxFlat.new()
	popup_style.bg_color = Color(0.13, 0.14, 0.17)
	popup_style.border_color = Color(0.35, 0.4, 0.5)
	popup_style.set_border_width_all(1)
	popup_style.set_corner_radius_all(4)
	popup_style.set_content_margin_all(10)
	_column_popup.add_theme_stylebox_override("panel", popup_style)
	add_child(_column_popup)
	_details = AcceptDialog.new()
	_details.title = "Calls"
	_details.ok_button_text = "Close"
	_details.min_size = Vector2i(900, 420)
	_details_tree = Tree.new()
	_details_tree.hide_root = true
	_details_tree.column_titles_visible = true
	_details.add_child(_details_tree)
	add_child(_details)

	_ensembl = EnsemblClient.new()
	add_child(_ensembl)
	_ensembl.progress.connect(func(done: int, total: int) -> void:
		status_message.emit("Ensembl: %d of %d looked up…" % [done, total]))
	_ensembl.finished.connect(func(results: Dictionary, error: String) -> void:
		if results:
			Workspace.add_live_annotations(results)
		var found := results.values().filter(func(r): return not r["consequence"].is_empty()).size()
		status_message.emit(("Ensembl: %d looked up, %d with consequences. " % [results.size(), found]) + error)
		_lookup.disabled = false)
	Workspace.changed.connect(func() -> void:
		_dirty = true
		if is_visible_in_tree():
			_sync.call_deferred())
	visibility_changed.connect(func() -> void:
		if is_visible_in_tree() and _dirty:
			_sync())
	_sync()


## Shows a Table node, e.g. from the graph's "View table" button.
func show_node(node_id: String) -> void:
	_node_id = node_id
	_dirty = true
	_sync()


## Re-reads the graph: the list of Table nodes and the shown table's data.
## Column filters, sort and page survive when only the data changed.
func _sync() -> void:
	_dirty = false
	var graph := Workspace.model.graph()
	var results := Workspace.graph_results()
	_node_ids = graph["nodes"].filter(func(n): return n["type"] == "table").map(func(n): return n["id"])
	if _node_id.is_empty() or not _node_id in _node_ids:
		_node_id = _node_ids[0] if _node_ids else ""
	_picker.clear()
	for i in _node_ids.size():
		var n := GraphModel.node(graph, _node_ids[i])
		var name: String = n["config"].get("name", "")
		if name.is_empty():
			name = GraphModel.custom_title(n)
		_picker.add_item(name if name else "Table %s" % n["id"].trim_prefix("n"))
		_picker.set_item_tooltip(i, str(n.get("description", "")))
		if _node_ids[i] == _node_id:
			_picker.select(i)
	var result: Dictionary = results.get(_node_id, {})
	var frame: MutationFrame = result.get("frame")
	if frame == null:
		_model = null
		_model_frame = null
		_tree.clear()
		_info.text = result.get("error", "Add a Table node on the Analysis tab and connect it to see mutations here.") if _node_id else \
			"Add a Table node on the Analysis tab and connect it to see mutations here."
		_filters_label.text = ""
		_clear_filters.visible = false
		_filter_bar.visible = false
		_page_label.text = ""
		_lookup.visible = false
		return
	if frame != _model_frame or Workspace.annotations_revision != _model_annotations_revision:
		_model_annotations_revision = Workspace.annotations_revision
		var previous := _model
		_model = TableViewModel.new(frame)
		_model_frame = frame
		if previous:
			_model.page_size = previous.page_size
			if frame.samples.size() == previous.frame.samples.size():
				_model.set_per_sample_columns(previous.per_sample_columns)
			for id in previous.filters:
				if not _model.column(id).is_empty():
					_model.filters[id] = previous.filters[id]
			if not _model.column(previous.sort["column"]).is_empty():
				_model.sort = previous.sort
			_model.refresh()
	_lookup.visible = "ensembl" in frame.annotation_sources
	_per_sample.set_pressed_no_signal(_model.per_sample_columns)
	_per_sample.disabled = frame.samples.size() > TableViewModel.PER_SAMPLE_MAX
	_per_sample.tooltip_text = "Not available above %d samples." % TableViewModel.PER_SAMPLE_MAX if _per_sample.disabled else ""
	_render()


func _render() -> void:
	if _model == null:
		return
	var cols := _model.columns
	_tree.clear()
	_tree.columns = cols.size()
	for i in cols.size():
		var title: String = cols[i]["title"]
		if cols[i]["kind"] == "sample_af" and title.length() > 14:
			title = title.left(13) + "…"
		if _model.filters.has(cols[i]["id"]):
			title = "[F] " + title
		if _model.sort["column"] == cols[i]["id"]:
			title += " (desc)" if _model.sort["descending"] else " (asc)"
		_tree.set_column_title(i, title)
		_tree.set_column_title_tooltip_text(i, cols[i]["title"])
		var kind: String = cols[i]["kind"]
		_tree.set_column_expand(i, cols[i]["id"] == "samples" or (cols[i]["id"] == "slippage" and not _model.per_sample_columns))
		_tree.set_column_clip_content(i, true)
		_tree.set_column_custom_minimum_width(i, 220 if cols[i]["id"] == "slippage" else {"sample_af": 90, "flag": 70, "number": 80}.get(kind, 120))
	var root := _tree.create_item()
	for row in _model.page_rows():
		var item := _tree.create_item(root)
		item.set_metadata(0, row)
		for i in cols.size():
			var text := _model.cell_text(row, cols[i])
			item.set_text(i, text)
			if text.length() > 18:
				item.set_tooltip_text(i, text.replace("; ", "\n"))
			if cols[i]["kind"] == "flag":
				item.set_custom_color(i, Color(0.7, 0.95, 0.75) if text == "yes" else Color(1, 0.6, 0.5))
	var total := _model.rows.size()
	var shown := _model.view.size()
	_info.text = "%d mutations across %d samples%s" % [
		shown, _model.frame.samples.size(), "" if shown == total else " (of %d before column filters)" % total]
	_filters_label.text = ("Column filters: " + _model.describe_filters()) if _model.filters else ""
	_clear_filters.visible = not _model.filters.is_empty()
	_filter_bar.visible = not _model.filters.is_empty()
	_page_label.text = "Page %d of %d" % [_model.page + 1, _model.page_count()]
	_prev.disabled = _model.page <= 0
	_next.disabled = _model.page >= _model.page_count() - 1


func _current_name() -> String:
	var n := GraphModel.node(Workspace.model.graph(), _node_id)
	return n.get("config", {}).get("name", "")


# --- Column popup: sort and filter -----------------------------------------

func _on_title_clicked(column_index: int, _mouse_button: int) -> void:
	if _model == null:
		return
	var col: Dictionary = _model.columns[column_index]
	for child in _column_popup.get_children():
		child.queue_free()
	var box := VBoxContainer.new()
	box.custom_minimum_size.x = 280
	_column_popup.add_child(box)
	var title := _label(col["title"])
	title.add_theme_color_override("font_color", Color(0.75, 0.85, 1.0))
	box.add_child(title)
	var sort_row := HBoxContainer.new()
	sort_row.add_child(_button("Sort ascending", func() -> void:
		_model.set_sort(col["id"], false)
		_column_popup.hide()
		_render()))
	sort_row.add_child(_button("Sort descending", func() -> void:
		_model.set_sort(col["id"], true)
		_column_popup.hide()
		_render()))
	box.add_child(sort_row)
	var spec: Dictionary = _model.filters.get(col["id"], {})
	var read_spec: Callable
	match col["kind"]:
		"number":
			var lo := _spin_text(spec.get("min"))
			var hi := _spin_text(spec.get("max"))
			box.add_child(_pair("At least", lo))
			box.add_child(_pair("At most", hi))
			var keep_missing := CheckBox.new()
			keep_missing.text = "Also keep rows with no value"
			keep_missing.tooltip_text = "e.g. variants not observed in gnomAD"
			keep_missing.button_pressed = spec.get("include_missing", false)
			box.add_child(keep_missing)
			read_spec = func() -> Dictionary:
				var out := {}
				if lo.text.strip_edges().is_valid_float():
					out["min"] = lo.text.to_float()
				if hi.text.strip_edges().is_valid_float():
					out["max"] = hi.text.to_float()
				if out and keep_missing.button_pressed:
					out["include_missing"] = true
				return out
		"text", "reason":
			var text := LineEdit.new()
			text.placeholder_text = "contains…"
			text.text = spec.get("contains", "")
			box.add_child(text)
			read_spec = func() -> Dictionary:
				return {"contains": text.text} if text.text.strip_edges() else {}
		"genes":
			var picker := CheckListButton.new("Any gene or region")
			var items: Array = Array(MitoGenes.symbols())
			items.push_front("(none)")
			picker.set_items(items)
			picker.set_selected(spec.get("genes", []))
			box.add_child(picker)
			box.add_child(_label("\"(none)\" = control region / intergenic"))
			read_spec = func() -> Dictionary:
				return {"genes": picker.get_selected()} if picker.get_selected() else {}
		"flag":
			var choice := _choice(["Any", "Only yes", "Only no"], 0 if spec.is_empty() else (1 if spec["value"] else 2))
			box.add_child(choice)
			read_spec = func() -> Dictionary:
				return {} if choice.selected == 0 else {"value": choice.selected == 1}
		"sample_af":
			var modes := ["any", "carried", "not_carried", "range"]
			var choice := _choice(["Any", "Carries the mutation", "Does not carry it", "AF in range"], modes.find(spec.get("mode", "any")))
			var lo := _spin_text(spec.get("min", 0.0))
			var hi := _spin_text(spec.get("max", 1.0))
			box.add_child(choice)
			box.add_child(_pair("AF from", lo))
			box.add_child(_pair("AF to", hi))
			read_spec = func() -> Dictionary:
				var mode: String = modes[choice.selected]
				if mode == "any":
					return {}
				var out := {"mode": mode}
				if mode == "range":
					out["min"] = lo.text.to_float() if lo.text.is_valid_float() else 0.0
					out["max"] = hi.text.to_float() if hi.text.is_valid_float() else 1.0
				return out
	if read_spec.is_valid():
		var actions := HBoxContainer.new()
		actions.add_child(_button("Apply filter", func() -> void:
			_model.set_filter(col["id"], read_spec.call())
			_column_popup.hide()
			_render()))
		actions.add_child(_button("Clear", func() -> void:
			_model.set_filter(col["id"], {})
			_column_popup.hide()
			_render()))
		box.add_child(actions)
	_column_popup.reset_size()
	_column_popup.position = Vector2i(get_viewport().get_mouse_position()) + Vector2i(0, 12)
	_column_popup.popup()


# --- Row details --------------------------------------------------------------

func _on_row_activated() -> void:
	var item := _tree.get_selected()
	if item == null or _model == null:
		return
	var row: Dictionary = item.get_metadata(0)
	var details := _model.row_details(row)
	_details.title = "m.%d %s>%s — %s" % [row["position"], row["event_ref"], row["alt"],
		", ".join(row["genes"]) if row["genes"] else row["region"].replace("_", " ")]
	var headers := ["Sample", "Subject", "AF", "FILTER", "As in the VCF", "Placement", "Multi-ALT row", "Repeats in ALT", "Reference context", "Raw value"]
	_details_tree.clear()
	_details_tree.columns = headers.size()
	var widths := [0, 110, 70, 90, 120, 100, 100, 110, 150, 0]
	for i in headers.size():
		_details_tree.set_column_title(i, headers[i])
		_details_tree.set_column_expand(i, i == 0 or i == 9)
		_details_tree.set_column_clip_content(i, true)
		_details_tree.set_column_custom_minimum_width(i, widths[i] if widths[i] else 160)
	var root := _details_tree.create_item()
	for d in details:
		var meta: Dictionary = d["metadata"]
		var repeats := TableViewModel.describe_repeats(meta)
		var context := "%s | %s" % [meta.get("REFERENCE_6_BEFORE", ""), meta.get("REFERENCE_6_AFTER", "")]
		var values := [d["sample"], d["subject"], "n/a" if d["af"] == null else String.num(d["af"], 4), d["filter"], d["vcf"],
			d["placement"], "yes" if d["multi_alt"] else "no", repeats, context, "%s = %s" % [d["format"], d["value"]]]
		var item2 := _details_tree.create_item(root)
		for i in values.size():
			item2.set_text(i, values[i])
			item2.set_tooltip_text(i, values[i])
	_details.popup_centered_ratio(0.8)
	_details_tree.custom_minimum_size = _details.size - Vector2i(20, 70)


# --- Small helpers --------------------------------------------------------------

static func _label(text: String) -> Label:
	var l := Label.new()
	l.text = text
	return l


static func _button(text: String, action: Callable) -> Button:
	var b := Button.new()
	b.text = text
	b.pressed.connect(action)
	return b


static func _spin_text(value) -> LineEdit:
	var l := LineEdit.new()
	l.text = "" if value == null else str(value)
	l.custom_minimum_size.x = 120
	return l


static func _pair(text: String, control: Control) -> HBoxContainer:
	var box := HBoxContainer.new()
	var l := _label(text)
	l.custom_minimum_size.x = 80
	box.add_child(l)
	control.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	box.add_child(control)
	return box


static func _choice(items: Array, selected: int) -> OptionButton:
	var o := OptionButton.new()
	for item in items:
		o.add_item(item)
	o.select(maxi(selected, 0))
	return o


func _lookup_missing() -> void:
	if _model == null or _ensembl.busy:
		return
	var alleles := _model.alleles_needing_lookup()
	if alleles.is_empty():
		status_message.emit("Every mutation in this table already has a consequence.")
		return
	_lookup.disabled = true
	status_message.emit("Asking Ensembl about %d mutations…" % alleles.size())
	_ensembl.lookup(alleles)
