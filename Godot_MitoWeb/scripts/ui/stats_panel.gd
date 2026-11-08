class_name StatsPanel
extends VBoxContainer
## The Statistics tab: a Statistics node's window counts as a chart and a
## table, with the measure and scaling chosen here (saved in the node).

signal status_message(text: String)

var _picker: OptionButton
var _measure: OptionButton
var _per_kb: CheckBox
var _per_sample: CheckBox
var _info: Label
var _chart: DensityChart
var _tree: Tree
var _node_ids: Array = []
var _node_id := ""
var _stats := {}
var _dirty := true


func _ready() -> void:
	var bar := HBoxContainer.new()
	add_child(bar)
	bar.add_child(_label("Statistics"))
	_picker = OptionButton.new()
	_picker.custom_minimum_size.x = 200
	_picker.item_selected.connect(func(i: int) -> void: show_node(_node_ids[i]))
	bar.add_child(_picker)
	bar.add_child(_label("  Count"))
	_measure = OptionButton.new()
	for m in WindowStats.MEASURES:
		_measure.add_item(WindowStats.MEASURE_LABELS[m])
	_measure.tooltip_text = "Distinct mutations: each mutation once. Carrier calls: once per sample carrying it. AF-weighted: the sum of carrier AF (heteroplasmy load)."
	_measure.item_selected.connect(func(_i) -> void: _option_changed())
	bar.add_child(_measure)
	_per_kb = CheckBox.new()
	_per_kb.text = "per kb"
	_per_kb.toggled.connect(func(_on) -> void: _option_changed())
	bar.add_child(_per_kb)
	_per_sample = CheckBox.new()
	_per_sample.text = "per sample"
	_per_sample.tooltip_text = "Divide by the number of samples in each group, so groups of different sizes are comparable."
	_per_sample.toggled.connect(func(_on) -> void: _option_changed())
	bar.add_child(_per_sample)
	_info = _label("")
	_info.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_info.text_overrun_behavior = TextServer.OVERRUN_TRIM_ELLIPSIS
	bar.add_child(_info)
	var download := Button.new()
	download.text = "Download CSV"
	download.pressed.connect(_download)
	bar.add_child(download)

	var split := VSplitContainer.new()
	split.size_flags_vertical = Control.SIZE_EXPAND_FILL
	add_child(split)
	_chart = DensityChart.new()
	_chart.size_flags_vertical = Control.SIZE_EXPAND_FILL
	_chart.custom_minimum_size.y = 300
	split.add_child(_chart)
	_tree = Tree.new()
	_tree.hide_root = true
	_tree.column_titles_visible = true
	_tree.custom_minimum_size.y = 160
	split.add_child(_tree)

	Workspace.changed.connect(func() -> void:
		_dirty = true
		if is_visible_in_tree():
			_sync.call_deferred())
	visibility_changed.connect(func() -> void:
		if is_visible_in_tree() and _dirty:
			_sync())
	_sync()


func show_node(node_id: String) -> void:
	_node_id = node_id
	_dirty = true
	_sync()


func _sync() -> void:
	_dirty = false
	var graph := Workspace.model.graph()
	var results := Workspace.graph_results()
	_node_ids = graph["nodes"].filter(func(n): return n["type"] == "statistics").map(func(n): return n["id"])
	if _node_id.is_empty() or not _node_id in _node_ids:
		_node_id = _node_ids[0] if _node_ids else ""
	_picker.clear()
	for i in _node_ids.size():
		var n := GraphModel.node(graph, _node_ids[i])
		var name: String = n["config"].get("name", "")
		_picker.add_item(name if name else "Statistics %s" % n["id"].trim_prefix("n"))
		if _node_ids[i] == _node_id:
			_picker.select(i)
	var result: Dictionary = results.get(_node_id, {})
	_stats = result.get("stats", {})
	if _stats.is_empty():
		_info.text = result.get("error", "Add a Statistics node on the Analysis tab and connect groups to it.") if _node_id else \
			"Add a Statistics node on the Analysis tab and connect groups to it."
		_chart.set_data([], [], false, "", 1, WindowStats.LENGTH)
		_tree.clear()
		return
	var config: Dictionary = GraphModel.node(graph, _node_id)["config"]
	_measure.select(maxi(WindowStats.MEASURES.find(config.get("measure", "calls")), 0))
	_per_kb.set_pressed_no_signal(config.get("per_kb", true))
	_per_sample.set_pressed_no_signal(config.get("per_sample", true))
	_render()


func _option_changed() -> void:
	var n := GraphModel.node(Workspace.model.graph(), _node_id)
	if n.is_empty():
		return
	n["config"]["measure"] = WindowStats.MEASURES[_measure.selected]
	n["config"]["per_kb"] = _per_kb.button_pressed
	n["config"]["per_sample"] = _per_sample.button_pressed
	Workspace.graph_edited()
	_render()


func _series() -> Array:
	var measure: String = WindowStats.MEASURES[_measure.selected]
	var windows: Array = _stats["windows"]
	var out: Array = []
	for g in _stats["groups"]:
		var raw: PackedFloat64Array = g["values"][measure]
		var values := PackedFloat64Array()
		for wi in windows.size():
			values.append(WindowStats.scaled(raw[wi], windows[wi], g["samples"], _per_kb.button_pressed, _per_sample.button_pressed))
		out.append({"label": "%s (n=%d)" % [g["label"], g["samples"]], "short": g["label"], "values": values})
	return out


func _y_title() -> String:
	var title: String = WindowStats.MEASURE_LABELS[WindowStats.MEASURES[_measure.selected]]
	if _per_sample.button_pressed:
		title += " per sample"
	if _per_kb.button_pressed:
		title += " per kb"
	return title


func _render() -> void:
	var config: Dictionary = GraphModel.node(Workspace.model.graph(), _node_id)["config"]
	var spec: Dictionary = config.get("window", {})
	var regions: bool = spec.get("mode") == "regions"
	var lo: int = 1 if regions else int(spec.get("range_start", 1))
	var hi: int = WindowStats.LENGTH if regions else int(spec.get("range_end", WindowStats.LENGTH))
	var series := _series()
	_chart.set_data(_stats["windows"], series, regions, _y_title(), lo, hi)
	var windows: Array = _stats["windows"]
	_info.text = "%d windows%s" % [windows.size(),
		"" if regions else " of %d bp, every %d bp" % [int(spec.get("size", 500)), int(spec.get("step", 500))]]
	_tree.clear()
	_tree.columns = 3 + series.size()
	var titles := ["Window", "Positions", "Length (bp)"]
	for s in series:
		titles.append(s["label"])
	for i in titles.size():
		_tree.set_column_title(i, titles[i])
		_tree.set_column_expand(i, i == 0 or i >= 3)
	var root := _tree.create_item()
	for wi in windows.size():
		var item := _tree.create_item(root)
		var w: Dictionary = windows[wi]
		item.set_text(0, w["label"])
		item.set_text(1, ", ".join(w["segments"].map(func(s): return "%d–%d" % [s[0], s[1]])))
		item.set_text(2, str(w["length"]))
		for si in series.size():
			item.set_text(3 + si, DensityChart._fmt(series[si]["values"][wi]))


func _download() -> void:
	if _stats.is_empty():
		return
	var series := _series()
	var windows: Array = _stats["windows"]
	var header := ["window", "positions", "length_bp"]
	for s in series:
		header.append("%s [%s]" % [s["label"], _y_title()])
	var lines := PackedStringArray([",".join(header.map(func(h): return TableViewModel._quote(h, ",")))])
	for wi in windows.size():
		var w: Dictionary = windows[wi]
		var row := [w["label"], " ".join(w["segments"].map(func(s): return "%d-%d" % [s[0], s[1]])), str(w["length"])]
		for s in series:
			row.append(str(s["values"][wi]))
		lines.append(",".join(row.map(func(c): return TableViewModel._quote(c, ","))))
	var n := GraphModel.node(Workspace.model.graph(), _node_id)
	var file_name: String = (n["config"].get("name", "") if n["config"].get("name", "") else "statistics").validate_filename().replace(" ", "_") + ".csv"
	WebFiles.save_file(file_name, ("\n".join(lines) + "\n").to_utf8_buffer(), "text/csv")
	status_message.emit("Downloaded %s." % file_name)


static func _label(text: String) -> Label:
	var l := Label.new()
	l.text = text
	return l
