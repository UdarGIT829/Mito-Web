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
var _cmp_group: OptionButton
var _cmp_a: OptionButton
var _cmp_b: OptionButton
var _cmp_unit: OptionButton
var _cmp_perms: SpinBox
var _cmp_seed: SpinBox
var _cmp_run: Button
var _cmp_scan: Button
var _cmp_cancel: Button
var _cmp_info: Label
var _test: PermutationTest
var _result := {}  ## last test: {"windows_key", "group", "tag_a", "tag_b", ...}
var _scan_dialog: AcceptDialog
var _scan_tree: Tree
var _scan_count: SpinBox
var _scan_candidates: Array = []
var _scan_order: Array = []
var _scan_result := {}
const BUDGET := 120_000_000  # labellings x units x windows


func _ready() -> void:
	var bar := HFlowContainer.new()  # wraps instead of widening the window
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
	_info.custom_minimum_size.x = 220
	_info.text_overrun_behavior = TextServer.OVERRUN_TRIM_ELLIPSIS
	bar.add_child(_info)
	var download := Button.new()
	download.text = "Download CSV"
	download.tooltip_text = "Window totals per group, as shown."
	download.pressed.connect(_download)
	bar.add_child(download)
	var matrix := Button.new()
	matrix.text = "Download matrix"
	matrix.tooltip_text = "Samples x windows: one row per sample (all groups), every tag as a 0/1 column, then the chosen count per window. Input for models outside this tool."
	matrix.pressed.connect(_download_matrix)
	bar.add_child(matrix)
	_build_compare_bar()

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
		if name.is_empty():
			name = GraphModel.custom_title(n)
		_picker.add_item(name if name else "Statistics %s" % n["id"].trim_prefix("n"))
		_picker.set_item_tooltip(i, str(n.get("description", "")))
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
	_fill_groups()
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
	# Regions and imported lists are drawn as one segment per window: joining
	# non-adjacent windows with lines would imply data between them.
	var regions: bool = spec.get("mode", "sliding") != "sliding"
	var lo: int = 1 if regions else int(spec.get("range_start", 1))
	var hi: int = WindowStats.LENGTH if regions else int(spec.get("range_end", WindowStats.LENGTH))
	var series := _series()
	_chart.set_data(_stats["windows"], series, regions, _y_title(), lo, hi)
	var windows: Array = _stats["windows"]
	match spec.get("mode", "sliding"):
		"custom": _info.text = "%d windows from a list" % windows.size()
		"regions": _info.text = "%d genes and regions" % windows.size()
		_: _info.text = "%d windows of %d bp, every %d bp" % [windows.size(), int(spec.get("size", 500)), int(spec.get("step", 500))]
	var tested := _result_applies()
	if tested:
		var marks := PackedByteArray()
		for p_adj in _result["p_adjusted"]:
			marks.append(1 if p_adj < 0.05 else 0)
		_chart.highlight = marks
		_chart.queue_redraw()
	_tree.clear()
	_tree.columns = 3 + series.size() + (5 if tested else 0)
	var titles := ["Window", "Positions", "Length (bp)"]
	for s in series:
		titles.append(s["label"])
	if tested:
		titles.append_array(["Mean %s" % _result["tag_a"], "Mean %s" % _result["tag_b"], "A minus B", "p", "Adjusted p"])
	for i in titles.size():
		_tree.set_column_title(i, titles[i])
		_tree.set_column_title_tooltip_text(i, titles[i])
		_tree.set_column_expand(i, false)
		_tree.set_column_clip_content(i, true)
		_tree.set_column_custom_minimum_width(i, [130, 130, 90][i] if i < 3 else 150)
	var root := _tree.create_item()
	for wi in windows.size():
		var item := _tree.create_item(root)
		var w: Dictionary = windows[wi]
		item.set_text(0, w["label"])
		item.set_text(1, ", ".join(w["segments"].map(func(s): return "%d–%d" % [s[0], s[1]])))
		item.set_text(2, str(w["length"]))
		for si in series.size():
			item.set_text(3 + si, DensityChart._fmt(series[si]["values"][wi]))
		if tested:
			var c := 3 + series.size()
			for key in ["mean_a", "mean_b", "difference", "p", "p_adjusted"]:
				item.set_text(c, DensityChart._fmt(_result[key][wi]))
				c += 1
			if _result["p_adjusted"][wi] < 0.05:
				item.set_custom_color(c - 1, Color(0.7, 0.95, 0.75))


func _download() -> void:
	if _stats.is_empty():
		return
	var series := _series()
	var windows: Array = _stats["windows"]
	var tested := _result_applies()
	var header := ["window", "positions", "length_bp"]
	for s in series:
		header.append("%s [%s]" % [s["label"], _y_title()])
	if tested:
		header.append_array(["mean %s" % _result["tag_a"], "mean %s" % _result["tag_b"], "A minus B", "p", "adjusted p (max-T)",
			"test: %s, %s; %s per %s in %s" % [_result["basis"], _result["unit"], _result["measure"], "kb" if _result["per_kb"] else "window", _result["group"]]])
	var lines := PackedStringArray([",".join(header.map(func(h): return TableViewModel._quote(h, ",")))])
	for wi in windows.size():
		var w: Dictionary = windows[wi]
		var row := [w["label"], " ".join(w["segments"].map(func(s): return "%d-%d" % [s[0], s[1]])), str(w["length"])]
		for s in series:
			row.append(str(s["values"][wi]))
		if tested:
			for key in ["mean_a", "mean_b", "difference", "p", "p_adjusted"]:
				row.append(str(_result[key][wi]))
			row.append("")
		lines.append(",".join(row.map(func(c): return TableViewModel._quote(c, ","))))
	var n := GraphModel.node(Workspace.model.graph(), _node_id)
	var file_name: String = (n["config"].get("name", "") if n["config"].get("name", "") else "statistics").validate_filename().replace(" ", "_") + ".csv"
	WebFiles.save_file(file_name, ("\n".join(lines) + "\n").to_utf8_buffer(), "text/csv")
	status_message.emit("Downloaded %s." % file_name)


static func _label(text: String) -> Label:
	var l := Label.new()
	l.text = text
	return l


# --- Matrix export ------------------------------------------------------------

func _matrix_measure() -> String:
	var m: String = WindowStats.MEASURES[_measure.selected]
	return "calls" if m == "mutations" else m


func _download_matrix() -> void:
	if _stats.is_empty():
		return
	var groups: Array = []
	for g in _stats["groups"]:
		groups.append({"label": g["label"], "matrix": SampleMatrix.build(g["frame"], _stats["windows"], _matrix_measure(), _per_kb.button_pressed)})
	var csv := SampleMatrix.to_csv(Workspace.model.store, groups, _stats["windows"], _matrix_measure(), _per_kb.button_pressed)
	var n := GraphModel.node(Workspace.model.graph(), _node_id)
	var base: String = n["config"].get("name", "") if n["config"].get("name", "") else "statistics"
	WebFiles.save_file(base.validate_filename().replace(" ", "_") + "_matrix.csv", csv.to_utf8_buffer(), "text/csv")
	var rows := 0
	for g in groups:
		rows += g["matrix"]["sample_ids"].size()
	status_message.emit("Downloaded a %d x %d matrix (%s%s)." % [rows, _stats["windows"].size(), _matrix_measure(), " per kb" if _per_kb.button_pressed else ""])


# --- Comparison (permutation test) ---------------------------------------------

func _build_compare_bar() -> void:
	var row := HFlowContainer.new()
	add_child(row)
	move_child(row, 1)
	row.add_child(_label("Compare within"))
	_cmp_group = OptionButton.new()
	_cmp_group.custom_minimum_size.x = 150
	_cmp_group.item_selected.connect(func(_i) -> void: _fill_tags())
	row.add_child(_cmp_group)
	row.add_child(_label("tag A"))
	_cmp_a = OptionButton.new()
	row.add_child(_cmp_a)
	row.add_child(_label("vs B"))
	_cmp_b = OptionButton.new()
	row.add_child(_cmp_b)
	_cmp_unit = OptionButton.new()
	_cmp_unit.add_item("by subject")
	_cmp_unit.add_item("by sample")
	_cmp_unit.tooltip_text = "Shuffle subjects when one subject has several samples (e.g. a cell line prepared several ways); those samples are not independent."
	row.add_child(_cmp_unit)
	row.add_child(_label("permutations"))
	_cmp_perms = SpinBox.new()
	_cmp_perms.min_value = 99
	_cmp_perms.max_value = 99999
	_cmp_perms.value = 999
	row.add_child(_cmp_perms)
	row.add_child(_label("seed"))
	_cmp_seed = SpinBox.new()
	_cmp_seed.min_value = 1
	_cmp_seed.max_value = 999999
	_cmp_seed.value = 1
	row.add_child(_cmp_seed)
	_cmp_run = Button.new()
	_cmp_run.text = "Run test"
	_cmp_run.tooltip_text = "Per-window difference in mean (A minus B) with a permutation p-value and a family-wise adjusted p (max-T)."
	_cmp_run.pressed.connect(func() -> void: _run_test(false))
	row.add_child(_cmp_run)
	_cmp_scan = Button.new()
	_cmp_scan.text = "Scan for windows…"
	_cmp_scan.tooltip_text = "Tests windows of 100, 200, 500, 1000 and 2000 bp at quarter-size steps around the genome, adjusting for the whole scan, and lets you keep the best non-overlapping ones."
	_cmp_scan.pressed.connect(func() -> void: _run_test(true))
	row.add_child(_cmp_scan)
	_cmp_cancel = Button.new()
	_cmp_cancel.text = "Cancel"
	_cmp_cancel.visible = false
	_cmp_cancel.pressed.connect(func() -> void:
		if _test:
			_test.cancelled = true)
	row.add_child(_cmp_cancel)
	_cmp_info = _label("")
	_cmp_info.custom_minimum_size.x = 300
	_cmp_info.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_cmp_info.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	row.add_child(_cmp_info)
	_build_scan_dialog()


func _fill_groups() -> void:
	var previous := _cmp_group.get_item_text(_cmp_group.selected) if _cmp_group.item_count else ""
	_cmp_group.clear()
	for g in _stats.get("groups", []):
		_cmp_group.add_item(g["label"])
		if g["label"] == previous:
			_cmp_group.select(_cmp_group.item_count - 1)
	_fill_tags()


func _fill_tags() -> void:
	var tags: Array = []
	var groups: Array = _stats.get("groups", [])
	if _cmp_group.selected >= 0 and _cmp_group.selected < groups.size():
		var frame: MutationFrame = groups[_cmp_group.selected]["frame"]
		for sample_id in frame.samples:
			for tag in frame.store.sample_by_id(sample_id).get("tags", PackedStringArray()):
				if not tag in tags:
					tags.append(tag)
	tags.sort()
	for option in [_cmp_a, _cmp_b]:
		var previous: String = option.get_item_text(option.selected) if option.item_count and option.selected >= 0 else ""
		option.clear()
		for tag in tags:
			option.add_item(tag)
			if tag == previous:
				option.select(option.item_count - 1)
	if _cmp_b.item_count > 1 and _cmp_b.selected == _cmp_a.selected:
		_cmp_b.select((_cmp_a.selected + 1) % _cmp_b.item_count)
	_cmp_run.disabled = tags.size() < 2
	_cmp_scan.disabled = tags.size() < 2


func _run_test(scan: bool) -> void:
	if _stats.is_empty() or _cmp_a.item_count < 2 or (_test and not _test.cancelled and _cmp_cancel.visible):
		return
	var group: Dictionary = _stats["groups"][_cmp_group.selected]
	var windows: Array = WindowStats.scan_windows([100, 200, 500, 1000, 2000]) if scan else _stats["windows"]
	var matrix := SampleMatrix.build(group["frame"], windows, _matrix_measure(), _per_kb.button_pressed)
	var tag_a := _cmp_a.get_item_text(_cmp_a.selected)
	var tag_b := _cmp_b.get_item_text(_cmp_b.selected)
	var by_subject := _cmp_unit.selected == 0
	var prepared := PermutationTest.prepare(Workspace.model.store, matrix, tag_a, tag_b, by_subject)
	if not prepared["ok"]:
		_cmp_info.text = prepared["error"]
		return
	var perms := int(_cmp_perms.value)
	var labellings := minf(PermutationTest.labelling_count(prepared["n_a"], prepared["n_b"]), perms + 1)
	var work: float = labellings * prepared["values"].size() * windows.size()
	if work > BUDGET:
		_cmp_info.text = "Too much work for the browser (%.0f million steps); use fewer permutations or windows." % (work / 1e6)
		return
	_test = PermutationTest.new()
	_test.progress.connect(func(done: int, total: int) -> void:
		_cmp_info.text = "Testing… %d of %d labellings" % [done, total])
	_cmp_cancel.visible = true
	_cmp_run.disabled = true
	_cmp_scan.disabled = true
	var result: Dictionary = await _test.run(prepared, perms, int(_cmp_seed.value))
	_cmp_cancel.visible = false
	_cmp_run.disabled = false
	_cmp_scan.disabled = false
	if not result["ok"]:
		_cmp_info.text = result["error"]
		return
	var unit := "subjects" if by_subject else "samples"
	var basis := "all %d labellings (exact)" % result["labellings"] if result["exact"] else "%d random labellings, seed %d" % [result["labellings"] - 1, int(_cmp_seed.value)]
	var floor_note := "; smallest possible p = %s" % DensityChart._fmt(result["min_p"]) if result["min_p"] > 0.01 else ""
	_cmp_info.text = "%s (%d %s) vs %s (%d): %s%s" % [tag_a, prepared["n_a"], unit, tag_b, prepared["n_b"], basis, floor_note]
	var meta := {"group": group["label"], "tag_a": tag_a, "tag_b": tag_b, "unit": unit, "basis": basis, "measure": _matrix_measure(), "per_kb": _per_kb.button_pressed}
	if scan:
		_show_scan(windows, result, meta)
	else:
		_result = result.duplicate()
		_result.merge(meta)
		_result["windows"] = windows
		_render()


func _result_applies() -> bool:
	return not _result.is_empty() and is_same(_result.get("windows"), _stats.get("windows")) \
		and _result.get("measure") == _matrix_measure() and _result.get("per_kb") == _per_kb.button_pressed


# --- Scan ------------------------------------------------------------------------

func _build_scan_dialog() -> void:
	_scan_dialog = AcceptDialog.new()
	_scan_dialog.title = "Window scan"
	_scan_dialog.ok_button_text = "Close"
	_scan_dialog.min_size = Vector2i(820, 480)
	var box := VBoxContainer.new()
	_scan_dialog.add_child(box)
	var note := _label("Every candidate window is tested; adjusted p accounts for the whole scan (max-T), so a window that survives it is unlikely to be a chance find.")
	note.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	note.custom_minimum_size.x = 780  # a wrapping label needs a width, or it sizes itself one word per line
	box.add_child(note)
	_scan_tree = Tree.new()
	_scan_tree.hide_root = true
	_scan_tree.column_titles_visible = true
	_scan_tree.size_flags_vertical = Control.SIZE_EXPAND_FILL
	_scan_tree.custom_minimum_size.y = 320
	box.add_child(_scan_tree)
	var row := HBoxContainer.new()
	box.add_child(row)
	row.add_child(_label("Keep the best"))
	_scan_count = SpinBox.new()
	_scan_count.min_value = 1
	_scan_count.max_value = 200
	_scan_count.value = 10
	row.add_child(_scan_count)
	row.add_child(_label("non-overlapping windows"))
	var use := Button.new()
	use.text = "Use as this node's windows"
	use.pressed.connect(_use_scan_windows)
	row.add_child(use)
	add_child(_scan_dialog)


func _show_scan(candidates: Array, result: Dictionary, meta: Dictionary) -> void:
	_scan_candidates = candidates
	_scan_result = result
	_scan_order = range(candidates.size())
	_scan_order.sort_custom(func(a: int, b: int) -> bool:
		if result["p_adjusted"][a] != result["p_adjusted"][b]:
			return result["p_adjusted"][a] < result["p_adjusted"][b]
		if result["p"][a] != result["p"][b]:
			return result["p"][a] < result["p"][b]
		return absf(result["difference"][a]) > absf(result["difference"][b]))
	_scan_tree.clear()
	var titles := ["Window", "Size", "Mean %s" % meta["tag_a"], "Mean %s" % meta["tag_b"], "A minus B", "p", "Adjusted p"]
	_scan_tree.columns = titles.size()
	for i in titles.size():
		_scan_tree.set_column_title(i, titles[i])
		_scan_tree.set_column_expand(i, i == 0)
		_scan_tree.set_column_clip_content(i, true)
		_scan_tree.set_column_custom_minimum_width(i, [150, 60, 110, 110, 100, 80, 100][i])
	var root := _scan_tree.create_item()
	for i in _scan_order.slice(0, 200):
		var w: Dictionary = candidates[i]
		var item := _scan_tree.create_item(root)
		var cells := [w["label"], str(w["length"]), DensityChart._fmt(result["mean_a"][i]), DensityChart._fmt(result["mean_b"][i]),
			DensityChart._fmt(result["difference"][i]), DensityChart._fmt(result["p"][i]), DensityChart._fmt(result["p_adjusted"][i])]
		for c in cells.size():
			item.set_text(c, cells[c])
	_scan_dialog.title = "Window scan: %s vs %s in %s (%d candidates; %s)" % [meta["tag_a"], meta["tag_b"], meta["group"], candidates.size(), meta["basis"]]
	_scan_dialog.popup_centered_clamped(Vector2i(900, 600), 0.9)


func _use_scan_windows() -> void:
	var picked := WindowStats.pick_non_overlapping(_scan_candidates, _scan_order, int(_scan_count.value))
	var intervals: Array = []
	for i in picked:
		var w: Dictionary = _scan_candidates[i]
		var segs: Array = w["segments"]
		intervals.append(["scan %s (p adj %s)" % [w["label"], DensityChart._fmt(_scan_result["p_adjusted"][i])], segs[0][0], segs[-1][1]])
	intervals.sort_custom(func(a, b): return a[1] < b[1])
	var n := GraphModel.node(Workspace.model.graph(), _node_id)
	n["config"]["window"]["mode"] = "custom"
	n["config"]["window"]["intervals"] = intervals
	_scan_dialog.hide()
	Workspace.graph_edited()
	Workspace.graph_rebuild_requested.emit()
	status_message.emit("This Statistics node now uses %d scanned windows (Imported list mode on the Analysis tab)." % intervals.size())
