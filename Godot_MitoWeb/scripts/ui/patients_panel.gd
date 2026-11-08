class_name PatientsPanel
extends VSplitContainer
## Patients (subjects) with their haplogroup, lineage and imported values,
## and the results of Compare Patients nodes as tables shaped like a paper's.

signal status_message(text: String)

const WARN := Color(1, 0.8, 0.4)

var _roster: Tree
var _roster_info: Label
var _picker: OptionButton
var _results: Tree
var _info: Label
var _node_ids: Array = []
var _node_id := ""
var _result: Dictionary = {}


func _ready() -> void:
	split_offset = 260
	var top := VBoxContainer.new()
	add_child(top)
	_roster_info = Label.new()
	_roster_info.text_overrun_behavior = TextServer.OVERRUN_TRIM_ELLIPSIS
	top.add_child(_roster_info)
	_roster = Tree.new()
	_roster.size_flags_vertical = Control.SIZE_EXPAND_FILL
	_roster.hide_root = true
	_roster.column_titles_visible = true
	top.add_child(_roster)

	var bottom := VBoxContainer.new()
	add_child(bottom)
	var bar := HBoxContainer.new()
	bottom.add_child(bar)
	var label := Label.new()
	label.text = "Results"
	label.add_theme_font_size_override("font_size", 16)
	bar.add_child(label)
	_picker = OptionButton.new()
	_picker.custom_minimum_size.x = 260
	_picker.item_selected.connect(func(i: int) -> void:
		_node_id = _node_ids[i]
		_show_result())
	bar.add_child(_picker)
	var export := Button.new()
	export.text = "Export CSV"
	export.pressed.connect(_export)
	bar.add_child(export)
	_info = Label.new()
	_info.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	_info.custom_minimum_size.x = 300
	bottom.add_child(_info)
	_results = Tree.new()
	_results.size_flags_vertical = Control.SIZE_EXPAND_FILL
	_results.hide_root = true
	_results.column_titles_visible = true
	bottom.add_child(_results)
	Workspace.changed.connect(refresh)
	visibility_changed.connect(func() -> void:
		if is_visible_in_tree():
			refresh())
	refresh()


func show_node(node_id: String) -> void:
	_node_id = node_id
	refresh()


func refresh() -> void:
	if Workspace.importing or not is_visible_in_tree():
		return
	_fill_roster()
	var graph := Workspace.model.graph()
	_node_ids = graph["nodes"].filter(func(n): return n["type"] == "patient_compare").map(func(n): return n["id"])
	if _node_id.is_empty() or not _node_id in _node_ids:
		_node_id = _node_ids[0] if _node_ids else ""
	_picker.clear()
	for i in _node_ids.size():
		var n := GraphModel.node(graph, _node_ids[i])
		var name: String = n["config"].get("name", "")
		if name.is_empty():
			name = GraphModel.custom_title(n)
		_picker.add_item(name if name else "Compare Patients %s" % n["id"].trim_prefix("n"))
		_picker.set_item_tooltip(i, str(n.get("description", "")))
		if _node_ids[i] == _node_id:
			_picker.select(i)
	_show_result()


func _fill_roster() -> void:
	var store := Workspace.model.store
	var by_subject := {}
	for sample in store.samples:
		by_subject.get_or_add(sample["subject"], []).append(sample)
	var attributes := PatientCompare.attribute_names(store, false)
	var titles := ["Patient", "Samples", "Haplogroup", "Lineage", "Mixture check"] + attributes
	_roster.clear()
	_roster.columns = titles.size()
	for i in titles.size():
		_roster.set_column_title(i, titles[i])
		_roster.set_column_expand(i, i == 0)
		_roster.set_column_custom_minimum_width(i, 140 if i < 5 else 90)
	var root := _roster.create_item()
	var subjects := by_subject.keys()
	subjects.sort_custom(func(a, b): return str(a).naturalnocasecmp_to(str(b)) < 0)
	var conflicts := 0
	for subject in subjects:
		var samples: Array = by_subject[subject]
		var item := _roster.create_item(root)
		item.set_text(0, subject)
		item.set_text(1, str(samples.size()))
		item.set_tooltip_text(1, "\n".join(samples.map(func(s): return s["name"])))
		var hgs := {}
		var manual := false
		for s in samples:
			var hg: Dictionary = s.get("haplogroup", {})
			if hg:
				hgs[hg["name"]] = true
				manual = manual or hg.get("manual", false)
		item.set_text(2, " / ".join(hgs.keys()) + (" (set by hand)" if manual else ""))
		var lineages := {}
		for s in samples:
			for t in s["tags"]:
				if t.begins_with(Lineages.PREFIX):
					lineages[t.trim_prefix(Lineages.PREFIX)] = true
		item.set_text(3, " / ".join(lineages.keys()))
		if hgs.size() > 1 or lineages.size() > 1:
			conflicts += 1
			for c in [2, 3]:
				item.set_custom_color(c, WARN)
			item.set_tooltip_text(2, "Its samples disagree. Set the haplogroup by hand on the Samples tab (select the samples, Set haplogroup…).")
		var mixed: Array = samples.filter(func(x): return x.get("haplogroup_check", {}).get("flag", false))
		var disagree: Array = samples.filter(func(x): return HaplogroupClassifier.manual_disagreement(x) != "")
		if mixed or disagree:
			var parts: Array = []
			if disagree:
				parts.append("%d file(s) look like another haplogroup" % disagree.size())
			if mixed:
				parts.append("%d mixed" % mixed.size())
			item.set_text(4, ", ".join(parts))
			item.set_custom_color(4, Color(1, 0.6, 0.35))
			var lines: Array = disagree.map(func(x): return "%s looks like %s (set to %s)" % [x["name"], HaplogroupClassifier.manual_disagreement(x), x["haplogroup"]["name"]])
			lines.append_array(mixed.map(func(x): return "%s: %s in %s" % [x["name"], x["haplogroup_check"]["minor"], x["haplogroup_check"]["major"]]))
			item.set_tooltip_text(4, "\n".join(lines))
		elif samples.any(func(x): return x.has("haplogroup_check")):
			item.set_text(4, "no mixture")
		for k in attributes.size():
			var value = null
			for s in samples:
				value = s.get("attributes", {}).get(attributes[k], value)
			if value != null:
				item.set_text(5 + k, PatientCompare._num(value) if typeof(value) == TYPE_FLOAT else str(value))
	_roster_info.text = ("%d patients, %d with values imported%s.%s" % [subjects.size(),
		subjects.filter(func(s): return by_subject[s].any(func(x): return not x.get("attributes", {}).is_empty())).size(),
		"; %d with disagreeing samples (highlighted)" % conflicts if conflicts else "",
		"" if attributes else " To add clinical values, upload the table (CSV) on the Files tab or use Samples > Import from a table."]) if subjects else \
		"No samples yet. Import VCFs on the Files tab, then values per patient with Samples > Import from a table."


func _show_result() -> void:
	_results.clear()
	_result = {}
	if _node_id.is_empty():
		_info.text = "Add a Compare Patients node on the Analysis tab (or insert the Cohort vs 1000 Genomes workflow) and connect groups to it."
		return
	var result: Dictionary = Workspace.graph_results().get(_node_id, {})
	if result.get("error"):
		_info.text = result["error"]
		return
	_result = result.get("patients", {})
	if _result.is_empty():
		_info.text = "Nothing to show yet."
		return
	var groups: Array = _result["groups"]
	var pair_names: Array = []
	for row in _result["rows"]:
		for pr in row["pairwise"]:
			var name := "%s vs %s" % [pr["a"], pr["b"]]
			if not name in pair_names:
				pair_names.append(name)
	var titles := ["Measure"]
	for g in groups:
		titles.append("%s (n = %d)" % [g["label"], g["n"]])
	titles.append_array(["Test", "p"])
	titles.append_array(pair_names)
	_results.columns = titles.size()
	for i in titles.size():
		_results.set_column_title(i, titles[i])
		_results.set_column_title_tooltip_text(i, titles[i])
		_results.set_column_expand(i, i == 0)
		_results.set_column_custom_minimum_width(i, 290 if i == 0 else 120)
	for i in groups.size():
		_results.set_column_title_tooltip_text(1 + i, "%s: %s" % [groups[i]["label"], ", ".join(groups[i]["subjects"])])
	var root := _results.create_item()
	for row in _result["rows"]:
		var item := _results.create_item(root)
		item.set_text(0, row["label"])
		item.set_tooltip_text(0, "%s%s" % [row["label"], ("\n" + row["note"]) if row["note"] else ""])
		for i in row["cells"].size():
			item.set_text(1 + i, row["cells"][i]["text"])
			item.set_tooltip_text(1 + i, "%d patients with a value" % row["cells"][i]["n"])
		var col := 1 + groups.size()
		item.set_text(col, row["test"])
		item.set_text(col + 1, PatientCompare.format_p(row["p"]))
		if row["p"] < 0.05:
			item.set_custom_color(col + 1, Color(0.6, 0.95, 0.7))
		for pr in row["pairwise"]:
			var at := col + 2 + pair_names.find("%s vs %s" % [pr["a"], pr["b"]])
			item.set_text(at, PatientCompare.format_p(pr["p"]))
			if pr["p"] < 0.05:
				item.set_custom_color(at, Color(0.6, 0.95, 0.7))
		if row["note"]:
			item.set_text(col, "%s (%s)" % [row["test"], row["note"]])
	var notes: Array = ["Patients are subjects: a patient's samples count once. Two-sided p-values; Wilcoxon uses the normal approximation with tie correction, as Stata's ranksum. p < 0.05 in green."]
	if _result["warnings"]:
		notes.append("Check: " + " ".join(_result["warnings"].slice(0, 4)) + (" (+%d more)" % (_result["warnings"].size() - 4) if _result["warnings"].size() > 4 else ""))
	_info.text = "\n".join(notes)
	_info.add_theme_color_override("font_color", Color(1, 0.8, 0.4) if _result.get("missing_attributes", []) else Color(0.85, 0.87, 0.9))


func _export() -> void:
	if _result.is_empty():
		status_message.emit("Nothing to export yet.")
		return
	var name := _picker.get_item_text(_picker.selected) if _picker.selected >= 0 else "patients"
	var file := name.validate_filename().replace(" ", "_") + ".csv"
	WebFiles.save_file(file, PatientCompare.to_csv(_result).to_utf8_buffer(), "text/csv")
	status_message.emit("Exported %s." % file)
