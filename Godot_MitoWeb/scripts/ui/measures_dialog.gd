class_name MeasuresDialog
extends ConfirmationDialog
## Edits a Compare Patients node's measures: one row per measure with its
## label, attribute, summary (median, mean, percent, correlation), the
## condition for percentages or the second attribute for correlations, and
## whether to add pairwise tests.

signal measures_changed(measures: Array)

var _rows_box: VBoxContainer
var _rows: Array = []  # [{label, attribute, summary, condition, with, pairwise, box}]
var _attributes: Array = []
var _hint: Label


func _init() -> void:
	title = "Measures"
	ok_button_text = "Save"
	min_size = Vector2i(900, 360)
	var root := VBoxContainer.new()
	root.add_theme_constant_override("separation", 8)
	add_child(root)
	var note := Label.new()
	note.text = "Each row is one line of the results table. Values come from columns imported on the Samples tab (Import from a table)."
	root.add_child(note)
	var head := HBoxContainer.new()
	for pair in [["Label", 170], ["Value (column)", 150], ["Summary and test", 210], ["Condition / with", 150], ["Pairwise", 80]]:
		var l := Label.new()
		l.text = pair[0]
		l.custom_minimum_size.x = pair[1]
		l.add_theme_color_override("font_color", Color(0.7, 0.75, 0.8))
		head.add_child(l)
	root.add_child(head)
	var scroll := ScrollContainer.new()
	scroll.custom_minimum_size.y = 220
	scroll.size_flags_vertical = Control.SIZE_EXPAND_FILL
	scroll.horizontal_scroll_mode = ScrollContainer.SCROLL_MODE_DISABLED
	root.add_child(scroll)
	_rows_box = VBoxContainer.new()
	_rows_box.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	scroll.add_child(_rows_box)
	var add := Button.new()
	add.text = "Add measure"
	add.size_flags_horizontal = Control.SIZE_SHRINK_BEGIN
	add.pressed.connect(func() -> void: _add_row({"summary": "median", "pairwise": true}))
	root.add_child(add)
	_hint = Label.new()
	_hint.add_theme_color_override("font_color", Color(1, 0.8, 0.4))
	root.add_child(_hint)
	confirmed.connect(func() -> void: measures_changed.emit(_collect()))


func open(measures: Array, attributes: Array) -> void:
	_attributes = attributes.duplicate()
	for m in measures:
		for key in ["attribute", "with"]:
			var a: String = m.get(key, "")
			if a and not a in _attributes:
				_attributes.append(a)
	for r in _rows:
		r["box"].queue_free()
	_rows.clear()
	for m in measures:
		_add_row(m)
	_hint.text = "" if attributes else "No values imported yet: on the Samples tab use Import from a table and choose value columns."
	popup_centered()


func _add_row(m: Dictionary) -> void:
	var box := HBoxContainer.new()
	var label := LineEdit.new()
	label.custom_minimum_size.x = 170
	label.placeholder_text = "e.g. SLEDAI-2K"
	label.text = m.get("label", "")
	box.add_child(label)
	var attribute := _attribute_option(m.get("attribute", ""))
	attribute.custom_minimum_size.x = 150
	box.add_child(attribute)
	var summary := OptionButton.new()
	summary.custom_minimum_size.x = 210
	for s in PatientCompare.SUMMARIES:
		summary.add_item("%s · %s" % [PatientCompare.SUMMARY_TITLES[s], PatientCompare.TESTS[s]])
	summary.select(maxi(0, PatientCompare.SUMMARIES.find(m.get("summary", "median"))))
	box.add_child(summary)
	var condition := LineEdit.new()
	condition.custom_minimum_size.x = 150
	condition.placeholder_text = "e.g. > 6, = F, = yes"
	condition.text = PatientCompare.condition_text(m.get("condition", {}))
	box.add_child(condition)
	var with := _attribute_option(m.get("with", ""))
	with.custom_minimum_size.x = 150
	box.add_child(with)
	var pairwise := CheckBox.new()
	pairwise.custom_minimum_size.x = 80
	pairwise.button_pressed = m.get("pairwise", false)
	pairwise.tooltip_text = "Also test every pair of groups (Wilcoxon rank-sum, or Fisher's exact for percentages)."
	box.add_child(pairwise)
	var remove := Button.new()
	remove.text = "Remove"
	remove.flat = true
	box.add_child(remove)
	var row := {"label": label, "attribute": attribute, "summary": summary, "condition": condition, "with": with, "pairwise": pairwise, "box": box}
	remove.pressed.connect(func() -> void:
		_rows.erase(row)
		box.queue_free())
	var show_extra := func(_i = 0) -> void:
		var s: String = PatientCompare.SUMMARIES[summary.selected]
		condition.visible = s == "percent"
		with.visible = s == "correlation"
		pairwise.disabled = s == "correlation"
	summary.item_selected.connect(show_extra)
	show_extra.call()
	_rows_box.add_child(box)
	_rows.append(row)


func _attribute_option(selected: String) -> OptionButton:
	var option := OptionButton.new()
	for a in _attributes:
		option.add_item(a)
	if _attributes.is_empty():
		option.add_item("(no values imported)")
		option.disabled = true
	option.select(maxi(0, _attributes.find(selected)))
	return option


func _collect() -> Array:
	var out: Array = []
	for r in _rows:
		if _attributes.is_empty():
			break
		var attribute: String = _attributes[r["attribute"].selected]
		var summary: String = PatientCompare.SUMMARIES[r["summary"].selected]
		var m := {"label": r["label"].text.strip_edges() if r["label"].text.strip_edges() else attribute,
			"attribute": attribute, "summary": summary, "pairwise": r["pairwise"].button_pressed and summary != "correlation"}
		if summary == "percent":
			m["condition"] = PatientCompare.parse_condition(r["condition"].text) if r["condition"].text.strip_edges() else {"op": "=", "value": "yes"}
		if summary == "correlation":
			m["with"] = _attributes[r["with"].selected]
		out.append(m)
	return out
