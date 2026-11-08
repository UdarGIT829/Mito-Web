class_name FilterDialog
extends ConfirmationDialog
## Edits a MutationFilters config. Emits `saved` with the new config when the
## settings are valid; invalid settings keep the dialog open with a message.

signal saved(config: Dictionary)

const META_FIELDS := ["repeat_count", "repeat_base", "repeat_2_count", "repeat_2_bases", "repeat_3_count",
	"repeat_3_bases", "reference_repeat", "reference_contains_alt", "polymorphism"]
const META_HINTS := {
	"repeat_count": "gte|3", "repeat_2_count": "gte|2", "repeat_3_count": "gte|2",
	"repeat_base": "C", "repeat_2_bases": "CA", "repeat_3_bases": "ACT",
	"reference_repeat": "before / after / one / both / none / either",
	"reference_contains_alt": "contains / not_contains", "polymorphism": "1 or 0",
}

var _pos_min: SpinBox
var _pos_max: SpinBox
var _genes: CheckListButton
var _refs: LineEdit
var _alts: LineEdit
var _presence: SpinBox
var _ranges: Array = []  # [enabled CheckBox, min SpinBox, max SpinBox, include max CheckBox]
var _af_percent: SpinBox
var _meta_rows: Array = []  # [OptionButton, LineEdit]
var _error: Label


func _init() -> void:
	title = "Edit filters"
	ok_button_text = "Save"
	dialog_hide_on_ok = false
	min_size = Vector2i(560, 0)
	var form := VBoxContainer.new()
	form.add_theme_constant_override("separation", 6)
	add_child(form)

	var grid := GridContainer.new()
	grid.columns = 2
	form.add_child(grid)
	_pos_min = _spin(0, MitoGenes.LENGTH, 1)
	_pos_max = _spin(0, MitoGenes.LENGTH, 1)
	_row(grid, "Position from / to (0 = open)", _pair(_pos_min, _pos_max))
	_genes = CheckListButton.new("Any gene")
	_genes.set_items(Array(MitoGenes.symbols()))
	_row(grid, "Genes (any of)", _genes)
	_refs = _line("e.g. C, T")
	_row(grid, "rCRS reference base", _refs)
	_alts = _line("e.g. G, A, AC")
	_row(grid, "ALT allele", _alts)
	_presence = _spin(0, 100, 0.01)
	_row(grid, "Carried by at least (% of samples)", _presence)

	form.add_child(_heading("Allele fraction (a carrier passes if its AF is in any enabled range)"))
	for i in 2:
		var enabled := CheckBox.new()
		enabled.text = "Range %d" % (i + 1)
		var lo := _spin(0, 1, 0.001)
		var hi := _spin(0, 1, 0.001)
		hi.value = 1
		var include_max := CheckBox.new()
		include_max.text = "include max"
		include_max.button_pressed = true
		var row := HBoxContainer.new()
		for c in [enabled, _label("from"), lo, _label("to"), hi, include_max]:
			row.add_child(c)
		form.add_child(row)
		_ranges.append([enabled, lo, hi, include_max])
	_af_percent = _spin(0, 100, 0.01)
	var af_row := HBoxContainer.new()
	af_row.add_child(_label("Carriers passing the AF ranges, at least (%)"))
	af_row.add_child(_af_percent)
	form.add_child(af_row)

	form.add_child(_heading("Slippage metadata (all must match on at least one carrier call)"))
	for i in 3:
		var field := OptionButton.new()
		field.add_item("(unused)")
		for f in META_FIELDS:
			field.add_item(f)
		var value := _line("")
		value.size_flags_horizontal = Control.SIZE_EXPAND_FILL
		field.item_selected.connect(func(index: int) -> void:
			value.placeholder_text = META_HINTS.get(field.get_item_text(index), ""))
		var row := HBoxContainer.new()
		row.add_child(field)
		row.add_child(value)
		form.add_child(row)
		_meta_rows.append([field, value])

	_error = _label("")
	_error.add_theme_color_override("font_color", Color(1, 0.5, 0.45))
	_error.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	form.add_child(_error)
	confirmed.connect(_on_save)


func edit(config: Dictionary) -> void:
	var position: Dictionary = config.get("position", {})
	_pos_min.value = position.get("minimum", position.get("exact", 0)) if position.get("minimum", position.get("exact")) != null else 0
	_pos_max.value = position.get("maximum", position.get("exact", 0)) if position.get("maximum", position.get("exact")) != null else 0
	_genes.set_selected(config.get("genes", []))
	_refs.text = ", ".join(config.get("reference_bases", []))
	_alts.text = ", ".join(config.get("alt_alleles", []))
	_presence.value = config.get("minimum_presence_percent", 0)
	var ranges: Array = config.get("af_ranges", [])
	for i in _ranges.size():
		var widgets: Array = _ranges[i]
		var r: Dictionary = ranges[i] if i < ranges.size() else {}
		widgets[0].button_pressed = not r.is_empty()
		widgets[1].value = r.get("minimum", 0) if r.get("minimum") != null else 0
		widgets[2].value = r.get("maximum", 1) if r.get("maximum") != null else 1
		widgets[3].button_pressed = r.get("include_maximum", true)
	_af_percent.value = config.get("minimum_af_percent", 0)
	var meta: Array = config.get("slippage_metadata_filters", [])
	for i in _meta_rows.size():
		var item: Dictionary = meta[i] if i < meta.size() else {}
		var index: int = META_FIELDS.find(item.get("field", "")) + 1
		_meta_rows[i][0].select(index)
		_meta_rows[i][1].text = str(item.get("value", ""))
	_error.text = ""
	popup_centered()


func _on_save() -> void:
	var config := {}
	var position := {}
	if _pos_min.value > 0:
		position["minimum"] = int(_pos_min.value)
	if _pos_max.value > 0:
		position["maximum"] = int(_pos_max.value)
	if position:
		config["position"] = position
	if _genes.get_selected():
		config["genes"] = _genes.get_selected()
	var refs := _split(_refs.text)
	if refs:
		config["reference_bases"] = refs
	var alts := _split(_alts.text)
	if alts:
		config["alt_alleles"] = alts
	if _presence.value > 0:
		config["minimum_presence_percent"] = _presence.value
	var ranges: Array = []
	for widgets in _ranges:
		if widgets[0].button_pressed:
			ranges.append({"minimum": widgets[1].value, "maximum": widgets[2].value,
				"include_minimum": true, "include_maximum": widgets[3].button_pressed})
	if ranges:
		config["af_ranges"] = ranges
		config["minimum_af_percent"] = _af_percent.value
	var meta: Array = []
	for row in _meta_rows:
		if row[0].selected > 0:
			meta.append({"field": row[0].get_item_text(row[0].selected), "value": row[1].text.strip_edges()})
	if meta:
		config["slippage_metadata_filters"] = meta
	var problem := MutationFilters.validate(config)
	if problem:
		_error.text = problem
		return
	saved.emit(config)
	hide()


## One-line description of a filter config for the node face.
static func describe(config: Dictionary) -> String:
	var parts: Array = []
	var position: Dictionary = config.get("position", {})
	if position:
		parts.append("pos %s–%s" % [position.get("minimum", "1"), position.get("maximum", "end")])
	if config.get("genes"):
		parts.append("genes " + ", ".join(config["genes"]))
	if config.get("reference_bases"):
		parts.append("ref " + ", ".join(config["reference_bases"]))
	if config.get("alt_alleles"):
		parts.append("alt " + ", ".join(config["alt_alleles"]))
	if config.get("minimum_presence_percent", 0) > 0:
		parts.append("carried by >= %s%%" % config["minimum_presence_percent"])
	for r in config.get("af_ranges", []):
		parts.append("AF %s–%s" % [r.get("minimum", 0), r.get("maximum", 1)])
	if config.get("af_ranges"):
		parts.append(">= %s%% of carriers in AF range" % config.get("minimum_af_percent", 0))
	for m in config.get("slippage_metadata_filters", []):
		parts.append("%s %s" % [m["field"], m["value"]])
	return "\n".join(parts) if parts else "No filters set (everything passes)"


static func _split(text: String) -> Array:
	return Array(text.to_upper().replace(";", ",").split(",", false)).map(func(s): return s.strip_edges()).filter(func(s): return s != "")


func _spin(lo: float, hi: float, step: float) -> SpinBox:
	var s := SpinBox.new()
	s.min_value = lo
	s.max_value = hi
	s.step = step
	s.custom_minimum_size.x = 100
	return s


func _line(placeholder: String) -> LineEdit:
	var l := LineEdit.new()
	l.placeholder_text = placeholder
	l.custom_minimum_size.x = 160
	return l


func _pair(a: Control, b: Control) -> HBoxContainer:
	var box := HBoxContainer.new()
	box.add_child(a)
	box.add_child(b)
	return box


func _row(grid: GridContainer, text: String, control: Control) -> void:
	grid.add_child(_label(text))
	grid.add_child(control)


func _heading(text: String) -> Label:
	var l := _label(text)
	l.add_theme_color_override("font_color", Color(0.75, 0.85, 1.0))
	return l


static func _label(text: String) -> Label:
	var l := Label.new()
	l.text = text
	return l
