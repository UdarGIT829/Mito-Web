class_name AnnotationFilterDialog
extends ConfirmationDialog
## Edits an AnnotationFilters config.

signal saved(config: Dictionary)

var _action: OptionButton
var _categories: CheckListButton
var _stars: OptionButton
var _af_min: LineEdit
var _af_max: LineEdit
var _keep_absent: CheckBox
var _marker: OptionButton
var _consequences: CheckListButton
var _trna: CheckBox
var _error: Label


func _init() -> void:
	title = "Annotation rules"
	ok_button_text = "Save"
	dialog_hide_on_ok = false
	min_size = Vector2i(560, 0)
	var form := VBoxContainer.new()
	form.add_theme_constant_override("separation", 6)
	add_child(form)
	var grid := GridContainer.new()
	grid.columns = 2
	form.add_child(grid)
	_action = _option(["Keep mutations matching every rule below", "Remove mutations matching every rule below"])
	_row(grid, "Action", _action)
	_categories = CheckListButton.new("Any (no ClinVar rule)")
	_categories.set_items(AnnotationFilters.CATEGORIES)
	_row(grid, "ClinVar category is one of", _categories)
	_stars = _option(["no minimum", "1 star (criteria provided)", "2 stars (multiple submitters, no conflicts)", "3 stars (expert panel)", "4 stars (practice guideline)"])
	_row(grid, "ClinVar review at least", _stars)
	_af_min = _line("e.g. 0")
	_af_max = _line("e.g. 0.001")
	var af_box := HBoxContainer.new()
	af_box.add_child(_af_min)
	af_box.add_child(FilterDialog._label("to"))
	af_box.add_child(_af_max)
	_row(grid, "gnomAD AF hom from / to", af_box)
	_keep_absent = CheckBox.new()
	_keep_absent.text = "Mutations not in gnomAD also match the frequency rule"
	_row(grid, "", _keep_absent)
	_marker = _option(["any", "yes (PhyloTree haplogroup marker)", "no (not a haplogroup marker)"])
	_row(grid, "Haplogroup marker", _marker)
	_consequences = CheckListButton.new("Any consequence")
	_consequences.set_items(AnnotationFilters.COMMON_CONSEQUENCES)
	_row(grid, "Consequence is any of", _consequences)
	_trna = CheckBox.new()
	_trna.text = "MitoTIP or PON-mt-tRNA predicts pathogenic"
	_row(grid, "tRNA prediction", _trna)
	var note := FilterDialog._label("Categories group ClinVar's wording: pathogenic = (likely) pathogenic, benign = (likely) benign,\nconflicting = conflicting or mixed, other = risk factor, affects, association…, absent = no ClinVar record.\nA mutation absent from a source is never treated as benign.")
	note.add_theme_color_override("font_color", Color(0.7, 0.7, 0.7))
	form.add_child(note)
	_error = FilterDialog._label("")
	_error.add_theme_color_override("font_color", Color(1, 0.5, 0.45))
	form.add_child(_error)
	confirmed.connect(_on_save)


func edit(config: Dictionary) -> void:
	_action.select(1 if config.get("action", "keep") == "remove" else 0)
	_categories.set_selected(config.get("clinvar_categories", []))
	_stars.select(int(config.get("clinvar_min_stars", 0)))
	_af_min.text = "" if config.get("gnomad_af_hom_min") == null else str(config["gnomad_af_hom_min"])
	_af_max.text = "" if config.get("gnomad_af_hom_max") == null else str(config["gnomad_af_hom_max"])
	_keep_absent.button_pressed = config.get("gnomad_keep_absent", false)
	_marker.select(["any", "yes", "no"].find(config.get("haplogroup_marker", "any")))
	_consequences.set_selected(config.get("consequences", []))
	_trna.button_pressed = config.get("trna_predicted_pathogenic", false)
	_error.text = ""
	popup_centered()


func _on_save() -> void:
	var config := {"action": "remove" if _action.selected == 1 else "keep"}
	if _categories.get_selected():
		config["clinvar_categories"] = _categories.get_selected()
	if _stars.selected > 0:
		config["clinvar_min_stars"] = _stars.selected
	for pair in [[_af_min, "gnomad_af_hom_min"], [_af_max, "gnomad_af_hom_max"]]:
		var text: String = pair[0].text.strip_edges()
		if text:
			if not text.is_valid_float():
				_error.text = "gnomAD frequencies must be numbers."
				return
			config[pair[1]] = text.to_float()
	if _keep_absent.button_pressed:
		config["gnomad_keep_absent"] = true
	if _marker.selected > 0:
		config["haplogroup_marker"] = ["any", "yes", "no"][_marker.selected]
	if _consequences.get_selected():
		config["consequences"] = _consequences.get_selected()
	if _trna.button_pressed:
		config["trna_predicted_pathogenic"] = true
	var problem := AnnotationFilters.validate(config)
	if problem:
		_error.text = problem
		return
	saved.emit(config)
	hide()


static func _option(items: Array) -> OptionButton:
	var o := OptionButton.new()
	for item in items:
		o.add_item(item)
	return o


static func _line(placeholder: String) -> LineEdit:
	var l := LineEdit.new()
	l.placeholder_text = placeholder
	l.custom_minimum_size.x = 110
	return l


static func _row(grid: GridContainer, text: String, control: Control) -> void:
	grid.add_child(FilterDialog._label(text))
	grid.add_child(control)
