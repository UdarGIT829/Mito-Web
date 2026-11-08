class_name SamplesPanel
extends VBoxContainer
## Imported samples grouped by subject, with call counts and how many calls
## could be placed on rCRS.

var _tree: Tree
var _summary: Label


var _tag_dialog: ConfirmationDialog
var _hg_dialog: ConfirmationDialog
var _hg_af: LineEdit
var _tag_table := {}
var _tag_match: OptionButton
var _tag_columns: CheckListButton
var _tag_prefix: CheckBox
var _tag_subject: OptionButton
var _tag_preview: Label
var _tag_values: CheckListButton
var _tag_haplogroup: OptionButton
var _set_dialog: ConfirmationDialog
var _set_name: LineEdit
var _set_preview: Label
var _set_ids: Array = []
var _lineage_dialog: ConfirmationDialog
var _lineage_preset: OptionButton
var _lineage_text: TextEdit
var _lineage_errors: Label


func _ready() -> void:
	var bar := HBoxContainer.new()
	add_child(bar)
	_summary = Label.new()
	_summary.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_summary.text_overrun_behavior = TextServer.OVERRUN_TRIM_ELLIPSIS  # long messages must not widen the page
	_summary.mouse_filter = Control.MOUSE_FILTER_PASS
	bar.add_child(_summary)
	var assign := Button.new()
	assign.text = "Assign haplogroups…"
	assign.tooltip_text = "Classifies each sample against PhyloTree 17 (HaploGrep's method) and adds haplogroup tags."
	assign.pressed.connect(func() -> void: _hg_dialog.popup_centered())
	bar.add_child(assign)
	var set_hg := Button.new()
	set_hg.text = "Set haplogroup…"
	set_hg.tooltip_text = "Sets the haplogroup of the selected samples (or patients) by hand. The classifier then leaves them alone."
	set_hg.pressed.connect(_open_set_dialog)
	bar.add_child(set_hg)
	var lineages := Button.new()
	lineages.text = "Lineages…"
	lineages.tooltip_text = "Groups samples into lineages (e.g. African, Amerindian) from their haplogroups, as lineage= tags."
	lineages.pressed.connect(_open_lineage_dialog)
	bar.add_child(lineages)
	var import_tags := Button.new()
	import_tags.text = "Import from a table…"
	import_tags.tooltip_text = "CSV or TSV with a header row: one column naming the sample (or subject); other columns become values (e.g. clinical scores), tags, or haplogroups set by hand."
	import_tags.pressed.connect(func() -> void: WebFiles.pick_files("tags", ".csv,.tsv,.txt", false))
	bar.add_child(import_tags)
	WebFiles.files_picked.connect(_on_files_picked)
	_build_tag_dialog()
	_build_haplogroup_dialog()
	_build_set_dialog()
	_build_lineage_dialog()
	Workspace.haplogroup_progress.connect(func(done: int, total: int) -> void:
		_summary.text = "Assigning haplogroups… %d of %d" % [done, total])
	_tree = Tree.new()
	_tree.size_flags_vertical = Control.SIZE_EXPAND_FILL
	_tree.columns = 7
	_tree.column_titles_visible = true
	_tree.hide_root = true
	_tree.select_mode = Tree.SELECT_MULTI
	var titles := ["Subject / sample", "Tags", "Haplogroup (score)", "Build", "Calls", "On rCRS", "Source file"]
	var widths := [0, 180, 160, 100, 70, 80, 0]
	for i in 7:
		_tree.set_column_title(i, titles[i])
		_tree.set_column_expand(i, i == 0 or i == 6)
		_tree.set_column_custom_minimum_width(i, widths[i])
	add_child(_tree)
	Workspace.changed.connect(refresh)
	refresh()


func refresh() -> void:
	if Workspace.importing:
		return
	var store := Workspace.model.store
	var summaries := store.sample_summaries()
	_tree.clear()
	var root := _tree.create_item()
	var by_subject := {}
	for sample in store.samples:
		by_subject.get_or_add(sample["subject"], []).append(sample)
	var subjects := by_subject.keys()
	subjects.sort_custom(func(a, b) -> bool: return str(a).naturalnocasecmp_to(str(b)) < 0)
	for subject in subjects:
		var members: Array = by_subject[subject]
		var parent := root
		if members.size() > 1:
			parent = _tree.create_item(root)
			parent.set_text(0, "%s  (%d samples)" % [subject, members.size()])
			parent.set_metadata(0, members.map(func(m): return m["id"]))
		for sample in members:
			_add_sample_row(parent, sample, summaries[sample["id"]], members.size() == 1)
	_summary.text = "%d samples from %d subjects; %d calls on %d distinct rCRS alleles." % [
		store.samples.size(), by_subject.size(), store.call_count(), store.allele_pos.size()]
	if store.samples.is_empty():
		_summary.text = "No samples imported yet. Upload VCF files on the Files tab and import them."


func _add_sample_row(parent: TreeItem, sample: Dictionary, summary: Dictionary, show_subject: bool) -> void:
	var item := _tree.create_item(parent)
	var label: String = sample["name"]
	if show_subject and sample["subject"] != sample["name"]:
		label = "%s  (subject %s)" % [sample["name"], sample["subject"]]
	item.set_text(0, label)
	item.set_metadata(0, [sample["id"]])
	var values: Dictionary = sample.get("attributes", {})
	item.set_tooltip_text(0, "Sample %s\nSubject %s%s" % [sample["name"], sample["subject"],
		("\n" + "\n".join(values.keys().map(func(k): return "%s: %s" % [k, PatientCompare._num(values[k]) if typeof(values[k]) == TYPE_FLOAT else str(values[k])]))) if values else ""])
	item.set_text(1, ", ".join(Array(sample["tags"]).filter(func(t): return not HaplogroupClassifier.TAG_PREFIXES.any(func(p): return t.begins_with(p)))))
	var hg: Dictionary = sample.get("haplogroup", {})
	if hg and hg.get("manual", false):
		item.set_text(2, "%s (set by hand)" % hg["name"])
		item.set_tooltip_text(2, "%s, set by hand%s. Clear it with Set haplogroup… to classify the sample again." % [hg["name"], " (written as %s)" % hg["given"] if hg.get("note") else ""])
		item.set_custom_color(2, Color(0.6, 0.85, 1.0))
	elif hg:
		item.set_text(2, "%s (%s)" % [hg["name"], String.num(hg["score"], 2)])
		item.set_tooltip_text(2, "%s, macro %s, score %.3f\nRunner-up: %s\nExpected variants found: %d\nMissing: %s\nExtra (in the tree, not expected): %s\nCalls used: AF >= %s or no AF\n%s" % [
			hg["name"], hg["macro"], hg["score"],
			", ".join(hg["top"].slice(1).map(func(x): return "%s %.3f" % [x[0], x[1]])),
			hg["found"], ", ".join(hg["missing"]) if hg["missing"] else "none",
			", ".join(hg["extra"]) if hg["extra"] else "none", hg["af_min"], hg["tree"]])
		if hg["score"] < 0.9:
			item.set_custom_color(2, Color(1, 0.8, 0.4))
	item.set_text(3, FilesPanel._short_build(sample["build"]))
	item.set_text(4, str(summary["calls"]))
	item.set_text(5, str(summary["placed"]))
	var lines: Array = []
	for status_name in summary["status"]:
		lines.append("%s: %d" % [status_name.capitalize(), summary["status"][status_name]])
	item.set_tooltip_text(5, "\n".join(lines))
	if summary["placed"] < summary["calls"]:
		item.set_custom_color(5, Color(1, 0.8, 0.4))
	item.set_text(6, Workspace.model.source_by_id(sample["source_id"]).get("name", ""))


func _build_tag_dialog() -> void:
	_tag_dialog = ConfirmationDialog.new()
	_tag_dialog.title = "Import from a table"
	_tag_dialog.ok_button_text = "Import"
	var form := VBoxContainer.new()
	_tag_dialog.add_child(form)
	var grid := GridContainer.new()
	grid.columns = 2
	form.add_child(grid)
	grid.add_child(_plain("Column naming the sample"))
	_tag_match = OptionButton.new()
	grid.add_child(_tag_match)
	grid.add_child(_plain("Columns that become values"))
	_tag_values = CheckListButton.new("Choose columns")
	_tag_values.tooltip_text = "Per-sample values such as clinical scores, age or sex, for Compare Patients nodes. Numbers stay numbers."
	grid.add_child(_tag_values)
	grid.add_child(_plain("Columns that become tags"))
	_tag_columns = CheckListButton.new("Choose columns")
	_tag_columns.tooltip_text = "Labels for grouping samples (tag filters, data providers)."
	grid.add_child(_tag_columns)
	grid.add_child(_plain("Haplogroup column (sets haplogroups by hand)"))
	_tag_haplogroup = OptionButton.new()
	_tag_haplogroup.tooltip_text = "Haplogroups decided outside the classifier, e.g. after resolving mixed-up samples. Names like C1c*4 resolve to the nearest PhyloTree 17 haplogroup (C1c)."
	grid.add_child(_tag_haplogroup)
	grid.add_child(_plain("Set subject from column"))
	_tag_subject = OptionButton.new()
	_tag_subject.tooltip_text = "e.g. a cell line or donor column, so a subject's samples can be analysed as one unit."
	grid.add_child(_tag_subject)
	_tag_prefix = CheckBox.new()
	_tag_prefix.text = "Prefix tags with the column name (e.g. haplogroup=H1)"
	form.add_child(_tag_prefix)
	_tag_preview = _plain("")
	_tag_preview.add_theme_color_override("font_color", Color(0.7, 0.7, 0.7))
	form.add_child(_tag_preview)
	_tag_dialog.confirmed.connect(func() -> void:
		var header: Array = _tag_table["header"]
		var columns: Array = _tag_columns.get_selected().map(func(name): return header.find(name))
		var value_columns: Array = _tag_values.get_selected().map(func(name): return header.find(name))
		var hg_column := _tag_haplogroup.selected - 1
		if columns.is_empty() and value_columns.is_empty() and hg_column < 0 and _tag_subject.selected == 0:
			return
		var result := Workspace.import_table(_tag_table, _tag_match.selected, columns, _tag_prefix.button_pressed, _tag_subject.selected - 1, value_columns, hg_column)
		var parts: Array = ["Matched %d sample(s)" % result["matched"]]
		if value_columns:
			parts.append("%d value(s) set" % result["attributes_set"])
		if columns:
			parts.append("%d new tag(s)" % result["tags_added"])
		if hg_column >= 0:
			parts.append("%d haplogroup(s) set by hand" % result["haplogroups_set"])
		if _tag_subject.selected > 0:
			parts.append("subjects set from %s" % _tag_subject.get_item_text(_tag_subject.selected))
		var message := ", ".join(parts) + "."
		if result["haplogroup_errors"]:
			message += " " + " ".join(result["haplogroup_errors"].slice(0, 3))
		if result["unmatched"]:
			message += " No sample matched %d row(s): %s" % [result["unmatched"].size(), ", ".join(result["unmatched"].slice(0, 5))]
		_summary.text = message)
	add_child(_tag_dialog)


func _on_files_picked(purpose: String, files: Array) -> void:
	if purpose != "tags" or files.is_empty():
		return
	var parsed := TagImport.parse(files[0]["bytes"].get_string_from_utf8())
	if not parsed["ok"]:
		_summary.text = parsed["error"]
		return
	_tag_table = parsed
	var header: Array = parsed["header"]
	_tag_match.clear()
	for name in header:
		_tag_match.add_item(name)
	_tag_subject.clear()
	_tag_subject.add_item("(keep subjects)")
	for name in header:
		_tag_subject.add_item(name)
	_tag_haplogroup.clear()
	_tag_haplogroup.add_item("(none)")
	var hg_guess := 0
	for i in header.size():
		_tag_haplogroup.add_item(header[i])
		if hg_guess == 0 and i > 0 and "haplogroup" in str(header[i]).to_lower():
			hg_guess = i + 1
	_tag_haplogroup.select(hg_guess)
	var rest: Array = header.slice(1).filter(func(h): return hg_guess == 0 or h != header[hg_guess - 1])
	var numeric := rest.filter(func(h): return _numeric_column(parsed, header.find(h)))
	_tag_columns.set_items(header.slice(1))
	_tag_values.set_items(header.slice(1))
	# Tables with numbers are usually values (clinical data); text-only tables are usually tags.
	if numeric.is_empty():
		_tag_columns.set_selected(rest)
		_tag_values.set_selected([])
	else:
		_tag_columns.set_selected([])
		_tag_values.set_selected(rest)
	_tag_preview.text = "%d rows. First row: %s" % [parsed["rows"].size(), ", ".join(parsed["rows"][0])]
	_tag_dialog.popup_centered()


## Whether every non-empty cell of a column is a number.
static func _numeric_column(table: Dictionary, col: int) -> bool:
	var seen := false
	for row in table["rows"]:
		var v: String = row[col].strip_edges()
		if v.is_empty():
			continue
		if not (v.is_valid_float() or v.is_valid_int()):
			return false
		seen = true
	return seen


func _selected_sample_ids() -> Array:
	var ids: Array = []
	var item := _tree.get_next_selected(null)
	while item:
		for id in item.get_metadata(0) if item.get_metadata(0) else []:
			if not id in ids:
				ids.append(id)
		item = _tree.get_next_selected(item)
	return ids


func _build_set_dialog() -> void:
	_set_dialog = ConfirmationDialog.new()
	_set_dialog.title = "Set haplogroup by hand"
	_set_dialog.ok_button_text = "Set"
	var box := VBoxContainer.new()
	_set_dialog.add_child(box)
	box.add_child(_plain("PhyloTree 17 haplogroup for the selected samples. Leave empty to clear\nthe override and let Assign haplogroups classify them again."))
	_set_name = LineEdit.new()
	_set_name.placeholder_text = "e.g. A2o, L3, C1c*4"
	_set_name.custom_minimum_size.x = 320
	_set_name.text_changed.connect(func(_t: String) -> void: _update_set_preview())
	box.add_child(_set_name)
	_set_dialog.register_text_enter(_set_name)
	_set_preview = _plain("")
	_set_preview.add_theme_color_override("font_color", Color(0.7, 0.75, 0.8))
	box.add_child(_set_preview)
	_set_dialog.confirmed.connect(func() -> void:
		_summary.text = Workspace.set_haplogroups(_set_ids, _set_name.text)
		_summary.tooltip_text = _summary.text)
	add_child(_set_dialog)


func _open_set_dialog() -> void:
	_set_ids = _selected_sample_ids()
	if _set_ids.is_empty():
		_summary.text = "Select samples (or a patient's row) first, then Set haplogroup…"
		return
	var current := {}
	for id in _set_ids:
		var hg: Dictionary = Workspace.model.store.sample_by_id(id).get("haplogroup", {})
		if hg.get("manual", false):
			current[hg["name"]] = true
	_set_name.text = current.keys()[0] if current.size() == 1 else ""
	_set_dialog.title = "Set haplogroup for %d sample(s)" % _set_ids.size()
	_update_set_preview()
	_set_dialog.popup_centered()
	_set_name.grab_focus.call_deferred()


func _update_set_preview() -> void:
	var text := _set_name.text.strip_edges()
	if text.is_empty():
		_set_preview.text = "Clears haplogroups set by hand."
		return
	var node := HaplogroupClassifier.resolve_name(text)
	if node < 0:
		_set_preview.text = "Not a PhyloTree 17 haplogroup."
		return
	var name := HaplogroupClassifier.node_name(node)
	var lineage := Lineages.group_for(node, Workspace.model.lineage_rules())
	_set_preview.text = "%s%s, in %s%s" % [name, " (from %s)" % text if name != text else "",
		", ".join(HaplogroupClassifier.clades(node)), "; lineage %s" % lineage if lineage else ""]


func _build_lineage_dialog() -> void:
	_lineage_dialog = ConfirmationDialog.new()
	_lineage_dialog.title = "Lineages"
	_lineage_dialog.ok_button_text = "Apply"
	var box := VBoxContainer.new()
	box.add_theme_constant_override("separation", 6)
	_lineage_dialog.add_child(box)
	box.add_child(_plain("One rule per line: clade = lineage. The most specific rule wins, so\n\"N = Asian\" with \"A2 = Amerindian\" puts A2 in Amerindian and other N in Asian.\nSamples get lineage= tags that data providers and tag filters can use."))
	var row := HBoxContainer.new()
	box.add_child(row)
	row.add_child(_plain("Start from"))
	_lineage_preset = OptionButton.new()
	_lineage_preset.add_item("(current rules)")
	for key in Lineages.PRESETS:
		_lineage_preset.add_item(Lineages.PRESETS[key]["title"])
		_lineage_preset.set_item_metadata(_lineage_preset.item_count - 1, key)
	_lineage_preset.item_selected.connect(func(i: int) -> void:
		if i == 0:
			_lineage_text.text = Lineages.to_text(Workspace.model.lineage_rules())
		else:
			_lineage_text.text = Lineages.to_text(Lineages.PRESETS[_lineage_preset.get_item_metadata(i)]["rules"])
		_check_lineage_text())
	row.add_child(_lineage_preset)
	_lineage_text = TextEdit.new()
	_lineage_text.custom_minimum_size = Vector2(460, 300)
	_lineage_text.text_changed.connect(_check_lineage_text)
	box.add_child(_lineage_text)
	_lineage_errors = _plain("")
	_lineage_errors.add_theme_color_override("font_color", Color(1, 0.8, 0.4))
	box.add_child(_lineage_errors)
	_lineage_dialog.confirmed.connect(func() -> void:
		var parsed := Lineages.parse_text(_lineage_text.text)
		if parsed["errors"]:
			_summary.text = "Lineages not applied: " + " ".join(parsed["errors"].slice(0, 3))
			return
		_summary.text = Workspace.set_lineage_rules(parsed["rules"])
		_summary.tooltip_text = _summary.text)
	add_child(_lineage_dialog)


func _open_lineage_dialog() -> void:
	var rules := Workspace.model.lineage_rules()
	if rules.is_empty():
		_lineage_preset.select(1)
		_lineage_text.text = Lineages.to_text(Lineages.PRESETS[_lineage_preset.get_item_metadata(1)]["rules"])
	else:
		_lineage_preset.select(0)
		_lineage_text.text = Lineages.to_text(rules)
	_check_lineage_text()
	_lineage_dialog.popup_centered()


func _check_lineage_text() -> void:
	var parsed := Lineages.parse_text(_lineage_text.text)
	_lineage_errors.text = " ".join(parsed["errors"].slice(0, 3)) if parsed["errors"] else "%d rules, %d lineages: %s" % [
		parsed["rules"].size(), Lineages.groups(parsed["rules"]).size(), ", ".join(Lineages.groups(parsed["rules"]))]


static func _plain(text: String) -> Label:
	var l := Label.new()
	l.text = text
	return l


func _build_haplogroup_dialog() -> void:
	_hg_dialog = ConfirmationDialog.new()
	_hg_dialog.title = "Assign haplogroups"
	_hg_dialog.ok_button_text = "Assign"
	var box := VBoxContainer.new()
	_hg_dialog.add_child(box)
	var note := _plain("Each sample is matched against PhyloTree Build 17 with HaploGrep's weighted\nKulczynski score. Results become tags (haplogroup=…, macro haplogroup=…, and haplogroup clade=…\nfor every major group it descends from) that tag filters can use. Calls without an AF (genotype-only data) always count.")
	box.add_child(note)
	var row := HBoxContainer.new()
	row.add_child(_plain("Use calls with AF at least"))
	_hg_af = LineEdit.new()
	_hg_af.text = "0.75"
	_hg_af.custom_minimum_size.x = 80
	_hg_af.tooltip_text = "0.75 = the germline cut-off. Lower values let heteroplasmic calls in; a second lineage at low AF can then show up (e.g. contamination)."
	row.add_child(_hg_af)
	box.add_child(row)
	_hg_dialog.confirmed.connect(func() -> void:
		var af := _hg_af.text.to_float() if _hg_af.text.is_valid_float() else 0.75
		_summary.text = "Assigning haplogroups…"
		_summary.text = await Workspace.assign_haplogroups(clampf(af, 0.0, 1.0))
		_summary.tooltip_text = _summary.text)
	add_child(_hg_dialog)
