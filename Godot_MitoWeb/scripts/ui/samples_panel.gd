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
	var import_tags := Button.new()
	import_tags.text = "Import tags from a table…"
	import_tags.tooltip_text = "CSV or TSV with a header row: one column naming the sample (or subject), other columns become tags."
	import_tags.pressed.connect(func() -> void: WebFiles.pick_files("tags", ".csv,.tsv,.txt", false))
	bar.add_child(import_tags)
	WebFiles.files_picked.connect(_on_files_picked)
	_build_tag_dialog()
	_build_haplogroup_dialog()
	Workspace.haplogroup_progress.connect(func(done: int, total: int) -> void:
		_summary.text = "Assigning haplogroups… %d of %d" % [done, total])
	_tree = Tree.new()
	_tree.size_flags_vertical = Control.SIZE_EXPAND_FILL
	_tree.columns = 7
	_tree.column_titles_visible = true
	_tree.hide_root = true
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
			parent.set_selectable(0, false)
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
	item.set_tooltip_text(0, "Sample %s\nSubject %s" % [sample["name"], sample["subject"]])
	item.set_text(1, ", ".join(Array(sample["tags"]).filter(func(t): return not HaplogroupClassifier.TAG_PREFIXES.any(func(p): return t.begins_with(p)))))
	var hg: Dictionary = sample.get("haplogroup", {})
	if hg:
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
	_tag_dialog.title = "Import tags from a table"
	_tag_dialog.ok_button_text = "Add tags"
	var form := VBoxContainer.new()
	_tag_dialog.add_child(form)
	var grid := GridContainer.new()
	grid.columns = 2
	form.add_child(grid)
	grid.add_child(_plain("Column naming the sample"))
	_tag_match = OptionButton.new()
	grid.add_child(_tag_match)
	grid.add_child(_plain("Columns that become tags"))
	_tag_columns = CheckListButton.new("Choose columns")
	grid.add_child(_tag_columns)
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
		if columns.is_empty() and _tag_subject.selected == 0:
			return
		var result := Workspace.import_tags(_tag_table, _tag_match.selected, columns, _tag_prefix.button_pressed, _tag_subject.selected - 1)
		var message := "Tagged %d sample(s) with %d new tag(s)%s." % [result["matched"], result["tags_added"],
			"; subjects set from %s" % _tag_subject.get_item_text(_tag_subject.selected) if _tag_subject.selected > 0 else ""]
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
	_tag_columns.set_items(header.slice(1))
	_tag_columns.set_selected(header.slice(1))
	_tag_preview.text = "%d rows. First row: %s" % [parsed["rows"].size(), ", ".join(parsed["rows"][0])]
	_tag_dialog.popup_centered()


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
	var note := _plain("Each sample is matched against PhyloTree Build 17 with HaploGrep's weighted\nKulczynski score. Results become tags (haplogroup=…, macro haplogroup=…)\nthat tag filters can use. Calls without an AF (genotype-only data) always count.")
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
