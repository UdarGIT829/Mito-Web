class_name SamplesPanel
extends VBoxContainer
## Imported samples grouped by subject, with call counts and how many calls
## could be placed on rCRS.

var _tree: Tree
var _summary: Label


var _tag_dialog: ConfirmationDialog
var _tag_table := {}
var _tag_match: OptionButton
var _tag_columns: CheckListButton
var _tag_prefix: CheckBox
var _tag_preview: Label


func _ready() -> void:
	var bar := HBoxContainer.new()
	add_child(bar)
	_summary = Label.new()
	_summary.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	bar.add_child(_summary)
	var import_tags := Button.new()
	import_tags.text = "Import tags from a table…"
	import_tags.tooltip_text = "CSV or TSV with a header row: one column naming the sample (or subject), other columns become tags."
	import_tags.pressed.connect(func() -> void: WebFiles.pick_files("tags", ".csv,.tsv,.txt", false))
	bar.add_child(import_tags)
	WebFiles.files_picked.connect(_on_files_picked)
	_build_tag_dialog()
	_tree = Tree.new()
	_tree.size_flags_vertical = Control.SIZE_EXPAND_FILL
	_tree.columns = 6
	_tree.column_titles_visible = true
	_tree.hide_root = true
	var titles := ["Subject / sample", "Tags", "Build", "Calls", "On rCRS", "Source file"]
	var widths := [0, 180, 100, 70, 80, 0]
	for i in 6:
		_tree.set_column_title(i, titles[i])
		_tree.set_column_expand(i, i == 0 or i == 5)
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
	item.set_text(1, ", ".join(sample["tags"]))
	item.set_text(2, FilesPanel._short_build(sample["build"]))
	item.set_text(3, str(summary["calls"]))
	item.set_text(4, str(summary["placed"]))
	var lines: Array = []
	for status_name in summary["status"]:
		lines.append("%s: %d" % [status_name.capitalize(), summary["status"][status_name]])
	item.set_tooltip_text(4, "\n".join(lines))
	if summary["placed"] < summary["calls"]:
		item.set_custom_color(4, Color(1, 0.8, 0.4))
	item.set_text(5, Workspace.model.source_by_id(sample["source_id"]).get("name", ""))


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
	_tag_prefix = CheckBox.new()
	_tag_prefix.text = "Prefix tags with the column name (e.g. haplogroup=H1)"
	form.add_child(_tag_prefix)
	_tag_preview = _plain("")
	_tag_preview.add_theme_color_override("font_color", Color(0.7, 0.7, 0.7))
	form.add_child(_tag_preview)
	_tag_dialog.confirmed.connect(func() -> void:
		var header: Array = _tag_table["header"]
		var columns: Array = _tag_columns.get_selected().map(func(name): return header.find(name))
		if columns.is_empty():
			return
		var result := Workspace.import_tags(_tag_table, _tag_match.selected, columns, _tag_prefix.button_pressed)
		var message := "Tagged %d sample(s) with %d new tag(s)." % [result["matched"], result["tags_added"]]
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
	_tag_columns.set_items(header.slice(1))
	_tag_columns.set_selected(header.slice(1))
	_tag_preview.text = "%d rows. First row: %s" % [parsed["rows"].size(), ", ".join(parsed["rows"][0])]
	_tag_dialog.popup_centered()


static func _plain(text: String) -> Label:
	var l := Label.new()
	l.text = text
	return l
