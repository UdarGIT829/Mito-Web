class_name SamplesPanel
extends VBoxContainer
## Imported samples grouped by subject, with call counts and how many calls
## could be placed on rCRS.

var _tree: Tree
var _summary: Label


func _ready() -> void:
	_summary = Label.new()
	add_child(_summary)
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
