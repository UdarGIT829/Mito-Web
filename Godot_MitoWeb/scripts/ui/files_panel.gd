class_name FilesPanel
extends HSplitContainer
## Uploaded files, their import plans, and the import controls.
## Left: file list (multi-select) with bulk actions for the selection.
## Right: details and the editable sample plan of the focused file.

signal status_message(text: String)

const SOURCE_ACCEPT := ".vcf,.gz,.bgz,.zip,.csv,.tsv"
const STATUS_LABELS := {"pending": "Ready", "imported": "Imported", "error": "Unreadable"}
const NO_BUILD := "Choose reference build…"

var _files: Tree
var _bulk_build: OptionButton
var _bulk_tags: LineEdit
var _bulk_regex: LineEdit
var _import_button: Button
var _undo_button: Button
var _remove_button: Button
var _details: VBoxContainer
var _focused_id := -1


func _ready() -> void:
	split_offset = 520
	_build_left()
	_build_right()
	Workspace.changed.connect(refresh)
	refresh()


func _build_left() -> void:
	var left := VBoxContainer.new()
	left.custom_minimum_size.x = 380
	add_child(left)

	var toolbar := HBoxContainer.new()
	left.add_child(toolbar)
	_add_button(toolbar, "Upload files…", func() -> void:
		WebFiles.pick_files("sources", SOURCE_ACCEPT, true))
	_import_button = _add_button(toolbar, "Import selected", func() -> void:
		Workspace.import_sources(_selected_ids()))
	_add_button(toolbar, "Import all ready", func() -> void:
		Workspace.import_sources(Workspace.model.sources().map(func(s): return s["id"])))

	_files = Tree.new()
	_files.size_flags_vertical = Control.SIZE_EXPAND_FILL
	_files.columns = 4
	_files.column_titles_visible = true
	_files.hide_root = true
	_files.select_mode = Tree.SELECT_MULTI
	var titles := ["File", "Status", "Build", "Samples"]
	var widths := [0, 120, 110, 70]
	for i in 4:
		_files.set_column_title(i, titles[i])
		_files.set_column_expand(i, i == 0)
		_files.set_column_custom_minimum_width(i, widths[i])
	_files.multi_selected.connect(func(_item, _col, _selected) -> void: _on_selection_changed.call_deferred())
	left.add_child(_files)

	var bulk := GridContainer.new()
	bulk.columns = 3
	left.add_child(bulk)
	bulk.add_child(_label("Build for selected"))
	_bulk_build = _build_option(true)
	_bulk_build.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	bulk.add_child(_bulk_build)
	_add_button(bulk, "Set", func() -> void:
		var build: String = _bulk_build.get_item_metadata(_bulk_build.selected)
		_for_selected_pending(func(source: Dictionary) -> void: source["plan"]["build"] = build)
		_after_bulk_edit("Set the reference build on the selected files."))

	bulk.add_child(_label("Add tags to selected"))
	_bulk_tags = LineEdit.new()
	_bulk_tags.placeholder_text = "e.g. Treated, EV"
	_bulk_tags.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	bulk.add_child(_bulk_tags)
	_add_button(bulk, "Add", func() -> void:
		var tags := VcfImporter.split_tags(_bulk_tags.text)
		if tags.is_empty():
			return
		_for_selected_pending(func(source: Dictionary) -> void:
			for sample in source["plan"]["samples"]:
				for tag in tags:
					if not tag in sample["tags"]:
						sample["tags"].append(tag))
		_after_bulk_edit("Added tags %s to the selected files." % ", ".join(tags)))

	bulk.add_child(_label("Subject from name"))
	_bulk_regex = LineEdit.new()
	_bulk_regex.placeholder_text = "regex, first group = subject, e.g. ^(C\\d+)_"
	_bulk_regex.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	bulk.add_child(_bulk_regex)
	_add_button(bulk, "Apply", func() -> void:
		var changed := 0
		var invalid := false
		for source in _selected_sources():
			if source["status"] != WorkspaceModel.STATUS_PENDING:
				continue
			var n := VcfImporter.apply_subject_regex(source["plan"], _bulk_regex.text)
			if n < 0:
				invalid = true
				break
			changed += n
		_after_bulk_edit("That regular expression is not valid." if invalid else "Changed %d subject(s)." % changed))

	var actions := HBoxContainer.new()
	left.add_child(actions)
	_undo_button = _add_button(actions, "Undo import", func() -> void:
		Workspace.undo_imports(_selected_ids())
		status_message.emit("Returned the selected files to ‘Ready to import’."))
	_remove_button = _add_button(actions, "Remove from project", func() -> void:
		var ids := _selected_ids()
		if ids:
			Workspace.remove_sources(ids)
			status_message.emit("Removed %d file(s) and their samples." % ids.size()))


func _build_right() -> void:
	var scroll := ScrollContainer.new()
	scroll.custom_minimum_size.x = 420
	scroll.horizontal_scroll_mode = ScrollContainer.SCROLL_MODE_DISABLED
	add_child(scroll)
	_details = VBoxContainer.new()
	_details.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_details.custom_minimum_size.x = 0
	_details.size_flags_vertical = Control.SIZE_EXPAND_FILL
	_details.add_theme_constant_override("separation", 8)
	scroll.add_child(_details)


func refresh() -> void:
	var selected := _selected_ids()
	_files.clear()
	var root := _files.create_item()
	for source in Workspace.model.sources():
		var item := _files.create_item(root)
		item.set_text(0, source["name"])
		item.set_tooltip_text(0, "%s\nSHA-256 %s" % [source["name"], source["sha256"]])
		item.set_metadata(0, source["id"])
		item.set_text(1, STATUS_LABELS.get(source["status"], source["status"]))
		if source["status"] == WorkspaceModel.STATUS_ERROR:
			item.set_custom_color(1, Color(1, 0.5, 0.45))
			item.set_tooltip_text(1, source["error"])
		elif source["status"] == WorkspaceModel.STATUS_PENDING and source["plan"]["build"].is_empty():
			item.set_text(1, "Needs build")
			item.set_custom_color(1, Color(1, 0.8, 0.4))
		var build: String = source["plan"].get("build", "")
		item.set_text(2, _short_build(build))
		item.set_tooltip_text(2, ReferenceSequence.label_for(build) if build else "Not chosen")
		item.set_text(3, str(source["plan"]["samples"].size()))
		if source["id"] in selected:
			item.select(0)
	if Workspace.model.sources().is_empty():
		var hint := _files.create_item(root)
		hint.set_text(0, "Upload VCF files (.vcf, .vcf.gz) or a .zip of them, optionally with a manifest.csv. They stay in this browser.")
		for i in 4:
			hint.set_selectable(i, false)
	var busy := Workspace.importing
	for button in [_import_button, _undo_button, _remove_button]:
		button.disabled = busy
	_show_details()


func _on_selection_changed() -> void:
	var ids := _selected_ids()
	_focused_id = ids[0] if ids.size() == 1 else -1
	_show_details()


func _show_details() -> void:
	for child in _details.get_children():
		child.queue_free()
	var source := Workspace.model.source_by_id(_focused_id)
	if source.is_empty():
		var count := _selected_ids().size()
		_details.add_child(_label(
			"%d files selected. Use the controls on the left to set their build, tags and subjects." % count
			if count > 1 else "Select a file to see its details and samples."))
		return
	var title := _label(source["name"])
	title.add_theme_font_size_override("font_size", 18)
	title.autowrap_mode = TextServer.AUTOWRAP_ARBITRARY
	_details.add_child(title)
	if source["status"] == WorkspaceModel.STATUS_ERROR:
		_details.add_child(_warning("This file could not be read: %s" % source["error"]))
		return
	_details.add_child(_facts(source))
	for warning in _warnings(source):
		_details.add_child(_warning(warning))

	var editable: bool = source["status"] == WorkspaceModel.STATUS_PENDING and not Workspace.importing
	if source["status"] == WorkspaceModel.STATUS_IMPORTED:
		_details.add_child(_label(_result_text(source["result"])))
		_details.add_child(_label("To change this file's settings, use “Undo import”, edit, and import again."))

	var build_row := HBoxContainer.new()
	build_row.add_child(_label("Reference build"))
	var build_option := _build_option(source["plan"]["build"].is_empty())
	build_option.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	build_option.disabled = not editable
	for i in build_option.item_count:
		if build_option.get_item_metadata(i) == source["plan"]["build"]:
			build_option.select(i)
	build_option.item_selected.connect(func(index: int) -> void:
		source["plan"]["build"] = build_option.get_item_metadata(index)
		Workspace.plan_edited())
	build_row.add_child(build_option)
	_details.add_child(build_row)

	_details.add_child(_label("Samples (double-click a name, subject or tags cell to edit)" if editable else "Samples"))
	var samples := Tree.new()
	samples.columns = 4
	samples.column_titles_visible = true
	samples.hide_root = true
	samples.custom_minimum_size.y = 260
	samples.size_flags_vertical = Control.SIZE_EXPAND_FILL
	var titles := ["Import", "Sample name", "Subject", "Tags"]
	for i in 4:
		samples.set_column_title(i, titles[i])
		samples.set_column_expand(i, i > 0)
	samples.set_column_custom_minimum_width(0, 60)
	var root := samples.create_item()
	for planned in source["plan"]["samples"]:
		var item := samples.create_item(root)
		item.set_cell_mode(0, TreeItem.CELL_MODE_CHECK)
		item.set_checked(0, planned.get("include", true))
		item.set_text(1, planned["name"])
		item.set_text(2, planned["subject"])
		item.set_text(3, ", ".join(planned["tags"]))
		item.set_tooltip_text(1, "VCF column: %s" % planned["column"] if planned["column"] else "Sites-only VCF")
		for col in 4:
			item.set_editable(col, editable)
		item.set_metadata(0, planned)
	samples.item_edited.connect(func() -> void:
		var item := samples.get_edited()
		var planned: Dictionary = item.get_metadata(0)
		planned["include"] = item.is_checked(0)
		planned["name"] = item.get_text(1).strip_edges()
		planned["subject"] = item.get_text(2).strip_edges()
		planned["tags"] = VcfImporter.split_tags(item.get_text(3))
		Workspace.plan_edited())
	_details.add_child(samples)


func _facts(source: Dictionary) -> Label:
	var info: Dictionary = source["info"]
	var lines: Array = []
	lines.append("%s records, %d sample column(s), %s" % [
		info["record_count"], info["sample_columns"].size(), String.humanize_size(source["size"])])
	if info["contig_id"]:
		lines.append("Contig: %s, length %d" % [info["contig_id"], info["contig_length"]])
	else:
		lines.append("Contig length: not stated in the file")
	if info["reference_line"]:
		lines.append("Reference: %s" % info["reference_line"].get_file())
	if info["source_line"]:
		lines.append("Caller: %s" % info["source_line"])
	if info["detected_build"]:
		lines.append("Detected build: %s (from the contig length)" % ReferenceSequence.label_for(info["detected_build"]))
	var agreement: Array = []
	for build in ReferenceSequence.build_keys():
		var counts: Array = info["ref_agreement"][build]
		agreement.append("%s %d/%d" % [_short_build(build), counts[0], counts[1]])
	lines.append("REF alleles matching each build: " + ", ".join(agreement))
	var label := _label("\n".join(lines))
	label.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	return label


func _warnings(source: Dictionary) -> Array:
	var info: Dictionary = source["info"]
	var out: Array = []
	if source["status"] == WorkspaceModel.STATUS_PENDING and source["plan"]["build"].is_empty():
		out.append("The file does not state its reference build. Choose it below. REF agreement is shown as evidence only; it cannot prove which reference was used.")
	if info.get("blank_contig_records", 0):
		out.append("%d record(s) have an empty CHROM field. They will be treated as chrM." % info["blank_contig_records"])
	if info.get("other_contig_records", 0):
		out.append("%d record(s) are on non-mitochondrial contigs and will be skipped." % info["other_contig_records"])
	return out


func _result_text(result: Dictionary) -> String:
	var counts: Dictionary = result.get("status_counts", {})
	var placed: int = counts.get("NATIVE", 0) + counts.get("MAPPED", 0)
	var text := "Imported %d calls; %d placed on rCRS." % [result.get("calls", 0), placed]
	var other: Array = []
	for key in counts:
		if key != "NATIVE" and key != "MAPPED":
			other.append("%s %d" % [key.capitalize(), counts[key]])
	if other:
		text += " Not placed: " + ", ".join(other) + "."
	return text


func _selected_ids() -> Array:
	var ids: Array = []
	if _files == null:
		return ids
	var item := _files.get_next_selected(null)
	while item:
		if item.get_metadata(0) != null:
			ids.append(item.get_metadata(0))
		item = _files.get_next_selected(item)
	return ids


func _selected_sources() -> Array:
	return _selected_ids().map(func(id: int) -> Dictionary: return Workspace.model.source_by_id(id))


func _for_selected_pending(action: Callable) -> void:
	for source in _selected_sources():
		if source["status"] == WorkspaceModel.STATUS_PENDING:
			action.call(source)


func _after_bulk_edit(message: String) -> void:
	if _selected_ids().is_empty():
		status_message.emit("Select one or more files first.")
		return
	Workspace.plan_edited()
	status_message.emit(message)


func _build_option(include_placeholder: bool) -> OptionButton:
	var option := OptionButton.new()
	option.fit_to_longest_item = false
	option.clip_text = true
	option.custom_minimum_size.x = 160
	if include_placeholder:
		option.add_item(NO_BUILD)
		option.set_item_metadata(0, "")
	for build in ReferenceSequence.build_keys():
		option.add_item(ReferenceSequence.label_for(build))
		option.set_item_metadata(option.item_count - 1, build)
	return option


static func _short_build(build: String) -> String:
	return {"rcrs": "rCRS/hg38", "hg19": "hg19", "hg19_plus58": "hg19+58", "": "—"}.get(build, build)


func _add_button(parent: Control, text: String, action: Callable) -> Button:
	var button := Button.new()
	button.text = text
	button.pressed.connect(action)
	parent.add_child(button)
	return button


static func _label(text: String) -> Label:
	var label := Label.new()
	label.text = text
	return label


static func _warning(text: String) -> Label:
	var label := _label(text)
	label.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	label.add_theme_color_override("font_color", Color(1, 0.8, 0.4))
	return label
