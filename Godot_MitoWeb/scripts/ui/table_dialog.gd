class_name TableDialog
extends AcceptDialog
## Shows a MutationFrame as a mutation table. Click a column title to sort.
## Copy or download the full table as tab-separated text.

const PAGE := 1000
const BASE_COLUMNS := ["Position", "Ref", "ALT", "Genes", "Carriers", "Presence %", "Samples (AF)"]
const SORT_FOR_COLUMN := {0: "position", 1: "reference", 2: "alt", 3: "gene", 4: "carrier_count", 5: "presence_percent"}

var _tree: Tree
var _info: Label
var _frame: MutationFrame
var _name := ""
var _sort := [{"field": "position", "direction": "asc"}]
var _rows: Array = []


func _init() -> void:
	title = "Table"
	ok_button_text = "Close"
	min_size = Vector2i(900, 560)
	var box := VBoxContainer.new()
	add_child(box)
	var bar := HBoxContainer.new()
	box.add_child(bar)
	_info = Label.new()
	_info.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	bar.add_child(_info)
	var copy := Button.new()
	copy.text = "Copy table"
	copy.pressed.connect(func() -> void:
		DisplayServer.clipboard_set(_tsv())
		_info.text = "Copied %d rows to the clipboard." % _rows.size())
	bar.add_child(copy)
	var download := Button.new()
	download.text = "Download .tsv"
	download.pressed.connect(func() -> void:
		var file_name := (_name if _name else "table").validate_filename().replace(" ", "_") + ".tsv"
		WebFiles.save_file(file_name, _tsv().to_utf8_buffer(), "text/tab-separated-values"))
	bar.add_child(download)
	_tree = Tree.new()
	_tree.size_flags_vertical = Control.SIZE_EXPAND_FILL
	_tree.hide_root = true
	_tree.column_titles_visible = true
	_tree.column_title_clicked.connect(_on_title_clicked)
	box.add_child(_tree)


func show_frame(frame: MutationFrame, table_name: String) -> void:
	_frame = frame
	_name = table_name
	title = table_name if table_name else "Table"
	_refresh()
	popup_centered_ratio(0.85)


func _refresh() -> void:
	var built := MutationTable.build(_frame, {"sort": _sort})
	_rows = built["rows"]
	var flag_names: Array = _frame.flags.keys()
	var columns: Array = BASE_COLUMNS.duplicate()
	columns.append_array(flag_names)
	_tree.clear()
	_tree.columns = columns.size()
	for i in columns.size():
		var arrow := ""
		if SORT_FOR_COLUMN.get(i) == _sort[0]["field"]:
			arrow = " (asc)" if _sort[0]["direction"] == "asc" else " (desc)"
		_tree.set_column_title(i, columns[i] + arrow)
		_tree.set_column_expand(i, i == 6)
		_tree.set_column_custom_minimum_width(i, 70 if i != 3 else 110)
	var root := _tree.create_item()
	for row in _rows.slice(0, PAGE):
		var item := _tree.create_item(root)
		var values := _values(row, flag_names)
		for i in values.size():
			item.set_text(i, values[i])
		item.set_tooltip_text(6, values[6].replace("; ", "\n"))
	_info.text = "%d mutations across %d samples%s" % [
		_rows.size(), _frame.samples.size(), " (showing the first %d)" % PAGE if _rows.size() > PAGE else ""]


func _values(row: Dictionary, flag_names: Array) -> Array:
	var s: Dictionary = row["sample_summary"]
	var samples: Array = row["sample_results"].map(func(r):
		return "%s %s" % [r["sample_label"], "—" if r["af"] == null else String.num(r["af"], 4)])
	var values := [
		str(row["position"]), row["event_ref"], row["alt"], ", ".join(row["genes"]) if row["genes"] else row["region"].replace("_", " "),
		"%d/%d" % [s["carrier_count"], s["selected_count"]], String.num(s["presence_percent"], 1), "; ".join(samples),
	]
	for name in flag_names:
		values.append("yes" if row.get(name, false) else "no")
	return values


func _tsv() -> String:
	var flag_names: Array = _frame.flags.keys()
	var header: Array = BASE_COLUMNS.duplicate()
	header.append_array(flag_names)
	var lines := PackedStringArray(["\t".join(header)])
	for row in _rows:
		lines.append("\t".join(_values(row, flag_names)))
	return "\n".join(lines) + "\n"


func _on_title_clicked(column: int, _mouse_button: int) -> void:
	var field: String = SORT_FOR_COLUMN.get(column, "")
	if field.is_empty():
		return
	var direction := "asc"
	if _sort[0]["field"] == field and _sort[0]["direction"] == "asc":
		direction = "desc"
	_sort = [{"field": field, "direction": direction}]
	_refresh()
