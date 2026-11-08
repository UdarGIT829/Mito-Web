class_name CheckListButton
extends MenuButton
## A dropdown of checkable items (tags, genes). The menu stays open while
## items are toggled.

signal selection_changed(selected: Array)

var empty_text := "None"
var _items: Array = []
var _selected: Array = []


func _init(placeholder: String = "None") -> void:
	empty_text = placeholder
	flat = false
	get_popup().hide_on_checkable_item_selection = false
	get_popup().id_pressed.connect(_on_id_pressed)
	_update_text()


func set_items(items: Array) -> void:
	_items = items.duplicate()
	var popup := get_popup()
	popup.clear()
	for i in _items.size():
		popup.add_check_item(str(_items[i]), i)
		popup.set_item_checked(i, _items[i] in _selected)
	if _items.is_empty():
		popup.add_item("(nothing to choose yet)", -1)
		popup.set_item_disabled(0, true)
	_update_text()


func set_selected(selected: Array) -> void:
	_selected = selected.duplicate()
	for i in _items.size():
		get_popup().set_item_checked(i, _items[i] in _selected)
	_update_text()


func get_selected() -> Array:
	return _selected.duplicate()


func _on_id_pressed(id: int) -> void:
	if id < 0 or id >= _items.size():
		return
	var item = _items[id]
	if item in _selected:
		_selected.erase(item)
	else:
		_selected.append(item)
	get_popup().set_item_checked(id, item in _selected)
	_update_text()
	selection_changed.emit(get_selected())


func _update_text() -> void:
	var shown := ", ".join(_selected.map(func(s): return str(s)))
	text = (shown if shown.length() <= 28 else shown.left(26) + "…") if _selected else empty_text
	tooltip_text = ", ".join(_selected.map(func(s): return str(s))) if _selected else empty_text
