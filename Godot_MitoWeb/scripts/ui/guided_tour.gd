class_name GuidedTour
extends CanvasLayer
## A short guided tour for first-time visitors (the web demo is mostly seen by
## recruiters). One step at a time, a box points at part of the screen: the
## target gets a gold outline and an arrow, the rest of the screen dims a
## little, and the box has Back / Next and a close button. Clicks pass through
## to the app, so visitors can try things while it is open. Keys: ← / → step,
## Esc closes.
##
## A "?" button (with a soft pulsing glow while the tour is closed) reopens it.
##
## Each step: {"title": String, "text": String (BBCode),
##             "target": Callable -> Rect2 in screen pixels (empty = no target),
##             "before": Callable (optional; e.g. switch tabs first)}

signal closed

const GOLD := Color("e6c98f")
const INK := Color("ece7dc")
const MUTED := Color("a3adb5")
const PANEL := Color(0.07, 0.105, 0.135, 0.97)
const SEEN_FILE := "user://guided_tour.cfg"
const BOX_WIDTH := 360.0

var steps: Array = []
var _index := -1
var _overlay: Control
var _bubble: PanelContainer
var _title: Label
var _body: RichTextLabel
var _count: Label
var _back: Button
var _next: Button
var _help: Button
var _help_style: StyleBoxFlat
var _help_anchor: Callable          # -> Vector2 top-left, when the button floats
var _time := 0.0


func _init() -> void:
	layer = 100


func _ready() -> void:
	_overlay = Control.new()
	_overlay.set_anchors_preset(Control.PRESET_FULL_RECT)
	_overlay.mouse_filter = Control.MOUSE_FILTER_IGNORE
	_overlay.draw.connect(_draw_overlay)
	add_child(_overlay)
	_build_bubble()


## Starts the tour, unless this browser has already seen (or closed) it.
func start_if_new() -> void:
	var cfg := ConfigFile.new()
	if cfg.load(SEEN_FILE) == OK and cfg.get_value("tour", "seen", false):
		return
	start()


func start(from := 0) -> void:
	if steps.is_empty():
		return
	_show(clampi(from, 0, steps.size() - 1))


func stop() -> void:
	if _index < 0:
		return
	_index = -1
	_bubble.visible = false
	_overlay.queue_redraw()
	var cfg := ConfigFile.new()
	cfg.set_value("tour", "seen", true)
	cfg.save(SEEN_FILE)
	closed.emit()


func is_running() -> bool:
	return _index >= 0


## A "?" button that reopens the tour. With a container it is added there;
## otherwise it floats, placed each frame at anchor.call() (screen pixels).
func add_help_button(container: Node = null, anchor: Callable = Callable()) -> Button:
	_help = Button.new()
	_help.text = "?"
	_help.tooltip_text = "Guided tour: what am I looking at?"
	_help.custom_minimum_size = Vector2(34, 34)
	_help.focus_mode = Control.FOCUS_NONE
	_help_style = _box_style(Color(GOLD, 0.16), GOLD, 17)
	_help_style.shadow_color = Color(GOLD, 0.45)
	for state in ["normal", "hover", "pressed", "focus"]:
		_help.add_theme_stylebox_override(state, _help_style)
	_help.add_theme_color_override("font_color", GOLD)
	_help.add_theme_color_override("font_hover_color", INK)
	_help.add_theme_font_size_override("font_size", 18)
	_help.pressed.connect(func() -> void: start())
	if container != null:
		container.add_child(_help)
	else:
		_help_anchor = anchor
		add_child(_help)
	return _help


func _process(delta: float) -> void:
	_time += delta
	if _help != null:
		# Pulse while the tour is closed, so the "?" is easy to notice.
		var pulse := 0.0 if is_running() else 0.5 + 0.5 * sin(_time * 3.0)
		_help_style.shadow_size = int(2 + 7 * pulse)
		_help_style.shadow_color = Color(GOLD, 0.15 + 0.4 * pulse)
		if _help_anchor.is_valid():
			_help.position = _help_anchor.call()
	if _index >= 0:
		_place_bubble()
		_overlay.queue_redraw()


func _unhandled_input(event: InputEvent) -> void:
	if _index < 0 or not (event is InputEventKey and event.pressed):
		return
	match event.keycode:
		KEY_ESCAPE:
			stop()
		KEY_RIGHT:
			if _index >= steps.size() - 1:
				stop()
			else:
				_show(_index + 1)
		KEY_LEFT:
			_show(max(0, _index - 1))
		_:
			return
	get_viewport().set_input_as_handled()


# ── Steps ─────────────────────────────────────────────────────────────────────

func _show(i: int) -> void:
	_index = i
	var step: Dictionary = steps[i]
	var before: Callable = step.get("before", Callable())
	if before.is_valid():
		before.call()
	_title.text = step.get("title", "")
	_body.text = step.get("text", "")
	_count.text = "%d of %d" % [i + 1, steps.size()]
	_back.disabled = i == 0
	_next.text = "Done" if i == steps.size() - 1 else "Next"
	_bubble.visible = true
	_bubble.reset_size()
	_place_bubble()


func _target_rect() -> Rect2:
	if _index < 0:
		return Rect2()
	var target: Callable = steps[_index].get("target", Callable())
	return target.call() if target.is_valid() else Rect2()


# ── Drawing ───────────────────────────────────────────────────────────────────

func _place_bubble() -> void:
	var view := _overlay.get_viewport_rect().size
	var size := _bubble.get_combined_minimum_size()
	size.x = max(size.x, BOX_WIDTH)
	_bubble.size = size
	var r := _target_rect()
	var pos: Vector2
	if not r.has_area():
		pos = Vector2((view.x - size.x) / 2.0, view.y - size.y - 40.0)
	elif r.end.y + 18.0 + size.y < view.y - 8.0:
		pos = Vector2(r.get_center().x - size.x / 2.0, r.end.y + 18.0)               # below
	elif r.position.y - 18.0 - size.y > 8.0:
		pos = Vector2(r.get_center().x - size.x / 2.0, r.position.y - 18.0 - size.y)  # above
	elif r.end.x + 18.0 + size.x < view.x - 8.0:
		pos = Vector2(r.end.x + 18.0, r.get_center().y - size.y / 2.0)               # right
	else:
		pos = Vector2(r.position.x - 18.0 - size.x, r.get_center().y - size.y / 2.0)  # left
	pos.x = clampf(pos.x, 8.0, max(8.0, view.x - size.x - 8.0))
	pos.y = clampf(pos.y, 8.0, max(8.0, view.y - size.y - 8.0))
	_bubble.position = pos


func _draw_overlay() -> void:
	if _index < 0:
		return
	var view := _overlay.get_viewport_rect().size
	var dim := Color(0, 0, 0, 0.32)
	var r := _target_rect()
	if not r.has_area():
		_overlay.draw_rect(Rect2(Vector2.ZERO, view), dim)
		return
	var hole := r.grow(6.0)
	# Dim everything but the target.
	_overlay.draw_rect(Rect2(0, 0, view.x, hole.position.y), dim)
	_overlay.draw_rect(Rect2(0, hole.end.y, view.x, view.y - hole.end.y), dim)
	_overlay.draw_rect(Rect2(0, hole.position.y, hole.position.x, hole.size.y), dim)
	_overlay.draw_rect(Rect2(hole.end.x, hole.position.y, view.x - hole.end.x, hole.size.y), dim)
	# Pulsing outline on the target.
	var glow := 0.55 + 0.45 * sin(_time * 4.0)
	_overlay.draw_rect(hole.grow(3.0), Color(GOLD, 0.25 * glow), false, 6.0)
	_overlay.draw_rect(hole, Color(GOLD, 0.9), false, 2.0)
	# Arrow from the box to the target.
	var box := Rect2(_bubble.position, _bubble.size)
	var from := _closest_on_rect(box, hole.get_center())
	var to := _closest_on_rect(hole, from)
	if from.distance_to(to) > 6.0:
		_overlay.draw_line(from, to, GOLD, 2.0, true)
		var dir := (to - from).normalized()
		var side := Vector2(-dir.y, dir.x)
		_overlay.draw_colored_polygon(PackedVector2Array([to, to - dir * 11.0 + side * 6.0, to - dir * 11.0 - side * 6.0]), GOLD)


func _closest_on_rect(rect: Rect2, point: Vector2) -> Vector2:
	var p := Vector2(clampf(point.x, rect.position.x, rect.end.x), clampf(point.y, rect.position.y, rect.end.y))
	if not rect.has_point(point):
		return p
	# Inside: snap to the nearest edge.
	var d := [p.x - rect.position.x, rect.end.x - p.x, p.y - rect.position.y, rect.end.y - p.y]
	var m: float = d.min()
	if m == d[0]: return Vector2(rect.position.x, p.y)
	if m == d[1]: return Vector2(rect.end.x, p.y)
	if m == d[2]: return Vector2(p.x, rect.position.y)
	return Vector2(p.x, rect.end.y)


# ── The box ───────────────────────────────────────────────────────────────────

func _build_bubble() -> void:
	_bubble = PanelContainer.new()
	_bubble.visible = false
	_bubble.add_theme_stylebox_override("panel", _box_style(PANEL, Color(GOLD, 0.8), 6, 16))
	add_child(_bubble)
	var box := VBoxContainer.new()
	box.add_theme_constant_override("separation", 10)
	_bubble.add_child(box)

	var top := HBoxContainer.new()
	box.add_child(top)
	_title = Label.new()
	_title.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_title.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	_title.add_theme_color_override("font_color", GOLD)
	_title.add_theme_font_size_override("font_size", 17)
	top.add_child(_title)
	var close := _button("×", func() -> void: stop(), true)
	close.tooltip_text = "Close the tour (the ? button reopens it)"
	top.add_child(close)

	_body = RichTextLabel.new()
	_body.bbcode_enabled = true
	_body.fit_content = true
	_body.scroll_active = false
	_body.custom_minimum_size = Vector2(BOX_WIDTH - 32.0, 0)
	_body.add_theme_color_override("default_color", INK)
	_body.add_theme_font_size_override("normal_font_size", 15)
	_body.add_theme_font_size_override("bold_font_size", 15)
	box.add_child(_body)

	var nav := HBoxContainer.new()
	nav.add_theme_constant_override("separation", 10)
	box.add_child(nav)
	_back = _button("Back", func() -> void: _show(max(0, _index - 1)))
	nav.add_child(_back)
	_count = Label.new()
	_count.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_count.horizontal_alignment = HORIZONTAL_ALIGNMENT_CENTER
	_count.add_theme_color_override("font_color", MUTED)
	_count.add_theme_font_size_override("font_size", 13)
	nav.add_child(_count)
	_next = _button("Next", func() -> void:
		if _index >= steps.size() - 1:
			stop()
		else:
			_show(_index + 1))
	nav.add_child(_next)


func _button(text: String, action: Callable, flat := false) -> Button:
	var b := Button.new()
	b.text = text
	b.focus_mode = Control.FOCUS_NONE
	var normal := _box_style(Color(0, 0, 0, 0) if flat else Color(GOLD, 0.12), Color(0, 0, 0, 0) if flat else Color(GOLD, 0.6), 4, 0)
	var hover := _box_style(Color(GOLD, 0.25), Color(GOLD, 0.9), 4, 0)
	for s in [normal, hover]:
		s.content_margin_left = 6.0 if flat else 14.0
		s.content_margin_right = s.content_margin_left
		s.content_margin_top = 4.0
		s.content_margin_bottom = 4.0
	b.add_theme_stylebox_override("normal", normal)
	b.add_theme_stylebox_override("hover", hover)
	b.add_theme_stylebox_override("pressed", hover)
	b.add_theme_stylebox_override("disabled", _box_style(Color(0, 0, 0, 0), Color(MUTED, 0.25), 4, 0))
	b.add_theme_color_override("font_color", INK)
	b.add_theme_color_override("font_hover_color", INK)
	b.add_theme_color_override("font_disabled_color", Color(MUTED, 0.5))
	b.add_theme_font_size_override("font_size", 20 if flat else 14)
	b.pressed.connect(action)
	return b


func _box_style(bg: Color, border: Color, radius: int, margin := 0) -> StyleBoxFlat:
	var s := StyleBoxFlat.new()
	s.bg_color = bg
	s.border_color = border
	s.set_border_width_all(1)
	s.set_corner_radius_all(radius)
	s.set_content_margin_all(margin)
	s.anti_aliasing = true
	return s


# ── Helpers for building steps ────────────────────────────────────────────────

## Screen rect of a control, or empty if it isn't visible.
static func rect_of(control: Control) -> Rect2:
	if control == null or not control.is_visible_in_tree():
		return Rect2()
	return control.get_global_rect()


## Screen rect of one tab's title in a TabContainer.
static func tab_rect(tabs: TabContainer, index: int) -> Rect2:
	if tabs == null or not tabs.is_visible_in_tree():
		return Rect2()
	var bar := tabs.get_tab_bar()
	var r := bar.get_tab_rect(index)
	return Rect2(bar.get_global_rect().position + r.position, r.size)
