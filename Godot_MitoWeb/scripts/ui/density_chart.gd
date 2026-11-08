class_name DensityChart
extends Control
## Mutation density along the genome, one 2px line per group (sliding
## windows through window midpoints; genes/regions as a segment spanning each
## region). One y-axis, recessive grid, a gene track under the x-axis, a
## legend for 2+ groups plus direct labels at the line ends for up to 4, and
## a crosshair tooltip on hover.

## Dark-mode categorical slots from the validated reference palette, in fixed
## order (checked with the dataviz validator against SURFACE).
const SERIES := [Color("#3987e5"), Color("#d95926"), Color("#199e70"), Color("#c98500"),
	Color("#d55181"), Color("#008300"), Color("#9085e9"), Color("#e66767")]
const SURFACE := Color("#1a1a19")
const TEXT_PRIMARY := Color("#f0efec")
const TEXT_SECONDARY := Color("#c3c2b7")
const GRID := Color("#33332f")
const AXIS := Color("#5a5955")
const TRACK := Color("#2a2a27")
const TRACK_ALT := Color("#34342f")
const MARGIN_LEFT := 64.0
const MARGIN_RIGHT := 150.0
const MARGIN_TOP := 46.0
const MARGIN_BOTTOM := 64.0

var windows: Array = []
var series: Array = []  ## [{"label": String, "values": PackedFloat64Array}]
var regions_mode := false
var y_title := ""
var range_start := 1
var range_end := WindowStats.LENGTH
var _hover := -1


func _init() -> void:
	custom_minimum_size = Vector2(600, 320)
	mouse_filter = Control.MOUSE_FILTER_STOP
	clip_contents = true


func set_data(new_windows: Array, new_series: Array, is_regions: bool, title: String, lo: int, hi: int) -> void:
	windows = new_windows
	series = new_series
	regions_mode = is_regions
	y_title = title
	range_start = lo
	range_end = hi
	_hover = -1
	queue_redraw()


func _plot_rect() -> Rect2:
	return Rect2(MARGIN_LEFT, MARGIN_TOP, maxf(size.x - MARGIN_LEFT - MARGIN_RIGHT, 50), maxf(size.y - MARGIN_TOP - MARGIN_BOTTOM, 50))


func _x_for(position: float, plot: Rect2) -> float:
	return plot.position.x + (position - range_start) / maxf(range_end - range_start, 1) * plot.size.x


func _y_max() -> float:
	var top := 0.0
	for s in series:
		for v in s["values"]:
			top = maxf(top, v)
	return _nice_ceiling(top)


func _midpoint(w: Dictionary) -> float:
	return minf(w["start"] + (w["length"] - 1) / 2.0, range_end)


func _draw() -> void:
	var font := get_theme_default_font()
	var fs := get_theme_default_font_size() - 2
	draw_rect(Rect2(Vector2.ZERO, size), SURFACE)
	var plot := _plot_rect()
	if windows.is_empty() or series.is_empty():
		draw_string(font, plot.get_center() - Vector2(120, 0), "Connect groups to a Statistics node.", HORIZONTAL_ALIGNMENT_LEFT, -1, fs, TEXT_SECONDARY)
		return
	var y_max := _y_max()
	# Grid and y-axis labels (recessive).
	for i in 5:
		var value := y_max * i / 4.0
		var y := plot.end.y - plot.size.y * i / 4.0
		draw_line(Vector2(plot.position.x, y), Vector2(plot.end.x, y), GRID if i > 0 else AXIS, 1.0)
		draw_string(font, Vector2(4, y + fs / 3.0), _fmt(value), HORIZONTAL_ALIGNMENT_RIGHT, MARGIN_LEFT - 10, fs, TEXT_SECONDARY)
	draw_string(font, Vector2(plot.position.x, MARGIN_TOP - 10), y_title, HORIZONTAL_ALIGNMENT_LEFT, -1, fs, TEXT_SECONDARY)
	# X ticks.
	var span := range_end - range_start
	var tick := 2000 if span > 8000 else (1000 if span > 3000 else (200 if span > 800 else 50))
	var t := ceili(float(range_start) / tick) * tick
	while t <= range_end:
		var x := _x_for(t, plot)
		draw_line(Vector2(x, plot.end.y), Vector2(x, plot.end.y + 4), AXIS, 1.0)
		draw_string(font, Vector2(x - 30, plot.end.y + 16), str(t), HORIZONTAL_ALIGNMENT_CENTER, 60, fs, TEXT_SECONDARY)
		t += tick
	# Gene track.
	var track_y := plot.end.y + 24
	var k := 0
	for gene in MitoGenes.GENES:
		if gene[3] < range_start or gene[2] > range_end:
			continue
		var x0 := _x_for(maxi(gene[2], range_start), plot)
		var x1 := _x_for(mini(gene[3], range_end), plot)
		draw_rect(Rect2(x0, track_y, maxf(x1 - x0, 1), 14), TRACK if k % 2 == 0 else TRACK_ALT)
		var name: String = gene[0].trim_prefix("MT-")
		if x1 - x0 > font.get_string_size(name, HORIZONTAL_ALIGNMENT_LEFT, -1, fs - 2).x + 4:
			draw_string(font, Vector2(x0 + 2, track_y + 11), name, HORIZONTAL_ALIGNMENT_LEFT, x1 - x0 - 2, fs - 2, TEXT_SECONDARY)
		k += 1
	draw_string(font, Vector2(4, track_y + 11), "genes", HORIZONTAL_ALIGNMENT_RIGHT, MARGIN_LEFT - 10, fs - 2, TEXT_SECONDARY)
	# Series.
	for si in series.size():
		var color: Color = SERIES[si % SERIES.size()]
		var values: PackedFloat64Array = series[si]["values"]
		if regions_mode:
			for wi in windows.size():
				var w: Dictionary = windows[wi]
				var y := plot.end.y - values[wi] / y_max * plot.size.y if y_max > 0 else plot.end.y
				for segment in w["segments"]:
					if segment[1] < range_start or segment[0] > range_end:
						continue
					draw_line(Vector2(_x_for(maxi(segment[0], range_start), plot), y), Vector2(_x_for(mini(segment[1], range_end), plot), y), color, 2.0, true)
		else:
			var points := PackedVector2Array()
			for wi in windows.size():
				var y := plot.end.y - values[wi] / y_max * plot.size.y if y_max > 0 else plot.end.y
				points.append(Vector2(_x_for(_midpoint(windows[wi]), plot), y))
			if points.size() >= 2:
				draw_polyline(points, color, 2.0, true)
			elif points.size() == 1:
				draw_circle(points[0], 4.0, color)
	# Legend (2+ groups) and direct labels at the line ends (up to 4).
	if series.size() >= 2:
		var lx := plot.position.x
		for si in series.size():
			var label: String = series[si]["label"]
			draw_rect(Rect2(lx, 8, 10, 10), SERIES[si % SERIES.size()])
			draw_string(font, Vector2(lx + 14, 17), label, HORIZONTAL_ALIGNMENT_LEFT, -1, fs, TEXT_PRIMARY)
			lx += 30 + font.get_string_size(label, HORIZONTAL_ALIGNMENT_LEFT, -1, fs).x
	if series.size() <= 4 and not regions_mode:
		var used: Array = []
		for si in series.size():
			var values: PackedFloat64Array = series[si]["values"]
			var y := plot.end.y - values[values.size() - 1] / y_max * plot.size.y if y_max > 0 else plot.end.y
			for u in used:
				if absf(u - y) < fs + 2:
					y = u + fs + 2
			used.append(y)
			var text: String = series[si].get("short", series[si]["label"])
			var room := MARGIN_RIGHT - 14
			while text.length() > 4 and font.get_string_size(text, HORIZONTAL_ALIGNMENT_LEFT, -1, fs).x > room:
				text = text.left(-2) + "…" if not text.ends_with("…") else text.left(-3) + "…"
			draw_string(font, Vector2(plot.end.x + 8, y + fs / 3.0), text, HORIZONTAL_ALIGNMENT_LEFT, -1, fs, TEXT_PRIMARY)
	# Hover crosshair and tooltip.
	if _hover >= 0 and _hover < windows.size():
		var w: Dictionary = windows[_hover]
		var hx := _x_for(_midpoint(w), plot)
		draw_line(Vector2(hx, plot.position.y), Vector2(hx, plot.end.y), TEXT_SECONDARY, 1.0)
		var lines: Array = ["%s  (%d bp)" % [w["label"], w["length"]]]
		for si in series.size():
			lines.append("%s: %s" % [series[si]["label"], _fmt(series[si]["values"][_hover])])
		var width := 0.0
		for line in lines:
			width = maxf(width, font.get_string_size(line, HORIZONTAL_ALIGNMENT_LEFT, -1, fs).x)
		var box := Rect2(hx + 12, plot.position.y + 8, width + 28, lines.size() * (fs + 6) + 10)
		if box.end.x > size.x - 4:
			box.position.x = hx - 12 - box.size.x
		draw_rect(box, Color("#262624"))
		draw_rect(box, AXIS, false, 1.0)
		for i in lines.size():
			var ty := box.position.y + 8 + (i + 1) * (fs + 6) - 6
			if i > 0:
				draw_rect(Rect2(box.position.x + 8, ty - fs + 3, 8, 8), SERIES[(i - 1) % SERIES.size()])
			draw_string(font, Vector2(box.position.x + (22 if i > 0 else 8), ty), lines[i], HORIZONTAL_ALIGNMENT_LEFT, -1, fs, TEXT_PRIMARY)


func _gui_input(event: InputEvent) -> void:
	if event is InputEventMouseMotion and not windows.is_empty():
		var plot := _plot_rect()
		var best := -1
		var best_d := INF
		if plot.grow(6).has_point(event.position):
			for wi in windows.size():
				var d := absf(_x_for(_midpoint(windows[wi]), plot) - event.position.x)
				if d < best_d:
					best_d = d
					best = wi
		if best != _hover:
			_hover = best
			queue_redraw()


func _notification(what: int) -> void:
	if what == NOTIFICATION_MOUSE_EXIT and _hover != -1:
		_hover = -1
		queue_redraw()
	elif what == NOTIFICATION_RESIZED:
		queue_redraw()


static func _nice_ceiling(value: float) -> float:
	if value <= 0:
		return 1.0
	var magnitude := pow(10.0, floor(log(value) / log(10.0)))
	for step in [1.0, 2.0, 2.5, 5.0, 10.0]:
		if value <= step * magnitude:
			return step * magnitude
	return 10.0 * magnitude


static func _fmt(value: float) -> String:
	if value == 0:
		return "0"
	if absf(value) >= 100:
		return String.num(value, 0)
	if absf(value) >= 1:
		return String.num(value, 2)
	return String.num(value, 4)
