class_name WindowImport
extends RefCounted
## Reads a window list from a table:
##   - CSV/TSV with a header naming "start" and "end" (1-based, inclusive),
##     and optionally "label" or "name"
##   - BED (no header; chrom, start, end, optional name; 0-based start)
## start > end means the window wraps past 16,569 to 1. Returns
## {"ok", "error", "intervals": [[label, start, end], ...]}.

const LENGTH := 16569


static func parse(text: String) -> Dictionary:
	text = text.trim_prefix("﻿").replace("\r\n", "\n").replace("\r", "\n")
	var lines := Array(text.split("\n")).filter(func(l): return l.strip_edges() != "" and not l.begins_with("#") and not l.begins_with("track") and not l.begins_with("browser"))
	if lines.is_empty():
		return _fail("The file has no windows.")
	var tab: bool = lines[0].count("\t") >= lines[0].count(",")
	var first: Array = Array(lines[0].split("\t" if tab else ",")).map(func(c): return c.strip_edges().to_lower())
	var intervals: Array = []
	if "start" in first and "end" in first:
		var s_col := first.find("start")
		var e_col := first.find("end")
		var l_col := first.find("label") if "label" in first else first.find("name")
		for i in range(1, lines.size()):
			var cells: Array = Array(lines[i].split("\t" if tab else ",")).map(func(c): return c.strip_edges())
			if cells.size() <= maxi(s_col, e_col) or not cells[s_col].is_valid_int() or not cells[e_col].is_valid_int():
				return _fail("Line %d: start and end must be whole numbers." % (i + 1))
			var label: String = cells[l_col] if l_col >= 0 and l_col < cells.size() and cells[l_col] else ""
			intervals.append([label, cells[s_col].to_int(), cells[e_col].to_int()])
	else:
		for i in lines.size():
			var cells: Array = Array(lines[i].split("\t" if tab else ",")).map(func(c): return c.strip_edges())
			if cells.size() < 3 or not cells[1].is_valid_int() or not cells[2].is_valid_int():
				return _fail("Use a header with start and end columns, or BED lines (chrom, start, end).")
			intervals.append([cells[3] if cells.size() > 3 else "", cells[1].to_int() + 1, cells[2].to_int()])
	for iv in intervals:
		if iv[1] < 1 or iv[1] > LENGTH or iv[2] < 1 or iv[2] > LENGTH:
			return _fail("Window %s–%s is outside 1–%d." % [iv[1], iv[2], LENGTH])
		if iv[0].is_empty():
			iv[0] = "%d–%d" % [iv[1], iv[2]]
	if intervals.size() > 20000:
		return _fail("Too many windows (%d)." % intervals.size())
	return {"ok": true, "error": "", "intervals": intervals}


static func _fail(message: String) -> Dictionary:
	return {"ok": false, "error": message, "intervals": []}
