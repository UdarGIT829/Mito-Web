class_name TagImport
extends RefCounted
## Assigns tags to imported samples from a table (CSV or TSV) such as a lab
## metadata sheet or a haplogroup list:
##   sample,cell_line,treatment,haplogroup
##   C10_EV_PCR...,C10,Treated,J1c
## One column identifies the sample (matched to sample name, then subject);
## each chosen column's value becomes a tag, optionally prefixed with the
## column name ("haplogroup=J1c").


## Parses CSV (with quotes) or TSV text. Returns {"ok", "error", "header",
## "rows"}; rows are Arrays of Strings padded to the header width.
static func parse(text: String) -> Dictionary:
	text = text.trim_prefix("﻿").replace("\r\n", "\n").replace("\r", "\n")
	var lines := Array(text.split("\n")).filter(func(l): return l.strip_edges() != "")
	if lines.size() < 2:
		return {"ok": false, "error": "The table needs a header row and at least one data row."}
	var tab: bool = lines[0].count("\t") > lines[0].count(",")
	var header := _split(lines[0], tab)
	header = header.map(func(h): return h.strip_edges())
	var rows: Array = []
	for i in range(1, lines.size()):
		var cells := _split(lines[i], tab).map(func(c): return c.strip_edges())
		while cells.size() < header.size():
			cells.append("")
		rows.append(cells)
	return {"ok": true, "error": "", "header": header, "rows": rows}


## Adds tags to matching samples (any dictionaries with "name", "subject"
## and "tags", e.g. imported samples or import-plan samples). Returns
## {"matched": int, "tags_added": int, "unmatched": [names]}. Tags already
## present are not duplicated.
## subject_column >= 0 also sets each matched sample's subject from that
## column (e.g. a cell line or donor), so analyses can treat a subject's
## samples as one unit.
##
## attribute_columns become per-sample values (sample["attributes"][column
## name]): numbers stay numbers (SLEDAI 8), other text stays text (sex F).
## haplogroup_column: the result lists [sample, text] pairs in
## "haplogroups" for the caller to apply (see Workspace.import_table).
static func apply(samples: Array, table: Dictionary, match_column: int, tag_columns: Array, prefix: bool, subject_column: int = -1,
		attribute_columns: Array = [], haplogroup_column: int = -1) -> Dictionary:
	var by_name := {}
	var by_subject := {}
	var by_lower := {}
	for sample in samples:
		by_name.get_or_add(sample["name"], []).append(sample)
		by_lower.get_or_add(sample["name"].to_lower(), []).append(sample)
		by_subject.get_or_add(sample["subject"], []).append(sample)
	var matched := 0
	var added := 0
	var attributes_set := 0
	var haplogroups: Array = []
	var unmatched: Array = []
	for row in table["rows"]:
		var id: String = row[match_column]
		if id.is_empty():
			continue
		var targets: Array = by_name.get(id, by_subject.get(id, by_lower.get(id.to_lower(), [])))
		if targets.is_empty():
			unmatched.append(id)
			continue
		matched += targets.size()
		for sample in targets:
			if subject_column >= 0 and row[subject_column] != "":
				sample["subject"] = row[subject_column]
			var tags := PackedStringArray(sample["tags"])
			for col in tag_columns:
				var value: String = row[col]
				if value.is_empty():
					continue
				var tag := "%s=%s" % [table["header"][col], value] if prefix else value
				if not tag in tags:
					tags.append(tag)
					added += 1
			sample["tags"] = tags
			for col in attribute_columns:
				var raw: String = row[col]
				if raw.is_empty():
					continue
				if not sample.has("attributes"):
					sample["attributes"] = {}
				sample["attributes"][table["header"][col]] = parse_value(raw)
				attributes_set += 1
			if haplogroup_column >= 0 and row[haplogroup_column] != "":
				haplogroups.append([sample, row[haplogroup_column]])
	return {"matched": matched, "tags_added": added, "attributes_set": attributes_set, "haplogroups": haplogroups, "unmatched": unmatched}


## A table cell as a value: a number when it reads as one, else trimmed text.
static func parse_value(raw: String) -> Variant:
	var text := raw.strip_edges()
	if text.is_valid_float() or text.is_valid_int():
		return text.to_float()
	return text


static func _split(line: String, tab: bool) -> Array:
	if tab:
		return Array(line.split("\t"))
	var out: Array = []
	var cell := ""
	var quoted := false
	var i := 0
	while i < line.length():
		var ch := line[i]
		if quoted:
			if ch == "\"":
				if i + 1 < line.length() and line[i + 1] == "\"":
					cell += "\""
					i += 1
				else:
					quoted = false
			else:
				cell += ch
		elif ch == "\"":
			quoted = true
		elif ch == ",":
			out.append(cell)
			cell = ""
		else:
			cell += ch
		i += 1
	out.append(cell)
	return out
