class_name TableViewModel
extends RefCounted
## Everything the Tables tab shows, independent of the UI: columns, column
## filters, sorting, paging, row details, and CSV/TSV export.
##
## View filters only narrow what is displayed; they never change the graph.
## Rows come from MutationTable.build() on the Table node's frame.

const PER_SAMPLE_LIMIT := 60  ## Above this, per-sample columns start switched off.
const PER_SAMPLE_MAX := 300  ## Above this, per-sample columns are not offered.
const PAGE_SIZES := [100, 250, 500, 1000]

## Column kinds and their filter specs:
##   number     {min, max}
##   text       {contains}               case-insensitive substring
##   genes      {genes: [...]}           any overlap; "(none)" = no gene
##   flag       {value: bool}
##   sample_af  {mode: "carried" | "not_carried" | "range", min, max}
var columns: Array = []  ## [{id, title, kind, sample_id?, flag?}]
var filters := {}  ## column id -> filter spec
var sort := {"column": "position", "descending": false}
var page := 0
var page_size := 250
var per_sample_columns := true

var frame: MutationFrame
var rows: Array = []  ## All rows from MutationTable.build, identity order.
var view: Array = []  ## Rows after filters and sorting.
var _af: Array = []  ## Per row: {sample_id: af or null}


func _init(source: MutationFrame) -> void:
	frame = source
	rows = MutationTable.build(frame, {})["rows"]
	for i in rows.size():
		var row: Dictionary = rows[i]
		row["_index"] = i
		var by_sample := {}
		for result in row["sample_results"]:
			by_sample[result["sample_id"]] = result["af"]
		_af.append(by_sample)
	per_sample_columns = frame.samples.size() <= PER_SAMPLE_LIMIT
	_build_columns()
	refresh()


func set_per_sample_columns(enabled: bool) -> void:
	per_sample_columns = enabled and frame.samples.size() <= PER_SAMPLE_MAX
	_build_columns()
	for id in filters.keys():
		if column(id).is_empty():
			filters.erase(id)
	refresh()


func column(id: String) -> Dictionary:
	for c in columns:
		if c["id"] == id:
			return c
	return {}


func set_filter(id: String, spec: Dictionary) -> void:
	if spec.is_empty():
		filters.erase(id)
	else:
		filters[id] = spec
	page = 0
	refresh()


func clear_filters() -> void:
	filters.clear()
	page = 0
	refresh()


func set_sort(id: String, descending: bool) -> void:
	sort = {"column": id, "descending": descending}
	refresh()


func page_count() -> int:
	return maxi(1, ceili(float(view.size()) / page_size))


func page_rows() -> Array:
	page = clampi(page, 0, page_count() - 1)
	return view.slice(page * page_size, (page + 1) * page_size)


## Re-applies filters and sorting.
func refresh() -> void:
	var indices: Array = []
	for i in rows.size():
		if _passes(i):
			indices.append(i)
	var col := column(sort["column"])
	if not col.is_empty():
		var keys := {}
		for i in indices:
			keys[i] = _sort_key(i, col)
		var desc: bool = sort["descending"]
		indices.sort_custom(func(a: int, b: int) -> bool:
			var ka = keys[a]
			var kb = keys[b]
			if ka != kb:
				# Missing values (null) always sort last.
				if ka == null:
					return false
				if kb == null:
					return true
				return ka > kb if desc else ka < kb
			return a < b)
	view = indices.map(func(i): return rows[i])
	page = clampi(page, 0, page_count() - 1)


const COMPACT_SAMPLES_SHOWN := 12

## `full` is used by exports; the screen shortens very long sample lists.
func cell_text(row: Dictionary, col: Dictionary, full := false) -> String:
	var s: Dictionary = row["sample_summary"]
	match col["id"]:
		"position": return str(row["position"])
		"event_ref": return row["event_ref"]
		"alt": return row["alt"]
		"genes": return ", ".join(row["genes"]) if row["genes"] else row["region"].replace("_", " ")
		"carriers": return "%d/%d" % [s["carrier_count"], s["selected_count"]]
		"presence": return String.num(s["presence_percent"], 1)
		"af_passing": return "" if s["af_passing_count"] == null else "%d" % s["af_passing_count"]
		"slippage": return row["slippage"]["footprint"] if full else compact_slippage(row["slippage"])
		"samples":
			var results: Array = row["sample_results"]
			var shown := results if full or results.size() <= COMPACT_SAMPLES_SHOWN else results.slice(0, COMPACT_SAMPLES_SHOWN)
			var text := "; ".join(shown.map(func(r):
				return "%s %s" % [r["sample_label"], "AF n/a" if r["af"] == null else String.num(r["af"], 4)]))
			if shown.size() < results.size():
				text += "; … +%d more" % (results.size() - shown.size())
			return text
	match col["kind"]:
		"flag":
			return "yes" if row.get(col["flag"], false) else "no"
		"reason":
			return row.get("flag_reasons", {}).get(col["flag"], "")
		"sample_af":
			var by_sample: Dictionary = _af_for(row)
			if not by_sample.has(col["sample_id"]):
				return ""
			var af = by_sample[col["sample_id"]]
			return "n/a" if af == null else String.num(af, 4)
	return ""


## Short on-screen form of a slippage summary: reference context around the
## position, then any repeated motifs in the ALT, e.g. "CACCAA|CCCCCC  CC x2".
static func compact_slippage(slippage: Dictionary) -> String:
	var signals: Dictionary = slippage.get("signals", {})
	if signals.is_empty():
		return ""
	var before := ",".join(signals.get("REFERENCE_6_BEFORE", PackedStringArray()))
	var after := ",".join(signals.get("REFERENCE_6_AFTER", PackedStringArray()))
	var text := "%s|%s" % [before, after]
	var repeats: Array = []
	for pair in [["REPEAT_1_BASE", "REPEAT_1_BASE_COUNT"], ["REPEAT_2_BASES", "REPEAT_2_BASES_COUNT"], ["REPEAT_3_BASES", "REPEAT_3_BASES_COUNT"]]:
		if signals.has(pair[0]):
			var counts := Array(signals.get(pair[1], PackedStringArray())).filter(func(c): return c.to_int() >= 2)
			if counts:
				repeats.append("%s x%s" % [",".join(signals[pair[0]]), ",".join(counts)])
	if repeats:
		text += "  " + "  ".join(repeats)
	return text


## Calls behind one row, for the details view.
func row_details(row: Dictionary) -> Array:
	var store := frame.store
	var features := CallFeatures.new(store)
	var out: Array = []
	for result in row["sample_results"]:
		var call: int = result["call_index"]
		var c := store.calls
		out.append({
			"sample": result["sample_label"],
			"subject": store.sample_by_id(result["sample_id"]).get("subject", ""),
			"af": result["af"],
			"filter": result["filter"],
			"vcf": "%d %s>%s" % [c["vcf_pos"][call], store.strings[c["vcf_ref"][call]], store.strings[c["vcf_alt"][call]]],
			"placement": VariantStore.CALL_STATUS_NAMES[c["status"][call]].capitalize(),
			"multi_alt": c["multi_alt"][call] == 1,
			"format": store.strings[c["format"][call]],
			"value": c["values"][call],
			"metadata": features.metadata(call),
		})
	return out


## The whole filtered, sorted table (all pages) as CSV or TSV text.
func export_text(separator: String = ",") -> String:
	var lines := PackedStringArray()
	lines.append(separator.join(columns.map(func(c):
		return _quote("Slippage footprint" if c["id"] == "slippage" else c["title"], separator))))
	for row in view:
		lines.append(separator.join(columns.map(func(c): return _quote(cell_text(row, c, true), separator))))
	return "\n".join(lines) + "\n"


func describe_filters() -> String:
	var parts: Array = []
	for id in filters:
		var col := column(id)
		var spec: Dictionary = filters[id]
		match col.get("kind"):
			"number":
				var lo = spec.get("min")
				var hi = spec.get("max")
				if lo != null and hi != null:
					parts.append("%s %s–%s" % [col["title"], _num(lo), _num(hi)])
				elif lo != null:
					parts.append("%s at least %s" % [col["title"], _num(lo)])
				elif hi != null:
					parts.append("%s at most %s" % [col["title"], _num(hi)])
			"text":
				parts.append("%s contains \"%s\"" % [col["title"], spec["contains"]])
			"genes":
				parts.append("%s: %s" % [col["title"], ", ".join(spec["genes"])])
			"flag":
				parts.append("%s = %s" % [col["title"], "yes" if spec["value"] else "no"])
			"sample_af":
				match spec["mode"]:
					"carried": parts.append("%s carries it" % col["title"])
					"not_carried": parts.append("%s does not carry it" % col["title"])
					_: parts.append("%s AF %s–%s" % [col["title"], spec.get("min", 0), spec.get("max", 1)])
	return " · ".join(parts)


func _build_columns() -> void:
	columns = [
		{"id": "position", "title": "Position", "kind": "number"},
		{"id": "event_ref", "title": "Ref", "kind": "text"},
		{"id": "alt", "title": "ALT", "kind": "text"},
		{"id": "genes", "title": "Gene / region", "kind": "genes"},
		{"id": "carriers", "title": "Carriers", "kind": "number"},
		{"id": "presence", "title": "Presence %", "kind": "number"},
	]
	if rows.any(func(r): return r["sample_summary"]["af_passing_count"] != null):
		columns.append({"id": "af_passing", "title": "AF passing", "kind": "number"})
	for name in frame.flags:
		columns.append({"id": "flag:" + name, "title": name, "kind": "flag", "flag": name})
		if frame.flag_reasons.has(name):
			columns.append({"id": "why:" + name, "title": name + ": why not", "kind": "reason", "flag": name})
	columns.append({"id": "slippage", "title": "Slippage (context | ALT repeats)", "kind": "text"})
	if per_sample_columns:
		for sample_id in frame.samples:
			var sample := frame.store.sample_by_id(sample_id)
			columns.append({"id": "sample:%d" % sample_id, "title": sample.get("name", str(sample_id)), "kind": "sample_af", "sample_id": sample_id})
	else:
		columns.append({"id": "samples", "title": "Samples (AF)", "kind": "text"})


func _passes(i: int) -> bool:
	var row: Dictionary = rows[i]
	for id in filters:
		var col := column(id)
		var spec: Dictionary = filters[id]
		match col.get("kind"):
			"number":
				var value = _sort_key(i, col)
				if value == null:
					return false
				if spec.get("min") != null and value < spec["min"]:
					return false
				if spec.get("max") != null and value > spec["max"]:
					return false
			"text", "reason":
				if not cell_text(row, col).to_lower().contains(str(spec.get("contains", "")).to_lower()):
					return false
			"genes":
				var wanted: Array = spec.get("genes", [])
				var genes: Array = Array(row["genes"])
				if genes.is_empty():
					genes = ["(none)"]
				if wanted and not wanted.any(func(g): return g in genes):
					return false
			"flag":
				if row.get(col["flag"], false) != spec.get("value", true):
					return false
			"sample_af":
				var by_sample: Dictionary = _af[i]
				var carried := by_sample.has(col["sample_id"])
				match spec.get("mode", "carried"):
					"carried":
						if not carried:
							return false
					"not_carried":
						if carried:
							return false
					_:
						var af = by_sample.get(col["sample_id"])
						if af == null or af < spec.get("min", 0.0) or af > spec.get("max", 1.0):
							return false
	return true


func _sort_key(i: int, col: Dictionary):
	var row: Dictionary = rows[i]
	var s: Dictionary = row["sample_summary"]
	match col["id"]:
		"position": return row["position"]
		"carriers": return s["carrier_count"]
		"presence": return s["presence_percent"]
		"af_passing": return s["af_passing_count"]
	match col["kind"]:
		"flag": return 1 if row.get(col["flag"], false) else 0
		"sample_af":
			return _af[i].get(col["sample_id"], null) if _af[i].has(col["sample_id"]) else null
	return cell_text(row, col).to_lower()


func _af_for(row: Dictionary) -> Dictionary:
	return _af[row["_index"]]


static func _num(value) -> String:
	return str(int(value)) if float(value) == floorf(float(value)) else str(value)


## Repeated motifs (count >= 2) in a call's ALT, e.g. "C x6, CA x2", or "none".
static func describe_repeats(meta: Dictionary) -> String:
	var parts: Array = []
	for pair in [["REPEAT_1_BASE", "REPEAT_1_BASE_COUNT"], ["REPEAT_2_BASES", "REPEAT_2_BASES_COUNT"], ["REPEAT_3_BASES", "REPEAT_3_BASES_COUNT"]]:
		if meta.has(pair[0]) and str(meta.get(pair[1], "0")).to_int() >= 2:
			parts.append("%s x%s" % [meta[pair[0]], meta[pair[1]]])
	return ", ".join(parts) if parts else "none"


static func _quote(text: String, separator: String) -> String:
	if separator == "\t":
		return text.replace("\t", " ").replace("\n", " ")
	if text.contains(",") or text.contains("\"") or text.contains("\n"):
		return "\"%s\"" % text.replace("\"", "\"\"")
	return text
