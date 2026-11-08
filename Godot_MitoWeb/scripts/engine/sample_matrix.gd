class_name SampleMatrix
extends RefCounted
## Per-sample values in each window: the samples x windows feature matrix.
##
##   calls    number of the sample's calls in the window
##   af_load  sum of the sample's call AF in the window (missing AF adds 0)
## Optionally per kb of window. Rows follow the frame's sample order.

const MEASURES := ["calls", "af_load"]


## Returns {"sample_ids": PackedInt32Array, "rows": Array[PackedFloat64Array]}.
static func build(frame: MutationFrame, window_list: Array, measure: String = "calls", per_kb: bool = false) -> Dictionary:
	# Frame rows are sorted by position, so appending calls row by row gives
	# each sample a sorted position list; window sums are then two binary
	# searches on it (plus a prefix sum for AF).
	var positions := {}
	var af_values := {}
	for sample_id in frame.samples:
		positions[sample_id] = PackedInt32Array()
		af_values[sample_id] = PackedFloat64Array()
	var af: PackedFloat64Array = frame.store.calls["af"]
	for r in frame.row_count():
		var pos: int = frame.store.allele_pos[frame.alleles[r]]
		for sample_id in frame.carriers[r]:
			if not positions.has(sample_id):
				continue
			positions[sample_id].append(pos)
			var v := af[frame.carriers[r][sample_id]]
			af_values[sample_id].append(0.0 if is_nan(v) else v)
	var rows: Array = []
	for sample_id in frame.samples:
		var pos_list: PackedInt32Array = positions[sample_id]
		var cumulative := PackedFloat64Array([0.0])
		for v in af_values[sample_id]:
			cumulative.append(cumulative[-1] + v)
		var row := PackedFloat64Array()
		row.resize(window_list.size())
		for wi in window_list.size():
			var w: Dictionary = window_list[wi]
			var total := 0.0
			for segment in w["segments"]:
				var lo := pos_list.bsearch(segment[0], true)
				var hi := pos_list.bsearch(segment[1], false)
				total += float(hi - lo) if measure == "calls" else cumulative[hi] - cumulative[lo]
			row[wi] = total * 1000.0 / maxi(w["length"], 1) if per_kb else total
		rows.append(row)
	return {"sample_ids": frame.samples, "rows": rows}


## CSV: one row per sample with its id columns, group, every tag as a 0/1
## column (labels for models), then one column per window.
static func to_csv(store: VariantStore, groups: Array, window_list: Array, measure: String, per_kb: bool) -> String:
	var all_tags: Array = []
	for g in groups:
		for sample_id in g["matrix"]["sample_ids"]:
			for tag in store.sample_by_id(sample_id).get("tags", PackedStringArray()):
				if not tag in all_tags:
					all_tags.append(tag)
	all_tags.sort()
	var header := ["sample", "subject", "group"]
	header.append_array(all_tags.map(func(t): return "tag:" + t))
	header.append_array(window_list.map(func(w): return "%s %s" % [w["label"], ("per kb" if per_kb else measure)]))
	var lines := PackedStringArray([",".join(header.map(func(h): return TableViewModel._quote(h, ",")))])
	for g in groups:
		var ids: PackedInt32Array = g["matrix"]["sample_ids"]
		for i in ids.size():
			var sample := store.sample_by_id(ids[i])
			var tags: PackedStringArray = sample.get("tags", PackedStringArray())
			var cells := [sample.get("name", str(ids[i])), sample.get("subject", ""), g["label"]]
			cells.append_array(all_tags.map(func(t): return "1" if t in tags else "0"))
			for v in g["matrix"]["rows"][i]:
				cells.append(str(v))
			lines.append(",".join(cells.map(func(c): return TableViewModel._quote(str(c), ","))))
	return "\n".join(lines) + "\n"
