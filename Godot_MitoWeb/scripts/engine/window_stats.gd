class_name WindowStats
extends RefCounted
## Mutation density along the genome, in windows.
##
## A window set is plain data, a list of intervals, so it can come from a
## sliding scan, gene/region boundaries, or (later) an imported or learned list:
##   [{"label", "segments": [[start, end], ...], "length"}]   1-based, inclusive
## Segments allow windows that wrap past 16,569 back to 1 (circular genome).
##
## Measures per window and group:
##   mutations  distinct mutations whose position is in the window
##   calls      carrier calls (a mutation carried by 3 samples counts 3)
##   af_load    sum of carrier AF (heteroplasmy load); calls without AF add 0
## Each also as per kb, and per sample (value / samples in the group).

const LENGTH := 16569
const MEASURES := ["mutations", "calls", "af_load"]
const MEASURE_LABELS := {"mutations": "Distinct mutations", "calls": "Carrier calls", "af_load": "AF-weighted calls"}


## Validates a window spec. Returns "" or a message.
## spec: {mode: "sliding" | "regions", size, step, range_start, range_end, wrap}
static func validate(spec: Dictionary) -> String:
	if spec.get("mode", "sliding") == "regions":
		return ""
	if spec.get("mode") == "custom":
		var intervals: Array = spec.get("intervals", [])
		if intervals.is_empty():
			return "Load a window list (CSV with start/end, or BED)."
		for iv in intervals:
			if typeof(iv) != TYPE_ARRAY or iv.size() < 3 or int(iv[1]) < 1 or int(iv[1]) > LENGTH or int(iv[2]) < 1 or int(iv[2]) > LENGTH:
				return "The window list has an invalid interval."
		return ""
	var size = spec.get("size", 500)
	var step = spec.get("step", size)
	var lo = spec.get("range_start", 1)
	var hi = spec.get("range_end", LENGTH)
	for v in [size, step, lo, hi]:
		if typeof(v) not in [TYPE_INT, TYPE_FLOAT] or int(v) != v:
			return "Window size, step and range must be whole numbers."
	if size < 1 or size > LENGTH:
		return "Window size must be 1 to %d bp." % LENGTH
	if step < 1:
		return "Window step must be at least 1 bp."
	if lo < 1 or hi > LENGTH or lo > hi:
		return "The range must lie within 1–%d and start before it ends." % LENGTH
	var count := ceili(float(hi - lo + 1) / step)
	if count > 20000:
		return "That makes %d windows; use a larger step." % count
	return ""


static func windows(spec: Dictionary) -> Array:
	if spec.get("mode", "sliding") == "regions":
		return region_windows()
	if spec.get("mode") == "custom":
		return custom_windows(spec.get("intervals", []))
	var size: int = int(spec.get("size", 500))
	var step: int = int(spec.get("step", size))
	var lo: int = int(spec.get("range_start", 1))
	var hi: int = int(spec.get("range_end", LENGTH))
	var wrap: bool = spec.get("wrap", false) and lo == 1 and hi == LENGTH
	var out: Array = []
	var start := lo
	while start <= hi:
		var end := start + size - 1
		var segments: Array
		if end <= hi:
			segments = [[start, end]]
		elif wrap:
			segments = [[start, LENGTH], [1, end - LENGTH]]
		else:
			segments = [[start, hi]]
		var label := "%d–%d" % [start, end if end <= hi or not wrap else end - LENGTH]
		out.append({"label": label, "segments": segments, "length": _length(segments), "start": start})
		if end >= hi and not wrap:
			break
		start += step
	return out


## Windows from [label, start, end] intervals; start > end wraps past the
## end of the genome to position 1.
static func custom_windows(intervals: Array) -> Array:
	var out: Array = []
	for iv in intervals:
		var s := int(iv[1])
		var e := int(iv[2])
		var segments: Array = [[s, e]] if s <= e else [[s, LENGTH], [1, e]]
		out.append({"label": str(iv[0]), "segments": segments, "length": _length(segments), "start": s})
	return out


## Candidate windows for a scan: every size, stepped by a quarter of the
## size, around the whole (circular) genome.
static func scan_windows(sizes: Array) -> Array:
	var out: Array = []
	for size in sizes:
		out.append_array(windows({"size": int(size), "step": maxi(int(size) / 4, 1), "wrap": true}))
	return out


## Greedy choice of non-overlapping windows, best first (by `order`).
static func pick_non_overlapping(candidates: Array, order: Array, limit: int) -> Array:
	var taken := PackedByteArray()
	taken.resize(LENGTH + 1)
	var chosen: Array = []
	for i in order:
		var w: Dictionary = candidates[i]
		var free := true
		for segment in w["segments"]:
			for p in range(segment[0], segment[1] + 1):
				if taken[p]:
					free = false
					break
		if not free:
			continue
		for segment in w["segments"]:
			for p in range(segment[0], segment[1] + 1):
				taken[p] = 1
		chosen.append(i)
		if chosen.size() >= limit:
			break
	return chosen


## The 37 NC_012920.1 genes, the control region (16024–576, wrapping the
## origin) and the intergenic stretches between them.
static func region_windows() -> Array:
	var covered := PackedByteArray()
	covered.resize(LENGTH + 1)
	var out: Array = []
	for gene in MitoGenes.GENES:
		out.append({"label": gene[0], "segments": [[gene[2], gene[3]]], "length": gene[3] - gene[2] + 1, "start": gene[2]})
		for p in range(gene[2], gene[3] + 1):
			covered[p] = 1
	var control := [[16024, LENGTH], [1, 576]]
	out.append({"label": "control region", "segments": control, "length": _length(control), "start": 16024})
	for segment in control:
		for p in range(segment[0], segment[1] + 1):
			covered[p] = 1
	var p := 1
	while p <= LENGTH:
		if covered[p]:
			p += 1
			continue
		var s := p
		while p <= LENGTH and not covered[p]:
			p += 1
		out.append({"label": "intergenic %d–%d" % [s, p - 1], "segments": [[s, p - 1]], "length": p - s, "start": s})
	out.sort_custom(func(a, b): return a["start"] < b["start"])
	return out


## Counts for one group over a window set. Returns
## {"samples": n, "values": {measure: PackedFloat64Array per window}}.
static func count(frame: MutationFrame, window_list: Array) -> Dictionary:
	# Per-position totals, then prefix sums: O(mutations + windows).
	var per := {}
	for m in MEASURES:
		var arr := PackedFloat64Array()
		arr.resize(LENGTH + 1)
		per[m] = arr
	var af: PackedFloat64Array = frame.store.calls["af"]
	for r in frame.row_count():
		var pos: int = frame.store.allele_pos[frame.alleles[r]]
		if pos < 1 or pos > LENGTH:
			continue
		per["mutations"][pos] += 1.0
		for sample_id in frame.carriers[r]:
			per["calls"][pos] += 1.0
			var value := af[frame.carriers[r][sample_id]]
			if not is_nan(value):
				per["af_load"][pos] += value
	var prefix := {}
	for m in MEASURES:
		var sums := PackedFloat64Array()
		sums.resize(LENGTH + 1)
		var running := 0.0
		for i in range(1, LENGTH + 1):
			running += per[m][i]
			sums[i] = running
		prefix[m] = sums
	var values := {}
	for m in MEASURES:
		var column := PackedFloat64Array()
		for w in window_list:
			var total := 0.0
			for segment in w["segments"]:
				total += prefix[m][segment[1]] - prefix[m][segment[0] - 1]
			column.append(total)
		values[m] = column
	return {"samples": frame.samples.size(), "values": values}


## Value shown for one window: raw, per kb of window, and/or per sample.
static func scaled(value: float, window: Dictionary, samples: int, per_kb: bool, per_sample: bool) -> float:
	var v := value
	if per_kb:
		v = v * 1000.0 / maxi(window["length"], 1)
	if per_sample:
		v = v / maxi(samples, 1)
	return v


static func _length(segments: Array) -> int:
	var total := 0
	for s in segments:
		total += s[1] - s[0] + 1
	return total
