class_name PermutationTest
extends RefCounted
## Two-group permutation test per window, on one group's samples split by a
## tag (A vs B). Statistic: difference in mean per-unit value (A minus B).
##
## Units are samples, or subjects (each subject's samples averaged first),
## which is the right choice when one subject has several samples, e.g. one
## cell line prepared four ways. Labels are shuffled between units, keeping
## the group sizes; when every distinct labelling fits in the budget they are
## enumerated exactly.
##
## Per window: p = share of labellings with |difference| at least the
## observed one. Adjusted p (family-wise, single-step max-T): share of
## labellings whose largest |difference| over all windows reaches the
## window's observed value. The observed labelling is always counted.

signal progress(done: int, total: int)

const FRAME_BUDGET_MSEC := 25
const EPS := 1e-12

var cancelled := false


## Prepares units from a samples x windows matrix. Returns {"ok", "error",
## "values": Array[PackedFloat64Array] per unit, "labels": PackedByteArray
## (1 = A), "units": [names], "n_a", "n_b"}.
static func prepare(store: VariantStore, matrix: Dictionary, tag_a: String, tag_b: String, by_subject: bool) -> Dictionary:
	if tag_a == tag_b:
		return _fail("Choose two different tags.")
	var groups := {}  # unit name -> {"label": 0/1, "rows": [...]}
	var ids: PackedInt32Array = matrix["sample_ids"]
	for i in ids.size():
		var sample := store.sample_by_id(ids[i])
		var tags: PackedStringArray = sample.get("tags", PackedStringArray())
		var in_a := tag_a in tags
		var in_b := tag_b in tags
		if in_a and in_b:
			return _fail("Sample %s has both tags." % sample.get("name", ""))
		if not in_a and not in_b:
			continue
		var unit: String = sample.get("subject", "") if by_subject else sample.get("name", str(ids[i]))
		var label := 1 if in_a else 0
		var entry: Dictionary = groups.get_or_add(unit, {"label": label, "rows": []})
		if entry["label"] != label:
			return _fail("Subject %s has samples in both groups; compare by sample, or fix the tags." % unit)
		entry["rows"].append(matrix["rows"][i])
	var values: Array = []
	var labels := PackedByteArray()
	var names: Array = groups.keys()
	names.sort()
	for unit in names:
		var rows: Array = groups[unit]["rows"]
		var mean: PackedFloat64Array = rows[0].duplicate()
		for k in range(1, rows.size()):
			for w in mean.size():
				mean[w] += rows[k][w]
		for w in mean.size():
			mean[w] /= rows.size()
		values.append(mean)
		labels.append(groups[unit]["label"])
	var n_a := labels.count(1)
	var n_b := labels.count(0)
	if n_a == 0 or n_b == 0:
		return _fail("Each group needs at least one %s (found %d with %s, %d with %s)." % ["subject" if by_subject else "sample", n_a, tag_a, n_b, tag_b])
	return {"ok": true, "error": "", "values": values, "labels": labels, "units": names, "n_a": n_a, "n_b": n_b}


## Number of distinct labellings: C(n_a + n_b, n_a).
static func labelling_count(n_a: int, n_b: int) -> float:
	var n := n_a + n_b
	var k := mini(n_a, n_b)
	var c := 1.0
	for i in k:
		c = c * (n - i) / (i + 1)
	return c


## Runs the test. Returns {"ok", "error", "cancelled", "exact", "labellings",
## "mean_a", "mean_b", "difference", "p", "p_adjusted" (PackedFloat64Arrays
## per window), "min_p"}.
func run(prepared: Dictionary, permutations: int = 999, seed: int = 1) -> Dictionary:
	var values: Array = prepared["values"]
	var labels: PackedByteArray = prepared["labels"]
	var n := labels.size()
	var windows: int = values[0].size()
	var n_a: int = prepared["n_a"]
	var possible := labelling_count(n_a, n - n_a)
	var exact := possible <= permutations + 1
	var totals := PackedFloat64Array()
	totals.resize(windows)
	for u in n:
		for w in windows:
			totals[w] += values[u][w]
	var observed := _differences(values, labels, totals, n_a, n)
	var mean_a := PackedFloat64Array()
	var mean_b := PackedFloat64Array()
	var sums_a := _sums(values, labels)
	for w in windows:
		mean_a.append(sums_a[w] / n_a)
		mean_b.append((totals[w] - sums_a[w]) / (n - n_a))
	var abs_obs := PackedFloat64Array()
	for d in observed:
		abs_obs.append(absf(d))
	var exceed := PackedInt32Array()
	exceed.resize(windows)
	var exceed_max := PackedInt32Array()
	exceed_max.resize(windows)
	var count := 0
	var labellings: Array = _all_labellings(n, n_a) if exact else []
	var total := labellings.size() if exact else permutations + 1
	var rng := RandomNumberGenerator.new()
	rng.seed = seed
	var shuffled := labels.duplicate()
	var last := Time.get_ticks_msec()
	for i in total:
		var current: PackedByteArray
		if exact:
			current = labellings[i]
		elif i == 0:
			current = labels  # the observed labelling counts once
		else:
			_shuffle(shuffled, rng)
			current = shuffled
		var diffs := _differences(values, current, totals, n_a, n)
		var biggest := 0.0
		for w in windows:
			var a := absf(diffs[w])
			if a >= abs_obs[w] - EPS:
				exceed[w] += 1
			biggest = maxf(biggest, a)
		for w in windows:
			if biggest >= abs_obs[w] - EPS:
				exceed_max[w] += 1
		count += 1
		if Time.get_ticks_msec() - last > FRAME_BUDGET_MSEC:
			progress.emit(count, total)
			await Engine.get_main_loop().process_frame
			last = Time.get_ticks_msec()
			if cancelled:
				return {"ok": false, "cancelled": true, "error": "Cancelled."}
	var p := PackedFloat64Array()
	var p_adj := PackedFloat64Array()
	for w in windows:
		p.append(float(exceed[w]) / count)
		p_adj.append(float(exceed_max[w]) / count)
	return {"ok": true, "error": "", "cancelled": false, "exact": exact, "labellings": count, "possible": possible,
		"mean_a": mean_a, "mean_b": mean_b, "difference": observed, "p": p, "p_adjusted": p_adj,
		"min_p": 1.0 / count}


static func _sums(values: Array, labels: PackedByteArray) -> PackedFloat64Array:
	var sums := PackedFloat64Array()
	sums.resize(values[0].size())
	for u in labels.size():
		if labels[u] == 1:
			var row: PackedFloat64Array = values[u]
			for w in sums.size():
				sums[w] += row[w]
	return sums


static func _differences(values: Array, labels: PackedByteArray, totals: PackedFloat64Array, n_a: int, n: int) -> PackedFloat64Array:
	var sums := _sums(values, labels)
	var out := PackedFloat64Array()
	out.resize(sums.size())
	for w in sums.size():
		out[w] = sums[w] / n_a - (totals[w] - sums[w]) / (n - n_a)
	return out


static func _shuffle(labels: PackedByteArray, rng: RandomNumberGenerator) -> void:
	for i in range(labels.size() - 1, 0, -1):
		var j := rng.randi_range(0, i)
		var t := labels[i]
		labels[i] = labels[j]
		labels[j] = t


## Every way to choose which n_a of n units are group A.
static func _all_labellings(n: int, n_a: int) -> Array:
	var out: Array = []
	var chosen := PackedInt32Array()
	for i in n_a:
		chosen.append(i)
	while true:
		var labels := PackedByteArray()
		labels.resize(n)
		for c in chosen:
			labels[c] = 1
		out.append(labels)
		var i := n_a - 1
		while i >= 0 and chosen[i] == n - n_a + i:
			i -= 1
		if i < 0:
			break
		chosen[i] += 1
		for j in range(i + 1, n_a):
			chosen[j] = chosen[j - 1] + 1
	return out


static func _fail(message: String) -> Dictionary:
	return {"ok": false, "error": message}
