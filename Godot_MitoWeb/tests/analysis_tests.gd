extends RefCounted
## Feature matrix, permutation test, window import and scan windows.
## Synthetic calls (engine tests VCF):
##   73 S1 .82 S2 .14 S4 .9 | 150 S1 .15 | 302AC S1 .3 | 302ACC S2 .4
##   750 S1-S4 1.0 | 8860 S1 .5 S2 n/a S3 .6
## Tags: S1 EV+Treated, S2 EV+Untreated, S3 Cell+Treated, S4 none.

const EngineTests := preload("res://tests/engine_tests.gd")

var t
var store: VariantStore
var frame: MutationFrame


func run(runner) -> void:
	t = runner
	store = VariantStore.new()
	var raw: PackedByteArray = EngineTests.VCF.to_utf8_buffer()
	var inspected := VcfImporter.inspect("synthetic.vcf", raw)
	var tags := [["EV", "Treated"], ["EV", "Untreated"], ["Cell", "Treated"], []]
	for i in 4:
		inspected["plan"]["samples"][i]["tags"] = PackedStringArray(tags[i])
	await VcfImporter.new().import_source({"id": 1, "name": "synthetic.vcf", "raw": raw, "plan": inspected["plan"]}, store)
	frame = MutationFrame.from_samples(store, SampleQuery.select(store))
	test_matrix()
	await test_permutation_on_data()
	await test_permutation_statistics()
	test_prepare_errors()
	test_window_import()
	test_scan_windows()


func check(condition: bool, label: String) -> void:
	t.check(condition, label)


func test_matrix() -> void:
	var tiles := WindowStats.windows({"size": 500, "step": 500})
	var m := SampleMatrix.build(frame, tiles, "calls")
	check(m["rows"].size() == 4 and m["rows"][0].size() == 34, "4 samples x 34 windows")
	check(Array(m["rows"].map(func(r): return r[0])) == [3.0, 2.0, 0.0, 1.0], "calls per sample in 1-500")
	check(Array(m["rows"].map(func(r): return r[17])) == [1.0, 1.0, 1.0, 0.0], "calls per sample in 8501-9000")
	var load := SampleMatrix.build(frame, tiles, "af_load")
	check(is_equal_approx(load["rows"][0][0], 1.27) and load["rows"][1][17] == 0.0, "AF load per sample; missing AF adds 0")
	var per_kb := SampleMatrix.build(frame, tiles, "calls", true)
	check(per_kb["rows"][0][0] == 6.0, "per kb (3 calls in 0.5 kb)")
	var wrapped := SampleMatrix.build(frame, WindowStats.windows({"size": 500, "step": 500, "wrap": true}), "calls")
	check(wrapped["rows"][0][33] == 3.0, "wrapping window 16501-431 counts calls after position 1")
	var csv := SampleMatrix.to_csv(store, [{"label": "all", "matrix": m}], tiles, "calls", false)
	var lines := csv.strip_edges().split("\n")
	check(lines.size() == 5, "matrix CSV: header + one row per sample")
	check(lines[0].begins_with("sample,subject,group,tag:Cell,tag:EV,tag:Treated,tag:Untreated,1–500 calls"), "tags become 0/1 label columns")
	check(lines[1].begins_with("S1,S1,all,0,1,1,0,3"), "sample row with labels and values")


func test_permutation_on_data() -> void:
	var tiles := WindowStats.windows({"size": 500, "step": 500})
	var m := SampleMatrix.build(frame, tiles, "calls")
	var prepared := PermutationTest.prepare(store, m, "Treated", "Untreated", false)
	check(prepared["ok"] and prepared["n_a"] == 2 and prepared["n_b"] == 1, "Treated (S1, S3) vs Untreated (S2); untagged S4 left out")
	check(PermutationTest.labelling_count(2, 1) == 3.0, "3 possible labellings")
	var result: Dictionary = await PermutationTest.new().run(prepared, 999, 7)
	check(result["exact"] and result["labellings"] == 3 and is_equal_approx(result["min_p"], 1.0 / 3.0), "few units: exact enumeration, smallest possible p = 1/3")
	check(result["difference"][1] == 0.0 and result["p"][1] == 1.0, "a window everyone shares has difference 0 and p = 1")
	check(is_equal_approx(result["mean_a"][0], 1.5) and result["mean_b"][0] == 2.0, "group means (1-500: A 3,0 -> 1.5; B 2)")


func test_permutation_statistics() -> void:
	# 10 units per group; window 0 separates them perfectly, window 1 is noise.
	var values: Array = []
	var labels := PackedByteArray()
	var rng := RandomNumberGenerator.new()
	rng.seed = 3
	for i in 20:
		values.append(PackedFloat64Array([5.0 if i < 10 else 0.0, rng.randf()]))
		labels.append(1 if i < 10 else 0)
	var prepared := {"values": values, "labels": labels, "n_a": 10, "n_b": 10}
	var test := PermutationTest.new()
	var r1: Dictionary = await test.run(prepared, 999, 11)
	var r2: Dictionary = await test.run(prepared, 999, 11)
	check(not r1["exact"] and r1["labellings"] == 1000, "many units: 999 random labellings + the observed one")
	check(r1["p"][0] <= 0.01 and r1["p"] == r2["p"], "a real difference gets a small p; the same seed reproduces it")
	check(r1["p"][1] > 0.05, "noise does not")
	check(r1["p_adjusted"][0] >= r1["p"][0] and r1["p_adjusted"][1] >= r1["p"][1], "adjusted p is never smaller than raw p")
	# Exact: 3 vs 3 with perfect separation; the labelling and its mirror
	# image both reach the observed |difference|, so p = 2/20.
	var small := {"values": [PackedFloat64Array([5.0]), PackedFloat64Array([5.0]), PackedFloat64Array([5.0]),
		PackedFloat64Array([0.0]), PackedFloat64Array([0.0]), PackedFloat64Array([0.0])],
		"labels": PackedByteArray([1, 1, 1, 0, 0, 0]), "n_a": 3, "n_b": 3}
	var r3: Dictionary = await test.run(small, 999, 1)
	check(r3["exact"] and r3["labellings"] == 20 and is_equal_approx(r3["p"][0], 0.1), "exact two-sided p for perfect separation of 3 vs 3 is 2/20")


func test_prepare_errors() -> void:
	var m := SampleMatrix.build(frame, WindowStats.windows({"size": 500, "step": 500}), "calls")
	check(not PermutationTest.prepare(store, m, "Treated", "Treated", false)["ok"], "same tag twice is rejected")
	check(not PermutationTest.prepare(store, m, "EV", "Treated", false)["ok"], "a sample with both tags is rejected")
	check(not PermutationTest.prepare(store, m, "Treated", "NoSuchTag", false)["ok"], "an empty group is rejected")
	for sample in store.samples:
		sample["subject"] = "line1"
	check(not PermutationTest.prepare(store, m, "Treated", "Untreated", true)["ok"], "a subject in both groups is rejected when comparing subjects")
	for sample in store.samples:
		sample["subject"] = sample["name"]


func test_window_import() -> void:
	var csv := WindowImport.parse("label,start,end\nD-loop,16024,576\nND1,3307,4262\n")
	check(csv["ok"] and csv["intervals"] == [["D-loop", 16024, 576], ["ND1", 3307, 4262]], "CSV with header")
	var bed := WindowImport.parse("track name=x\nchrM\t3306\t4262\tND1\nchrM\t0\t500\n")
	check(bed["ok"] and bed["intervals"][0] == ["ND1", 3307, 4262] and bed["intervals"][1] == ["1–500", 1, 500], "BED is 0-based; missing names are generated")
	check(not WindowImport.parse("start,end\n1,20000\n")["ok"], "positions past the genome are rejected")
	check(not WindowImport.parse("a,b\nx,y\n")["ok"], "unrecognised tables are rejected")
	var windows := WindowStats.custom_windows(csv["intervals"])
	check(windows[0]["segments"] == [[16024, 16569], [1, 576]] and windows[0]["length"] == 1122, "start > end wraps the origin")
	check(WindowStats.validate({"mode": "custom", "intervals": csv["intervals"]}) == "" and WindowStats.validate({"mode": "custom"}) != "", "custom spec validation")


func test_scan_windows() -> void:
	var scan := WindowStats.scan_windows([500, 1000])
	check(scan.size() == 133 + 67, "scan: every size, stepped by a quarter, around the circle (%d)" % scan.size())
	check(scan.all(func(w): return w["length"] == 500 or w["length"] == 1000), "scan windows keep their full size by wrapping")
	var order: Array = range(scan.size())
	var picked := WindowStats.pick_non_overlapping(scan, order, 5)
	check(picked == [0, 4, 8, 12, 16], "greedy non-overlapping choice (%s)" % [picked])
