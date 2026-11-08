extends RefCounted
## Window statistics. Synthetic calls (engine tests VCF):
##   73 S1 .82 S2 .14 S4 .9 | 150 S1 .15 | 302AC S1 .3 | 302ACC S2 .4
##   750 S1-S4 1.0 | 8860 S1 .5 S2 n/a S3 .6

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
	test_windows()
	test_regions()
	test_counts()
	test_statistics_node()


func check(condition: bool, label: String) -> void:
	t.check(condition, label)


func test_windows() -> void:
	var tiles := WindowStats.windows({"size": 500, "step": 500})
	check(tiles.size() == 34 and tiles[0]["segments"] == [[1, 500]], "500 bp tiles: 34 windows")
	check(tiles[-1]["segments"] == [[16501, 16569]] and tiles[-1]["length"] == 69, "last tile is clipped to the genome end")
	var sliding := WindowStats.windows({"size": 500, "step": 100})
	check(sliding[1]["segments"] == [[101, 600]], "step smaller than size overlaps")
	check(sliding[-1]["segments"][0][1] == 16569 and sliding.size() == 162, "sliding windows stop at the end (%d)" % sliding.size())
	var wrapped := WindowStats.windows({"size": 500, "step": 500, "wrap": true})
	check(wrapped.size() == 34 and wrapped[-1]["segments"] == [[16501, 16569], [1, 431]] and wrapped[-1]["length"] == 500, "wrap continues at position 1")
	var gene := WindowStats.windows({"size": 200, "step": 50, "range_start": 3307, "range_end": 4262})
	check(gene[0]["segments"] == [[3307, 3506]] and gene[-1]["segments"][0][1] == 4262, "range limits the scan (MT-ND1)")
	check(WindowStats.windows({"size": 200, "step": 50, "range_start": 3307, "range_end": 4262, "wrap": true})[-1]["segments"].size() == 1, "wrap only applies to the whole genome")
	for bad in [{"size": 0}, {"size": 500, "step": 0}, {"size": 500, "range_start": 900, "range_end": 100}, {"size": 500, "range_end": 20000}, {"size": 2.5}]:
		check(WindowStats.validate(bad) != "", "invalid window spec rejected: %s" % bad)
	check(WindowStats.validate({"mode": "regions"}) == "", "regions need no settings")


func test_regions() -> void:
	var regions := WindowStats.region_windows()
	var labels: Array = regions.map(func(w): return w["label"])
	check("MT-ND1" in labels and "control region" in labels and labels.filter(func(l): return l.begins_with("intergenic")).size() > 0, "genes, control region and intergenic stretches")
	var control: Dictionary = regions.filter(func(w): return w["label"] == "control region")[0]
	check(control["segments"] == [[16024, 16569], [1, 576]] and control["length"] == 1122, "control region wraps the origin (1,122 bp)")
	var covered := {}
	for w in regions:
		for s in w["segments"]:
			for p in range(s[0], s[1] + 1):
				covered[p] = true
	check(covered.size() == 16569, "regions cover the whole genome")


func test_counts() -> void:
	var tiles := WindowStats.windows({"size": 500, "step": 500})
	var c := WindowStats.count(frame, tiles)
	check(c["samples"] == 4, "group size")
	check(c["values"]["mutations"][0] == 4 and c["values"]["calls"][0] == 6, "1-500: 4 mutations, 6 carrier calls")
	check(is_equal_approx(c["values"]["af_load"][0], 2.71), "1-500: AF load 2.71")
	check(c["values"]["calls"][1] == 4 and is_equal_approx(c["values"]["af_load"][1], 4.0), "501-1000: 750 in all four samples")
	check(c["values"]["calls"][17] == 3 and is_equal_approx(c["values"]["af_load"][17], 1.1), "8501-9000: missing AF adds nothing to the load")
	var total := 0.0
	for v in c["values"]["calls"]:
		total += v
	check(total == 13, "every call counted once in non-overlapping tiles (3+1+1+1+4+3)")
	check(is_equal_approx(WindowStats.scaled(6, tiles[0], 4, true, true), 3.0), "per kb and per sample: 6 / 0.5 kb / 4")
	var regions := WindowStats.region_windows()
	var rc := WindowStats.count(frame, regions)
	var by_label := {}
	for i in regions.size():
		by_label[regions[i]["label"]] = rc["values"]["mutations"][i]
	check(by_label["control region"] == 4 and by_label["MT-RNR1"] == 1 and by_label["MT-ATP6"] == 1 and by_label["MT-ATP8"] == 0, "region counts")


func test_statistics_node() -> void:
	var g := GraphModel.new_graph()
	var provider := GraphModel.add_node(g, "data_provider")
	var treated := GraphModel.add_node(g, "tag_filter", Vector2.ZERO, {"tags": ["Treated"], "restrict_samples": true})
	var untreated := GraphModel.add_node(g, "tag_filter", Vector2.ZERO, {"tags": ["Untreated"], "restrict_samples": true})
	var stats := GraphModel.add_node(g, "statistics")
	for f in [treated, untreated]:
		GraphModel.connect_nodes(g, provider["id"], 0, f["id"], 0)
		GraphModel.connect_nodes(g, f["id"], 0, stats["id"], 0)
	GraphModel.connect_nodes(g, provider["id"], 0, stats["id"], 0)
	var result: Dictionary = GraphEvaluator.new().evaluate(g, store, 1)[stats["id"]]
	check(result["error"] == "" and result["summary"] == "34 windows · 3 groups", "statistics node evaluates (%s)" % result.get("summary"))
	var groups: Array = result["stats"]["groups"]
	check(groups.map(func(x): return x["label"]) == ["Tags: Treated", "Tags: Untreated", "Samples: all"], "groups are labelled by where they come from")
	check(groups[0]["samples"] == 2 and groups[1]["samples"] == 1, "group sizes follow restricted tag filters")
	stats["config"]["window"] = {"mode": "sliding", "size": 0}
	check(GraphEvaluator.new().evaluate(g, store, 1)[stats["id"]]["error"] != "", "invalid windows reported on the node")
	stats["config"]["window"] = {"mode": "regions"}
	check(GraphEvaluator.new().evaluate(g, store, 1)[stats["id"]]["stats"]["windows"].size() > 38, "regions mode")
