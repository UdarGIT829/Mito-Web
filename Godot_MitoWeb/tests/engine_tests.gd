extends RefCounted
## Analysis engine tests. Cases mirror the Python unit tests
## (tests/test_table_phase_3.py, test_reference_features.py and the
## comparison characterization) on a small synthetic rCRS VCF.

## Four samples. REF bases match rCRS. 302 is a multi-ALT row; 8860 fails FILTER.
const VCF := """##fileformat=VCFv4.2
##contig=<ID=chrM,length=16569>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
chrM	73	.	A	G	.	PASS	.	GT:AF	1:0.82	1:0.14	0:.	1:0.9
chrM	150	.	C	T	.	PASS	.	GT:AF	1:0.15	0:.	0:.	0:.
chrM	302	.	A	AC,ACC	.	PASS	.	GT:AF	1:0.3,.	2:.,0.4	0:.	0:.
chrM	750	.	A	G	.	PASS	.	GT:AF	1:1.0	1:1.0	1:1.0	1:1.0
chrM	8860	.	A	G	.	lowq	.	GT:AF	1:0.5	1:.	1:0.6	0:.
"""

var t  # the test runner, for check()
var store: VariantStore
var ids := {}  # name -> sample id


func run(runner) -> void:
	t = runner
	store = VariantStore.new()
	var raw := VCF.to_utf8_buffer()
	var inspected := VcfImporter.inspect("synthetic.vcf", raw)
	var tags := [["EV", "Treated"], ["EV", "Untreated"], ["Cell", "Treated"], []]
	for i in 4:
		inspected["plan"]["samples"][i]["tags"] = PackedStringArray(tags[i])
	await VcfImporter.new().import_source({"id": 1, "name": "synthetic.vcf", "raw": raw, "plan": inspected["plan"]}, store)
	for sample in store.samples:
		ids[sample["name"]] = sample["id"]
	test_percentages()
	test_af_ranges()
	test_genes()
	test_sample_query()
	test_frame()
	test_filters()
	test_flag_mode_and_sort()
	test_comparison()
	test_set_operations()
	test_metadata_filters()
	test_af_float32_migration()


func check(condition: bool, label: String) -> void:
	t.check(condition, label)


func _ids(names: Array) -> PackedInt32Array:
	return PackedInt32Array(names.map(func(n): return ids[n]))


func _positions(frame: MutationFrame) -> Array:
	return Array(frame.alleles).map(func(a): return "%d%s>%s" % [store.allele_pos[a], store.allele_ref[a], store.allele_alt[a]])


func test_percentages() -> void:
	check(MutationFilters.percent_passes(6, 10, 60), "presence threshold is inclusive (6/10 >= 60%)")
	check(not MutationFilters.percent_passes(6, 10, 60.0000000001), "presence just above the boundary fails")
	check(MutationFilters.percent_passes(4, 6, 66.6666666666), "4/6 passes 66.6666666666%")
	check(not MutationFilters.percent_passes(4, 6, 66.6666666667), "4/6 fails 66.6666666667%")
	check(MutationFilters.percent_passes(0, 10, 0) and MutationFilters.percent_passes(10, 10, 100), "0% and 100% edges")


func test_af_ranges() -> void:
	var ranges := [{"minimum": 0.1, "maximum": 0.2, "include_maximum": false}, {"minimum": 0.8, "maximum": 1.0}]
	check(MutationFilters.af_in_any_range(0.15, ranges), "AF inside a range")
	check(not MutationFilters.af_in_any_range(0.2, ranges), "excluded maximum is respected")
	check(MutationFilters.af_in_any_range(0.8, ranges), "included minimum is respected")
	check(not MutationFilters.af_in_any_range(NAN, ranges), "missing AF never matches")
	check(not MutationFilters.af_in_any_range(1.5, ranges), "AF outside 0..1 never matches")


func test_genes() -> void:
	check(MitoGenes.GENES.size() == 37, "37 NCBI gene records")
	check(Array(MitoGenes.locus(9208, "T")["genes"]) == ["MT-CO3"], "position inside one gene")
	check(Array(MitoGenes.locus(9207, "T")["genes"]) == ["MT-ATP6", "MT-CO3"], "overlapping genes are all returned")
	check(Array(MitoGenes.locus(10760, "A")["genes"]) == ["MT-ND4L", "MT-ND4"], "ND4L/ND4 overlap")
	var control := MitoGenes.locus(73, "A")
	check(control["region"] == "control_region" and control["genes"].is_empty(), "control region claims no genes")
	check(Array(MitoGenes.locus(9206, "TT")["genes"]) == ["MT-ATP6", "MT-CO3"], "reference spans are used")
	check(MitoGenes.locus(16569, "AA")["region"] == "control_region", "circular origin wraps")


func test_sample_query() -> void:
	check(SampleQuery.select(store, {"tags_all": ["EV", "Treated"]}) == _ids(["S1"]), "tags_all requires every tag")
	check(SampleQuery.select(store, {"tags_any": ["Treated"]}) == _ids(["S1", "S3"]), "tags_any")
	check(SampleQuery.select(store, {"tags_all": ["<NONE>"]}) == _ids(["S4"]), "<NONE> selects untagged samples")
	check(SampleQuery.select(store, {"tags_none": ["EV"]}) == _ids(["S3", "S4"]), "tags_none")
	check(SampleQuery.tag_counts(store)[0] == ["Cell", 1], "tag counts are listed")


func test_frame() -> void:
	var frame := MutationFrame.from_samples(store, _ids(["S1", "S2", "S3", "S4"]))
	check(_positions(frame) == ["73A>G", "150C>T", "302A>AC", "302A>ACC", "750A>G", "8860A>G"], "frame rows sorted by identity (got %s)" % [_positions(frame)])
	check(frame.carrier_count(0) == 3, "73 has three carriers")
	var passing := MutationFrame.from_samples(store, _ids(["S1", "S2", "S3", "S4"]), {"filter_values": ["PASS"]})
	check(not "8860A>G" in _positions(passing), "FILTER pre-filter drops non-PASS calls")
	var nothing := MutationFrame.from_samples(store, _ids(["S1"]), {"filter_values": ["NO_SUCH_FILTER"]})
	check(nothing.row_count() == 0 and store.string_id("NO_SUCH_FILTER") == -1, "unknown FILTER keeps nothing and adds no strings")


func test_filters() -> void:
	var frame := MutationFrame.from_samples(store, _ids(["S1", "S2", "S3", "S4"]))
	var table := MutationTable.build(frame, {"filters": {"minimum_presence_percent": 75}})
	check(table["rows"].map(func(r): return r["position"]) == [73, 750, 8860], "presence >= 75%")
	var af := MutationTable.build(frame, {"filters": {"af_ranges": [{"minimum": 0.8, "maximum": 1.0}], "minimum_af_percent": 66.6666666666, "position": {"exact": 73}}})
	check(af["matched_count"] == 1 and af["rows"][0]["sample_summary"]["af_passing_count"] == 2, "AF carrier percentage counts carriers (2 of 3)")
	var af_fail := MutationTable.build(frame, {"filters": {"af_ranges": [{"minimum": 0.8, "maximum": 1.0}], "minimum_af_percent": 66.6666666667, "position": {"exact": 73}}})
	check(af_fail["matched_count"] == 0, "AF carrier percentage boundary is exact")
	var missing := MutationTable.build(frame, {"filters": {"af_ranges": [{"minimum": 0.0, "maximum": 1.0}], "minimum_af_percent": 100, "position": {"exact": 8860}}})
	check(missing["matched_count"] == 0, "a carrier with missing AF fails the AF condition")
	var row: Dictionary = table["rows"][0]
	check(row["sample_results"].map(func(s): return s["sample_label"]) == ["S1", "S2", "S4"], "carriers in selected-sample order")
	check(row["reference"] == "A" and row["region"] == "control_region", "row has rCRS base and region")
	check(MutationTable.build(frame, {"filters": {"genes": ["MT-RNR1"]}})["rows"].map(func(r): return r["position"]) == [750], "gene filter")
	check(MutationTable.build(frame, {"filters": {"alt_alleles": ["T"]}})["matched_count"] == 1, "ALT filter")
	check(MutationTable.build(frame, {"filters": {"reference_bases": ["c"]}})["rows"][0]["position"] == 150, "reference base filter is case-insensitive")
	check(MutationTable.build(frame, {"filters": {"position": {"minimum": 300, "maximum": 800}}})["matched_count"] == 3, "position range")
	for bad in [{"position": {"exact": 5, "minimum": 1}}, {"minimum_af_percent": 10}, {"af_ranges": [{"minimum": 0.9, "maximum": 0.1}]},
			{"minimum_presence_percent": 101}, {"slippage_metadata_filters": [{"field": "repeat_count", "value": "gt"}]}]:
		check(MutationFilters.validate(bad) != "", "invalid config rejected: %s" % bad)


func test_flag_mode_and_sort() -> void:
	var frame := MutationFrame.from_samples(store, _ids(["S1", "S2", "S3", "S4"]))
	var flagged := MutationFilters.apply(frame, {"minimum_presence_percent": 75}, "flag", "common")
	check(flagged.row_count() == frame.row_count(), "flag mode keeps every row")
	check(Array(flagged.flags["common"]) == [1, 0, 0, 0, 1, 1], "flag column marks passing rows")
	var rows: Array = MutationTable.build(flagged, {"sort": [{"field": "carrier_count", "direction": "desc"}]})["rows"]
	check(rows.map(func(r): return r["position"]) == [750, 73, 8860, 150, 302, 302], "sort by carriers desc, ties by identity")
	check(rows[0]["common"] == true and rows[3]["common"] == false, "flag values appear on table rows")
	check(MutationTable.build(frame, {"sort": [{"field": "nonsense"}]})["error"] != "", "unknown sort field rejected")


func test_comparison() -> void:
	var both := _ids(["S1", "S2"])
	var frame := MutationFrame.from_samples(store, _ids(["S1", "S2", "S3", "S4"]))
	var default := Comparison.compare(frame, both)
	check(_positions(default["frame"]) == ["73A>G", "750A>G", "8860A>G"], "default rules require presence in every compared sample")
	check(Array(default["statuses"]) == ["common", "common", "common"], "statuses")
	var only_s1 := Comparison.compare(frame, both, Comparison.DEFAULT_STATUSES, {ids["S1"]: ["present"], ids["S2"]: ["not_in"]})
	check(_positions(only_s1["frame"]) == ["150C>T", "302A>AC"], "present + not_in")
	var uniques := Comparison.compare(frame, both, ["unique"], {ids["S1"]: ["unique"], ids["S2"]: ["unique"]})
	check(_positions(uniques["frame"]) == ["150C>T", "302A>AC", "302A>ACC"], "unique in either sample")
	check(Comparison.compare(frame, _ids(["S1"]))["frame"].row_count() == 0, "fewer than two samples compares nothing")
	check(uniques["frame"].samples == both, "compared samples become the universe")


func test_set_operations() -> void:
	var treated := MutationFrame.from_samples(store, SampleQuery.select(store, {"tags_any": ["Treated"]}))
	var untreated := MutationFrame.from_samples(store, SampleQuery.select(store, {"tags_any": ["Untreated"]}))
	var in_both := treated.intersect(untreated)
	check(_positions(in_both) == ["73A>G", "750A>G", "8860A>G"], "intersect")
	check(_positions(treated.subtract(in_both)) == ["150C>T", "302A>AC"], "treated-exclusive = treated minus both")
	check(_positions(untreated.subtract(treated)) == ["302A>ACC"], "untreated-exclusive")
	var union := treated.union(untreated)
	check(union.row_count() == 6 and union.samples.size() == 3, "union merges rows and universes")


func test_metadata_filters() -> void:
	var frame := MutationFrame.from_samples(store, _ids(["S1", "S2", "S3", "S4"]))
	var poly := MutationTable.build(frame, {"filters": {"slippage_metadata_filters": [{"field": "polymorphism", "value": "1"}]}})
	check(poly["rows"].map(func(r): return r["alt"]) == ["AC", "ACC"], "polymorphism filter keeps multi-ALT rows")
	var repeat := MutationTable.build(frame, {"filters": {"slippage_metadata_filters": [{"field": "repeat_count", "value": "gte|2"}]}})
	check(repeat["rows"].map(func(r): return r["alt"]) == ["ACC"], "repeat count uses the VCF ALT (ACC has CC)")
	var context := MutationTable.build(frame, {"filters": {"slippage_metadata_filters": [{"field": "reference_repeat", "value": "both"}]}})
	# 150 CCTCAT|CTATTA, 302 CACCAA|CCCCCC, 750 AAAAGG|ACAAGC; 73 and 8860 only before.
	check(context["rows"].map(func(r): return r["position"]) == [150, 302, 302, 750], "reference repeat on both sides")
	var features := CallFeatures.new(store)
	var meta := features.metadata(frame.carriers[0][ids["S1"]])
	check(meta["REFERENCE_6_BEFORE"] == "GGGGGT" and meta["REFERENCE_6_AFTER"] == "TGCACG", "reference context around 73")
	check(not meta.has("REPEAT_1_BASE") and meta["REPEAT_1_BASE_COUNT"] == "0", "empty repeat motif omitted, count kept")
	check(CallFeatures.slippage([]) == {"risk": "unavailable", "footprint": "", "signals": {}}, "no signals means unavailable")


func test_af_float32_migration() -> void:
	var data := store.to_dict()
	var calls: Dictionary = data["calls"].duplicate()
	calls["af"] = PackedFloat32Array([0.8, 0.006738, NAN, 1.0])
	for name in calls:
		if name != "af":
			calls[name] = calls[name].slice(0, 4)
	data["calls"] = calls
	var restored := VariantStore.from_dict(data)
	check(restored["ok"], "store with float32 AF loads")
	var af: PackedFloat64Array = restored["store"].calls["af"]
	check(af[0] == 0.8 and af[1] == 0.006738 and is_nan(af[2]) and af[3] == 1.0, "float32 AF restored to exact decimals (got %s)" % af)
