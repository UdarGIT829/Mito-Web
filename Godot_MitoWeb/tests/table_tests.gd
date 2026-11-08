extends RefCounted
## TableViewModel tests on the engine tests' synthetic VCF.
## Rows: 73A>G, 150C>T, 302A>AC, 302A>ACC, 750A>G, 8860A>G.
## AF: 73 S1 .82 S2 .14 S4 .9 | 150 S1 .15 | 302AC S1 .3 | 302ACC S2 .4
##     750 all 1.0 | 8860 S1 .5 S2 n/a S3 .6

const EngineTests := preload("res://tests/engine_tests.gd")

var t
var store: VariantStore
var frame: MutationFrame


func run(runner) -> void:
	t = runner
	store = VariantStore.new()
	var raw: PackedByteArray = EngineTests.VCF.to_utf8_buffer()
	var inspected := VcfImporter.inspect("synthetic.vcf", raw)
	await VcfImporter.new().import_source({"id": 1, "name": "synthetic.vcf", "raw": raw, "plan": inspected["plan"]}, store)
	frame = MutationFrame.from_samples(store, SampleQuery.select(store))
	test_columns()
	test_filters()
	test_sorting_and_paging()
	test_flags_and_reasons()
	test_export_and_details()


func check(condition: bool, label: String) -> void:
	t.check(condition, label)


func _ids(model: TableViewModel) -> Array:
	return model.view.map(func(r): return "%d%s" % [r["position"], r["alt"]])


func _sample_col(model: TableViewModel, name: String) -> String:
	for c in model.columns:
		if c["kind"] == "sample_af" and c["title"] == name:
			return c["id"]
	return ""


func test_columns() -> void:
	var model := TableViewModel.new(frame)
	var titles: Array = model.columns.map(func(c): return c["title"])
	check(titles.slice(0, 6) == ["Position", "Ref", "ALT", "Gene / region", "Carriers", "Presence %"], "base columns")
	check(titles.slice(-4) == ["S1", "S2", "S3", "S4"], "one AF column per sample when few samples")
	var row: Dictionary = model.view[0]
	check(model.cell_text(row, model.column(_sample_col(model, "S1"))) == "0.82", "AF cell")
	check(model.cell_text(row, model.column(_sample_col(model, "S3"))) == "", "non-carrier cell is empty")
	check(model.cell_text(model.view[5], model.column(_sample_col(model, "S2"))) == "n/a", "carrier without AF shows n/a")
	check(model.cell_text(row, model.column("genes")) == "control region", "region shown when no gene")
	model.set_per_sample_columns(false)
	check(model.columns[-1]["id"] == "samples" and model.cell_text(row, model.columns[-1]) == "S1 0.82; S2 0.14; S4 0.9", "compact samples column")


func test_filters() -> void:
	var model := TableViewModel.new(frame)
	model.set_filter("position", {"min": 100, "max": 800})
	check(_ids(model) == ["150T", "302AC", "302ACC", "750G"], "position range filter")
	model.clear_filters()
	model.set_filter("alt", {"contains": "c"})
	check(_ids(model) == ["302AC", "302ACC"], "text filter is case-insensitive")
	model.clear_filters()
	model.set_filter("genes", {"genes": ["(none)"]})
	check(_ids(model) == ["73G", "150T", "302AC", "302ACC"], "(none) selects control region rows")
	model.set_filter("genes", {"genes": ["MT-RNR1", "MT-ATP6"]})
	check(_ids(model) == ["750G", "8860G"], "gene filter")
	model.clear_filters()
	var s2 := _sample_col(model, "S2")
	model.set_filter(s2, {"mode": "carried"})
	check(_ids(model) == ["73G", "302ACC", "750G", "8860G"], "sample carries it")
	model.set_filter(s2, {"mode": "not_carried"})
	check(_ids(model) == ["150T", "302AC"], "sample does not carry it")
	model.clear_filters()
	model.set_filter(_sample_col(model, "S1"), {"mode": "range", "min": 0.5, "max": 1.0})
	check(_ids(model) == ["73G", "750G", "8860G"], "sample AF range (inclusive)")
	model.set_filter("carriers", {"min": 4})
	check(_ids(model) == ["750G"], "filters combine")
	check(model.describe_filters().contains("S1 AF 0.5–1") and model.describe_filters().contains("Carriers at least 4"), "filters are described")
	model.set_filter("carriers", {})
	check(model.filters.size() == 1, "an empty spec clears that filter")


func test_sorting_and_paging() -> void:
	var model := TableViewModel.new(frame)
	model.set_sort("carriers", true)
	check(_ids(model) == ["750G", "73G", "8860G", "150T", "302AC", "302ACC"], "sort by carriers desc, ties keep identity order")
	model.set_sort(_sample_col(model, "S1"), false)
	check(_ids(model) == ["150T", "302AC", "8860G", "73G", "750G", "302ACC"], "sort by a sample's AF; non-carriers last")
	model.set_sort(_sample_col(model, "S1"), true)
	check(_ids(model).back() == "302ACC", "missing values stay last when descending")
	model.page_size = 4
	check(model.page_count() == 2 and model.page_rows().size() == 4, "page 1")
	model.page = 1
	check(model.page_rows().size() == 2, "page 2")
	model.page = 9
	check(model.page_rows().size() == 2 and model.page == 1, "page is clamped")


func test_flags_and_reasons() -> void:
	var flagged := MutationFilters.apply(frame, {"minimum_presence_percent": 75, "alt_alleles": ["G", "T"]}, "flag", "keep")
	var model := TableViewModel.new(flagged)
	check(model.column("flag:keep")["kind"] == "flag" and not model.column("why:keep").is_empty(), "flag and reason columns")
	var reasons: Array = model.view.map(func(r): return model.cell_text(r, model.column("why:keep")))
	check(reasons == ["", "too few carriers", "ALT, too few carriers", "ALT, too few carriers", "", ""], "reasons explain failures (got %s)" % [reasons])
	model.set_filter("flag:keep", {"value": false})
	check(_ids(model) == ["150T", "302AC", "302ACC"], "filter on flag")
	model.set_filter("why:keep", {"contains": "alt"})
	check(_ids(model) == ["302AC", "302ACC"], "filter on reason text")


func test_export_and_details() -> void:
	var model := TableViewModel.new(frame)
	model.set_filter("position", {"max": 200})
	var csv := model.export_text(",")
	var lines := csv.strip_edges().split("\n")
	check(lines.size() == 3, "export has a header plus every filtered row (all pages)")
	check(lines[0].begins_with("Position,Ref,ALT,Gene / region,Carriers"), "CSV header")
	check(TableViewModel._quote("a,b \"c\"", ",") == "\"a,b \"\"c\"\"\"", "CSV quoting")
	check(not model.export_text("\t").contains(","), "TSV export has no commas here")
	var details := model.row_details(model.view[0])
	check(details.size() == 3 and details[0]["vcf"] == "73 A>G" and details[0]["placement"] == "Native", "row details list each call")
	check(details[0]["value"] == "1:0.82" and details[0]["metadata"].has("REFERENCE_6_BEFORE"), "details include raw value and metadata")
	check(TableViewModel.describe_repeats(details[0]["metadata"]) == "none", "no repeats described as none")
	check(TableViewModel.describe_repeats({"REPEAT_1_BASE": "C", "REPEAT_1_BASE_COUNT": "6", "REPEAT_2_BASES_COUNT": "0"}) == "C x6", "repeats described")
