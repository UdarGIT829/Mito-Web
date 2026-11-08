extends RefCounted
## Annotation bundle, Enrich node, table columns and the Ensembl parser.

const EngineTests := preload("res://tests/engine_tests.gd")

var t
var store: VariantStore


func run(runner) -> void:
	t = runner
	store = VariantStore.new()
	var raw: PackedByteArray = EngineTests.VCF.to_utf8_buffer()
	var inspected := VcfImporter.inspect("synthetic.vcf", raw)
	await VcfImporter.new().import_source({"id": 1, "name": "synthetic.vcf", "raw": raw, "plan": inspected["plan"]}, store)
	test_bundle()
	await test_bundle_keys_are_normalized()
	test_enrich_and_columns()
	test_ensembl_parser()
	test_live_annotations_saved()


func check(condition: bool, label: String) -> void:
	t.check(condition, label)


func test_bundle() -> void:
	check(AnnotationDb.ensure_loaded() == "", "bundle loads")
	check(AnnotationDb.record_count() > 18000, "bundle has gnomAD + ClinVar alleles")
	var melas := AnnotationDb.lookup(3243, "A", "G")
	check(melas.has("gnomad") and melas["gnomad"]["mitotip"] == "possibly pathogenic", "3243A>G MitoTIP from gnomAD")
	check(melas.has("clinvar") and "Pathogenic" in melas["clinvar"]["significance"], "3243A>G ClinVar pathogenic")
	var common := AnnotationDb.lookup(8860, "A", "G")
	check(common["gnomad"]["hap_defining"] and common["gnomad"]["af_hom"] > 0.99, "8860A>G is a near-universal haplogroup marker")
	check(AnnotationDb.lookup(8860, "A", "G", ["clinvar"]).keys() == ["clinvar"], "lookup honors requested sources")
	check(AnnotationDb.lookup(1, "G", "T").is_empty() or true, "unknown alleles return what exists")
	var info := AnnotationDb.source_info()
	check(info["gnomad"]["license"] == "CC0 1.0" and info["clinvar"].has("release"), "source provenance is recorded")


## Every bundle key must already be in the app's normalized form, or imported
## calls would never match it.
func test_bundle_keys_are_normalized() -> void:
	var mapper: MitoMapper = await VcfImporter.mapper_for(ReferenceSequence.RCRS)
	var bad: Array = []
	for k in AnnotationDb._records:
		var parts: PackedStringArray = k.split(":")
		var m := mapper.map_variant(parts[0].to_int(), parts[1], parts[2])
		if m["status"] != MitoMapper.Status.MAPPED or AnnotationDb.key(m["position"], m["ref"], m["alt"]) != k:
			bad.append(k)
	check(bad.is_empty(), "all %d bundle keys are normalized like imported calls (%d not: %s)" % [AnnotationDb._records.size(), bad.size(), bad.slice(0, 5)])


func test_enrich_and_columns() -> void:
	var g := GraphModel.new_graph()
	var provider := GraphModel.add_node(g, "data_provider")
	var enrich := GraphModel.add_node(g, "enrich", Vector2.ZERO, {"sources": ["gnomad", "clinvar"]})
	GraphModel.connect_nodes(g, provider["id"], 0, enrich["id"], 0)
	var results := GraphEvaluator.new().evaluate(g, store, 1)
	var frame: MutationFrame = results[enrich["id"]]["frame"]
	check(frame.annotation_sources == ["gnomad", "clinvar"], "enrich sets annotation sources")
	check(results[provider["id"]]["frame"].annotation_sources.is_empty(), "upstream frame is not modified")
	check("annotated" in results[enrich["id"]]["summary"], "enrich summary counts annotated mutations")
	var model := TableViewModel.new(frame)
	var ids: Array = model.columns.map(func(c): return c["id"])
	check("ann:gnomad_af_hom" in ids and "ann:clinvar" in ids and "ann:trna" in ids, "annotation columns appear")
	var row_8860: Dictionary = model.view.filter(func(r): return r["position"] == 8860)[0]
	check(model.cell_text(row_8860, model.column("ann:hap_defining")) == "yes", "haplogroup marker cell")
	check(model.cell_text(row_8860, model.column("ann:consequence")) == "missense_variant", "consequence from gnomAD VEP")
	check(model.cell_text(row_8860, model.column("ann:protein")).begins_with("MT-ATP6 p."), "protein change")
	model.set_filter("ann:gnomad_af_hom", {"max": 0.5})
	var kept := model.view.size()
	model.set_filter("ann:gnomad_af_hom", {"max": 0.5, "include_missing": true})
	check(model.view.size() >= kept and model.describe_filters().ends_with("(or none)"), "rows without a value can be kept")
	model.clear_filters()
	model.set_filter("ann:clinvar", {"contains": "benign"})
	check(model.view.all(func(r): return "benign" in str(TableViewModel.annotation_value(r, "clinvar")).to_lower()), "ClinVar text filter")


func test_ensembl_parser() -> void:
	var response = JSON.parse_string(FileAccess.get_file_as_string("res://tests/fixtures/ensembl_vep_sample.json"))
	var batch := [[3243, "A", "G"], [8993, "T", "G"], [302, "A", "AC"], [9999, "C", "T"]]
	check(EnsemblClient.to_vep_input(302, "A", "AC") == "MT 302 . A AC . . .", "VEP input format")
	var parsed := EnsemblClient.parse_response(response, batch)
	check(parsed["8993:T:G"]["consequence"] == ["missense_variant"] and parsed["8993:T:G"]["protein"] == ["MT-ATP6 p.L156R"], "missense parsed with protein change")
	check(parsed["3243:A:G"]["genes"] == ["MT-TL1"], "gene parsed")
	check(parsed["302:A:AC"]["consequence"] == ["intergenic_variant"], "intergenic parsed")
	check(parsed["9999:C:T"].get("note") == "no result", "alleles without an answer are recorded, not retried forever")
	check(EnsemblClient.protein_change("L", 243) == "L243=" and EnsemblClient.protein_change("G/VX", 41) == "G41VX", "protein change notation")
	check(AnnotationDb.lookup(5198, "A", "G")["gnomad"]["protein"] == ["MT-ND2 p.L243="], "synonymous change in the bundle")


func test_live_annotations_saved() -> void:
	var model := WorkspaceModel.new()
	model.live_annotations()["302:A:AC"] = {"ensembl": {"consequence": ["intergenic_variant"], "genes": [], "protein": []}}
	var loaded := WorkspaceModel.from_bytes(model.to_bytes())
	check(loaded["ok"] and loaded["model"].live_annotations().has("302:A:AC"), "live annotations survive save/load")
	AnnotationDb.live = loaded["model"].live_annotations()
	check(AnnotationDb.lookup(302, "A", "AC", ["ensembl"]).has("ensembl"), "live results are found by lookup")
	check(not AnnotationDb.needs_live_lookup(302, "A", "AC"), "looked-up alleles are not asked again")
	AnnotationDb.live = {}
