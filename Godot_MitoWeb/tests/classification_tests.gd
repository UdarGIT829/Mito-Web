extends RefCounted
## Germline / somatic AF classes, annotation rules and tag import.
## Synthetic calls (engine tests VCF):
##   73   S1 .82 G | S2 .14 S | S4 .9 G
##   150  S1 .15 S          302AC S1 .3 S       302ACC S2 .4 S
##   750  S1-S4 1.0 G       8860 S1 .5 I | S2 n/a ? | S3 .6 I

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
	test_af_classes()
	test_germline_somatic_node()
	test_clinvar_categories()
	test_annotation_rules()
	test_tag_import()


func check(condition: bool, label: String) -> void:
	t.check(condition, label)


func _labels(f: MutationFrame) -> Array:
	return Array(f.alleles).map(func(a): return "%d%s" % [store.allele_pos[a], store.allele_alt[a]])


func test_af_classes() -> void:
	var cases := [[1.0, "germline"], [0.75, "germline"], [0.7499, "indeterminate"], [0.5, "indeterminate"],
		[0.4999, "somatic"], [0.0, "somatic"], [NAN, "unknown"], [1.2, "unknown"]]
	for c in cases:
		check(AfClasses.classify(c[0]) == c[1], "AF %s is %s" % c)
	check(AfClasses.classify(0.6, {"germline_min": 0.9, "somatic_below": 0.6}) == "indeterminate", "custom cut-offs")
	check(AfClasses.validate({"germline_min": 0.4, "somatic_below": 0.5}) != "", "cut-offs out of order are rejected")


func test_germline_somatic_node() -> void:
	var somatic := AfClasses.apply(frame, AfClasses.DEFAULT_THRESHOLDS, ["somatic"])
	check(_labels(somatic) == ["73G", "150T", "302AC", "302ACC"], "somatic calls only (got %s)" % [_labels(somatic)])
	check(somatic.carriers[0].keys() == [store.samples[1]["id"]], "a mutation keeps only its somatic carriers (73: S2)")
	check(somatic.samples.size() == 4, "sample universe is unchanged")
	var germline := AfClasses.apply(frame, AfClasses.DEFAULT_THRESHOLDS, ["germline"])
	check(_labels(germline) == ["73G", "750G"] and germline.carriers[0].size() == 2, "germline calls only")
	var middle := AfClasses.apply(frame, AfClasses.DEFAULT_THRESHOLDS, ["indeterminate", "unknown"])
	check(_labels(middle) == ["8860G"] and middle.carriers[0].size() == 3, "indeterminate and unknown calls")
	var labelled := AfClasses.apply(frame, AfClasses.DEFAULT_THRESHOLDS, [], "label")
	check(labelled.row_count() == frame.row_count(), "label mode keeps every call")
	var model := TableViewModel.new(labelled)
	var row_8860: Dictionary = model.view[5]
	check(model.cell_text(row_8860, model.column("af_classes")) == "2 indeterminate, 1 unknown", "call classes column")
	var s1_col: Dictionary = model.columns.filter(func(c): return c.get("title") == "S1")[0]
	check(model.cell_text(model.view[0], s1_col) == "0.82 G", "per-sample cell shows the class")
	var g := GraphModel.new_graph()
	var provider := GraphModel.add_node(g, "data_provider")
	var node := GraphModel.add_node(g, "germline_somatic")
	GraphModel.connect_nodes(g, provider["id"], 0, node["id"], 0)
	var results := GraphEvaluator.new().evaluate(g, store, 1)
	check(results[node["id"]]["summary"].contains("calls: 4 somatic"), "node defaults to somatic calls and reports counts (%s)" % results[node["id"]]["summary"])
	node["config"]["keep"] = []
	check(GraphEvaluator.new().evaluate(g, store, 1)[node["id"]]["error"] != "", "keeping no classes is an error")


func test_clinvar_categories() -> void:
	var cases := [
		[{"clinvar": {"significance": "Pathogenic"}}, "pathogenic"],
		[{"clinvar": {"significance": "Likely pathogenic"}}, "pathogenic"],
		[{"clinvar": {"significance": "Benign/Likely benign"}}, "benign"],
		[{"clinvar": {"significance": "Conflicting classifications of pathogenicity"}}, "conflicting"],
		[{"clinvar": {"significance": "Pathogenic / Benign"}}, "conflicting"],
		[{"clinvar": {"significance": "Uncertain significance"}}, "uncertain"],
		[{"clinvar": {"significance": "Affects"}}, "other"],
		[{}, "absent"],
	]
	for c in cases:
		check(AnnotationFilters.clinvar_category(c[0]) == c[1], "ClinVar '%s' -> %s" % [c[0].get("clinvar", {}).get("significance", "(none)"), c[1]])
	check(AnnotationFilters.clinvar_category(AnnotationDb.lookup(3243, "A", "G")) == "pathogenic", "3243A>G is pathogenic")
	check(AnnotationFilters.clinvar_category(AnnotationDb.lookup(8860, "A", "G")) == "benign", "8860A>G is benign")
	check(AnnotationFilters.clinvar_stars({"clinvar": {"review": "criteria provided, multiple submitters, no conflicts"}}) == 2, "two stars")
	check(AnnotationFilters.clinvar_stars({"clinvar": {"review": "reviewed by expert panel"}}) == 3, "three stars")
	check(AnnotationFilters.clinvar_stars({"clinvar": {"review": "no assertion criteria provided"}}) == 0, "zero stars")


func test_annotation_rules() -> void:
	var category := func(f: MutationFrame) -> Array:
		return Array(f.alleles).map(func(a): return AnnotationFilters.clinvar_category(AnnotationDb.lookup(store.allele_pos[a], store.allele_ref[a], store.allele_alt[a])))
	var no_benign := AnnotationFilters.apply(frame, {"clinvar_categories": ["benign"], "action": "remove"})
	check(not "benign" in category.call(no_benign) and no_benign.row_count() < frame.row_count(), "remove known benign")
	check("absent" in category.call(no_benign), "mutations without ClinVar are kept when removing benign")
	var markers := AnnotationFilters.apply(frame, {"haplogroup_marker": "yes"})
	check(markers.row_count() > 0 and Array(markers.alleles).all(func(a): return AnnotationDb.lookup(store.allele_pos[a], store.allele_ref[a], store.allele_alt[a]).get("gnomad", {}).get("hap_defining", false)), "keep haplogroup markers")
	var rare := AnnotationFilters.apply(frame, {"gnomad_af_hom_max": 0.001})
	var rare_or_absent := AnnotationFilters.apply(frame, {"gnomad_af_hom_max": 0.001, "gnomad_keep_absent": true})
	check(rare_or_absent.row_count() >= rare.row_count(), "absent-from-gnomAD option only adds rows")
	var flagged := AnnotationFilters.apply(frame, {"clinvar_categories": ["pathogenic"]}, "flag", "pathogenic")
	check(flagged.row_count() == frame.row_count() and flagged.flag_reasons["pathogenic"][0] == "ClinVar category", "flag mode explains failures")
	check(AnnotationFilters.validate({"action": "keep"}) != "", "a rule needs at least one condition")
	check(AnnotationFilters.validate({"clinvar_categories": ["nonsense"]}) != "", "unknown category rejected")
	check(AnnotationFilters.describe({"clinvar_categories": ["benign"], "action": "remove"}).begins_with("Remove mutations matching"), "description")


func test_tag_import() -> void:
	var csv := "﻿sample,cell line,treatment,\"note, quoted\"\r\nS1,C10,Treated,\"a, b\"\r\ns2,AP4,Untreated,\r\nNOPE,X,Y,\r\n"
	var parsed := TagImport.parse(csv)
	check(parsed["ok"] and parsed["header"] == ["sample", "cell line", "treatment", "note, quoted"], "CSV header with BOM and quotes")
	check(parsed["rows"][0][3] == "a, b", "quoted cell with comma")
	var tsv := TagImport.parse("id\thaplogroup\nS3\tH1\n")
	check(tsv["header"] == ["id", "haplogroup"] and tsv["rows"][0] == ["S3", "H1"], "TSV detected")
	var samples: Array = [
		{"name": "S1", "subject": "S1", "tags": PackedStringArray(["EV"])},
		{"name": "S2", "subject": "S2", "tags": PackedStringArray()},
	]
	var result := TagImport.apply(samples, parsed, 0, [1, 2], false)
	check(result["matched"] == 2 and result["unmatched"] == ["NOPE"], "matches by name, case-insensitively; reports unmatched rows")
	check(Array(samples[0]["tags"]) == ["EV", "C10", "Treated"], "tags added after existing ones")
	var again := TagImport.apply(samples, parsed, 0, [1, 2], false)
	check(again["tags_added"] == 0, "importing twice adds nothing")
	TagImport.apply(samples, tsv, 0, [1], true)
	var by_subject: Array = [{"name": "PL-20-19_NCRef", "subject": "19", "tags": PackedStringArray()},
		{"name": "PL20-19-A", "subject": "19", "tags": PackedStringArray()}]
	var subject_table := TagImport.parse("subject,group\n19,Lupus\n")
	check(TagImport.apply(by_subject, subject_table, 0, [1], true)["matched"] == 2 and Array(by_subject[1]["tags"]) == ["group=Lupus"], "a subject id tags all its samples, with prefix")
	check(not TagImport.parse("only a header\n")["ok"], "a table without data rows is rejected")
	var lines: Array = [{"name": "C10_EV", "subject": "C10_EV", "tags": PackedStringArray()}, {"name": "C10_cell", "subject": "C10_cell", "tags": PackedStringArray()}]
	TagImport.apply(lines, TagImport.parse("sample,line\nC10_EV,C10\nC10_cell,C10\n"), 0, [], false, 1)
	check(lines.all(func(x): return x["subject"] == "C10"), "subject set from a column")
