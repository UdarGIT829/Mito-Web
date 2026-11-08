extends RefCounted
## Duplicate detection: same content in another form, the same sample in
## another call set, and cleanup.

const EngineTests := preload("res://tests/engine_tests.gd")

var t


func run(runner) -> void:
	t = runner
	test_same_content()
	test_same_sample()
	await test_imported_original_and_removal()


func check(condition: bool, label: String) -> void:
	t.check(condition, label)


static func _vcf(column: String, records: Array) -> String:
	var lines := ["##fileformat=VCFv4.2", "##contig=<ID=chrM,length=16569>", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + column]
	for r in records:
		lines.append("chrM\t%d\t.\t%s\t%s\t.\tPASS\t.\tGT:AF\t0/1:%s" % r)
	return "\n".join(lines) + "\n"


const STRICT := [[73, "A", "G", "0.99"], [150, "C", "T", "0.98"], [750, "A", "G", "1.0"]]
const FILTERED := [[73, "A", "G", "0.99"], [150, "C", "T", "0.98"], [750, "A", "G", "1.0"], [2706, "A", "G", "0.02"], [8860, "A", "G", "0.01"]]


func test_same_content() -> void:
	var model := WorkspaceModel.new()
	var text := _vcf("AP4_EV_PCR", STRICT)
	check(model.add_source("a.vcf", text.to_utf8_buffer())[1], "first copy added")
	check(not model.add_source("a.vcf.gz", text.to_utf8_buffer().compress(FileAccess.COMPRESSION_GZIP))[1], "gzipped copy of the same content is not added")
	check(not model.add_source("a_crlf.vcf", text.replace("\n", "\r\n").to_utf8_buffer())[1], "CRLF copy is not added")
	check(model.sources().size() == 1, "one source")


func test_same_sample() -> void:
	var model := WorkspaceModel.new()
	var strict: Dictionary = model.add_source("AP4.strict.vcf", _vcf("AP4_EV_PCR", STRICT).to_utf8_buffer())[0]
	var filtered: Dictionary = model.add_source("AP4.filtered.vcf", _vcf("AP4_EV_PCR", FILTERED).to_utf8_buffer())[0]
	var other: Dictionary = model.add_source("C10.vcf", _vcf("C10_EV_PCR", FILTERED).to_utf8_buffer())[0]
	var report := model.duplicate_report()
	check(report.keys() == [filtered["id"]], "same sample column in another call set is flagged; other samples are not")
	check(report[filtered["id"]]["of"] == strict["id"] and "100% of calls shared" in report[filtered["id"]]["reason"], "reason names the original and the overlap (%s)" % report[filtered["id"]]["reason"])
	# Generic column names (Illumina sample numbers) need shared calls too.
	var generic := WorkspaceModel.new()
	generic.add_source("PL20-19-A_S16.vcf", _vcf("16", STRICT).to_utf8_buffer())
	generic.add_source("PL20-53-A_S16.vcf", _vcf("16", [[3243, "A", "G", "0.5"], [9999, "C", "T", "0.4"]]).to_utf8_buffer())
	check(generic.duplicate_report().is_empty(), "a generic column name alone (\"16\") is not a duplicate")
	generic.add_source("PL20-19-A_copy.vcf", _vcf("16", FILTERED).to_utf8_buffer())
	check(generic.duplicate_report().size() == 1, "generic name plus shared calls is a duplicate")
	var multi := WorkspaceModel.new()
	multi.add_source("cohort.vcf", EngineTests.VCF.to_utf8_buffer())
	multi.add_source("cohort_v2.vcf", EngineTests.VCF.replace("PASS\t.\tGT:AF\t1:0.82", "PASS\t.\tGT:AF\t1:0.80").to_utf8_buffer())
	check(multi.duplicate_report().size() == 1, "a multi-sample file sharing most sample columns is flagged")


func test_imported_original_and_removal() -> void:
	var model := WorkspaceModel.new()
	var first: Dictionary = model.add_source("AP4.filtered.vcf", _vcf("AP4_EV_PCR", FILTERED).to_utf8_buffer())[0]
	var second: Dictionary = model.add_source("AP4.strict.vcf", _vcf("AP4_EV_PCR", STRICT).to_utf8_buffer())[0]
	for source in [first, second]:
		var result: Dictionary = await VcfImporter.new().import_source(source, model.store)
		model.mark_imported(source["id"], result)
	check(model.store.samples.size() == 2, "both copies imported (the accident)")
	check(model.duplicate_report().keys() == [second["id"]], "the later copy is the duplicate")
	check(model.remove_duplicates() == 1 and model.sources().size() == 1 and model.sources()[0]["id"] == first["id"], "remove duplicates keeps the original")
	check(model.store.samples.size() == 1 and model.store.call_count() == 5, "its imported samples and calls are removed too")
	var mixed := WorkspaceModel.new()
	var pending: Dictionary = mixed.add_source("x.strict.vcf", _vcf("S", STRICT).to_utf8_buffer())[0]
	var imported: Dictionary = mixed.add_source("x.filtered.vcf", _vcf("S", FILTERED).to_utf8_buffer())[0]
	var r: Dictionary = await VcfImporter.new().import_source(imported, mixed.store)
	mixed.mark_imported(imported["id"], r)
	check(mixed.duplicate_report().keys() == [pending["id"]], "an imported copy counts as the original over a pending one")
