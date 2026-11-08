extends RefCounted
## Zip uploads and manifests.

const EngineTests := preload("res://tests/engine_tests.gd")

## A single-sample file without a contig line (build must come from the manifest).
const NO_CONTIG := "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tX\nchrM\t150\t.\tC\tT\t.\tPASS\t.\tGT:VF\t0/1:0.2\n"

var t


func run(runner) -> void:
	t = runner
	var zip_bytes := _make_zip({
		"bundle/one.vcf": NO_CONTIG,
		"bundle/cohort.vcf": EngineTests.VCF,
		"bundle/manifest.csv": "# comment line\nfile,sample,subject,build,tags\none.vcf,,line1,hg38,EV;Treated\ncohort.vcf,S1,donorA,,Treated\ncohort.vcf,S2,donorA,,Untreated\nmissing.vcf,,,,x\n",
		"bundle/notes.txt": "ignored",
		"__MACOSX/bundle/._one.vcf": "junk",
	})
	var expanded := UploadBundle.expand([{"name": "bundle.zip", "bytes": zip_bytes}])
	check(expanded["files"].map(func(f): return f["name"]) == ["one.vcf", "cohort.vcf"], "VCFs extracted from the zip (got %s)" % [expanded["files"].map(func(f): return f["name"])])
	check(expanded["manifests"].size() == 1 and expanded["skipped"] == ["bundle/notes.txt"], "manifest found; other files skipped; __MACOSX ignored")
	var model := WorkspaceModel.new()
	for f in expanded["files"]:
		model.add_source(f["name"], f["bytes"])
	var one: Dictionary = model.sources()[0]
	check(one["plan"]["build"] == "", "without the manifest the build is unknown")
	var result := UploadBundle.apply_manifests(model.sources(), expanded["manifests"])
	check(result["applied"] == 3 and result["unmatched"] == ["missing.vcf"], "manifest rows applied; unknown files reported")
	check(one["plan"]["build"] == "rcrs", "build from the manifest (hg38 alias)")
	var planned: Dictionary = one["plan"]["samples"][0]
	check(planned["subject"] == "line1" and Array(planned["tags"]) == ["EV", "Treated"], "subject and tags for a single-sample file")
	var cohort: Array = model.sources()[1]["plan"]["samples"]
	check(cohort[0]["subject"] == "donorA" and Array(cohort[0]["tags"]) == ["Treated"], "multi-sample rows match by sample column")
	check(Array(cohort[1]["tags"]) == ["Untreated"] and cohort[2]["subject"] == "S3" and cohort[2]["tags"].is_empty(), "other samples untouched")
	UploadBundle.apply_manifests(model.sources(), expanded["manifests"])
	check(Array(cohort[0]["tags"]) == ["Treated"], "applying twice does not duplicate tags")
	var bad := UploadBundle.expand([{"name": "x.zip", "bytes": "not a zip".to_utf8_buffer()}])
	check(bad["errors"].size() == 1 and bad["files"].is_empty(), "a broken zip is reported")
	var loose := UploadBundle.expand([{"name": "a.vcf", "bytes": NO_CONTIG.to_utf8_buffer()},
		{"name": "my_manifest.csv", "bytes": "file,build\na.vcf,hg19+58\n".to_utf8_buffer()}])
	check(loose["files"].size() == 1 and loose["manifests"].size() == 1, "a manifest can be uploaded next to plain VCFs")
	var m2 := WorkspaceModel.new()
	m2.add_source("a.vcf", NO_CONTIG.to_utf8_buffer())
	UploadBundle.apply_manifests(m2.sources(), loose["manifests"])
	check(m2.sources()[0]["plan"]["build"] == "hg19_plus58", "build alias hg19+58")
	var unknown := UploadBundle.apply_manifests(m2.sources(), UploadBundle.expand([{"name": "manifest.csv", "bytes": "file,build\na.vcf,mm10\n".to_utf8_buffer()}])["manifests"])
	check(unknown["errors"].size() == 1, "unknown builds are reported")
	# Values tables (not manifests) in a zip or uploaded alone.
	var with_table := UploadBundle.expand([{"name": "b.zip", "bytes": _make_zip({"b/one.vcf": NO_CONTIG, "b/clinical.csv": "Subject,Haplogroup,Score,Sex\nline1,A2o,8,F\n"})},
		{"name": "more.tsv", "bytes": "Subject\tGroup\nline1\tcase\n".to_utf8_buffer()}])
	check(with_table["tables"].size() == 2 and with_table["files"].size() == 1, "values tables are recognised in zips and on their own")
	var roles := TagImport.default_roles(with_table["tables"][0]["table"])
	check(roles == {"match": 0, "values": [2, 3], "tags": [], "haplogroup": 1}, "default roles: haplogroup column, numbers make the rest values (got %s)" % roles)
	check(TagImport.default_roles(with_table["tables"][1]["table"])["tags"] == [1], "text-only tables become tags")
	# Values and hand-set haplogroups reach samples imported later.
	var m3 := WorkspaceModel.new()
	m3.add_source("one.vcf", NO_CONTIG.to_utf8_buffer())
	var src: Dictionary = m3.sources()[0]
	src["plan"]["build"] = "rcrs"
	src["plan"]["samples"][0]["subject"] = "line1"
	var pending := TagImport.apply(src["plan"]["samples"], with_table["tables"][0]["table"], 0, [], false, -1, [2, 3], 1)
	for pair in pending["haplogroups"]:
		pair[0]["haplogroup_manual"] = pair[1]
	await VcfImporter.new().import_source(src, m3.store)
	var imported: Dictionary = m3.store.samples[0]
	check(imported["attributes"] == {"Score": 8.0, "Sex": "F"} and imported.get("haplogroup_manual") == "A2o", "pending samples carry values and the hand-set haplogroup into the import")


func check(condition: bool, label: String) -> void:
	t.check(condition, label)


func _make_zip(entries: Dictionary) -> PackedByteArray:
	var path := "user://_test_bundle.zip"
	var packer := ZIPPacker.new()
	packer.open(path)
	for name in entries:
		packer.start_file(name)
		packer.write_file(str(entries[name]).to_utf8_buffer())
		packer.close_file()
	packer.close()
	var bytes := FileAccess.get_file_as_bytes(path)
	DirAccess.remove_absolute(path)
	return bytes
