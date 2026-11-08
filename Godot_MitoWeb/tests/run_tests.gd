extends SceneTree
## Headless tests. Run from Godot_MitoWeb/:
##   godot --headless --script res://tests/run_tests.gd
## Exits with status 1 if any check fails.

var _failures := 0
var _checks := 0

const SINGLE_SAMPLE_VCF := """##fileformat=VCFv4.2
##contig=<ID=chrM,length=16569>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chrM	73	.	A	G	.	PASS	DP=10	GT:AF	0/1:0.99
chrM	302	.	A	AC,C	.	strand_bias	DP=10	GT:AF	0/1/2:0.2,0.05
chrM	3107	.	N	C	.	PASS	DP=10	GT:AF	0/1:0.5
chr1	100	.	A	G	.	PASS	DP=10	GT:AF	0/1:0.5
"""

## Mimics the Illumina NCRef exports: no ID column, quoted sample values,
## CRLF line endings, VF instead of AF.
const QUOTED_VCF := "##fileformat=VCFv4.1\r\n#CHROM\tPOS\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tX.bam\r\nchrM\t150\tC\tT\t100\tPASS\t5\tGT:VF\t\"0/1:0.078\"\r\n"

const MULTI_SAMPLE_VCF := """##fileformat=VCFv4.1
##contig=<ID=MT,length=16569>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B	C
MT	150	.	C	T	100	PASS	.	GT	1	0	.
MT	263	.	A	G	100	PASS	.	GT	1	1	1
MT	310	.	T	C,TC	100	PASS	.	GT	2	0	1
"""


func _init() -> void:
	_run.call_deferred()


func _run() -> void:
	test_workspace_round_trip()
	test_workspace_rejects_bad_input()
	test_legacy_json_migration()
	test_gzip_and_bgzf()
	test_vcf_reader()
	test_sequence_features()
	await test_mapper_matches_python()
	await test_import_single_sample()
	await test_import_quoted_illumina()
	await test_import_multi_sample()
	await test_remove_and_undo()
	await preload("res://tests/engine_tests.gd").new().run(self)
	await preload("res://tests/graph_tests.gd").new().run(self)
	await preload("res://tests/table_tests.gd").new().run(self)
	await preload("res://tests/annotation_tests.gd").new().run(self)
	await preload("res://tests/classification_tests.gd").new().run(self)
	print("%d checks, %d failures" % [_checks, _failures])
	quit(1 if _failures else 0)


func check(condition: bool, label: String) -> void:
	_checks += 1
	if not condition:
		_failures += 1
		printerr("FAIL: " + label)


func test_workspace_round_trip() -> void:
	var model := WorkspaceModel.new()
	model.set_name("Treated vs untreated")
	var a := model.add_source("a.vcf", SINGLE_SAMPLE_VCF.to_utf8_buffer())
	var b := model.add_source("renamed.vcf", SINGLE_SAMPLE_VCF.to_utf8_buffer())
	check(a[1] and not b[1], "identical content is not added twice")
	check(a[0]["status"] == "pending", "a valid VCF is ready to import")
	check(a[0]["plan"]["build"] == "rcrs", "build detected from contig length 16569")
	check(a[0]["plan"]["samples"][0]["name"] == "a", "single-sample name comes from the file name")
	var bad := model.add_source("notes.txt", "hello".to_utf8_buffer())
	check(bad[0]["status"] == "error", "a non-VCF file is marked unreadable")
	model.store.intern_allele(73, "A", "G")
	var result := WorkspaceModel.from_bytes(model.to_bytes())
	check(result["ok"], "binary workspace parses: %s" % result.get("error", ""))
	if result["ok"]:
		var copy: WorkspaceModel = result["model"]
		check(copy.get_name() == "Treated vs untreated", "name survives round trip")
		check(copy.sources().size() == 2, "sources survive round trip")
		check(copy.sources()[0]["raw"] == SINGLE_SAMPLE_VCF.to_utf8_buffer(), "raw bytes survive round trip")
		check(copy.store.allele_pos.size() == 1, "store survives round trip")
		check(copy.add_source("c.vcf", "x".to_utf8_buffer())[0]["id"] == 3, "ids continue after load")


func test_workspace_rejects_bad_input() -> void:
	check(not WorkspaceModel.from_bytes("not a workspace".to_utf8_buffer())["ok"], "rejects arbitrary bytes")
	check(not WorkspaceModel.from_bytes("{\"format\": \"other\"}".to_utf8_buffer())["ok"], "rejects other JSON")
	var future := WorkspaceModel.new().to_bytes()
	future.encode_u32(WorkspaceModel.MAGIC.length(), 99)
	check(not WorkspaceModel.from_bytes(future)["ok"], "rejects newer format versions")
	var truncated := WorkspaceModel.new().to_bytes()
	truncated.resize(truncated.size() - 5)
	check(not WorkspaceModel.from_bytes(truncated)["ok"], "rejects damaged files")


func test_legacy_json_migration() -> void:
	var legacy := {
		"format": "mito-web-workspace", "format_version": 1, "name": "Old", "created": "2026-11-07T00:00:00Z",
		"modified": "", "next_id": 2,
		"sources": [{"id": 1, "name": "x.vcf", "size": 1, "sha256": "", "added": "", "text": MULTI_SAMPLE_VCF}],
	}
	var result := WorkspaceModel.from_bytes(JSON.stringify(legacy).to_utf8_buffer())
	check(result["ok"] and result["migrated"], "format 1 workspaces are migrated")
	if result["ok"]:
		var model: WorkspaceModel = result["model"]
		check(model.get_name() == "Old" and model.sources()[0]["status"] == "pending", "migrated sources become pending")


func test_gzip_and_bgzf() -> void:
	var text := MULTI_SAMPLE_VCF.to_utf8_buffer()
	var gz := text.compress(FileAccess.COMPRESSION_GZIP)
	var plain := GzipReader.decompress(gz)
	check(plain["ok"] and plain["bytes"] == text, "plain gzip decompresses")
	# Two hand-built BGZF blocks plus the empty EOF block.
	var half := text.size() / 2
	var bgzf := _bgzf_block(text.slice(0, half))
	bgzf.append_array(_bgzf_block(text.slice(half)))
	bgzf.append_array(_bgzf_block(PackedByteArray()))
	var blocked := GzipReader.decompress(bgzf)
	check(blocked["ok"] and blocked["bytes"] == text, "multi-block BGZF decompresses fully")
	check(GzipReader.decompress(text)["bytes"] == text, "uncompressed input passes through")


func test_vcf_reader() -> void:
	var reader := VcfReader.parse(SINGLE_SAMPLE_VCF)
	check(reader.ok, "single-sample VCF parses")
	check(reader.contig_length == 16569 and reader.record_count == 4, "header facts")
	check(reader.other_contig_records == 1, "non-mito records are counted")
	var quoted := VcfReader.parse(QUOTED_VCF)
	check(quoted.ok and not quoted.column_index.has("ID"), "VCF without ID column parses")
	check(Array(quoted.sample_columns) == ["X.bam"], "sample columns after FORMAT (CRLF stripped)")
	check(VcfReader.unquote("\"0/1:0.1\"") == "0/1:0.1", "quoted values are unquoted")
	check(not VcfReader.parse("no header here").ok, "rejects text without #CHROM")


func test_sequence_features() -> void:
	check(SequenceFeatures.longest_repeat("ACCCCCCT", 1) == ["C", 6], "1-base repeat")
	check(SequenceFeatures.longest_repeat("ACACACAG", 2) == ["AC", 3], "2-base repeat")
	check(SequenceFeatures.longest_repeat("AT", 1) == ["", 0], "no repeat")
	check(not SequenceFeatures.is_minimal_motif("AA"), "AA is not a minimal motif")
	var rcrs := ReferenceSequence.get_build(ReferenceSequence.RCRS).sequence
	check(SequenceFeatures.reference_context(rcrs, 1, 3) == [rcrs.substr(16566, 3), rcrs.substr(1, 3)], "context wraps around the origin")
	check(VcfImporter.default_sample_name("PL20-19-A_S16.vcf.gz") == "PL20-19-A", "Illumina _S suffix is stripped")


func test_mapper_matches_python() -> void:
	var cases = JSON.parse_string(FileAccess.get_file_as_string("res://tests/fixtures/mapper_cases.json"))
	var mismatches := 0
	for case in cases:
		var mapper: MitoMapper = await VcfImporter.mapper_for(case["build"])
		var got := mapper.map_variant(int(case["position"]), case["ref"], case["alt"])
		var status: String = MitoMapper.STATUS_NAMES[got["status"]]
		var converted = [got["position"], got["ref"], got["alt"]] if status == "MAPPED" else null
		var expected = case["converted"]
		if expected != null:
			expected = [int(expected[0]), expected[1], expected[2]]
		if status != case["status"] or converted != expected:
			mismatches += 1
			if mismatches <= 5:
				printerr("  mapper mismatch %s: got %s %s" % [case, status, converted])
	check(mismatches == 0, "mapper matches Python on %d fixture cases (%d mismatches)" % [cases.size(), mismatches])


func test_import_single_sample() -> void:
	var store := VariantStore.new()
	var result := await _import("s1.vcf", SINGLE_SAMPLE_VCF, store)
	check(result["ok"], "single-sample import succeeds")
	check(result["calls"] == 4, "one call per ALT on mito rows (got %d)" % result["calls"])
	check(result["skipped_rows"] == 1, "non-mito row skipped")
	var c := store.calls
	check(c["af"][0] == 0.99 or is_equal_approx(c["af"][0], 0.99), "AF read from FORMAT AF")
	check(c["multi_alt"][1] == 1 and c["multi_alt"][0] == 0, "multi-ALT rows are flagged")
	check(store.strings[c["filter"][1]] == "strand_bias", "FILTER kept")
	check(c["status"][3] == VariantStore.CallStatus.UNSUPPORTED and c["allele"][3] == -1, "N REF is kept but not placed")
	check(store.sample_summaries()[store.samples[0]["id"]]["placed"] == 3, "placed count excludes unplaced calls")


func test_import_quoted_illumina() -> void:
	var store := VariantStore.new()
	var result := await _import("PL-20-01_NCRef.vcf", QUOTED_VCF, store, ReferenceSequence.RCRS)
	check(result["ok"] and result["calls"] == 1, "quoted Illumina VCF imports")
	check(is_equal_approx(store.calls["af"][0], 0.078), "AF read from VF")
	check(store.calls["values"][0] == "0/1:0.078", "quotes removed from the stored value")


func test_import_multi_sample() -> void:
	var store := VariantStore.new()
	var result := await _import("cohort.vcf.gz", MULTI_SAMPLE_VCF, store)
	check(result["ok"] and result["samples"] == 3, "multi-sample import creates one sample per column")
	var per_sample := {}
	for sample in store.samples:
		per_sample[sample["name"]] = store.sample_summaries()[sample["id"]]["calls"]
	check(per_sample == {"A": 3, "B": 1, "C": 2}, "calls follow genotypes (got %s)" % per_sample)


func test_remove_and_undo() -> void:
	var model := WorkspaceModel.new()
	var first: Dictionary = model.add_source("a.vcf", SINGLE_SAMPLE_VCF.to_utf8_buffer())[0]
	var second: Dictionary = model.add_source("cohort.vcf", MULTI_SAMPLE_VCF.to_utf8_buffer())[0]
	for source in [first, second]:
		var result: Dictionary = await VcfImporter.new().import_source(source, model.store)
		model.mark_imported(source["id"], result)
	check(model.store.samples.size() == 4 and model.store.call_count() == 10, "two files imported")
	model.undo_import(first["id"])
	check(first["status"] == "pending" and model.store.samples.size() == 3, "undo import removes its samples")
	check(model.store.call_count() == 6, "undo import removes its calls")
	model.remove_source(second["id"])
	check(model.store.samples.is_empty() and model.store.call_count() == 0, "removing a file removes its samples")


func _import(file_name: String, text: String, store: VariantStore, build: String = "") -> Dictionary:
	var raw := text.to_utf8_buffer()
	var inspected := VcfImporter.inspect(file_name, raw)
	if build:
		inspected["plan"]["build"] = build
	return await VcfImporter.new().import_source({"id": 1, "name": file_name, "raw": raw, "plan": inspected["plan"]}, store)


static func _bgzf_block(data: PackedByteArray) -> PackedByteArray:
	var deflated := data.compress(FileAccess.COMPRESSION_DEFLATE)
	# zlib-wrapped output: drop the 2-byte header and 4-byte Adler-32 trailer.
	var raw_deflate := deflated.slice(2, deflated.size() - 4)
	var block := PackedByteArray([0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, 66, 67, 2, 0, 0, 0])
	block.append_array(raw_deflate)
	var trailer := PackedByteArray()
	trailer.resize(8)
	trailer.encode_u32(0, _crc32(data))
	trailer.encode_u32(4, data.size())
	block.append_array(trailer)
	block.encode_u16(16, block.size() - 1)
	return block


static func _crc32(data: PackedByteArray) -> int:
	var crc := 0xFFFFFFFF
	for byte in data:
		crc ^= byte
		for _i in 8:
			crc = (crc >> 1) ^ (0xEDB88320 if crc & 1 else 0)
	return crc ^ 0xFFFFFFFF
