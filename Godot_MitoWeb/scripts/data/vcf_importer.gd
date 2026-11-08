class_name VcfImporter
extends RefCounted
## Inspects uploaded VCFs (header facts, detected build, default sample plan)
## and imports them into a VariantStore, projecting every ALT onto rCRS.
##
## Import is chunked: it yields to the main loop regularly so the browser tab
## stays responsive, and samples and calls are only added once the whole file
## has been read, so a cancelled or failed import adds no data (at most a few
## unused lookup-table entries).

signal progress(fraction: float, message: String)

const AF_KEYS := ["AF", "VF", "HF"]  # Mutect2, Illumina/Pisces, mtDNA-Server
const FRAME_BUDGET_MSEC := 25

static var _mappers := {}

var cancelled := false


## Reads a file's header and proposes an import plan.
## Returns {"ok", "error", "info", "plan"}.
static func inspect(file_name: String, raw: PackedByteArray) -> Dictionary:
	var unzipped := GzipReader.decompress(raw)
	if not unzipped["ok"]:
		return {"ok": false, "error": unzipped["error"]}
	var reader := VcfReader.parse(unzipped["bytes"].get_string_from_utf8())
	if not reader.ok:
		return {"ok": false, "error": reader.error}
	var detected := ReferenceSequence.detect_from_contig_length(reader.contig_length)
	var agreement := {}
	for build in ReferenceSequence.build_keys():
		agreement[build] = ref_agreement(reader, build)
	var info := {
		"content_sha256": content_fingerprint(unzipped["bytes"]),
		"record_keys": record_keys(reader),
		"contig_id": reader.contig_id,
		"contig_length": reader.contig_length,
		"reference_line": reader.reference_line,
		"source_line": reader.source_line,
		"record_count": reader.record_count,
		"blank_contig_records": reader.blank_contig_records,
		"other_contig_records": reader.other_contig_records,
		"sample_columns": reader.sample_columns,
		"detected_build": detected,
		"ref_agreement": agreement,
	}
	var samples: Array = []
	var single := reader.sample_columns.size() <= 1
	var stem := default_sample_name(file_name)
	if reader.sample_columns.is_empty():
		samples.append(_plan_sample("", stem))
	for column in reader.sample_columns:
		samples.append(_plan_sample(column, stem if single else column))
	return {"ok": true, "error": "", "info": info, "plan": {"build": detected, "samples": samples}}


## SHA-256 of the decompressed text with line endings normalized, so a file
## and its gzipped copy (or a CRLF copy) are recognized as the same content.
static func content_fingerprint(text_bytes: PackedByteArray) -> String:
	var text := text_bytes.get_string_from_utf8().replace("\r\n", "\n").strip_edges()
	return WorkspaceModel.sha256_hex(text.to_utf8_buffer())


## The file's variant records as "POS:REF:ALT", sorted and unique; used to
## tell whether two files describe the same sample's calls.
static func record_keys(reader: VcfReader) -> PackedStringArray:
	var pos_col: int = reader.column_index["POS"]
	var ref_col: int = reader.column_index["REF"]
	var alt_col: int = reader.column_index["ALT"]
	var seen := {}
	for i in range(reader.first_record_line, reader.lines.size()):
		var line := reader.lines[i]
		if line.is_empty() or line.begins_with("#") or line == "\r":
			continue
		var cells := line.trim_suffix("\r").split("\t")
		if cells.size() > alt_col:
			seen["%s:%s:%s" % [cells[pos_col], cells[ref_col].to_upper(), cells[alt_col].to_upper()]] = true
	var keys := PackedStringArray(seen.keys())
	keys.sort()
	return keys


## File name without .vcf/.vcf.gz and without an Illumina "_S<n>" suffix,
## matching main.py's default_sample_id.
static func default_sample_name(file_name: String) -> String:
	var stem := file_name.get_file()
	for suffix in [".vcf.gz", ".vcf.bgz", ".vcf", ".gz"]:
		if stem.to_lower().ends_with(suffix):
			stem = stem.left(-suffix.length())
			break
	var regex := RegEx.create_from_string("_S\\d+$")
	return regex.sub(stem, "")


## [matching, total] for record REF alleles against a build (evidence only).
static func ref_agreement(reader: VcfReader, build: String) -> Array:
	var reference := ReferenceSequence.get_build(build)
	var pos_col: int = reader.column_index["POS"]
	var ref_col: int = reader.column_index["REF"]
	var matching := 0
	var total := 0
	for i in range(reader.first_record_line, reader.lines.size()):
		var line := reader.lines[i]
		if line.is_empty() or line.begins_with("#") or line == "\r":
			continue
		total += 1
		var ref := line.get_slice("\t", ref_col).to_upper()
		var position := line.get_slice("\t", pos_col).to_int()
		if position >= 1 and reference.supplied_bases(position, ref.length()) == ref:
			matching += 1
	return [matching, total]


## Applies a regex to sample names to set subjects. The first capture group
## is used if present, otherwise the whole match; non-matching names keep
## their current subject. Returns the number of subjects changed, or -1 if
## the pattern is invalid.
static func apply_subject_regex(plan: Dictionary, pattern: String) -> int:
	var regex := RegEx.create_from_string(pattern)
	if not regex.is_valid():
		return -1
	var changed := 0
	for sample in plan["samples"]:
		var match := regex.search(sample["name"])
		if match == null:
			continue
		var subject := (match.get_string(1) if match.get_group_count() >= 1 else match.get_string()).strip_edges()
		if subject and subject != sample["subject"]:
			sample["subject"] = subject
			changed += 1
	return changed


static func split_tags(text: String) -> PackedStringArray:
	var tags := PackedStringArray()
	for part in RegEx.create_from_string("[,;|]").sub(text, ",", true).split(","):
		var tag := part.strip_edges()
		if tag and not tag in tags:
			tags.append(tag)
	return tags


## Returns a ready mapper from `build` onto rCRS, building it on first use.
static func mapper_for(build: String) -> MitoMapper:
	if _mappers.has(build):
		var cached: MitoMapper = _mappers[build]
		while not cached.is_ready:
			await Engine.get_main_loop().process_frame
		return cached
	var mapper := MitoMapper.new(ReferenceSequence.get_build(build), ReferenceSequence.get_build(ReferenceSequence.RCRS))
	_mappers[build] = mapper
	var error := await mapper.build(200)
	if error:
		_mappers.erase(build)
		push_error(error)
		return null
	return mapper


## Imports one source using its plan. Returns {"ok", "error", "cancelled",
## "samples", "calls", "status_counts", "skipped_rows"}.
func import_source(source: Dictionary, store: VariantStore) -> Dictionary:
	var plan: Dictionary = source["plan"]
	var build: String = plan.get("build", "")
	if build.is_empty():
		return _fail("Choose the reference build this file was aligned to.")
	progress.emit(0.0, "Reading %s…" % source["name"])
	await Engine.get_main_loop().process_frame
	var unzipped := GzipReader.decompress(source["raw"])
	if not unzipped["ok"]:
		return _fail(unzipped["error"])
	var reader := VcfReader.parse(unzipped["bytes"].get_string_from_utf8())
	if not reader.ok:
		return _fail(reader.error)
	if build != ReferenceSequence.RCRS:
		progress.emit(0.0, "Aligning %s to rCRS…" % ReferenceSequence.label_for(build))
	var mapper: MitoMapper = await mapper_for(build)
	if mapper == null:
		return _fail("Could not align %s to rCRS." % ReferenceSequence.label_for(build))

	var ci: Dictionary = reader.column_index
	var chrom_col: int = ci["CHROM"]
	var pos_col: int = ci["POS"]
	var ref_col: int = ci["REF"]
	var alt_col: int = ci["ALT"]
	var filter_col: int = ci.get("FILTER", -1)
	var format_col: int = ci.get("FORMAT", -1)

	var new_samples: Array = []
	var targets: Array = []  # [column index or -1, sample id]
	for planned in plan["samples"]:
		if not planned.get("include", true):
			continue
		var sample := {
			"id": store.new_sample_id(),
			"source_id": source["id"],
			"column": planned["column"],
			"name": planned["name"],
			"subject": planned["subject"],
			"tags": PackedStringArray(planned["tags"]),
			"build": build,
		}
		if planned.get("attributes"):
			sample["attributes"] = planned["attributes"].duplicate()
		new_samples.append(sample)
		targets.append([ci.get(planned["column"], -1) if planned["column"] else -1, sample["id"]])
	if targets.is_empty():
		return _fail("No samples are selected for import in this file.")

	var batch := VariantStore.new_call_columns()
	var status_counts := {}
	var skipped_rows := 0
	var blank_contig_rows := 0
	var native := build == ReferenceSequence.RCRS
	var last_yield := Time.get_ticks_msec()
	var total_lines := reader.lines.size() - reader.first_record_line
	var lines := reader.lines

	for li in range(reader.first_record_line, lines.size()):
		var line := lines[li]
		if line.is_empty() or line.begins_with("#") or line == "\r":
			continue
		var fields := line.trim_suffix("\r").split("\t")
		if fields.size() <= alt_col:
			skipped_rows += 1
			continue
		var chrom := fields[chrom_col]
		if chrom.is_empty():
			# Some exports leave CHROM blank in mito-only files; the Python
			# pipeline accepted these too (allow_blank_contig).
			blank_contig_rows += 1
		elif not VcfReader.is_mito_contig(chrom):
			skipped_rows += 1
			continue
		var position := fields[pos_col].to_int()
		var ref := fields[ref_col]
		var alts := fields[alt_col].split(",")
		var multi_alt := 1 if alts.size() > 1 else 0
		var filter_id := store.intern(fields[filter_col] if filter_col >= 0 and filter_col < fields.size() else "")
		var format := fields[format_col] if format_col >= 0 and format_col < fields.size() else ""
		var format_id := store.intern(format)
		var keys := format.split(":")
		var gt_index := keys.find("GT")
		var af_index := -1
		for key in AF_KEYS:
			af_index = keys.find(key)
			if af_index >= 0:
				break
		var gt_only := format == "GT"
		var ref_id := store.intern(ref)

		for target in targets:
			var column: int = target[0]
			var value := VcfReader.unquote(fields[column]) if column >= 0 and column < fields.size() else ""
			if gt_only and (value == "0" or value == "." or value.is_empty()):
				continue  # fast path for large haploid genotype matrices
			var parts := value.split(":")
			var gt := parts[gt_index] if gt_index >= 0 and gt_index < parts.size() else ""
			var gt_alleles := gt.replace("|", "/").split("/")
			var af_parts := parts[af_index].split(",") if af_index >= 0 and af_index < parts.size() else PackedStringArray()
			for ai in alts.size():
				var alt := alts[ai]
				if alt.is_empty() or alt == ".":
					continue
				if gt and not str(ai + 1) in gt_alleles:
					continue  # this sample does not carry this ALT
				var mapped := mapper.map_variant(position, ref, alt)
				var status: int
				var allele := VariantStore.NO_ALLELE
				if mapped["status"] == MitoMapper.Status.MAPPED:
					status = VariantStore.CallStatus.NATIVE if native else VariantStore.CallStatus.MAPPED
					allele = store.intern_allele(mapped["position"], mapped["ref"], mapped["alt"])
				else:
					status = mapped["status"] + 1
				var af := NAN
				if ai < af_parts.size() and af_parts[ai].is_valid_float():
					af = af_parts[ai].to_float()
				batch["sample"].append(target[1])
				batch["allele"].append(allele)
				batch["status"].append(status)
				batch["af"].append(af)
				batch["multi_alt"].append(multi_alt)
				batch["filter"].append(filter_id)
				batch["vcf_pos"].append(position)
				batch["vcf_ref"].append(ref_id)
				batch["vcf_alt"].append(store.intern(alt))
				batch["format"].append(format_id)
				batch["values"].append(value)
				var status_name: String = VariantStore.CALL_STATUS_NAMES[status]
				status_counts[status_name] = status_counts.get(status_name, 0) + 1

		if Time.get_ticks_msec() - last_yield > FRAME_BUDGET_MSEC:
			progress.emit(float(li - reader.first_record_line) / maxi(total_lines, 1),
				"Importing %s…" % source["name"])
			await Engine.get_main_loop().process_frame
			last_yield = Time.get_ticks_msec()
			if cancelled:
				return {"ok": false, "cancelled": true, "error": "Import cancelled."}

	store.append(new_samples, batch)
	progress.emit(1.0, "Imported %s." % source["name"])
	return {
		"ok": true, "error": "", "cancelled": false,
		"samples": new_samples.size(), "calls": batch["sample"].size(),
		"status_counts": status_counts, "skipped_rows": skipped_rows,
		"blank_contig_rows": blank_contig_rows,
	}


static func _plan_sample(column: String, name: String) -> Dictionary:
	return {"column": column, "name": name, "subject": name, "tags": PackedStringArray(), "include": true}


static func _fail(message: String) -> Dictionary:
	return {"ok": false, "cancelled": false, "error": message}
