class_name UploadBundle
extends RefCounted
## Zip uploads and manifests.
##
## A .zip is expanded: every VCF inside (.vcf, .vcf.gz, .bgz) becomes an
## upload, and a CSV/TSV whose name contains "manifest" supplies each file's
## import settings. A manifest can also be uploaded next to plain VCFs.
##
## Manifest columns (header row; all but "file" optional):
##   file     VCF file name (as inside the zip)
##   sample   VCF sample column, for multi-sample files (blank = every sample)
##   subject  subject (e.g. cell line, donor, or the sample itself)
##   build    rcrs (also hg38, GRCh37, rCRS), hg19, or hg19_plus58 (hg19+58)
##   tags     tags separated by ; or |
## Lines starting with # are comments.

const TEMP_ZIP := "user://_upload.zip"
const BUILD_ALIASES := {
	"rcrs": ReferenceSequence.RCRS, "hg38": ReferenceSequence.RCRS, "grch38": ReferenceSequence.RCRS,
	"grch37": ReferenceSequence.RCRS, "b37": ReferenceSequence.RCRS, "nc_012920.1": ReferenceSequence.RCRS,
	"hg19": ReferenceSequence.HG19, "hg19_plus58": ReferenceSequence.HG19_PLUS58, "hg19+58": ReferenceSequence.HG19_PLUS58,
}


static func is_vcf_name(name: String) -> bool:
	var lower := name.to_lower()
	return lower.ends_with(".vcf") or lower.ends_with(".vcf.gz") or lower.ends_with(".vcf.bgz") or lower.ends_with(".bgz")


static func is_manifest_name(name: String) -> bool:
	var lower := name.get_file().to_lower()
	return "manifest" in lower and (lower.ends_with(".csv") or lower.ends_with(".tsv") or lower.ends_with(".txt"))


## Expands uploads. Returns {"files": [{name, bytes}], "manifests": [table],
## "errors": [String], "skipped": [String]}.
static func expand(uploads: Array) -> Dictionary:
	var out := {"files": [], "manifests": [], "errors": [], "skipped": []}
	for upload in uploads:
		var name: String = upload["name"]
		if name.to_lower().ends_with(".zip"):
			_expand_zip(upload, out)
		elif is_manifest_name(name):
			_add_manifest(name, upload["bytes"], out)
		else:
			out["files"].append(upload)
	return out


static func _expand_zip(upload: Dictionary, out: Dictionary) -> void:
	var file := FileAccess.open(TEMP_ZIP, FileAccess.WRITE)
	if file == null:
		out["errors"].append("Could not open %s (browser storage unavailable)." % upload["name"])
		return
	file.store_buffer(upload["bytes"])
	file.close()
	var zip := ZIPReader.new()
	if zip.open(TEMP_ZIP) != OK:
		out["errors"].append("%s is not a readable zip file." % upload["name"])
		DirAccess.remove_absolute(TEMP_ZIP)
		return
	for path in zip.get_files():
		if path.ends_with("/") or path.begins_with("__MACOSX") or path.get_file().begins_with("."):
			continue
		if is_vcf_name(path):
			out["files"].append({"name": path.get_file(), "bytes": zip.read_file(path)})
		elif is_manifest_name(path):
			_add_manifest(path.get_file(), zip.read_file(path), out)
		else:
			out["skipped"].append(path)
	zip.close()
	DirAccess.remove_absolute(TEMP_ZIP)


static func _add_manifest(name: String, bytes: PackedByteArray, out: Dictionary) -> void:
	var text := "\n".join(Array(bytes.get_string_from_utf8().split("\n")).filter(func(l): return not l.strip_edges().begins_with("#")))
	var table := TagImport.parse(text)
	if not table["ok"]:
		out["errors"].append("%s: %s" % [name, table["error"]])
		return
	var header: Array = table["header"].map(func(h): return h.to_lower())
	if not "file" in header:
		out["errors"].append("%s has no \"file\" column." % name)
		return
	table["columns"] = {}
	for key in ["file", "sample", "subject", "build", "tags"]:
		table["columns"][key] = header.find(key)
	table["name"] = name
	out["manifests"].append(table)


## Applies manifests to sources' import plans (pending sources only).
## Returns {"applied": rows used, "unmatched": [file names], "errors": [..]}.
static func apply_manifests(sources: Array, manifests: Array) -> Dictionary:
	var by_file := {}
	for source in sources:
		by_file[source["name"]] = source
	var applied := 0
	var unmatched: Array = []
	var errors: Array = []
	for table in manifests:
		var col: Dictionary = table["columns"]
		for row in table["rows"]:
			var file_name: String = row[col["file"]]
			if file_name.is_empty():
				continue
			var source: Dictionary = by_file.get(file_name, by_file.get(file_name.get_file(), {}))
			if source.is_empty():
				if not file_name in unmatched:
					unmatched.append(file_name)
				continue
			if source.get("status") != "pending":
				continue
			var plan: Dictionary = source["plan"]
			if col["build"] >= 0 and row[col["build"]] != "":
				var build: String = BUILD_ALIASES.get(row[col["build"]].strip_edges().to_lower(), "")
				if build.is_empty():
					errors.append("%s: unknown build \"%s\"." % [file_name, row[col["build"]]])
				else:
					plan["build"] = build
			var column: String = row[col["sample"]] if col["sample"] >= 0 else ""
			for planned in plan["samples"]:
				if column != "" and planned["column"] != column:
					continue
				if col["subject"] >= 0 and row[col["subject"]] != "":
					planned["subject"] = row[col["subject"]]
				if col["tags"] >= 0 and row[col["tags"]] != "":
					var tags := PackedStringArray(planned["tags"])
					for tag in VcfImporter.split_tags(row[col["tags"]]):
						if not tag in tags:
							tags.append(tag)
					planned["tags"] = tags
			applied += 1
	return {"applied": applied, "unmatched": unmatched, "errors": errors}
