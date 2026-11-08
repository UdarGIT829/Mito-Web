class_name VcfReader
extends RefCounted
## Reads VCF text: header facts, column layout, and record lines.
## Columns are found by name (some callers omit ID). Sample columns are the
## columns after FORMAT. Values wrapped in double quotes are unquoted.

const MITO_CONTIGS := ["chrM", "MT", "M", "chrMT", "NC_012920.1", "NC_012920"]

var ok := false
var error := ""
var columns := PackedStringArray()
var column_index := {}
var sample_columns := PackedStringArray()
var contig_id := ""
var contig_length := 0
var reference_line := ""
var source_line := ""
var lines := PackedStringArray()  ## All text lines.
var first_record_line := 0  ## Index into `lines` of the first data line.
var record_count := 0
var blank_contig_records := 0  ## Records with an empty CHROM field.
var other_contig_records := 0  ## Records on a non-mitochondrial contig.


static func parse(text: String) -> VcfReader:
	var reader := VcfReader.new()
	reader._parse(text)
	return reader


func _parse(text: String) -> void:
	lines = text.split("\n")
	var header_found := false
	var contigs := {}
	for i in lines.size():
		var line := lines[i].trim_suffix("\r")
		if line.begins_with("##"):
			if line.begins_with("##contig=<"):
				var id := _header_field(line, "ID")
				if id:
					contigs[id] = _header_field(line, "length").to_int()
			elif line.begins_with("##reference="):
				reference_line = line.substr(12)
			elif line.begins_with("##source="):
				source_line = line.substr(9)
			continue
		if line.begins_with("#CHROM"):
			columns = line.substr(1).split("\t")
			for c in columns.size():
				column_index[columns[c]] = c
			header_found = true
			first_record_line = i + 1
			break
		if line.strip_edges():
			error = "Data appears before the #CHROM header line."
			return
	if not header_found:
		error = "No #CHROM header line was found. Is this a VCF file?"
		return
	for required in ["CHROM", "POS", "REF", "ALT"]:
		if not column_index.has(required):
			error = "The VCF header has no %s column." % required
			return
	if column_index.has("FORMAT"):
		sample_columns = columns.slice(column_index["FORMAT"] + 1)
	for id in MITO_CONTIGS:
		if contigs.has(id):
			contig_id = id
			contig_length = contigs[id]
			break
	var chrom_col: int = column_index["CHROM"]
	for i in range(first_record_line, lines.size()):
		var line := lines[i]
		if line and not line.begins_with("#") and line != "\r":
			record_count += 1
			var chrom := line.get_slice("\t", chrom_col)
			if chrom.is_empty():
				blank_contig_records += 1
			elif not is_mito_contig(chrom):
				other_contig_records += 1
	ok = true


static func is_mito_contig(chrom: String) -> bool:
	return chrom in MITO_CONTIGS


static func unquote(value: String) -> String:
	if value.length() >= 2 and value.begins_with("\"") and value.ends_with("\""):
		return value.substr(1, value.length() - 2)
	return value.trim_prefix("\"").trim_suffix("\"")


## Value of KEY=value inside a ##meta=<...> header line.
static func _header_field(line: String, key: String) -> String:
	var body := line.substr(line.find("<") + 1).trim_suffix(">")
	for part in body.split(","):
		var eq := part.find("=")
		if eq > 0 and part.substr(0, eq).strip_edges() == key:
			return part.substr(eq + 1).strip_edges()
	return ""
