class_name ReferenceSequence
extends RefCounted
## Mitochondrial reference builds bundled with the app. rCRS (NC_012920.1,
## identical to hg38/GRCh37 chrM/MT) is the analysis coordinate system; every
## imported call is projected onto it.

const RCRS := "rcrs"
const HG19 := "hg19"
const HG19_PLUS58 := "hg19_plus58"

const _SPECS := {
	RCRS: {
		"label": "rCRS / hg38 / GRCh37 (NC_012920.1, 16,569 bp)",
		"fasta": "res://reference/rcrs_chrM.fa",
		"inferred": false,
		"circular_extension": 0,
		"provenance": "UCSC hg38 chrM / NC_012920.1",
		"sha256": "f156ff3f65bbcc80c7ebb9936dceb96b1477b4f8f535c4e1dbe7baea225cbc66",
	},
	HG19: {
		"label": "UCSC hg19 chrM (Yoruba, 16,571 bp)",
		"fasta": "res://reference/hg19_chrM.fa",
		"inferred": false,
		"circular_extension": 0,
		"provenance": "UCSC hg19 chrM",
		"sha256": "8264337d37c28c396e70204867ede43ba223090ff5d27b46625c1752d6d8bd96",
	},
	HG19_PLUS58: {
		"label": "UCSC hg19 chrM + 58-base circular extension (16,629 bp)",
		"fasta": "res://reference/hg19_chrM.fa",
		"inferred": true,
		"circular_extension": 58,
		"provenance": "UCSC hg19 chrM plus inferred first 58 bases appended (Illumina 'withExtra58' references)",
		"sha256": "8264337d37c28c396e70204867ede43ba223090ff5d27b46625c1752d6d8bd96",
	},
}

static var _loaded := {}

var key: String
var label: String
var sequence: String  ## Unextended core sequence, uppercase.
var inferred: bool
var circular_extension: int
var provenance: String


static func build_keys() -> Array:
	return [RCRS, HG19, HG19_PLUS58]


static func label_for(build_key: String) -> String:
	return _SPECS[build_key]["label"] if _SPECS.has(build_key) else build_key


static func get_build(build_key: String) -> ReferenceSequence:
	if not _loaded.has(build_key):
		var spec: Dictionary = _SPECS[build_key]
		var ref := ReferenceSequence.new()
		ref.key = build_key
		ref.label = spec["label"]
		ref.sequence = _read_fasta(spec["fasta"])
		ref.inferred = spec["inferred"]
		ref.circular_extension = spec["circular_extension"]
		ref.provenance = spec["provenance"]
		if spec["sha256"] and ref.sha256() != spec["sha256"]:
			push_error("Reference %s does not match its pinned fingerprint." % build_key)
		_loaded[build_key] = ref
	return _loaded[build_key]


## Picks a build from the VCF ##contig length. Returns "" when unknown; the
## user must then choose explicitly. REF agreement alone is never used to
## choose a build.
static func detect_from_contig_length(length: int) -> String:
	match length:
		16569:
			return RCRS
		16571:
			return HG19
		16629:
			return HG19_PLUS58
	return ""


func length() -> int:
	return sequence.length()


## Length of the coordinate space a VCF may use, including the extension.
func supplied_length() -> int:
	return sequence.length() + circular_extension


## Reference bases at a 1-based position in the supplied coordinate space.
func supplied_bases(position: int, count: int) -> String:
	var supplied := sequence + sequence.substr(0, circular_extension)
	return supplied.substr(position - 1, count)


func sha256() -> String:
	var ctx := HashingContext.new()
	ctx.start(HashingContext.HASH_SHA256)
	ctx.update(sequence.to_ascii_buffer())
	return ctx.finish().hex_encode()


static func _read_fasta(path: String) -> String:
	var parts := PackedStringArray()
	for line in FileAccess.get_file_as_string(path).split("\n"):
		line = line.strip_edges()
		if line and not line.begins_with(">"):
			parts.append(line.to_upper())
	return "".join(parts)
