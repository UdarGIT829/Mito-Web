class_name AnnotationDb
extends RefCounted
## Public annotations keyed by rCRS-normalized allele "pos:REF:ALT".
##
## The bundled file (tools/build_annotations.py) holds gnomAD v3.1 mtDNA and
## ClinVar. `live` holds results fetched in this workspace from Ensembl VEP
## for alleles the bundle does not cover; it is saved with the workspace.

const BUNDLE_PATH := "res://annotations/mito_annotations.json.gz"
const SOURCES := ["gnomad", "clinvar", "ensembl"]

static var _bundle := {}
static var _records := {}
static var _load_error := ""
static var live := {}  ## key -> {"ensembl": {...}}; owned by the workspace data.


static func key(position: int, ref: String, alt: String) -> String:
	return "%d:%s:%s" % [position, ref, alt]


## Loads the bundle on first use. Returns "" or an error message.
static func ensure_loaded() -> String:
	if not _bundle.is_empty() or _load_error:
		return _load_error
	var raw := FileAccess.get_file_as_bytes(BUNDLE_PATH)
	if raw.is_empty():
		_load_error = "The annotation bundle is missing from this build."
		return _load_error
	var unzipped := GzipReader.decompress(raw)
	var parsed = JSON.parse_string(unzipped["bytes"].get_string_from_utf8()) if unzipped["ok"] else null
	if typeof(parsed) != TYPE_DICTIONARY or parsed.get("format") != "mito-web-annotations":
		_load_error = "The annotation bundle could not be read."
		return _load_error
	_bundle = parsed
	_records = parsed["records"]
	return ""


static func source_info() -> Dictionary:
	ensure_loaded()
	return _bundle.get("sources", {})


## Annotations for one allele from the requested sources ("gnomad",
## "clinvar", "ensembl" = live Ensembl VEP results). Missing sources are
## simply absent; absence is not a classification.
static func lookup(position: int, ref: String, alt: String, sources: Array = SOURCES) -> Dictionary:
	ensure_loaded()
	var k := key(position, ref, alt)
	var out := {}
	var bundled: Dictionary = _records.get(k, {})
	for source in sources:
		if source == "ensembl":
			if live.has(k) and live[k].has("ensembl"):
				out["ensembl"] = live[k]["ensembl"]
		elif bundled.has(source):
			out[source] = bundled[source]
	return out


## Consequence terms from gnomAD's VEP, else from a live Ensembl lookup.
static func consequences(annotations: Dictionary) -> Array:
	if annotations.has("gnomad") and annotations["gnomad"].get("consequence"):
		return annotations["gnomad"]["consequence"]
	if annotations.has("ensembl"):
		return annotations["ensembl"].get("consequence", [])
	return []


static func needs_live_lookup(position: int, ref: String, alt: String) -> bool:
	ensure_loaded()
	var k := key(position, ref, alt)
	var bundled: Dictionary = _records.get(k, {})
	return not (bundled.has("gnomad") and bundled["gnomad"].get("consequence")) and not live.has(k)


static func record_count() -> int:
	ensure_loaded()
	return _records.size()
