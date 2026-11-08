class_name CallFeatures
extends RefCounted
## Per-call sequencing-slippage metadata, with the same keys and meaning as
## the Python importer and table pipeline:
##   POLYMORPHISM            "1" when the VCF row had several ALTs
##   REPEAT_n_BASE(S)/_COUNT longest 1/2/3-base repeat in the VCF ALT allele
##   REFERENCE_6_BEFORE/AFTER six reference bases either side of the VCF
##                           position, on the reference the file was called
##                           against (identical to Python for rCRS files)

const SIGNAL_KEYS := [
	"REFERENCE_6_BEFORE", "REFERENCE_6_AFTER",
	"REPEAT_1_BASE", "REPEAT_1_BASE_COUNT",
	"REPEAT_2_BASES", "REPEAT_2_BASES_COUNT",
	"REPEAT_3_BASES", "REPEAT_3_BASES_COUNT",
]

var _store: VariantStore
var _sample_build := {}
var _cache := {}


func _init(store: VariantStore) -> void:
	_store = store
	for sample in store.samples:
		_sample_build[sample["id"]] = sample["build"]


func metadata(call_index: int) -> Dictionary:
	var c := _store.calls
	var build: String = _sample_build.get(c["sample"][call_index], ReferenceSequence.RCRS)
	var position: int = c["vcf_pos"][call_index]
	var alt: String = _store.strings[c["vcf_alt"][call_index]]
	var multi: int = c["multi_alt"][call_index]
	var key := "%s:%d:%s:%d" % [build, position, alt, multi]
	if _cache.has(key):
		return _cache[key]
	var profile := SequenceFeatures.repeat_profile(alt)
	var meta := {
		"POLYMORPHISM": str(multi),
		"REPEAT_1_BASE_COUNT": str(profile[1][1]),
		"REPEAT_2_BASES_COUNT": str(profile[2][1]),
		"REPEAT_3_BASES_COUNT": str(profile[3][1]),
	}
	# Empty motifs are omitted, as in the Python repository layer.
	if profile[1][0]:
		meta["REPEAT_1_BASE"] = profile[1][0]
	if profile[2][0]:
		meta["REPEAT_2_BASES"] = profile[2][0]
	if profile[3][0]:
		meta["REPEAT_3_BASES"] = profile[3][0]
	var reference := ReferenceSequence.get_build(build)
	if position >= 1 and position <= reference.supplied_length():
		var wrapped := posmod(position - 1, reference.length()) + 1
		var context := SequenceFeatures.reference_context(reference.sequence, wrapped)
		meta["REFERENCE_6_BEFORE"] = context[0]
		meta["REFERENCE_6_AFTER"] = context[1]
	_cache[key] = meta
	return meta


## Summarizes carrier-call signals like the Python _slippage_result.
## Risk is "unavailable" or "unclassified" until a classification is defined.
static func slippage(metadata_rows: Array) -> Dictionary:
	var signals := {}
	for key in SIGNAL_KEYS:
		var values := PackedStringArray()
		for meta in metadata_rows:
			var value := str(meta.get(key, "")).strip_edges()
			if value and not value in values:
				values.append(value)
		if values:
			signals[key] = values
	if signals.is_empty():
		return {"risk": "unavailable", "footprint": "", "signals": {}}
	var parts := PackedStringArray()
	for key in signals:
		parts.append("%s=%s" % [key, ",".join(signals[key])])
	return {"risk": "unclassified", "footprint": "; ".join(parts), "signals": signals}
