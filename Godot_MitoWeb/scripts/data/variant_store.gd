class_name VariantStore
extends RefCounted
## Columnar store of imported samples and their calls.
##
## Every call keeps its original VCF identity (position, REF, ALT) and points
## at an rCRS-normalized allele when it could be placed on rCRS. Calls that
## could not be placed keep allele = NO_ALLELE and a status explaining why;
## they are retained for transparency and excluded from analysis.

const NO_ALLELE := -1

## Call placement status. Values after MAPPED mirror MitoMapper.Status + 1.
enum CallStatus { NATIVE, MAPPED, TARGET_REFERENCE, AMBIGUOUS, REFERENCE_MISMATCH, UNSUPPORTED, INVALID_INPUT }
const CALL_STATUS_NAMES := ["NATIVE", "MAPPED", "TARGET_REFERENCE", "AMBIGUOUS", "REFERENCE_MISMATCH", "UNSUPPORTED", "INVALID_INPUT"]

const CALL_COLUMNS := {
	"sample": TYPE_PACKED_INT32_ARRAY,  # sample id
	"allele": TYPE_PACKED_INT32_ARRAY,  # index into allele_* or NO_ALLELE
	"status": TYPE_PACKED_BYTE_ARRAY,  # CallStatus
	"af": TYPE_PACKED_FLOAT64_ARRAY,  # allele fraction; NAN when absent
	"multi_alt": TYPE_PACKED_BYTE_ARRAY,  # VCF row had more than one ALT
	"filter": TYPE_PACKED_INT32_ARRAY,  # string id
	"vcf_pos": TYPE_PACKED_INT32_ARRAY,
	"vcf_ref": TYPE_PACKED_INT32_ARRAY,  # string id
	"vcf_alt": TYPE_PACKED_INT32_ARRAY,  # string id
	"format": TYPE_PACKED_INT32_ARRAY,  # string id of the FORMAT keys
	"values": TYPE_PACKED_STRING_ARRAY,  # raw sample column value
}

var strings := PackedStringArray()
var allele_pos := PackedInt32Array()
var allele_ref := PackedStringArray()
var allele_alt := PackedStringArray()
var samples: Array = []  ## {id, source_id, column, name, subject, tags, build}
var next_sample_id := 1
var calls := {}

var _string_index := {}
var _allele_index := {}
var _sample_index := {}  ## sample id -> sample dictionary


func _init() -> void:
	calls = new_call_columns()


static func new_call_columns() -> Dictionary:
	var columns := {}
	for name in CALL_COLUMNS:
		columns[name] = type_convert([], CALL_COLUMNS[name])
	return columns


func call_count() -> int:
	return calls["sample"].size()


func intern(text: String) -> int:
	var id = _string_index.get(text)
	if id == null:
		id = strings.size()
		strings.append(text)
		_string_index[text] = id
	return id


## Id of an already-interned string, or -1.
func string_id(text: String) -> int:
	return _string_index.get(text, -1)


func intern_allele(position: int, ref: String, alt: String) -> int:
	var key := "%d:%s:%s" % [position, ref, alt]
	var id = _allele_index.get(key)
	if id == null:
		id = allele_pos.size()
		allele_pos.append(position)
		allele_ref.append(ref)
		allele_alt.append(alt)
		_allele_index[key] = id
	return id


func new_sample_id() -> int:
	next_sample_id += 1
	return next_sample_id - 1


## Adds samples and a batch of calls (same columns as `calls`) in one step.
func append(new_samples: Array, batch: Dictionary) -> void:
	samples.append_array(new_samples)
	_reindex_samples()
	for name in CALL_COLUMNS:
		calls[name].append_array(batch[name])


func sample_by_id(sample_id: int) -> Dictionary:
	if _sample_index.size() != samples.size():
		_reindex_samples()
	return _sample_index.get(sample_id, {})


func _reindex_samples() -> void:
	_sample_index.clear()
	for sample in samples:
		_sample_index[sample["id"]] = sample


func samples_for_source(source_id: int) -> Array:
	return samples.filter(func(s: Dictionary) -> bool: return s["source_id"] == source_id)


## Removes a source's samples and their calls.
func remove_source(source_id: int) -> void:
	var removed := {}
	for sample in samples_for_source(source_id):
		removed[sample["id"]] = true
	if removed.is_empty():
		return
	samples = samples.filter(func(s: Dictionary) -> bool: return not removed.has(s["id"]))
	_reindex_samples()
	var kept := new_call_columns()
	var sample_col: PackedInt32Array = calls["sample"]
	for i in sample_col.size():
		if removed.has(sample_col[i]):
			continue
		for name in CALL_COLUMNS:
			kept[name].append(calls[name][i])
	calls = kept


## Per-sample totals: {sample_id: {"calls": int, "placed": int, "status": {name: count}}}.
func sample_summaries() -> Dictionary:
	var out := {}
	for sample in samples:
		out[sample["id"]] = {"calls": 0, "placed": 0, "status": {}}
	var sample_col: PackedInt32Array = calls["sample"]
	var status_col: PackedByteArray = calls["status"]
	var allele_col: PackedInt32Array = calls["allele"]
	for i in sample_col.size():
		var summary: Dictionary = out.get(sample_col[i], {})
		if summary.is_empty():
			continue
		summary["calls"] += 1
		if allele_col[i] != NO_ALLELE:
			summary["placed"] += 1
		var status_name: String = CALL_STATUS_NAMES[status_col[i]]
		summary["status"][status_name] = summary["status"].get(status_name, 0) + 1
	return out


func to_dict() -> Dictionary:
	return {
		"strings": strings,
		"allele_pos": allele_pos,
		"allele_ref": allele_ref,
		"allele_alt": allele_alt,
		"samples": samples,
		"next_sample_id": next_sample_id,
		"calls": calls,
	}


## Returns {"ok": true, "store": VariantStore} or {"ok": false, "error": String}.
static func from_dict(data: Dictionary) -> Dictionary:
	var store := VariantStore.new()
	var expected := {
		"strings": TYPE_PACKED_STRING_ARRAY, "allele_pos": TYPE_PACKED_INT32_ARRAY,
		"allele_ref": TYPE_PACKED_STRING_ARRAY, "allele_alt": TYPE_PACKED_STRING_ARRAY,
		"samples": TYPE_ARRAY, "next_sample_id": TYPE_INT, "calls": TYPE_DICTIONARY,
	}
	for key in expected:
		if typeof(data.get(key)) != expected[key]:
			return {"ok": false, "error": "The variant store is missing or has a malformed \"%s\"." % key}
	var call_data: Dictionary = data["calls"]
	if typeof(call_data.get("af")) == TYPE_PACKED_FLOAT32_ARRAY:
		call_data["af"] = _af_from_float32(call_data["af"])
	var rows := -1
	for name in CALL_COLUMNS:
		if typeof(call_data.get(name)) != CALL_COLUMNS[name]:
			return {"ok": false, "error": "The variant store has a malformed call column \"%s\"." % name}
		if rows >= 0 and call_data[name].size() != rows:
			return {"ok": false, "error": "The variant store's call columns have different lengths."}
		rows = call_data[name].size()
	store.strings = data["strings"]
	store.allele_pos = data["allele_pos"]
	store.allele_ref = data["allele_ref"]
	store.allele_alt = data["allele_alt"]
	store.samples = data["samples"]
	store.next_sample_id = data["next_sample_id"]
	store.calls = call_data
	for i in store.strings.size():
		store._string_index[store.strings[i]] = i
	for i in store.allele_pos.size():
		store._allele_index["%d:%s:%s" % [store.allele_pos[i], store.allele_ref[i], store.allele_alt[i]]] = i
	return {"ok": true, "store": store}


## Early format-2 files stored AF as 32-bit floats (0.8 became 0.800000012).
## Rounding to the 7 significant digits a float32 carries restores the
## decimal values VCFs write, so AF boundary comparisons stay exact.
static func _af_from_float32(values: PackedFloat32Array) -> PackedFloat64Array:
	var out := PackedFloat64Array()
	out.resize(values.size())
	for i in values.size():
		var x: float = values[i]
		if is_nan(x) or x == 0.0:
			out[i] = x
			continue
		var digits := 6 - int(floor(log(absf(x)) / log(10.0)))
		var scale := pow(10.0, digits)
		out[i] = roundf(x * scale) / scale
	return out
