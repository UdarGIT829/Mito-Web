class_name MutationFrame
extends RefCounted
## The table that flows between analysis steps: one row per rCRS mutation
## (normalized allele) carried by at least one sample in the frame's sample
## universe, plus which call each carrier sample contributes.
##
## Only calls placed on rCRS take part. If one sample has several calls that
## normalize to the same allele, the first is used and the rest are counted
## in `duplicate_calls`.
##
## `flags` holds named per-row true/false columns written by filters running
## in flag mode, so a table can show why a row would be excluded.

var store: VariantStore
var samples := PackedInt32Array()  ## Sample universe (denominator), in order.
var alleles := PackedInt32Array()  ## Row allele ids, sorted by position/REF/ALT.
var carriers: Array = []  ## Per row: {sample_id: call_index}.
var flags := {}  ## name -> PackedByteArray parallel to rows.
var duplicate_calls := 0


func row_count() -> int:
	return alleles.size()


## Builds a frame from the calls of the given samples.
## Options: "filter_values" (Array of VCF FILTER values to keep; empty keeps
## all), "positions" ([min, max] to pre-limit rows).
static func from_samples(source: VariantStore, sample_ids: PackedInt32Array, options: Dictionary = {}) -> MutationFrame:
	var frame := MutationFrame.new()
	frame.store = source
	frame.samples = sample_ids
	var wanted := {}
	for id in sample_ids:
		wanted[id] = true
	var allowed_filters := {}
	var filter_values: Array = Array(options.get("filter_values", []))
	for value in filter_values:
		allowed_filters[source.string_id(value)] = true
	var by_allele := {}
	var c := source.calls
	var sample_col: PackedInt32Array = c["sample"]
	var allele_col: PackedInt32Array = c["allele"]
	var filter_col: PackedInt32Array = c["filter"]
	for i in sample_col.size():
		var allele := allele_col[i]
		if allele == VariantStore.NO_ALLELE or not wanted.has(sample_col[i]):
			continue
		if filter_values and not allowed_filters.has(filter_col[i]):
			continue
		var row: Dictionary = by_allele.get_or_add(allele, {})
		if row.has(sample_col[i]):
			frame.duplicate_calls += 1
		else:
			row[sample_col[i]] = i
	var keys := by_allele.keys()
	_sort_alleles(source, keys)
	for allele in keys:
		frame.alleles.append(allele)
		frame.carriers.append(by_allele[allele])
	return frame


## A frame with the same universe and a subset of rows (indices in order).
func take_rows(rows: Array) -> MutationFrame:
	var out := _empty_like(samples)
	for r in rows:
		out.alleles.append(alleles[r])
		out.carriers.append(carriers[r])
		for name in flags:
			out.flags[name].append(flags[name][r])
	return out


func carrier_count(row: int) -> int:
	return carriers[row].size()


func row_index(allele: int) -> int:
	return alleles.find(allele)


## Rows of this frame whose mutation is in `other` (sample universe and
## carriers are this frame's).
func intersect(other: MutationFrame) -> MutationFrame:
	var keep := _allele_set(other)
	return take_rows(range(row_count()).filter(func(r: int) -> bool: return keep.has(alleles[r])))


## Rows of this frame whose mutation is not in `other`.
func subtract(other: MutationFrame) -> MutationFrame:
	var drop := _allele_set(other)
	return take_rows(range(row_count()).filter(func(r: int) -> bool: return not drop.has(alleles[r])))


## Mutations in either frame; the universe and carriers are combined.
func union(other: MutationFrame) -> MutationFrame:
	var universe := samples.duplicate()
	for id in other.samples:
		if not id in universe:
			universe.append(id)
	var out := _empty_like(universe)
	out.flags = {}
	var merged := {}
	for frame in [self, other]:
		for r in frame.row_count():
			var row: Dictionary = merged.get_or_add(frame.alleles[r], {})
			for sample_id in frame.carriers[r]:
				if not row.has(sample_id):
					row[sample_id] = frame.carriers[r][sample_id]
	var keys := merged.keys()
	_sort_alleles(store, keys)
	for allele in keys:
		out.alleles.append(allele)
		out.carriers.append(merged[allele])
	return out


func set_flag(name: String, values: PackedByteArray) -> void:
	flags[name] = values


func _empty_like(universe: PackedInt32Array) -> MutationFrame:
	var out := MutationFrame.new()
	out.store = store
	out.samples = universe
	out.duplicate_calls = duplicate_calls
	for name in flags:
		out.flags[name] = PackedByteArray()
	return out


static func _allele_set(frame: MutationFrame) -> Dictionary:
	var out := {}
	for allele in frame.alleles:
		out[allele] = true
	return out


static func _sort_alleles(source: VariantStore, keys: Array) -> void:
	keys.sort_custom(func(a: int, b: int) -> bool:
		if source.allele_pos[a] != source.allele_pos[b]:
			return source.allele_pos[a] < source.allele_pos[b]
		if source.allele_ref[a] != source.allele_ref[b]:
			return source.allele_ref[a] < source.allele_ref[b]
		return source.allele_alt[a] < source.allele_alt[b])
