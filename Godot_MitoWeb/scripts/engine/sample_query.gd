class_name SampleQuery
extends RefCounted
## Selects samples from a VariantStore. All given criteria must hold.
##
## Criteria (all optional):
##   tags_all     sample has every tag (Python viewer semantics)
##   tags_any     sample has at least one of the tags
##   tags_none    sample has none of the tags
##   untagged     true: sample has no tags (the viewer's "<NONE>")
##   subjects     sample subject is one of these
##   source_ids   sample came from one of these uploaded files
##   sample_ids   explicit sample ids
## Results keep store (import) order.

const UNTAGGED := "<NONE>"


static func select(store: VariantStore, criteria: Dictionary = {}) -> PackedInt32Array:
	var tags_all: Array = Array(criteria.get("tags_all", []))
	var untagged: bool = criteria.get("untagged", false) or UNTAGGED in tags_all
	tags_all.erase(UNTAGGED)
	var tags_any: Array = Array(criteria.get("tags_any", []))
	var tags_none: Array = Array(criteria.get("tags_none", []))
	var subjects: Array = Array(criteria.get("subjects", []))
	var source_ids: Array = Array(criteria.get("source_ids", []))
	var sample_ids: Array = Array(criteria.get("sample_ids", []))
	var out := PackedInt32Array()
	for sample in store.samples:
		var tags: PackedStringArray = sample["tags"]
		if untagged and not tags.is_empty():
			continue
		if not tags_all.all(func(t) -> bool: return t in tags):
			continue
		if tags_any and not tags_any.any(func(t) -> bool: return t in tags):
			continue
		if tags_none.any(func(t) -> bool: return t in tags):
			continue
		if subjects and not sample["subject"] in subjects:
			continue
		if source_ids and not sample["source_id"] in source_ids:
			continue
		if sample_ids and not sample["id"] in sample_ids:
			continue
		out.append(sample["id"])
	return out


## Every tag in use, sorted, with sample counts: [[tag, count], ...].
static func tag_counts(store: VariantStore) -> Array:
	var counts := {}
	var untagged := 0
	for sample in store.samples:
		if sample["tags"].is_empty():
			untagged += 1
		for tag in sample["tags"]:
			counts[tag] = counts.get(tag, 0) + 1
	var tags := counts.keys()
	tags.sort_custom(func(a, b) -> bool: return a.naturalnocasecmp_to(b) < 0)
	var out: Array = tags.map(func(t) -> Array: return [t, counts[t]])
	if untagged:
		out.append([UNTAGGED, untagged])
	return out
