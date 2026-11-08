class_name Comparison
extends RefCounted
## Sample comparison, ported from mito_viewer/domain/comparison.py.
##
## Each mutation in the frame gets a status from how many of the compared
## samples carry it: "common" (all), "unique" (exactly one) or "partial".
## Per-sample rules then constrain which mutations pass:
##   "present"  the sample must carry it
##   "not_in"   the sample must not carry it
##   "unique"   it passes if this sample is its only carrier
## A sample without explicit rules defaults to "present" (Python behavior);
## an empty rule set means the sample is unconstrained.

const DEFAULT_STATUSES := ["common", "partial", "unique"]


static func status_for(present_count: int, total_count: int) -> String:
	if present_count == total_count:
		return "common"
	if present_count == 1:
		return "unique"
	return "partial"


static func sample_rules_match(compare_ids: PackedInt32Array, sample_rules: Dictionary, present: Dictionary) -> bool:
	var present_required: Array = []
	var absent_required: Array = []
	var unique_allowed: Array = []
	for sample_id in compare_ids:
		var allowed: Array = sample_rules.get(sample_id, ["present"])
		if allowed.is_empty():
			continue
		if "not_in" in allowed:
			absent_required.append(sample_id)
		if "present" in allowed:
			present_required.append(sample_id)
		if "unique" in allowed:
			unique_allowed.append(sample_id)
	if absent_required.any(func(id) -> bool: return present.has(id)):
		return false
	var present_branch := not present_required.is_empty() and present_required.all(func(id) -> bool: return present.has(id))
	var unique_branch := not unique_allowed.is_empty() and present.size() == 1 \
		and unique_allowed.any(func(id) -> bool: return present.has(id))
	if present_required or unique_allowed:
		return present_branch or unique_branch
	return true


## Returns {"frame": MutationFrame, "statuses": PackedStringArray} with the
## passing rows. The frame's samples are the compared samples.
## `sample_rules`: {sample_id: ["present" | "not_in" | "unique", ...]}.
static func compare(frame: MutationFrame, compare_ids: PackedInt32Array, statuses: Array = DEFAULT_STATUSES, sample_rules: Dictionary = {}) -> Dictionary:
	var empty := {"frame": frame.take_rows([]), "statuses": PackedStringArray()}
	empty["frame"].samples = compare_ids
	if compare_ids.size() < 2:
		return empty
	var wanted := {}
	for id in compare_ids:
		wanted[id] = true
	var rows: Array = []
	var row_statuses := PackedStringArray()
	var restricted_carriers: Array = []
	for r in frame.row_count():
		var present := {}
		for sample_id in frame.carriers[r]:
			if wanted.has(sample_id):
				present[sample_id] = frame.carriers[r][sample_id]
		if present.is_empty():
			continue
		var status := status_for(present.size(), compare_ids.size())
		if not status in statuses:
			continue
		if not sample_rules_match(compare_ids, sample_rules, present):
			continue
		rows.append(r)
		row_statuses.append(status)
		restricted_carriers.append(present)
	var out := frame.take_rows(rows)
	out.samples = compare_ids
	out.carriers = restricted_carriers
	return {"frame": out, "statuses": row_statuses}
