class_name MutationTable
extends RefCounted
## Builds display rows for a MutationFrame: one row per mutation with
## mutation-level fields once and AF per carrier sample, ported from the
## Python mutation-centered Table (MutationResult.to_json).

const SORT_FIELDS := ["position", "reference", "alt", "gene", "carrier_count", "presence_percent", "af_passing_percent"]


## Returns {"rows": Array, "matched_count": int, "error": String}.
## options: filters (MutationFilters config), sort ([{field, direction}]),
## offset, limit (default all).
static func build(frame: MutationFrame, options: Dictionary = {}) -> Dictionary:
	var config: Dictionary = options.get("filters", {})
	var problem := MutationFilters.validate(config)
	if problem:
		return {"rows": [], "matched_count": 0, "error": problem}
	var sort: Array = options.get("sort", [{"field": "position", "direction": "asc"}])
	for item in sort:
		if not item.get("field") in SORT_FIELDS or not item.get("direction", "asc") in ["asc", "desc"]:
			return {"rows": [], "matched_count": 0, "error": "Unsupported sort: %s." % item}
	var features := CallFeatures.new(frame.store)
	var evaluations := MutationFilters.evaluate(frame, config, features)
	var rows: Array = []
	for r in frame.row_count():
		if evaluations[r]["passes"]:
			rows.append(row_for(frame, r, evaluations[r]["summary"], features))
	_sort(rows, sort)
	var offset: int = options.get("offset", 0)
	var limit: int = options.get("limit", rows.size())
	return {"rows": rows.slice(offset, offset + limit), "matched_count": rows.size(), "error": ""}


static func row_for(frame: MutationFrame, r: int, summary: Dictionary, features: CallFeatures) -> Dictionary:
	var store := frame.store
	var allele := frame.alleles[r]
	var pos := store.allele_pos[allele]
	var locus := MitoGenes.locus(pos, store.allele_ref[allele])
	var carriers: Dictionary = frame.carriers[r]
	var sample_results: Array = []
	var metadata_rows: Array = []
	for sample_id in frame.samples:  # selected-sample order
		if not carriers.has(sample_id):
			continue
		var call: int = carriers[sample_id]
		var af: float = store.calls["af"][call]
		var valid_af := not is_nan(af) and af >= 0.0 and af <= 1.0
		sample_results.append({
			"sample_id": sample_id,
			"sample_label": store.sample_by_id(sample_id).get("name", str(sample_id)),
			"af": af if valid_af else null,
			"filter": store.strings[store.calls["filter"][call]],
			"call_index": call,
		})
		metadata_rows.append(features.metadata(call))
	var row := {
		"allele": allele,
		"position": pos,
		"reference": ReferenceSequence.get_build(ReferenceSequence.RCRS).sequence[pos - 1],
		"event_ref": store.allele_ref[allele],
		"alt": store.allele_alt[allele],
		"genes": locus["genes"],
		"region": locus["region"],
		"slippage": CallFeatures.slippage(metadata_rows),
		"sample_summary": summary,
		"sample_results": sample_results,
	}
	for name in frame.flags:
		row[name] = frame.flags[name][r] == 1
	return row


static func _sort(rows: Array, sort: Array) -> void:
	rows.sort_custom(func(a: Dictionary, b: Dictionary) -> bool:
		for item in sort:
			var va = _sort_value(a, item["field"])
			var vb = _sort_value(b, item["field"])
			if va != vb:
				return va > vb if item.get("direction", "asc") == "desc" else va < vb
		# Ties fall back to mutation identity, as in the Python stable sort.
		if a["position"] != b["position"]:
			return a["position"] < b["position"]
		if a["event_ref"] != b["event_ref"]:
			return a["event_ref"] < b["event_ref"]
		return a["alt"] < b["alt"])


static func _sort_value(row: Dictionary, field: String):
	match field:
		"position": return row["position"]
		"reference": return row["reference"]
		"alt": return row["alt"]
		"gene": return row["genes"][0] if row["genes"] else ""
		"carrier_count": return row["sample_summary"]["carrier_count"]
		"presence_percent":
			var s: Dictionary = row["sample_summary"]
			return float(s["carrier_count"]) / s["selected_count"] if s["selected_count"] else 0.0
		"af_passing_percent":
			var s: Dictionary = row["sample_summary"]
			return -1.0 if s["af_passing_count"] == null else float(s["af_passing_count"]) / maxi(s["carrier_count"], 1)
	return null
