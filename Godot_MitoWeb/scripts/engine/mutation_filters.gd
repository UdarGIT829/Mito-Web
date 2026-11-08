class_name MutationFilters
extends RefCounted
## Mutation-level filters, ported from the Python Table pipeline
## (mito_viewer/domain/table_requests.py and catalog/mutation_table.py).
##
## Config dictionary (every key optional):
##   position: {exact, minimum, maximum}
##   genes: [symbol, ...]                       any overlap passes
##   reference_bases: [base, ...]               rCRS base at the position
##   alt_alleles: [alt, ...]                    normalized ALT
##   slippage_risks: [risk, ...]
##   slippage_metadata_filters: [{field, value}, ...]   all must match on
##                                              at least one carrier call
##   minimum_presence_percent: 0..100           carriers / selected samples
##   af_ranges: [{minimum, maximum, include_minimum, include_maximum}, ...]
##   minimum_af_percent: 0..100                 carriers with AF in any range
##
## Percentages are inclusive and compared on counts (n*100 >= p*d), never on
## rounded percentages.

const AF_OPERATORS := ["gt", "gte", "lt", "lte", "eq", "neq"]
const METADATA_FIELDS := [
	"polymorphism", "reference_contains_alt", "reference_context", "reference_repeat",
	"repeat_base", "repeat_count", "repeat_2_bases", "repeat_2_count", "repeat_3_bases", "repeat_3_count",
]
const COUNT_FIELDS := ["repeat_count", "repeat_2_count", "repeat_3_count"]
const REFERENCE_REPEAT_VALUES := ["before", "after", "one", "both", "none", "either"]
const METADATA_KEYS := {
	"repeat_base": "REPEAT_1_BASE", "repeat_count": "REPEAT_1_BASE_COUNT",
	"repeat_2_bases": "REPEAT_2_BASES", "repeat_2_count": "REPEAT_2_BASES_COUNT",
	"repeat_3_bases": "REPEAT_3_BASES", "repeat_3_count": "REPEAT_3_BASES_COUNT",
	"polymorphism": "POLYMORPHISM",
}


## Returns "" when the config is valid, otherwise a message.
static func validate(config: Dictionary) -> String:
	var position: Dictionary = config.get("position", {})
	var exact = position.get("exact")
	var minimum = position.get("minimum")
	var maximum = position.get("maximum")
	for value in [exact, minimum, maximum]:
		if value != null and (typeof(value) not in [TYPE_INT, TYPE_FLOAT] or int(value) < 1 or int(value) != value):
			return "Positions must be positive whole numbers."
	if exact != null and (minimum != null or maximum != null):
		return "Use an exact position or a position range, not both."
	if minimum != null and maximum != null and minimum > maximum:
		return "Minimum position cannot exceed maximum position."
	for key in ["minimum_presence_percent", "minimum_af_percent"]:
		var p = config.get(key, 0)
		if typeof(p) not in [TYPE_INT, TYPE_FLOAT] or p < 0 or p > 100:
			return "Percentages must be between 0 and 100."
	var ranges: Array = config.get("af_ranges", [])
	if ranges.is_empty() and config.get("minimum_af_percent", 0) != 0:
		return "An AF carrier percentage needs at least one AF range."
	for r in ranges:
		var lo = r.get("minimum")
		var hi = r.get("maximum")
		if lo == null and hi == null:
			return "An AF range needs a minimum or maximum."
		for bound in [lo, hi]:
			if bound != null and (typeof(bound) not in [TYPE_INT, TYPE_FLOAT] or bound < 0 or bound > 1):
				return "AF range bounds must be between 0 and 1."
		if lo != null and hi != null and lo > hi:
			return "An AF range minimum cannot exceed its maximum."
	for item in config.get("slippage_metadata_filters", []):
		var problem := _validate_metadata_filter(str(item.get("field", "")), str(item.get("value", "")))
		if problem:
			return problem
	return ""


static func _validate_metadata_filter(field: String, value: String) -> String:
	if not field in METADATA_FIELDS:
		return "Unknown metadata filter field: %s." % field
	if value.strip_edges().is_empty():
		return "Metadata filter value cannot be empty."
	if field == "polymorphism" and not value in ["0", "1"]:
		return "Polymorphism filter must be 0 or 1."
	if field == "reference_contains_alt" and not value in ["contains", "not_contains"]:
		return "Reference ALT filter must be contains or not_contains."
	if field == "reference_repeat" and not value in REFERENCE_REPEAT_VALUES:
		return "Unknown reference-repeat filter value."
	if field in COUNT_FIELDS:
		var parts := value.split("|")
		if parts.size() != 2 or not parts[0] in AF_OPERATORS or not parts[1].is_valid_int():
			return "Repeat-count filter must use operator|integer."
	return ""


## Evaluates every row. Returns per-row results:
##   {"passes": bool, "reasons": [String], "summary": {...}}
## `features` may be shared across calls to reuse slippage metadata.
static func evaluate(frame: MutationFrame, config: Dictionary, features: CallFeatures = null) -> Array:
	if features == null:
		features = CallFeatures.new(frame.store)
	var store := frame.store
	var rcrs := ReferenceSequence.get_build(ReferenceSequence.RCRS).sequence
	var position: Dictionary = config.get("position", {})
	var genes_wanted: Array = Array(config.get("genes", []))
	var refs_wanted: Array = Array(config.get("reference_bases", [])).map(func(v): return str(v).to_upper())
	var alts_wanted: Array = Array(config.get("alt_alleles", [])).map(func(v): return str(v).to_upper())
	var risks_wanted: Array = Array(config.get("slippage_risks", []))
	var meta_filters: Array = config.get("slippage_metadata_filters", [])
	var presence_min: float = config.get("minimum_presence_percent", 0)
	var ranges: Array = config.get("af_ranges", [])
	var af_min: float = config.get("minimum_af_percent", 0)
	var selected := frame.samples.size()
	var out: Array = []
	for r in frame.row_count():
		var allele := frame.alleles[r]
		var pos := store.allele_pos[allele]
		var alt := store.allele_alt[allele]
		var reasons: Array = []
		if not _position_matches(position, pos):
			reasons.append("position")
		if genes_wanted:
			var genes: PackedStringArray = MitoGenes.locus(pos, store.allele_ref[allele])["genes"]
			if not genes_wanted.any(func(g) -> bool: return g in genes):
				reasons.append("gene")
		if refs_wanted and not rcrs[pos - 1] in refs_wanted:
			reasons.append("reference")
		if alts_wanted and not alt in alts_wanted:
			reasons.append("alt")
		var carrier_calls: Dictionary = frame.carriers[r]
		if risks_wanted or meta_filters:
			var metadata_rows: Array = []
			for sample_id in carrier_calls:
				metadata_rows.append(features.metadata(carrier_calls[sample_id]))
			if risks_wanted and not CallFeatures.slippage(metadata_rows)["risk"] in risks_wanted:
				reasons.append("slippage_risk")
			if meta_filters and not metadata_rows.any(func(m) -> bool: return metadata_matches(m, meta_filters, alt)):
				reasons.append("slippage_metadata")
		var carriers := carrier_calls.size()
		var summary := {
			"selected_count": selected,
			"carrier_count": carriers,
			"presence_percent": 100.0 * carriers / selected if selected else 0.0,
			"af_passing_count": null,
			"af_passing_percent": null,
		}
		if not percent_passes(carriers, selected, presence_min):
			reasons.append("presence")
		if ranges:
			var passing := 0
			for sample_id in carrier_calls:
				if af_in_any_range(store.calls["af"][carrier_calls[sample_id]], ranges):
					passing += 1
			summary["af_passing_count"] = passing
			summary["af_passing_percent"] = 100.0 * passing / carriers if carriers else 0.0
			if not percent_passes(passing, carriers, af_min):
				reasons.append("af")
		out.append({"passes": reasons.is_empty(), "reasons": reasons, "summary": summary})
	return out


## Applies filters. mode "drop" removes failing rows; mode "flag" keeps every
## row and writes a `flag_name` column (1 = passes) instead.
static func apply(frame: MutationFrame, config: Dictionary, mode: String = "drop", flag_name: String = "passes") -> MutationFrame:
	var results := evaluate(frame, config)
	if mode == "flag":
		var out := frame.take_rows(range(frame.row_count()))
		var column := PackedByteArray()
		for result in results:
			column.append(1 if result["passes"] else 0)
		out.set_flag(flag_name, column)
		return out
	return frame.take_rows(range(frame.row_count()).filter(func(r: int) -> bool: return results[r]["passes"]))


static func percent_passes(numerator: int, denominator: int, minimum_percent: float) -> bool:
	if denominator < 1:
		return minimum_percent <= 0
	return float(numerator) * 100.0 >= minimum_percent * float(denominator)


static func af_in_any_range(af: float, ranges: Array) -> bool:
	if is_nan(af) or af < 0.0 or af > 1.0:
		return false
	for r in ranges:
		var lo = r.get("minimum")
		var hi = r.get("maximum")
		if lo != null and (af < lo or (af == lo and not r.get("include_minimum", true))):
			continue
		if hi != null and (af > hi or (af == hi and not r.get("include_maximum", true))):
			continue
		return true
	return false


## All filters must match this one call's metadata.
static func metadata_matches(meta: Dictionary, filters: Array, alt: String) -> bool:
	for item in filters:
		var field: String = item["field"]
		var value: String = str(item["value"])
		match field:
			"reference_contains_alt":
				var context := (str(meta.get("REFERENCE_6_BEFORE", "")) + str(meta.get("REFERENCE_6_AFTER", ""))).to_upper()
				if (value == "contains") != (alt.to_upper() in context):
					return false
			"reference_repeat":
				var before := _single_base_repeat(meta.get("REFERENCE_6_BEFORE", ""))
				var after := _single_base_repeat(meta.get("REFERENCE_6_AFTER", ""))
				var ok: bool = {
					"before": before and not after, "after": after and not before,
					"one": before != after, "both": before and after,
					"none": not before and not after, "either": before or after,
				}[value]
				if not ok:
					return false
			_:
				var key: String = METADATA_KEYS.get(field, "")
				if key.is_empty():
					continue  # e.g. reference_context: accepted, never applied (as in Python)
				var observed := str(meta.get(key, "")).strip_edges()
				if field in COUNT_FIELDS:
					var parts := value.split("|")
					if not observed.is_valid_int() or not _count_matches(observed.to_int(), parts[0], parts[1].to_int()):
						return false
				elif observed.to_upper() != value.to_upper():
					return false
	return true


static func _position_matches(position: Dictionary, pos: int) -> bool:
	if position.get("exact") != null:
		return pos == int(position["exact"])
	if position.get("minimum") != null and pos < int(position["minimum"]):
		return false
	if position.get("maximum") != null and pos > int(position["maximum"]):
		return false
	return true


static func _single_base_repeat(sequence) -> bool:
	var text := str(sequence).to_upper()
	for base in ["AA", "CC", "GG", "TT"]:
		if base in text:
			return true
	return false


static func _count_matches(value: int, operator: String, threshold: int) -> bool:
	match operator:
		"gt": return value > threshold
		"gte": return value >= threshold
		"lt": return value < threshold
		"lte": return value <= threshold
		"eq": return value == threshold
		"neq": return value != threshold
	return false
