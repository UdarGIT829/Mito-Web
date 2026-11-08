class_name AfClasses
extends RefCounted
## Germline / somatic classification of each call from its allele fraction,
## using fixed quarters of the AF scale (rule "af-quarters-v1"). Mitochondrial
## heteroplasmy is bimodal (one hump near 1, one near 0), so the middle band is
## nearly empty and calls there are flagged rather than labelled.
##
##   AF >= germline_min (0.75)                 germline
##   somatic_below (0.50) <= AF < germline_min indeterminate
##   AF < somatic_below                        somatic
##   no AF                                     unknown (never guessed)

const RULE := "af-quarters-v1"
const GERMLINE := "germline"
const INDETERMINATE := "indeterminate"
const SOMATIC := "somatic"
const UNKNOWN := "unknown"
const CLASSES := [GERMLINE, INDETERMINATE, SOMATIC, UNKNOWN]
const SHORT := {GERMLINE: "G", INDETERMINATE: "I", SOMATIC: "S", UNKNOWN: "?"}
const DEFAULT_THRESHOLDS := {"germline_min": 0.75, "somatic_below": 0.5}


static func classify(af: float, thresholds: Dictionary = DEFAULT_THRESHOLDS) -> String:
	if is_nan(af) or af < 0.0 or af > 1.0:
		return UNKNOWN
	if af >= thresholds.get("germline_min", 0.75):
		return GERMLINE
	if af >= thresholds.get("somatic_below", 0.5):
		return INDETERMINATE
	return SOMATIC


static func validate(thresholds: Dictionary) -> String:
	var g = thresholds.get("germline_min")
	var s = thresholds.get("somatic_below")
	if typeof(g) not in [TYPE_FLOAT, TYPE_INT] or typeof(s) not in [TYPE_FLOAT, TYPE_INT]:
		return "Set both AF cut-offs."
	if not (0.0 < s and s <= g and g <= 1.0):
		return "AF cut-offs must satisfy 0 < somatic-below <= germline-from <= 1."
	return ""


## mode "keep": carriers are reduced to calls whose class is in `keep`;
## mutations left with no carriers are removed. The sample universe (the
## presence denominator) is unchanged.
## mode "label": every call is kept; tables show each call's class.
## Both modes record the thresholds on the frame so tables can show classes.
static func apply(frame: MutationFrame, thresholds: Dictionary, keep: Array, mode: String = "keep") -> MutationFrame:
	var af: PackedFloat64Array = frame.store.calls["af"]
	var out := frame.take_rows(range(frame.row_count()))
	out.af_thresholds = thresholds.duplicate()
	if mode != "keep":
		return out
	var rows: Array = []
	var carriers: Array = []
	for r in out.row_count():
		var kept := {}
		for sample_id in out.carriers[r]:
			var call: int = out.carriers[r][sample_id]
			if classify(af[call], thresholds) in keep:
				kept[sample_id] = call
		if not kept.is_empty():
			rows.append(r)
			carriers.append(kept)
	var result := out.take_rows(rows)
	result.carriers = carriers
	return result


## Count of carrier calls per class for one row: {"germline": n, ...}.
static func row_counts(frame: MutationFrame, r: int) -> Dictionary:
	var af: PackedFloat64Array = frame.store.calls["af"]
	var counts := {}
	for sample_id in frame.carriers[r]:
		var c := classify(af[frame.carriers[r][sample_id]], frame.af_thresholds)
		counts[c] = counts.get(c, 0) + 1
	return counts


static func describe_counts(counts: Dictionary) -> String:
	var parts: Array = []
	for c in CLASSES:
		if counts.get(c, 0):
			parts.append("%d %s" % [counts[c], c])
	return ", ".join(parts)
