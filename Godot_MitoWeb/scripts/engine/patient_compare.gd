class_name PatientCompare
extends RefCounted
## Compares groups of patients (subjects) on per-sample attributes such as
## clinical scores, producing a table shaped like a paper's results table:
## one row per measure, one column per group, then the omnibus p and the
## pairwise p-values.
##
## A patient is a subject: all of its samples count once. Its value for an
## attribute comes from its samples in the group (they should agree).
##
## Measure: {"label", "attribute", "summary", "condition", "with", "pairwise"}
##   summary "median"       median [25th, 75th]; Kruskal-Wallis, pairwise Wilcoxon rank-sum
##           "mean"         mean (SD); one-way ANOVA, pairwise Wilcoxon rank-sum
##           "percent"      n (%) meeting `condition`; Fisher's exact (r x c), pairwise Fisher
##           "correlation"  Spearman rho of attribute vs `with`, per group; overall rho and p
##   condition {"op": ">", ">=", "<", "<=", "=", "!=", "value": number or text}
##           For "percent"; "= yes" also matches 1, true and yes.
##   pairwise true/false

const SUMMARIES := ["median", "mean", "percent", "correlation"]
const SUMMARY_TITLES := {"median": "Median [25th, 75th]", "mean": "Mean (SD)", "percent": "n (%) where…", "correlation": "Spearman correlation with…"}
const TESTS := {"median": "Kruskal-Wallis", "mean": "ANOVA", "percent": "Fisher's exact", "correlation": "Spearman"}
const OPS := [">", ">=", "<", "<=", "=", "!="]

## Values computed from the data rather than imported. Mutation counts come
## from the group's own mutations (so filters upstream apply), per sample by
## AF class (Germline / Somatic thresholds of the frame, else 0.75 / 0.50);
## a patient with several samples gets their mean.
## Value: [AF class ("" = any), minimum AF]. The 5% floor drops the
## low-AF noise that amplicon data carries in large numbers.
const COMPUTED_MUTATIONS := {
	"[Mutations: all]": ["", 0.0],
	"[Mutations: somatic, AF < 0.50]": [AfClasses.SOMATIC, 0.0],
	"[Mutations: somatic, 0.05 <= AF < 0.50]": [AfClasses.SOMATIC, 0.05],
	"[Mutations: indeterminate, 0.50-0.75]": [AfClasses.INDETERMINATE, 0.0],
	"[Mutations: germline, AF >= 0.75]": [AfClasses.GERMLINE, 0.0],
}
const COMPUTED_TAGS := {
	"[Haplogroup]": "haplogroup=",
	"[Major haplogroup]": "macro haplogroup=",
	"[Lineage]": "lineage=",
}


## Every attribute name present on any sample, sorted, then the computed ones.
static func attribute_names(store: VariantStore, with_computed: bool = true) -> Array:
	var names := {}
	for sample in store.samples:
		for key in sample.get("attributes", {}):
			names[key] = true
	var out := names.keys()
	out.sort_custom(func(a, b): return str(a).naturalnocasecmp_to(str(b)) < 0)
	if with_computed:
		out.append_array(COMPUTED_MUTATIONS.keys())
		out.append_array(COMPUTED_TAGS.keys())
	return out


static func is_computed(attribute: String) -> bool:
	return COMPUTED_MUTATIONS.has(attribute) or COMPUTED_TAGS.has(attribute)


## Per-sample mutation counts in a frame: {sample_id: {computed name: n}}.
static func mutation_counts(frame: MutationFrame) -> Dictionary:
	var thresholds: Dictionary = frame.af_thresholds if not frame.af_thresholds.is_empty() else AfClasses.DEFAULT_THRESHOLDS
	var af_col: PackedFloat64Array = frame.store.calls["af"]
	var out := {}
	for id in frame.samples:
		out[id] = {}
		for key in COMPUTED_MUTATIONS:
			out[id][key] = 0
	for r in frame.row_count():
		var row: Dictionary = frame.carriers[r]
		for id in row:
			var counts: Dictionary = out[id]
			var af := af_col[row[id]]
			var cls := AfClasses.classify(af, thresholds)
			for key in COMPUTED_MUTATIONS:
				var spec: Array = COMPUTED_MUTATIONS[key]
				if (spec[0] == "" or spec[0] == cls) and (spec[1] <= 0.0 or (not is_nan(af) and af >= spec[1])):
					counts[key] += 1
	return out


## groups: [{"label": String, "frame": MutationFrame}]. Returns
## {"groups": [{"label", "n", "subjects"}], "rows": [...], "warnings": [...]}.
static func run(groups: Array, store: VariantStore, measures: Array) -> Dictionary:
	var warnings: Array = []
	var group_subjects: Array = []  # per group: {subject: [sample dicts]}
	var contexts: Array = []  # per group: {"frame", "counts" (computed lazily)}
	var where := {}  # subject -> [group labels]
	for g in groups:
		var by_subject := {}
		var frame: MutationFrame = g["frame"]
		for id in frame.samples:
			var sample := store.sample_by_id(id)
			by_subject.get_or_add(sample["subject"], []).append(sample)
		group_subjects.append(by_subject)
		contexts.append({"frame": frame, "counts": null})
		for subject in by_subject:
			where.get_or_add(subject, []).append(g["label"])
	var shared: Array = where.keys().filter(func(subject): return where[subject].size() > 1)
	if shared:
		shared.sort_custom(func(a, b): return str(a).naturalnocasecmp_to(str(b)) < 0)
		warnings.append("%d patient(s) are in more than one group, so the groups are not independent: %s%s. Set their haplogroup or tags so each is in one group." % [
			shared.size(), ", ".join(shared.slice(0, 5).map(func(x): return "%s (%s)" % [x, ", ".join(where[x])])), ", …" if shared.size() > 5 else ""])
	var out_groups: Array = []
	for i in groups.size():
		var subjects: Array = group_subjects[i].keys()
		subjects.sort_custom(func(a, b): return str(a).naturalnocasecmp_to(str(b)) < 0)
		out_groups.append({"label": groups[i]["label"], "n": subjects.size(), "subjects": subjects})
	var rows: Array = []
	for m in measures:
		rows.append(_measure_row(m, group_subjects, contexts, out_groups, warnings))
	return {"groups": out_groups, "rows": rows, "warnings": warnings}


## A subject's value for an attribute (first sample that has it), noting
## disagreements between its samples. Computed values use the group context.
static func _value(samples: Array, attribute: String, warnings: Array, ctx: Dictionary = {}) -> Variant:
	if COMPUTED_MUTATIONS.has(attribute):
		if ctx["counts"] == null:
			ctx["counts"] = mutation_counts(ctx["frame"])
		var total := 0.0
		for sample in samples:
			total += ctx["counts"].get(sample["id"], {}).get(attribute, 0)
		return total / samples.size()
	if COMPUTED_TAGS.has(attribute):
		var prefix: String = COMPUTED_TAGS[attribute]
		for sample in samples:
			for tag in sample["tags"]:
				if tag.begins_with(prefix):
					return tag.trim_prefix(prefix)
		return null
	var value = null
	for sample in samples:
		var v = sample.get("attributes", {}).get(attribute)
		if v == null:
			continue
		if value == null:
			value = v
		elif v != value:
			var message := "%s has different %s values in its samples (%s, %s); using %s." % [samples[0]["subject"], attribute, value, v, value]
			if not message in warnings:
				warnings.append(message)
	return value


static func _measure_row(m: Dictionary, group_subjects: Array, contexts: Array, groups: Array, warnings: Array) -> Dictionary:
	var summary: String = m.get("summary", "median")
	var attribute: String = m.get("attribute", "")
	var row := {"label": m.get("label", attribute), "summary": summary, "cells": [], "test": TESTS.get(summary, ""),
		"p": NAN, "pairwise": [], "missing": 0, "note": ""}
	var values: Array = []  # per group: Array of numbers (or 0/1 for percent)
	var pairs: Array = []   # per group, for correlation: [[x, y], ...]
	for i in group_subjects.size():
		var by_subject: Dictionary = group_subjects[i]
		var vals: Array = []
		var xy: Array = []
		for subject in by_subject:
			var v = _value(by_subject[subject], attribute, warnings, contexts[i])
			if summary == "correlation":
				var w = _value(by_subject[subject], m.get("with", ""), warnings, contexts[i])
				if typeof(v) == TYPE_FLOAT and typeof(w) == TYPE_FLOAT:
					xy.append([v, w])
				else:
					row["missing"] += 1
				continue
			if summary == "percent":
				if v == null:
					row["missing"] += 1
				else:
					vals.append(1 if matches(v, m.get("condition", {})) else 0)
				continue
			if typeof(v) == TYPE_FLOAT:
				vals.append(v)
			else:
				row["missing"] += 1
		values.append(vals)
		pairs.append(xy)
	match summary:
		"median":
			for vals in values:
				var q := GroupTests.median_iqr(vals)
				row["cells"].append({"text": "%s [%s, %s]" % [_num(q[0]), _num(q[1]), _num(q[2])] if vals else "–", "n": vals.size()})
			row["p"] = GroupTests.kruskal_wallis(values)["p"]
		"mean":
			for vals in values:
				var ms := GroupTests.mean_sd(vals)
				row["cells"].append({"text": "%s (%s)" % [_num(ms[0], 1), _num(ms[1], 1)] if vals else "–", "n": vals.size()})
			row["p"] = GroupTests.anova(values)["p"]
		"percent":
			var table: Array = []
			for vals in values:
				var yes: int = vals.reduce(func(s, x): return s + x, 0)
				table.append([yes, vals.size() - yes])
				row["cells"].append({"text": "%d (%s%%)" % [yes, _num(100.0 * yes / vals.size(), 0)] if vals else "–", "n": vals.size()})
			var fisher := GroupTests.fisher_exact(table)
			row["p"] = fisher["p"]
			if fisher.get("method") == "Monte Carlo":
				row["test"] = "Fisher's exact (Monte Carlo, 20,000 tables)"
		"correlation":
			var all_x: Array = []
			var all_y: Array = []
			for xy in pairs:
				var r := GroupTests.spearman(xy.map(func(p): return p[0]), xy.map(func(p): return p[1]))
				row["cells"].append({"text": ("rho %s" % _num(r["rho"], 2)) if not is_nan(r["rho"]) else "–", "n": xy.size()})
				for p in xy:
					all_x.append(p[0])
					all_y.append(p[1])
			var overall := GroupTests.spearman(all_x, all_y)
			row["p"] = overall["p"]
			row["note"] = "overall rho %s, n %d" % [_num(overall["rho"], 2), overall["n"]]
	if m.get("pairwise", false) and summary != "correlation":
		for i in values.size():
			for j in range(i + 1, values.size()):
				var p := NAN
				if summary == "percent":
					var a: Array = values[i]
					var b: Array = values[j]
					if a.size() and b.size():
						var ya: int = a.reduce(func(s, x): return s + x, 0)
						var yb: int = b.reduce(func(s, x): return s + x, 0)
						if ya + yb > 0 and ya + yb < a.size() + b.size():
							p = GroupTests.fisher_exact([[ya, a.size() - ya], [yb, b.size() - yb]])["p"]
				else:
					p = GroupTests.rank_sum(values[i], values[j])["p"]
				row["pairwise"].append({"a": groups[i]["label"], "b": groups[j]["label"], "p": p})
	return row


## Whether a value meets a condition {"op", "value"}.
static func matches(v: Variant, condition: Dictionary) -> bool:
	var op: String = condition.get("op", "=")
	var target = condition.get("value", "")
	if typeof(target) == TYPE_STRING and (target.is_valid_float() or target.is_valid_int()):
		target = target.to_float()
	if typeof(target) == TYPE_FLOAT and typeof(v) == TYPE_FLOAT:
		match op:
			">": return v > target
			">=": return v >= target
			"<": return v < target
			"<=": return v <= target
			"=": return is_equal_approx(v, target)
			"!=": return not is_equal_approx(v, target)
		return false
	var a := str(v).strip_edges().to_lower()
	var b := str(target).strip_edges().to_lower()
	var truthy := ["yes", "true", "y", "1", "1.0"]
	var equal := a == b or (b in truthy and a in truthy) or (b in ["no", "false", "n", "0"] and a in ["no", "false", "n", "0", "0.0"])
	match op:
		"=": return equal
		"!=": return not equal
	return false


static func condition_text(condition: Dictionary) -> String:
	if condition.is_empty():
		return ""
	var value = condition.get("value", "")
	return "%s %s" % [condition.get("op", "="), _num(value) if typeof(value) == TYPE_FLOAT else str(value)]


## Parses "> 6", ">= 6", "= F", "yes" into a condition.
static func parse_condition(text: String) -> Dictionary:
	var t := text.strip_edges()
	if t.is_empty():
		return {}
	for op in [">=", "<=", "!=", ">", "<", "="]:
		if t.begins_with(op):
			var rest := t.substr(op.length()).strip_edges()
			return {"op": op, "value": rest.to_float() if rest.is_valid_float() else rest}
	return {"op": "=", "value": t.to_float() if t.is_valid_float() else t}


static func format_p(p: float) -> String:
	if is_nan(p):
		return "NE"
	if p < 0.001:
		return "<0.001"
	return "%.3f" % p


static func _num(x: float, decimals: int = -1) -> String:
	if is_nan(x):
		return "–"
	if decimals < 0:
		return str(snappedf(x, 0.01)).trim_suffix(".0")
	return String.num(x, decimals) if decimals > 0 else str(int(round(x)))


## CSV text of a result (for export).
static func to_csv(result: Dictionary) -> String:
	var groups: Array = result["groups"]
	var header := ["Measure"]
	for g in groups:
		header.append("%s (n = %d)" % [g["label"], g["n"]])
	header.append_array(["Test", "p"])
	var pair_names: Array = []
	for row in result["rows"]:
		for pr in row["pairwise"]:
			var name := "%s vs %s" % [pr["a"], pr["b"]]
			if not name in pair_names:
				pair_names.append(name)
	header.append_array(pair_names)
	var lines := [",".join(header.map(_csv))]
	for row in result["rows"]:
		var cells := [row["label"]]
		cells.append_array(row["cells"].map(func(c): return c["text"]))
		cells.append_array([row["test"], format_p(row["p"])])
		for name in pair_names:
			var hit: Array = row["pairwise"].filter(func(pr): return "%s vs %s" % [pr["a"], pr["b"]] == name)
			cells.append(format_p(hit[0]["p"]) if hit else "")
		lines.append(",".join(cells.map(_csv)))
	return "\n".join(lines) + "\n"


static func _csv(value) -> String:
	var s := str(value)
	if "," in s or "\"" in s or "\n" in s:
		return "\"" + s.replace("\"", "\"\"") + "\""
	return s
