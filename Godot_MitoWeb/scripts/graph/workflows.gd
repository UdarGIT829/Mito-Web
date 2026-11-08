class_name Workflows
extends RefCounted
## Preset analyses: named, parameterised recipes that build an analysis
## graph. Built-in workflows live here; users can also save the current
## graph as a workflow (stored in the workspace, exportable as a file).
##
## A workflow: {"id", "title", "description", "needs", "params": [param],
## "build": Callable(params) -> graph}  (built-ins)
## or {"id", "title", "description", "graph"}  (saved).
## A param: {"id", "label", "type": "tag" | "number" | "haplogroup",
## "default", "min", "max", "help"}.

const FILE_FORMAT := "mito-web-workflow"
const FILE_VERSION := 1


static func builtin() -> Array:
	return [
		{
			"id": "treated_vs_untreated",
			"title": "Treated vs untreated",
			"description": "Mutations only in group A, only in group B, and in both; the exclusive sets are annotated, and density statistics compare them with all somatic calls (where Treated vs Untreated can be tested by subject).",
			"needs": "Samples tagged with the two groups (Files tab, or Import tags from a table).",
			"params": [
				{"id": "tag_a", "label": "Group A tag", "type": "tag", "default": "Treated"},
				{"id": "tag_b", "label": "Group B tag", "type": "tag", "default": "Untreated"},
			],
			"build": func(p: Dictionary) -> Dictionary: return GraphModel.example_graph(p["tag_a"], p["tag_b"]),
		},
		{
			"id": "haplogroup_markers",
			"title": "Recreate haplogroup markers",
			"description": "Finds the variants that define a haplogroup from the data alone, then checks them against PhyloTree. Gained: carried by nearly every sample of the haplogroup and almost no one else. Absent: carried by nearly everyone else but almost no sample of the haplogroup (how H, the rCRS haplogroup, is defined: by lacking 2706G and 7028T). Each set gets a PhyloTree marker column. By default the haplogroup is its whole clade (R includes H, J, T, U, K...), so ancestral haplogroups like R, N, M and L3 work too. Markers that recur in other lineages are excluded.",
			"needs": "Haplogroup tags: Samples tab, Assign haplogroups. Works best with many samples (e.g. 1000 Genomes).",
			"params": [
				{"id": "haplogroup", "label": "Haplogroup", "type": "haplogroup", "default": "H"},
				{"id": "scope", "label": "Samples in the group", "type": "choice", "default": "clade",
					"options": [["clade", "Whole clade (with descendant haplogroups)"], ["macro", "This major haplogroup only"]]},
				{"id": "common", "label": "Common means carried by at least (%)", "type": "number", "default": 90, "min": 50, "max": 100},
				{"id": "rare", "label": "Rare means carried by under (%)", "type": "number", "default": 5, "min": 0, "max": 50},
			],
			"build": func(p: Dictionary) -> Dictionary: return haplogroup_markers_graph(str(p["haplogroup"]), float(p["common"]), float(p["rare"]), str(p["scope"])),
		},
		_population_entry(),
	]


static func _population_entry() -> Dictionary:
	return {
		"id": "cohort_vs_population",
		"title": "Cohort vs 1000 Genomes",
		"description": "Compares the lineage make-up of a patient cohort with 1000 Genomes reference populations: share of African, Amerindian, Caucasian and Asian lineages in the cohort, in Mexican Ancestry from Los Angeles (MXL), in the other admixed American populations and in the rest of 1000 Genomes, with Fisher's exact tests. Shows whether the cohort's lineages differ from the population it was drawn from.",
		"needs": "The 1000 Genomes bundle (tags 1000G and population codes) and a patient cohort, all with haplogroups. Lineage rules are set to the four-lineage preset if none exist.",
		"needs_haplogroups": true,
		"lineage_preset": "four_lineages",
		"params": [
			{"id": "tag", "label": "Cohort tag", "type": "tag", "default": "patients"},
			{"id": "population", "label": "Local population tag", "type": "tag", "default": "MXL"},
		],
		"build": func(p: Dictionary) -> Dictionary: return population_graph(str(p["tag"]), str(p["population"])),
	}


## A cohort's lineage shares against 1000 Genomes populations.
static func population_graph(tag: String, population: String) -> Dictionary:
	var g := GraphModel.new_graph()
	# Non-overlapping groups, so the across-group test is valid.
	var other_amr: Array = ["MXL", "PUR", "CLM", "PEL"].filter(func(code): return code != population)
	var specs := [
		[tag, {"tags_all": [tag]}, "The patient cohort (samples tagged %s)." % tag],
		["1000G %s" % population, {"tags_all": ["1000G", population]}, "1000 Genomes samples of population %s (MXL: Mexican Ancestry from Los Angeles)." % population],
		["1000G other AMR", {"tags_any": other_amr}, "The other admixed American populations in 1000 Genomes (%s)." % ", ".join(other_amr)],
		["1000G rest", {"tags_any": ["AFR", "EAS", "EUR", "SAS"]}, "All other 1000 Genomes samples (African, East Asian, European, South Asian)."],
	]
	var groups: Array = []
	var y := 0.0
	for spec in specs:
		var n := GraphModel.add_node(g, "data_provider", Vector2(0, y), {"criteria": spec[1]})
		GraphModel.describe(n, spec[0], spec[2])
		groups.append(n)
		y += 300.0
	var measures: Array = []
	for lineage in ["African", "Amerindian", "Caucasian", "Asian"]:
		measures.append({"label": "%s lineage, n (%%)" % lineage, "attribute": "[Lineage]", "summary": "percent",
			"condition": {"op": "=", "value": lineage}, "pairwise": true})
	var compare := GraphModel.add_node(g, "patient_compare", Vector2(560, 0), {"name": "Lineage shares", "measures": measures})
	GraphModel.describe(compare, "Lineage shares", "Share of each lineage per group, Fisher's exact test across groups and for each pair. The cohort vs %s column asks whether the patients' maternal lineages differ from the local population." % population)
	for group in groups:
		GraphModel.connect_nodes(g, group["id"], 0, compare["id"], 0)
	return g


static func find_builtin(id: String) -> Dictionary:
	for w in builtin():
		if w["id"] == id:
			return w
	return {}


## Parameter values with defaults filled in.
static func resolve_params(workflow: Dictionary, values: Dictionary) -> Dictionary:
	var out := {}
	for p in workflow.get("params", []):
		out[p["id"]] = values.get(p["id"], p["default"])
	return out


## Builds a workflow's graph. Built-ins use their recipe; saved workflows
## return a copy of their stored graph.
static func build(workflow: Dictionary, values: Dictionary = {}) -> Dictionary:
	if workflow.has("build"):
		return workflow["build"].call(resolve_params(workflow, values))
	return workflow["graph"].duplicate(true)


## The haplogroup marker recreation graph for haplogroup X:
##   In X / Everyone else (tag "haplogroup clade=X" for the whole clade,
##   or "macro haplogroup=X" for samples whose nearest major group is X)
##   Gained in X  = common in X  minus  seen (>= rare%) in everyone else
##   Absent in X  = common in everyone else  minus  seen (>= rare%) in X
## Both annotated (gnomAD) and flagged against PhyloTree markers.
static func haplogroup_markers_graph(haplogroup: String, common: float, rare: float, scope: String = "clade") -> Dictionary:
	var clade := scope != "macro"
	var tag := (HaplogroupClassifier.CLADE_PREFIX if clade else HaplogroupClassifier.MACRO_PREFIX) + haplogroup
	var x := haplogroup
	var members := ("in clade %s (%s and every haplogroup descended from it)" % [x, x]) if clade else \
		("whose nearest major haplogroup is %s (descendant major haplogroups not included)" % x)
	var c := _pct(common)
	var r := _pct(rare)
	var g := GraphModel.new_graph()
	var all := GraphModel.describe(GraphModel.add_node(g, "data_provider", Vector2(0, 40)), "All samples",
		"Every imported sample, before splitting by haplogroup.")
	var others := GraphModel.describe(GraphModel.add_node(g, "data_provider", Vector2(0, 640), {"criteria": {"tags_none": [tag]}}), "Everyone else",
		"Samples without the \"%s\" tag (not %s): the comparison group." % [tag, ("in clade " + x) if clade else x])
	var in_x := GraphModel.describe(GraphModel.add_node(g, "tag_filter", Vector2(330, 40), {"tags": [tag], "restrict_samples": true}), "In %s" % x,
		"Samples %s, tagged \"%s\" by the haplogroup classifier. Only these samples count toward presence below." % [members, tag])
	var common_x := _presence(g, Vector2(660, 0), common, "Common in %s" % x,
		"Mutations carried by at least %s%% of %s samples: candidate markers %s has gained." % [c, x, x])
	var seen_x := _presence(g, Vector2(660, 320), rare, "Seen in %s" % x,
		"Mutations carried by at least %s%% of %s samples. Anything here is not absent from %s." % [r, x, x])
	var common_o := _presence(g, Vector2(660, 640), common, "Common in everyone else",
		"Mutations carried by at least %s%% of the other samples: candidates for what %s lacks." % [c, x])
	var seen_o := _presence(g, Vector2(660, 960), rare, "Seen in everyone else",
		"Mutations carried by at least %s%% of the other samples. Anything here is not specific to %s." % [r, x])
	var gained := GraphModel.describe(GraphModel.add_node(g, "comparison", Vector2(1000, 120)), "Gained in %s" % x,
		"Common in %s minus Seen in everyone else: carried by >= %s%% of %s and < %s%% of everyone else." % [x, c, x, r])
	var absent := GraphModel.describe(GraphModel.add_node(g, "comparison", Vector2(1000, 760)), "Absent in %s" % x,
		"Common in everyone else minus Seen in %s: carried by >= %s%% of the others and < %s%% of %s.%s" % [x, c, r, x, _absent_note(x)])
	GraphModel.connect_nodes(g, all["id"], 0, in_x["id"], 0)
	for pair in [[in_x, common_x], [in_x, seen_x], [others, common_o], [others, seen_o]]:
		GraphModel.connect_nodes(g, pair[0]["id"], 0, pair[1]["id"], 0)
	GraphModel.connect_nodes(g, common_x["id"], 0, gained["id"], 0)
	GraphModel.connect_nodes(g, seen_o["id"], 0, gained["id"], 1)
	GraphModel.connect_nodes(g, common_o["id"], 0, absent["id"], 0)
	GraphModel.connect_nodes(g, seen_x["id"], 0, absent["id"], 1)
	var y := 60
	for pair in [[gained, "%s: gained markers" % x, "gained"], [absent, "%s: absent markers" % x, "absent"]]:
		var enrich := GraphModel.describe(GraphModel.add_node(g, "enrich", Vector2(1340, y), {"sources": ["gnomad", "clinvar"]}),
			"Annotate %s markers" % pair[2],
			"Adds gnomAD frequencies, including the haplogroups each variant is common in, and ClinVar significance.")
		var marker := GraphModel.describe(GraphModel.add_node(g, "annotation_filter", Vector2(1680, y), {
			"rules": {"haplogroup_marker": "yes", "action": "keep"}, "mode": "flag", "flag_name": "PhyloTree marker"}),
			"Check against PhyloTree",
			"Keeps every mutation and adds a \"PhyloTree marker\" column: yes if gnomAD lists it as haplogroup-defining.")
		var table := GraphModel.describe(GraphModel.add_node(g, "table", Vector2(2020, y), {"name": pair[1]}), "Table: %s markers" % pair[2],
			"The %s markers found from the data, with the PhyloTree check. Open on the Tables tab." % pair[2])
		GraphModel.connect_nodes(g, pair[0]["id"], 0, enrich["id"], 0)
		GraphModel.connect_nodes(g, enrich["id"], 0, marker["id"], 0)
		GraphModel.connect_nodes(g, marker["id"], 0, table["id"], 0)
		y += 700
	return g


static func _presence(g: Dictionary, position: Vector2, percent: float, title: String, description: String) -> Dictionary:
	var n := GraphModel.add_node(g, "logical_filter", position, {"filters": {"minimum_presence_percent": percent}})
	return GraphModel.describe(n, title, description)


## Haplogroups the reference (rCRS, in H2a2a1) belongs to.
const RCRS_CLADES := ["H", "HV", "R", "N", "L3"]


## Why "absent" markers matter for haplogroups that contain the reference.
static func _absent_note(haplogroup: String) -> String:
	if haplogroup == "H":
		return " This is how H shows up: rCRS is itself H, so H is defined by lacking 2706G and 7028T, which everyone else carries."
	if haplogroup in RCRS_CLADES:
		return " rCRS is itself in %s, so %s's own defining variants show up here: everyone else carries the older base." % [haplogroup, haplogroup]
	return ""


## 95.0 -> "95", 2.5 -> "2.5".
static func _pct(value: float) -> String:
	return str(snappedf(value, 0.1)).trim_suffix(".0")


## Haplogroups present in the data, most common first: [[name, count]].
## Counts whole clades by default, or samples per major haplogroup with
## prefix MACRO_PREFIX.
static func haplogroups_in(store: VariantStore, prefix: String = HaplogroupClassifier.CLADE_PREFIX) -> Array:
	var counts := {}
	for sample in store.samples:
		for tag in sample["tags"]:
			if tag.begins_with(prefix):
				var name: String = tag.trim_prefix(prefix)
				counts[name] = counts.get(name, 0) + 1
	var names := counts.keys()
	names.sort_custom(func(a, b): return counts[a] > counts[b] or (counts[a] == counts[b] and a < b))
	return names.map(func(n): return [n, counts[n]])


# --- Saved workflows and files -------------------------------------------------

## A saved workflow from the current graph.
static func from_graph(graph: Dictionary, title: String, description: String) -> Dictionary:
	return {"id": "saved-%d" % Time.get_unix_time_from_system(), "title": title, "description": description,
		"graph": graph.duplicate(true)}


## JSON text for sharing a workflow (positions as [x, y]).
static func to_file_text(workflow: Dictionary) -> String:
	var graph: Dictionary = workflow["graph"].duplicate(true)
	for n in graph["nodes"]:
		n["position"] = [n["position"].x, n["position"].y]
	var data := {"format": FILE_FORMAT, "version": FILE_VERSION, "title": workflow["title"],
		"description": workflow.get("description", ""), "graph": graph}
	if workflow.get("lineage_rules"):
		data["lineage_rules"] = workflow["lineage_rules"]  # applied on insert when a project has none
	return JSON.stringify(data, "  ")


## Parses a workflow file. Returns {"ok", "error", "workflow"}.
static func from_file_text(text: String) -> Dictionary:
	var json := JSON.new()
	var parsed = json.data if json.parse(text) == OK else null
	if typeof(parsed) != TYPE_DICTIONARY or parsed.get("format") != FILE_FORMAT:
		return {"ok": false, "error": "This is not a Mito Web workflow file."}
	if int(parsed.get("version", 0)) > FILE_VERSION:
		return {"ok": false, "error": "This workflow was saved by a newer version of Mito Web."}
	var graph = parsed.get("graph")
	var problem := GraphModel.validate(graph)
	if problem:
		return {"ok": false, "error": problem}
	graph["next_id"] = int(graph.get("next_id", 1))
	for n in graph["nodes"]:
		var p = n.get("position", [0, 0])
		n["position"] = Vector2(float(p[0]), float(p[1])) if typeof(p) == TYPE_ARRAY and p.size() == 2 else Vector2.ZERO
	for e in graph["edges"]:
		e["from_port"] = int(e["from_port"])
		e["to_port"] = int(e["to_port"])
	var workflow := from_graph(graph, str(parsed.get("title", "Imported workflow")), str(parsed.get("description", "")))
	var rules = parsed.get("lineage_rules", [])
	if typeof(rules) == TYPE_ARRAY and not rules.is_empty():
		var clean: Array = rules.filter(func(r): return typeof(r) == TYPE_ARRAY and r.size() == 2).map(func(r): return [str(r[0]), str(r[1])])
		if Lineages.validate(clean).is_empty():
			workflow["lineage_rules"] = clean
	return {"ok": true, "error": "", "workflow": workflow}
