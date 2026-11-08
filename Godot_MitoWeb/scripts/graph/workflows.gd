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
			"description": "Finds the variants that define a haplogroup from the data alone, then checks them against PhyloTree. Gained: carried by nearly every sample of the haplogroup and almost no one else. Absent: carried by nearly everyone else but almost no sample of the haplogroup (how H, the rCRS haplogroup, is defined: by lacking 2706G and 7028T). Each set gets a PhyloTree marker column. Markers that recur in other lineages are excluded, and a haplogroup whose descendants are counted as everyone else (e.g. L3, ancestral to M and N) finds few.",
			"needs": "Haplogroup tags: Samples tab, Assign haplogroups. Works best with many samples (e.g. 1000 Genomes).",
			"params": [
				{"id": "haplogroup", "label": "Major haplogroup", "type": "haplogroup", "default": "H"},
				{"id": "common", "label": "Common means carried by at least (%)", "type": "number", "default": 95, "min": 50, "max": 100},
				{"id": "rare", "label": "Rare means carried by under (%)", "type": "number", "default": 5, "min": 0, "max": 50},
			],
			"build": func(p: Dictionary) -> Dictionary: return haplogroup_markers_graph(str(p["haplogroup"]), float(p["common"]), float(p["rare"])),
		},
	]


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


## The haplogroup marker recreation graph for major haplogroup X:
##   In X / Everyone else (tag "macro haplogroup=X")
##   Gained in X  = common in X  minus  seen (>= rare%) in everyone else
##   Absent in X  = common in everyone else  minus  seen (>= rare%) in X
## Both annotated (gnomAD) and flagged against PhyloTree markers.
static func haplogroup_markers_graph(haplogroup: String, common: float, rare: float) -> Dictionary:
	var tag := "macro haplogroup=" + haplogroup
	var g := GraphModel.new_graph()
	var all := GraphModel.add_node(g, "data_provider", Vector2(0, 40))
	all["title"] = "All samples"
	var others := GraphModel.add_node(g, "data_provider", Vector2(0, 420), {"criteria": {"tags_none": [tag]}})
	others["title"] = "Everyone else"
	var in_x := GraphModel.add_node(g, "tag_filter", Vector2(330, 40), {"tags": [tag], "restrict_samples": true})
	in_x["title"] = "In %s" % haplogroup
	var common_x := _presence(g, Vector2(660, 0), common, "Common in %s" % haplogroup)
	var seen_x := _presence(g, Vector2(660, 300), rare, "Seen in %s" % haplogroup)
	var common_o := _presence(g, Vector2(660, 560), common, "Common in everyone else")
	var seen_o := _presence(g, Vector2(660, 860), rare, "Seen in everyone else")
	var gained := GraphModel.add_node(g, "comparison", Vector2(1000, 120))
	gained["title"] = "Gained in %s" % haplogroup
	var absent := GraphModel.add_node(g, "comparison", Vector2(1000, 640))
	absent["title"] = "Absent in %s" % haplogroup
	GraphModel.connect_nodes(g, all["id"], 0, in_x["id"], 0)
	for pair in [[in_x, common_x], [in_x, seen_x], [others, common_o], [others, seen_o]]:
		GraphModel.connect_nodes(g, pair[0]["id"], 0, pair[1]["id"], 0)
	GraphModel.connect_nodes(g, common_x["id"], 0, gained["id"], 0)
	GraphModel.connect_nodes(g, seen_o["id"], 0, gained["id"], 1)
	GraphModel.connect_nodes(g, common_o["id"], 0, absent["id"], 0)
	GraphModel.connect_nodes(g, seen_x["id"], 0, absent["id"], 1)
	var y := 60
	for pair in [[gained, "%s: gained markers" % haplogroup], [absent, "%s: absent markers" % haplogroup]]:
		var enrich := GraphModel.add_node(g, "enrich", Vector2(1340, y), {"sources": ["gnomad", "clinvar"]})
		var marker := GraphModel.add_node(g, "annotation_filter", Vector2(1680, y), {
			"rules": {"haplogroup_marker": "yes", "action": "keep"}, "mode": "flag", "flag_name": "PhyloTree marker"})
		marker["title"] = "Check against PhyloTree"
		var table := GraphModel.add_node(g, "table", Vector2(2020, y), {"name": pair[1]})
		GraphModel.connect_nodes(g, pair[0]["id"], 0, enrich["id"], 0)
		GraphModel.connect_nodes(g, enrich["id"], 0, marker["id"], 0)
		GraphModel.connect_nodes(g, marker["id"], 0, table["id"], 0)
		y += 520
	return g


static func _presence(g: Dictionary, position: Vector2, percent: float, title: String) -> Dictionary:
	var n := GraphModel.add_node(g, "logical_filter", position, {"filters": {"minimum_presence_percent": percent}})
	n["title"] = title
	return n


## Major haplogroups present in the data, most common first: [[name, count]].
static func haplogroups_in(store: VariantStore) -> Array:
	var counts := {}
	for sample in store.samples:
		for tag in sample["tags"]:
			if tag.begins_with("macro haplogroup="):
				var name: String = tag.trim_prefix("macro haplogroup=")
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
	return JSON.stringify({"format": FILE_FORMAT, "version": FILE_VERSION, "title": workflow["title"],
		"description": workflow.get("description", ""), "graph": graph}, "  ")


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
	return {"ok": true, "error": "", "workflow": workflow}
