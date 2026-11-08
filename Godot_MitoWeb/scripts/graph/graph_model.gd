class_name GraphModel
extends RefCounted
## The analysis graph as plain data, stored in the workspace:
##   {"next_id": int, "nodes": [{id, type, title, position: Vector2, config}],
##    "edges": [{from, from_port, to, to_port}]}
## Every port carries a MutationFrame. Node types are described in TYPES.

const TYPES := {
	"data_provider": {
		"title": "Data Provider",
		"help": "Selects samples and outputs their mutations.",
		"inputs": [],
		"outputs": 1,
		"config": {"criteria": {}, "pass_only": false},
	},
	"tag_filter": {
		"title": "Filter on Tag",
		"help": "Keeps mutations carried by samples with the chosen tags.",
		"inputs": [{"name": "Mutations", "multi": false}],
		"outputs": 1,
		"config": {"tags": [], "match": "any", "restrict_samples": false},
	},
	"comparison": {
		"title": "Logical Comparison",
		"help": "Combines Include inputs, then removes mutations found in any Exclude input.",
		"inputs": [{"name": "Include", "multi": true}, {"name": "Exclude", "multi": true}],
		"outputs": 1,
		"config": {"combine": "intersect"},
	},
	"logical_filter": {
		"title": "Logical Filter",
		"help": "Filters mutations, or flags them in a column instead of removing them.",
		"inputs": [{"name": "Mutations", "multi": false}],
		"outputs": 1,
		"config": {"filters": {}, "mode": "drop", "flag_name": "passes"},
	},
	"germline_somatic": {
		"title": "Germline / Somatic",
		"help": "Classifies each sample's call by AF in fixed quarters: >= 0.75 germline, 0.50-0.75 indeterminate, < 0.50 somatic, no AF unknown. Keeps only the calls of the chosen classes, or labels them.",
		"inputs": [{"name": "Mutations", "multi": false}],
		"outputs": 1,
		"config": {"thresholds": {"germline_min": 0.75, "somatic_below": 0.5}, "keep": ["somatic"], "mode": "keep"},
	},
	"annotation_filter": {
		"title": "Annotation Filter",
		"help": "Keeps or removes mutations by ClinVar category and review, gnomAD frequency, haplogroup-marker status, consequence or tRNA prediction.",
		"inputs": [{"name": "Mutations", "multi": false}],
		"outputs": 1,
		"config": {"rules": {"clinvar_categories": ["benign"], "action": "remove"}, "mode": "drop", "flag_name": "passes"},
	},
	"enrich": {
		"title": "Enrich with Annotations",
		"help": "Adds public annotations: gnomAD mtDNA frequencies, haplogroup markers and consequences, ClinVar classifications, and Ensembl VEP lookups you run from the Tables tab.",
		"inputs": [{"name": "Mutations", "multi": false}],
		"outputs": 1,
		"config": {"sources": ["gnomad", "clinvar", "ensembl"]},
	},
	"statistics": {
		"title": "Statistics",
		"help": "Mutation density along the genome for each connected group: sliding windows (size, step, range, wrap) or genes and regions. View it on the Statistics tab.",
		"inputs": [{"name": "Groups", "multi": true}],
		"outputs": 0,
		"config": {
			"window": {"mode": "sliding", "size": 500, "step": 500, "range_start": 1, "range_end": 16569, "wrap": false},
			"measure": "calls", "per_kb": true, "per_sample": true, "name": "",
		},
	},
	"patient_compare": {
		"title": "Compare Patients",
		"help": "Compares the connected groups patient by patient (a subject's samples count once) on values imported per sample, such as clinical scores: medians with Kruskal-Wallis and Wilcoxon, means with ANOVA, percentages with Fisher's exact test, or Spearman correlations. View it on the Patients tab.",
		"inputs": [{"name": "Groups", "multi": true}],
		"outputs": 0,
		"config": {"name": "", "measures": []},
	},
	"table": {
		"title": "Table",
		"help": "Shows the mutations arriving at this point.",
		"inputs": [{"name": "Mutations", "multi": false}],
		"outputs": 0,
		"config": {"name": ""},
	},
}


static func new_graph() -> Dictionary:
	return {"next_id": 1, "nodes": [], "edges": []}


static func add_node(graph: Dictionary, type: String, position: Vector2 = Vector2.ZERO, config: Dictionary = {}) -> Dictionary:
	var spec: Dictionary = TYPES[type]
	var merged: Dictionary = spec["config"].duplicate(true)
	merged.merge(config.duplicate(true), true)
	var node := {
		"id": "n%d" % graph["next_id"],
		"type": type,
		"title": spec["title"],
		"description": "",
		"position": position,
		"config": merged,
	}
	graph["next_id"] += 1
	graph["nodes"].append(node)
	return node


static func node(graph: Dictionary, node_id: String) -> Dictionary:
	for n in graph["nodes"]:
		if n["id"] == node_id:
			return n
	return {}


static func remove_node(graph: Dictionary, node_id: String) -> void:
	graph["nodes"] = graph["nodes"].filter(func(n): return n["id"] != node_id)
	graph["edges"] = graph["edges"].filter(func(e): return e["from"] != node_id and e["to"] != node_id)


## Adds an edge. Single-input ports replace their existing edge. Returns ""
## or a reason the connection is not allowed.
static func connect_nodes(graph: Dictionary, from_id: String, from_port: int, to_id: String, to_port: int) -> String:
	var from_node := node(graph, from_id)
	var to_node := node(graph, to_id)
	if from_node.is_empty() or to_node.is_empty():
		return "Unknown node."
	if from_id == to_id:
		return "A node cannot feed itself."
	var inputs: Array = TYPES[to_node["type"]]["inputs"]
	if to_port < 0 or to_port >= inputs.size() or from_port < 0 or from_port >= TYPES[from_node["type"]]["outputs"]:
		return "Unknown port."
	if is_upstream(graph, to_id, from_id):
		return "That connection would create a loop."
	for e in graph["edges"]:
		if e["from"] == from_id and e["from_port"] == from_port and e["to"] == to_id and e["to_port"] == to_port:
			return ""
	if not inputs[to_port]["multi"]:
		graph["edges"] = graph["edges"].filter(func(e): return not (e["to"] == to_id and e["to_port"] == to_port))
	graph["edges"].append({"from": from_id, "from_port": from_port, "to": to_id, "to_port": to_port})
	return ""


static func disconnect_nodes(graph: Dictionary, from_id: String, from_port: int, to_id: String, to_port: int) -> void:
	graph["edges"] = graph["edges"].filter(func(e): return not (
		e["from"] == from_id and e["from_port"] == from_port and e["to"] == to_id and e["to_port"] == to_port))


## Node ids feeding a node's input port, in connection order.
static func inputs_of(graph: Dictionary, node_id: String, port: int) -> Array:
	return graph["edges"].filter(func(e): return e["to"] == node_id and e["to_port"] == port).map(func(e): return e["from"])


## True if `candidate` is `node_id` or feeds it, directly or indirectly.
static func is_upstream(graph: Dictionary, candidate: String, node_id: String) -> bool:
	if candidate == node_id:
		return true
	var stack := [node_id]
	var seen := {}
	while stack:
		var current: String = stack.pop_back()
		for e in graph["edges"]:
			if e["to"] == current and not seen.has(e["from"]):
				if e["from"] == candidate:
					return true
				seen[e["from"]] = true
				stack.append(e["from"])
	return false


## Node ids in dependency order (inputs before the nodes that use them).
static func topological_order(graph: Dictionary) -> Array:
	var order: Array = []
	var visited := {}
	for n in graph["nodes"]:
		_visit(graph, n["id"], visited, order)
	return order


static func _visit(graph: Dictionary, node_id: String, visited: Dictionary, order: Array) -> void:
	if visited.has(node_id):
		return
	visited[node_id] = true
	for e in graph["edges"]:
		if e["to"] == node_id:
			_visit(graph, e["from"], visited, order)
	order.append(node_id)


## Sets a node's name and description; returns the node.
static func describe(n: Dictionary, title: String, description: String) -> Dictionary:
	n["title"] = title
	n["description"] = description
	return n


## The name the user gave a node, or "" if it still has its type's name.
static func custom_title(n: Dictionary) -> String:
	var title := str(n.get("title", "")).strip_edges()
	return "" if title.is_empty() or title == TYPES.get(n.get("type"), {}).get("title", "") else title


## Copies every node and edge of `source` into `target` (new ids), placed
## below whatever `target` already holds. Returns the new node ids.
static func merge(target: Dictionary, source: Dictionary) -> Array:
	var offset := Vector2.ZERO
	for n in target["nodes"]:
		offset.y = maxf(offset.y, n["position"].y + 320)
	var renamed := {}
	for n in source["nodes"]:
		var copy := add_node(target, n["type"], n["position"] + offset, n["config"])
		copy["title"] = n["title"]
		copy["description"] = str(n.get("description", ""))
		renamed[n["id"]] = copy["id"]
	for e in source["edges"]:
		connect_nodes(target, renamed[e["from"]], e["from_port"], renamed[e["to"]], e["to_port"])
	return renamed.values()


## Validates graph data loaded from a file. Returns "" or a message.
static func validate(graph) -> String:
	if typeof(graph) != TYPE_DICTIONARY or typeof(graph.get("nodes")) != TYPE_ARRAY or typeof(graph.get("edges")) != TYPE_ARRAY:
		return "The analysis graph is malformed."
	var ids := {}
	for n in graph["nodes"]:
		if typeof(n) != TYPE_DICTIONARY or not TYPES.has(n.get("type")) or typeof(n.get("id")) != TYPE_STRING \
				or typeof(n.get("config")) != TYPE_DICTIONARY:
			return "The analysis graph contains an unknown or malformed node."
		ids[n["id"]] = true
	for e in graph["edges"]:
		if typeof(e) != TYPE_DICTIONARY or not ids.has(e.get("from")) or not ids.has(e.get("to")):
			return "The analysis graph contains a connection to a missing node."
	return ""


## The worked example from the project design: treated vs untreated
## exclusive mutations (annotated), and those in both, each with a table,
## plus density statistics for the two exclusive sets and for all somatic
## calls (the group to test Treated vs Untreated in, by subject).
static func example_graph(tag_a: String = "Treated", tag_b: String = "Untreated") -> Dictionary:
	var g := new_graph()
	var source := describe(add_node(g, "data_provider", Vector2(0, 330)), "All samples",
		"Every imported sample and its calls: the starting point for both groups and for the somatic baseline.")
	# "Keep only these samples" makes each group's sample count its own
	# samples, so per-sample statistics divide by the right number.
	var treated := describe(add_node(g, "tag_filter", Vector2(340, 0), {"tags": [tag_a], "restrict_samples": true}), "%s samples" % tag_a,
		"Mutations carried by at least one %s sample. Only %s samples are kept, so counts and per-sample rates refer to this group alone." % [tag_a, tag_a])
	var both := describe(add_node(g, "tag_filter", Vector2(340, 340), {"tags": [tag_a, tag_b], "match": "all"}), "In both groups",
		"Mutations carried by at least one %s and at least one %s sample. Subtracted from each group to leave its exclusive mutations." % [tag_a, tag_b])
	var untreated := describe(add_node(g, "tag_filter", Vector2(340, 680), {"tags": [tag_b], "restrict_samples": true}), "%s samples" % tag_b,
		"Mutations carried by at least one %s sample. Only %s samples are kept, so counts and per-sample rates refer to this group alone." % [tag_b, tag_b])
	var treated_only := describe(add_node(g, "comparison", Vector2(700, 0)), "%s-exclusive" % tag_a,
		"%s samples minus In both: mutations seen in %s samples and in no %s sample." % [tag_a, tag_a, tag_b])
	var untreated_only := describe(add_node(g, "comparison", Vector2(700, 680)), "%s-exclusive" % tag_b,
		"%s samples minus In both: mutations seen in %s samples and in no %s sample." % [tag_b, tag_b, tag_a])
	var annotate := "Adds gnomAD population frequency, ClinVar significance and Ensembl consequence, to tell known polymorphisms from candidate new mutations."
	var e1 := describe(add_node(g, "enrich", Vector2(1040, 0)), "Annotate %s-exclusive" % tag_a, annotate)
	var e3 := describe(add_node(g, "enrich", Vector2(1040, 680)), "Annotate %s-exclusive" % tag_b, annotate)
	var view := "Open on the Tables tab; double-click a row to see each sample's call and AF."
	var t1 := describe(add_node(g, "table", Vector2(1420, 0), {"name": "%s-exclusive" % tag_a}), "Table: %s-exclusive" % tag_a, view)
	var t2 := describe(add_node(g, "table", Vector2(1040, 340), {"name": "In both"}), "Table: in both",
		"Mutations shared by the two groups (often germline variants of the same lineage). " + view)
	var t3 := describe(add_node(g, "table", Vector2(1420, 780), {"name": "%s-exclusive" % tag_b}), "Table: %s-exclusive" % tag_b, view)
	for f in [treated, both, untreated]:
		connect_nodes(g, source["id"], 0, f["id"], 0)
	connect_nodes(g, treated["id"], 0, treated_only["id"], 0)
	connect_nodes(g, both["id"], 0, treated_only["id"], 1)
	connect_nodes(g, untreated["id"], 0, untreated_only["id"], 0)
	connect_nodes(g, both["id"], 0, untreated_only["id"], 1)
	connect_nodes(g, treated_only["id"], 0, e1["id"], 0)
	connect_nodes(g, e1["id"], 0, t1["id"], 0)
	connect_nodes(g, both["id"], 0, t2["id"], 0)
	connect_nodes(g, untreated_only["id"], 0, e3["id"], 0)
	connect_nodes(g, e3["id"], 0, t3["id"], 0)
	var somatic := describe(add_node(g, "germline_somatic", Vector2(340, 1020)), "Somatic calls",
		"Every sample's calls with AF below 0.5 (germline is AF >= 0.75; 0.5-0.75 is indeterminate). The baseline for density, and where %s vs %s can be tested by subject on the Statistics tab." % [tag_a, tag_b])
	connect_nodes(g, source["id"], 0, somatic["id"], 0)
	var stats := describe(add_node(g, "statistics", Vector2(1420, 340), {"name": "Mutation density"}), "Mutation density",
		"Mutations per window along the genome for each connected group. The Statistics tab shows the chart, the feature matrix and permutation tests.")
	connect_nodes(g, treated_only["id"], 0, stats["id"], 0)
	connect_nodes(g, untreated_only["id"], 0, stats["id"], 0)
	connect_nodes(g, somatic["id"], 0, stats["id"], 0)
	return g
