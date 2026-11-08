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
	"enrich": {
		"title": "Enrich with Annotations",
		"help": "Adds public annotations: gnomAD mtDNA frequencies, haplogroup markers and consequences, ClinVar classifications, and Ensembl VEP lookups you run from the Tables tab.",
		"inputs": [{"name": "Mutations", "multi": false}],
		"outputs": 1,
		"config": {"sources": ["gnomad", "clinvar", "ensembl"]},
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
## exclusive mutations (annotated), and those in both, each with a table.
static func example_graph() -> Dictionary:
	var g := new_graph()
	var source := add_node(g, "data_provider", Vector2(0, 160))
	var treated := add_node(g, "tag_filter", Vector2(340, 0), {"tags": ["Treated"]})
	var both := add_node(g, "tag_filter", Vector2(340, 200), {"tags": ["Treated", "Untreated"], "match": "all"})
	var untreated := add_node(g, "tag_filter", Vector2(340, 400), {"tags": ["Untreated"]})
	var treated_only := add_node(g, "comparison", Vector2(700, 40))
	treated_only["title"] = "Treated-exclusive"
	var untreated_only := add_node(g, "comparison", Vector2(700, 380))
	untreated_only["title"] = "Untreated-exclusive"
	var e1 := add_node(g, "enrich", Vector2(1040, 0))
	var e3 := add_node(g, "enrich", Vector2(1040, 400))
	var t1 := add_node(g, "table", Vector2(1420, 0), {"name": "Treated-exclusive"})
	var t2 := add_node(g, "table", Vector2(1040, 230), {"name": "In both"})
	var t3 := add_node(g, "table", Vector2(1420, 420), {"name": "Untreated-exclusive"})
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
	return g
