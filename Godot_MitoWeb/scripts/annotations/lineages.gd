class_name Lineages
extends RefCounted
## Groups samples into named lineages (e.g. African, Amerindian) from their
## haplogroups with rules "clade -> group". The most specific rule wins: a
## sample takes the group of the nearest haplogroup on its path to the root
## that has a rule, so "N = Asian" and "A2 = Amerindian" together put A2
## samples in Amerindian and other N samples in Asian. Results are tags
## "lineage=<group>".

const PREFIX := "lineage="

## Built-in rule sets. Rules are [clade, group].
const PRESETS := {
	"four_lineages": {
		"title": "Four maternal lineages: African, Amerindian, Caucasian, Asian",
		"description": "African = L lineages. Amerindian = the Native American founder clades of A, B, C and D (A2, B2, C1b, C1c, C1d, C4c, D1, D4h3a) and X2a. Caucasian = H, HV, J, T, V, U, K, W, X, I. Asian = everything else in M and N (e.g. B4, D4, D5).",
		"rules": [
			["L0", "African"], ["L1", "African"], ["L2", "African"], ["L3", "African"], ["L4", "African"], ["L5", "African"], ["L6", "African"],
			["M", "Asian"], ["N", "Asian"],
			["A2", "Amerindian"], ["B2", "Amerindian"], ["C1b", "Amerindian"], ["C1c", "Amerindian"], ["C1d", "Amerindian"],
			["C4c", "Amerindian"], ["D1", "Amerindian"], ["D4h3a", "Amerindian"], ["X2a", "Amerindian"],
			["H", "Caucasian"], ["HV", "Caucasian"], ["J", "Caucasian"], ["T", "Caucasian"], ["V", "Caucasian"],
			["U", "Caucasian"], ["K", "Caucasian"], ["W", "Caucasian"], ["X", "Caucasian"], ["I", "Caucasian"],
		],
	},
}


## Rules as editable text, one "clade = group" per line.
static func to_text(rules: Array) -> String:
	return "\n".join(rules.map(func(r): return "%s = %s" % r))


## Parses "clade = group" lines (also "clade, group" or tab-separated; #
## starts a comment). Returns {"rules": [[clade, group]], "errors": [String]}.
static func parse_text(text: String) -> Dictionary:
	var rules: Array = []
	var errors: Array = []
	var lines := text.split("\n")
	for i in lines.size():
		var line := lines[i].get_slice("#", 0).strip_edges()
		if line.is_empty():
			continue
		var parts: PackedStringArray
		for sep in ["=", "\t", ","]:
			if sep in line:
				parts = line.split(sep, false, 1)
				break
		if parts.size() != 2 or parts[0].strip_edges().is_empty() or parts[1].strip_edges().is_empty():
			errors.append("Line %d: write it as \"clade = group\"." % (i + 1))
			continue
		rules.append([parts[0].strip_edges(), parts[1].strip_edges()])
	errors.append_array(validate(rules))
	return {"rules": rules, "errors": errors}


## Messages for rules whose clade is not a PhyloTree 17 haplogroup or cluster.
static func validate(rules: Array) -> Array:
	HaplogroupClassifier.ensure_loaded()
	var errors: Array = []
	var seen := {}
	for r in rules:
		var clade: String = r[0]
		if HaplogroupClassifier.node_index(clade) < 0 and not HaplogroupClassifier._clusters.has(clade):
			errors.append("%s is not a PhyloTree 17 haplogroup." % clade)
		if seen.has(clade) and seen[clade] != r[1]:
			errors.append("%s is given two groups (%s and %s)." % [clade, seen[clade], r[1]])
		seen[clade] = r[1]
	return errors


## The group for a tree node, or "" when no rule covers it.
static func group_for(node: int, rules: Array) -> String:
	if node < 0:
		return ""
	var by_clade := {}
	for r in rules:
		by_clade[r[0]] = r[1]
	for n in HaplogroupClassifier.true_ancestors(node):
		var name := HaplogroupClassifier.node_name(n)
		if by_clade.has(name):
			return by_clade[name]
		var cluster: String = HaplogroupClassifier._cluster_nodes.get(name, "")
		if cluster and by_clade.has(cluster):
			return by_clade[cluster]
	return ""


## Sets "lineage=" tags on samples with a haplogroup. Returns {group: count},
## with "" counting samples no rule covers.
static func apply(samples: Array, rules: Array) -> Dictionary:
	HaplogroupClassifier.ensure_loaded()
	var counts := {}
	for sample in samples:
		var tags := PackedStringArray(Array(sample.get("tags", [])).filter(func(t): return not t.begins_with(PREFIX)))
		var hg: Dictionary = sample.get("haplogroup", {})
		if hg and not rules.is_empty():
			var group := group_for(HaplogroupClassifier.node_index(hg["name"]), rules)
			counts[group] = counts.get(group, 0) + 1
			if group:
				tags.append(PREFIX + group)
		sample["tags"] = tags
	return counts


## Groups used by rules, in first-appearance order.
static func groups(rules: Array) -> Array:
	var out: Array = []
	for r in rules:
		if not r[1] in out:
			out.append(r[1])
	return out


## Colours for lineages, in order of first appearance in the rules
## (readable on the dark theme; validated categorical palette).
const PALETTE := ["#3987e5", "#d95926", "#199e70", "#c98500", "#9d5ce0", "#d1477a", "#2aa6b8", "#8b9a46"]
const NO_GROUP_COLOR := "#7a8591"


static func color_for(group: String, rules: Array) -> Color:
	var i := groups(rules).find(group)
	return Color(PALETTE[i % PALETTE.size()]) if i >= 0 else Color(NO_GROUP_COLOR)


## The rules as a tree for display: every clade named by a rule plus every
## major haplogroup the samples fall in, each under its nearest such
## ancestor. Returns {"nodes": {name: {"parent", "children", "group",
## "own_rule", "split", "total", "here"}}, "roots": [names]}.
## total = samples in the clade; here = samples whose nearest shown clade is
## this one. split = the node's lineage differs from its parent's.
static func tree(rules: Array, samples: Array) -> Dictionary:
	HaplogroupClassifier.ensure_loaded()
	var by_clade := {}
	for r in rules:
		by_clade[r[0]] = r[1]
	var shown := {}
	for r in rules:
		shown[r[0]] = true
	var sample_nodes: Array = []
	for sample in samples:
		var hg: Dictionary = sample.get("haplogroup", {})
		if hg.is_empty():
			continue
		var n := HaplogroupClassifier.node_index(hg["name"])
		if n < 0:
			continue
		sample_nodes.append(n)
		for clade in HaplogroupClassifier.clades(n):
			shown[clade] = true
	var nodes := {}
	for name in shown:
		var rep := _representative(name)
		var parent := ""
		if rep >= 0:
			var path := HaplogroupClassifier.true_ancestors(rep)
			for k in range(1, path.size()):
				var hit := _shown_name(path[k], shown)
				if hit and hit != name:
					parent = hit
					break
		var group: String = group_for(rep, rules) if rep >= 0 else by_clade.get(name, "")
		nodes[name] = {"parent": parent, "children": [], "group": group, "own_rule": by_clade.has(name),
			"split": false, "total": 0, "here": 0, "node": rep}
	var roots: Array = []
	for name in nodes:
		var parent: String = nodes[name]["parent"]
		if parent:
			nodes[parent]["children"].append(name)
		else:
			roots.append(name)
	for name in nodes:
		var parent: String = nodes[name]["parent"]
		nodes[name]["split"] = parent != "" and nodes[parent]["group"] != nodes[name]["group"]
		nodes[name]["children"].sort_custom(_name_order)
	roots.sort_custom(_name_order)
	for n in sample_nodes:
		var counted := {}  # a cluster such as B can match several nodes on one path
		for a in HaplogroupClassifier.true_ancestors(n):
			var hit := _shown_name(a, shown)
			if hit == "" or counted.has(hit):
				continue
			if counted.is_empty():
				nodes[hit]["here"] += 1
			counted[hit] = true
			nodes[hit]["total"] += 1
	return {"nodes": nodes, "roots": roots}


static func _name_order(a, b) -> bool:
	return str(a).naturalnocasecmp_to(str(b)) < 0


## The name a tree node is shown under ("" when not shown): its own name,
## or the cluster it belongs to (B4'5 -> B).
static func _shown_name(node: int, shown: Dictionary) -> String:
	var name := HaplogroupClassifier.node_name(node)
	if shown.has(name):
		return name
	var cluster: String = HaplogroupClassifier._cluster_nodes.get(name, "")
	if cluster and shown.has(cluster) and HaplogroupClassifier.node_index(cluster) < 0:
		return cluster
	return ""


## A tree node standing for a clade name: the node itself, or for a
## cluster without its own node (B) the member nearest the root.
static func _representative(name: String) -> int:
	var n := HaplogroupClassifier.node_index(name)
	if n >= 0:
		return n
	var best := -1
	var best_depth := 1 << 30
	for member in HaplogroupClassifier._cluster_nodes:
		if HaplogroupClassifier._cluster_nodes[member] != name:
			continue
		var m := HaplogroupClassifier.node_index(member)
		var depth := HaplogroupClassifier.true_ancestors(m).size()
		if m >= 0 and depth < best_depth:
			best = m
			best_depth = depth
	return best
