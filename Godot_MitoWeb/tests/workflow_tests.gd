extends RefCounted
## Workflows: built-in recipes, the haplogroup marker graph, merging into a
## graph, saving in the workspace, and workflow files.
## Synthetic calls (engine tests VCF):
##   73 S1 S2 S4 | 150 S1 | 302AC S1 | 302ACC S2 | 750 all | 8860 S1 S2 S3

const EngineTests := preload("res://tests/engine_tests.gd")

var t
var store: VariantStore


func run(runner) -> void:
	t = runner
	store = VariantStore.new()
	var raw: PackedByteArray = EngineTests.VCF.to_utf8_buffer()
	var inspected := VcfImporter.inspect("synthetic.vcf", raw)
	var tags := [["macro haplogroup=J", "Treated"], ["macro haplogroup=H", "Untreated"], ["macro haplogroup=J"], []]
	for i in 4:
		inspected["plan"]["samples"][i]["tags"] = PackedStringArray(tags[i])
	await VcfImporter.new().import_source({"id": 1, "name": "synthetic.vcf", "raw": raw, "plan": inspected["plan"]}, store)
	test_builtins()
	test_haplogroup_graph()
	test_merge()
	test_files_and_saving()
	test_node_names()


func check(condition: bool, label: String) -> void:
	t.check(condition, label)


func _table(g: Dictionary, results: Dictionary, name: String) -> MutationFrame:
	for n in g["nodes"]:
		if n["type"] == "table" and n["config"]["name"] == name:
			return results[n["id"]]["frame"]
	return null


func _labels(f: MutationFrame) -> Array:
	return Array(f.alleles).map(func(a): return "%d%s" % [store.allele_pos[a], store.allele_alt[a]])


func test_builtins() -> void:
	var ids: Array = Workflows.builtin().map(func(w): return w["id"])
	check(ids == ["treated_vs_untreated", "haplogroup_markers"], "built-in workflows")
	var wf := Workflows.find_builtin("treated_vs_untreated")
	check(Workflows.resolve_params(wf, {}) == {"tag_a": "Treated", "tag_b": "Untreated"}, "defaults fill missing settings")
	var g := Workflows.build(wf, {"tag_a": "EV", "tag_b": "Cell"})
	var titles: Array = g["nodes"].map(func(n): return n["title"])
	check("EV-exclusive" in titles and "Cell-exclusive" in titles, "settings change the graph (tags EV / Cell)")
	check(GraphModel.validate(g) == "", "built graph is valid")
	check(Workflows.haplogroups_in(store) == [["J", 2], ["H", 1]], "haplogroups present, most common first")


func test_haplogroup_graph() -> void:
	# In J = S1, S3; everyone else = S2, S4 (S4 untagged counts as "else").
	var wf := Workflows.find_builtin("haplogroup_markers")
	var g := Workflows.build(wf, {"haplogroup": "J", "common": 95, "rare": 60})
	var results := GraphEvaluator.new().evaluate(g, store, 1)
	check(results.values().all(func(r): return r["error"] == ""), "haplogroup workflow evaluates without errors")
	var gained := _table(g, results, "J: gained markers")
	var absent := _table(g, results, "J: absent markers")
	check(_labels(gained) == ["8860G"], "gained: common in J, rare elsewhere (got %s)" % [_labels(gained)])
	check(_labels(absent) == ["73G"], "absent: common elsewhere, rare in J (got %s)" % [_labels(absent)])
	check(gained.flags.has("PhyloTree marker") and gained.flags["PhyloTree marker"][0] == 1, "rows are flagged against PhyloTree markers (8860G is one)")
	var titles: Array = g["nodes"].map(func(n): return n["title"])
	check("In J" in titles and "Everyone else" in titles and "Gained in J" in titles, "nodes are named for the haplogroup")


func test_merge() -> void:
	var target := GraphModel.new_graph()
	GraphModel.add_node(target, "data_provider", Vector2(0, 100))
	var source := Workflows.build(Workflows.find_builtin("treated_vs_untreated"))
	var added := GraphModel.merge(target, source)
	check(added.size() == source["nodes"].size() and target["nodes"].size() == source["nodes"].size() + 1, "merge adds every node with new ids")
	check(target["edges"].size() == source["edges"].size(), "merge reconnects every edge")
	var lowest: float = target["nodes"].slice(1).map(func(n): return n["position"].y).min()
	check(lowest >= 420, "merged nodes are placed below the existing graph")
	check(GraphModel.validate(target) == "", "merged graph is valid")


func test_files_and_saving() -> void:
	var graph := Workflows.build(Workflows.find_builtin("haplogroup_markers"), {"haplogroup": "U"})
	var saved := Workflows.from_graph(graph, "U markers", "Recreates U")
	var text := Workflows.to_file_text(saved)
	var back := Workflows.from_file_text(text)
	check(back["ok"], "workflow file round trip")
	var g: Dictionary = back["workflow"]["graph"]
	check(g["nodes"].size() == graph["nodes"].size() and g["nodes"][0]["position"] is Vector2, "positions restored as vectors")
	check(g["edges"].all(func(e): return typeof(e["to_port"]) == TYPE_INT), "ports restored as integers")
	check(back["workflow"]["title"] == "U markers" and back["workflow"]["description"] == "Recreates U", "title and description kept")
	var results := GraphEvaluator.new().evaluate(g, store, 1)
	check(results.values().all(func(r): return r["error"] == "" or "No samples" in r["error"]), "an imported workflow evaluates")
	check(not Workflows.from_file_text("{\"format\": \"other\"}")["ok"] and not Workflows.from_file_text("nonsense")["ok"], "other files are rejected")
	var model := WorkspaceModel.new()
	model.saved_workflows().append(saved)
	var loaded := WorkspaceModel.from_bytes(model.to_bytes())
	check(loaded["ok"] and loaded["model"].saved_workflows().size() == 1 and loaded["model"].saved_workflows()[0]["title"] == "U markers", "saved workflows survive save/load")


func test_node_names() -> void:
	var g := GraphModel.new_graph()
	var source := GraphModel.add_node(g, "data_provider")
	var stats := GraphModel.add_node(g, "statistics")
	GraphModel.connect_nodes(g, source["id"], 0, stats["id"], 0)
	check(source["description"] == "" and GraphModel.custom_title(source) == "", "new nodes have their type's name and no description")
	var evaluator := GraphEvaluator.new()
	var before: Dictionary = evaluator.evaluate(g, store, 1)[stats["id"]]
	check(before["stats"]["groups"][0]["label"] == "Samples: all", "an unnamed group is labelled from its settings")
	source["title"] = "Cohort"
	source["description"] = "Every sample, PASS or not."
	var after: Dictionary = evaluator.evaluate(g, store, 1)[stats["id"]]
	check(GraphModel.custom_title(source) == "Cohort" and after["stats"]["groups"][0]["label"] == "Cohort", "renaming a node relabels its group (cache refreshed)")
	var target := GraphModel.new_graph()
	GraphModel.merge(target, g)
	check(target["nodes"][0]["title"] == "Cohort" and target["nodes"][0]["description"] == "Every sample, PASS or not.", "merge keeps names and descriptions")
	var back := Workflows.from_file_text(Workflows.to_file_text(Workflows.from_graph(g, "x", "")))
	check(back["workflow"]["graph"]["nodes"][0]["description"] == "Every sample, PASS or not.", "workflow files keep descriptions")
	var model := WorkspaceModel.new()
	model.graph()["nodes"] = g["nodes"]
	model.graph()["edges"] = g["edges"]
	model.graph()["next_id"] = g["next_id"]
	var loaded := WorkspaceModel.from_bytes(model.to_bytes())
	check(loaded["ok"] and loaded["model"].graph()["nodes"][0]["description"] == "Every sample, PASS or not.", "descriptions survive save/load")
