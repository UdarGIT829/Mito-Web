extends RefCounted
## Graph model and evaluator tests, using the engine tests' synthetic VCF
## (S1 EV+Treated, S2 EV+Untreated, S3 Cell+Treated, S4 untagged).

const EngineTests := preload("res://tests/engine_tests.gd")

var t
var store: VariantStore


func run(runner) -> void:
	t = runner
	store = VariantStore.new()
	var raw: PackedByteArray = EngineTests.VCF.to_utf8_buffer()
	var inspected := VcfImporter.inspect("synthetic.vcf", raw)
	var tags := [["EV", "Treated"], ["EV", "Untreated"], ["Cell", "Treated"], []]
	for i in 4:
		inspected["plan"]["samples"][i]["tags"] = PackedStringArray(tags[i])
	await VcfImporter.new().import_source({"id": 1, "name": "synthetic.vcf", "raw": raw, "plan": inspected["plan"]}, store)
	test_model_rules()
	test_example_graph()
	test_caching_and_errors()
	test_filter_and_tag_nodes()
	test_workspace_round_trip()


func check(condition: bool, label: String) -> void:
	t.check(condition, label)


func _labels(frame: MutationFrame) -> Array:
	return Array(frame.alleles).map(func(a): return "%d%s>%s" % [store.allele_pos[a], store.allele_ref[a], store.allele_alt[a]])


func _by_name(graph: Dictionary, results: Dictionary, table_name: String) -> Dictionary:
	for n in graph["nodes"]:
		if n["type"] == "table" and n["config"]["name"] == table_name:
			return results[n["id"]]
	return {}


func test_model_rules() -> void:
	var g := GraphModel.new_graph()
	var a := GraphModel.add_node(g, "data_provider")
	var f := GraphModel.add_node(g, "logical_filter")
	var b := GraphModel.add_node(g, "data_provider")
	var c := GraphModel.add_node(g, "comparison")
	check(GraphModel.connect_nodes(g, a["id"], 0, f["id"], 0) == "", "connect provider to filter")
	check(GraphModel.connect_nodes(g, b["id"], 0, f["id"], 0) == "" and GraphModel.inputs_of(g, f["id"], 0) == [b["id"]], "single input is replaced")
	GraphModel.connect_nodes(g, a["id"], 0, c["id"], 0)
	GraphModel.connect_nodes(g, b["id"], 0, c["id"], 0)
	check(GraphModel.inputs_of(g, c["id"], 0).size() == 2, "multi input keeps both")
	GraphModel.connect_nodes(g, f["id"], 0, c["id"], 1)
	check(GraphModel.connect_nodes(g, c["id"], 0, f["id"], 0) != "", "loops are refused")
	check(GraphModel.connect_nodes(g, f["id"], 0, f["id"], 0) != "", "self connection refused")
	check(GraphModel.connect_nodes(g, a["id"], 0, c["id"], 5) != "", "unknown port refused")
	var order := GraphModel.topological_order(g)
	check(order.find(b["id"]) < order.find(f["id"]) and order.find(f["id"]) < order.find(c["id"]), "topological order")
	GraphModel.remove_node(g, f["id"])
	check(g["edges"].all(func(e): return e["from"] != f["id"] and e["to"] != f["id"]), "removing a node removes its edges")
	check(GraphModel.validate(g) == "" and GraphModel.validate({"nodes": [{"type": "nope"}], "edges": []}) != "", "graph validation")


func test_example_graph() -> void:
	var g := GraphModel.example_graph()
	var results := GraphEvaluator.new().evaluate(g, store, 1)
	check(_labels(_by_name(g, results, "Treated-exclusive")["frame"]) == ["150C>T", "302A>AC"], "example: treated-exclusive")
	check(_labels(_by_name(g, results, "In both")["frame"]) == ["73A>G", "750A>G", "8860A>G"], "example: in both")
	check(_labels(_by_name(g, results, "Untreated-exclusive")["frame"]) == ["302A>ACC"], "example: untreated-exclusive")
	check(results.values().all(func(r): return r["error"] == ""), "example graph has no errors")


func test_caching_and_errors() -> void:
	var g := GraphModel.example_graph()
	var evaluator := GraphEvaluator.new()
	var first := evaluator.evaluate(g, store, 1)
	var second := evaluator.evaluate(g, store, 1)
	check(first.keys().all(func(id): return first[id]["frame"] == second[id]["frame"]), "unchanged graph reuses cached frames")
	var provider: Dictionary = g["nodes"][0]
	var untreated_filter: Dictionary = g["nodes"][3]
	untreated_filter["config"]["tags"] = ["Untreated", "Cell"]
	var third := evaluator.evaluate(g, store, 1)
	check(third[provider["id"]]["frame"] == first[provider["id"]]["frame"], "upstream of an edit stays cached")
	check(third[untreated_filter["id"]]["frame"] != first[untreated_filter["id"]]["frame"], "edited node recomputes")
	var fourth := evaluator.evaluate(g, store, 2)
	check(fourth[provider["id"]]["frame"] != third[provider["id"]]["frame"], "new data revision recomputes everything")
	provider["config"]["criteria"] = {"tags_all": ["NoSuchTag"]}
	var broken := evaluator.evaluate(g, store, 2)
	check(broken[provider["id"]]["error"] != "", "provider with no samples reports an error")
	check(broken[g["nodes"][-1]["id"]]["error"] != "", "errors propagate downstream")
	var lonely := GraphModel.new_graph()
	var table := GraphModel.add_node(lonely, "table")
	check(GraphEvaluator.new().evaluate(lonely, store, 1)[table["id"]]["error"] == "Connect an input.", "unconnected node asks for input")


func test_filter_and_tag_nodes() -> void:
	var g := GraphModel.new_graph()
	var provider := GraphModel.add_node(g, "data_provider", Vector2.ZERO, {"pass_only": true})
	var filter := GraphModel.add_node(g, "logical_filter", Vector2.ZERO, {"filters": {"minimum_presence_percent": 75}, "mode": "flag", "flag_name": "common"})
	var bad := GraphModel.add_node(g, "logical_filter", Vector2.ZERO, {"filters": {"minimum_af_percent": 50}})
	var ev := GraphModel.add_node(g, "tag_filter", Vector2.ZERO, {"tags": ["EV"], "restrict_samples": true})
	GraphModel.connect_nodes(g, provider["id"], 0, filter["id"], 0)
	GraphModel.connect_nodes(g, provider["id"], 0, bad["id"], 0)
	GraphModel.connect_nodes(g, provider["id"], 0, ev["id"], 0)
	var results := GraphEvaluator.new().evaluate(g, store, 1)
	check(not "8860A>G" in _labels(results[provider["id"]]["frame"]), "provider PASS-only option")
	var flagged: MutationFrame = results[filter["id"]]["frame"]
	check(flagged.row_count() == 5 and Array(flagged.flags["common"]) == [1, 0, 0, 0, 1], "flag mode keeps rows and marks them")
	check("common: 2" in results[filter["id"]]["summary"], "summary counts flagged rows")
	check(results[bad["id"]]["error"] != "", "invalid filter settings are reported on the node")
	var ev_frame: MutationFrame = results[ev["id"]]["frame"]
	check(ev_frame.samples.size() == 2 and _labels(ev_frame) == ["73A>G", "150C>T", "302A>AC", "302A>ACC", "750A>G"], "restricted tag filter narrows samples")
	check(ev_frame.carriers[0].size() == 2, "restricted carriers exclude other samples")


func test_workspace_round_trip() -> void:
	var model := WorkspaceModel.new()
	model.data["graph"] = GraphModel.example_graph()
	var loaded := WorkspaceModel.from_bytes(model.to_bytes())
	check(loaded["ok"], "workspace with graph loads")
	if loaded["ok"]:
		var g: Dictionary = loaded["model"].graph()
		check(g["nodes"].size() == 11 and g["edges"].size() == 12, "graph survives save/load")
		check(g["nodes"][0]["position"] is Vector2, "node positions keep their type")
	var legacy := WorkspaceModel.new()
	legacy.data.erase("graph")
	var old := WorkspaceModel.from_bytes(legacy.to_bytes())
	check(old["ok"] and old["model"].graph()["nodes"].is_empty(), "files without a graph get an empty one")
