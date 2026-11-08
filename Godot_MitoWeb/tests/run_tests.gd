extends SceneTree
## Headless tests. Run from Godot_MitoWeb/:
##   godot --headless --script res://tests/run_tests.gd

const WorkspaceModelScript := preload("res://scripts/workspace_model.gd")

var _failures := 0
var _checks := 0


func _init() -> void:
	test_new_workspace()
	test_add_source_dedupes_by_content()
	test_round_trip()
	test_rejects_bad_input()
	print("%d checks, %d failures" % [_checks, _failures])
	quit(1 if _failures else 0)


func check(condition: bool, label: String) -> void:
	_checks += 1
	if not condition:
		_failures += 1
		printerr("FAIL: " + label)


func test_new_workspace() -> void:
	var model = WorkspaceModelScript.new()
	check(model.data["format"] == "mito-web-workspace", "new workspace has format tag")
	check(model.sources().is_empty(), "new workspace has no sources")


func test_add_source_dedupes_by_content() -> void:
	var model = WorkspaceModelScript.new()
	var a = model.add_source("a.vcf", "##fileformat=VCFv4.2\n".to_utf8_buffer())
	var b = model.add_source("renamed.vcf", "##fileformat=VCFv4.2\n".to_utf8_buffer())
	model.add_source("c.vcf", "other".to_utf8_buffer())
	check(a["id"] == b["id"], "identical content returns the existing source")
	check(model.sources().size() == 2, "duplicate content is not stored twice")
	check(model.remove_source(a["id"]), "remove_source finds the source")
	check(model.sources().size() == 1, "remove_source removes it")


func test_round_trip() -> void:
	var model = WorkspaceModelScript.new()
	model.set_name("Treated vs untreated")
	model.add_source("x.vcf", "chrM\t3243\t.\tA\tG".to_utf8_buffer())
	var result = WorkspaceModelScript.from_json(model.to_json())
	check(result["ok"], "serialized workspace parses")
	var copy = result["model"]
	check(copy.get_name() == "Treated vs untreated", "name survives round trip")
	check(copy.sources()[0]["text"] == "chrM\t3243\t.\tA\tG", "source text survives round trip")
	check(typeof(copy.sources()[0]["id"]) == TYPE_INT, "ids are restored as ints")
	var next = copy.add_source("y.vcf", "new".to_utf8_buffer())
	check(next["id"] == 2, "ids continue after import")


func test_rejects_bad_input() -> void:
	check(not WorkspaceModelScript.from_json("not json")["ok"], "rejects non-JSON")
	check(not WorkspaceModelScript.from_json("{\"format\": \"other\"}")["ok"], "rejects other formats")
	var future = WorkspaceModelScript.new_workspace_data()
	future["format_version"] = 99
	check(not WorkspaceModelScript.from_json(JSON.stringify(future))["ok"], "rejects newer format versions")
	var broken = WorkspaceModelScript.new_workspace_data()
	broken["sources"] = [{"id": 1}]
	check(not WorkspaceModelScript.from_json(JSON.stringify(broken))["ok"], "rejects malformed sources")
