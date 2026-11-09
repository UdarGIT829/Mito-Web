extends SceneTree

## Screenshots for the project page, run by site/tools/take_screenshots.sh in
## a disposable copy of the app (its own user-data folder). Starts a new
## project, loads only the public 1000 Genomes example data, and saves the
## window at each stage. Nothing is sent to the desktop as input.

var out_dir := ""
var main: Control
var tabs: TabContainer
var workspace: Node


func _initialize() -> void:
	call_deferred("run")


func run() -> void:
	out_dir = OS.get_environment("OUT")
	DirAccess.make_dir_recursive_absolute(out_dir)
	workspace = root.get_node("/root/Workspace")
	main = load("res://scenes/main.tscn").instantiate()
	root.add_child(main)
	await settle(10)
	workspace.new_project()                    # a fresh, empty project in the copy's own user folder
	tabs = main.get("_tabs")
	await settle(5)

	# Files: the example data, imported.
	tabs.get_child(0)._load_example()
	print("importing the example data…")
	print(await workspace.import_finished)
	await shot("files", 0)

	# Samples: PhyloTree 17 haplogroups for all 2,534 people.
	print(await workspace.assign_haplogroups(0.75))
	await shot("samples", 1)

	# Analysis: the "Recreate haplogroup markers" workflow for haplogroup J.
	var analysis = tabs.get_child(2)
	var graph: Dictionary = Workflows.haplogroup_markers_graph("J", 90.0, 5.0, "clade")   # the workflow defaults: common >= 90 %, rare < 5 %
	analysis.insert_graph(graph, true, "Recreate haplogroup markers")
	tabs.current_tab = 2
	await settle(30)
	zoom_graph(analysis, 0.8)
	await shot("graph", 2, 20)

	# Tables: the markers it found.
	for n in workspace.model.graph()["nodes"]:
		if n["type"] == "table":
			analysis.open_table.emit(n["id"])
			break
	await shot("table", -1, 40)

	# Statistics: mutation density along the genome, African vs European ancestry,
	# from the "Treated vs untreated" workflow with population tags as the groups.
	analysis.insert_graph(GraphModel.example_graph("AFR", "EUR"), true, "Treated vs untreated")
	await settle(30)
	for n in workspace.model.graph()["nodes"]:
		if n["type"] == "statistics":
			analysis.open_stats.emit(n["id"])
			break
	await shot("stats", -1, 60)
	quit(0)


## Zooms the graph in on its top-left nodes (the fitted view is too small to read).
func zoom_graph(analysis, zoom: float) -> void:
	var ge: GraphEdit = analysis.get("_graph_edit")
	var corner := Vector2(INF, INF)
	for child in ge.get_children():
		if child is GraphNode:
			corner = Vector2(minf(corner.x, child.position_offset.x), minf(corner.y, child.position_offset.y))
	ge.zoom = zoom
	ge.scroll_offset = corner * zoom - Vector2(30, 30)


func settle(frames := 20) -> void:
	for i in frames:
		await process_frame


func shot(name: String, tab: int, frames := 20) -> void:
	if tab >= 0:
		tabs.current_tab = tab
	await settle(frames)
	await RenderingServer.frame_post_draw
	var image := root.get_texture().get_image()
	if image.get_width() > 1600:
		image.resize(1600, roundi(float(image.get_height()) * 1600.0 / float(image.get_width())), Image.INTERPOLATE_LANCZOS)
	image.save_jpg("%s/%s.jpg" % [out_dir, name], 0.86)
	print("saved ", name)
