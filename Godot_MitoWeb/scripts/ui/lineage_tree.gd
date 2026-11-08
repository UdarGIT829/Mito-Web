class_name LineageTree
extends VBoxContainer
## Lineage rules drawn as a haplogroup tree: each clade a rule names, and
## each major haplogroup the samples fall in, under its nearest shown
## ancestor. Colour = lineage; a solid band marks the clade where a rule sets
## it, a faint band a clade that inherits it; "split" marks where a clade
## leaves its parent's lineage. Counts are samples in this project.

var _legend: HFlowContainer
var _tree: Tree


func _init() -> void:
	add_theme_constant_override("separation", 6)
	_legend = HFlowContainer.new()
	_legend.add_theme_constant_override("h_separation", 14)
	add_child(_legend)
	_tree = Tree.new()
	_tree.size_flags_vertical = Control.SIZE_EXPAND_FILL
	_tree.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_tree.custom_minimum_size = Vector2(620, 300)
	_tree.hide_root = true
	_tree.columns = 3
	_tree.column_titles_visible = true
	for pair in [[0, "Haplogroup"], [1, "Lineage"], [2, "Samples"]]:
		_tree.set_column_title(pair[0], pair[1])
	_tree.set_column_expand(0, true)
	_tree.set_column_custom_minimum_width(0, 200)
	_tree.set_column_custom_minimum_width(1, 280)
	_tree.set_column_custom_minimum_width(2, 110)
	_tree.set_column_expand(1, false)
	_tree.set_column_expand(2, false)
	add_child(_tree)
	var key := Label.new()
	key.text = "Solid band: a rule sets the lineage at this clade. Faint band: inherited from above. \"split from\" marks where a clade leaves its parent's lineage."
	key.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	key.custom_minimum_size.x = 600
	key.add_theme_color_override("font_color", Color(0.7, 0.75, 0.8))
	add_child(key)


func show_rules(rules: Array, samples: Array) -> void:
	var data := Lineages.tree(rules, samples)
	var nodes: Dictionary = data["nodes"]
	# Legend: one chip per lineage with its sample count (from top-level clades).
	for child in _legend.get_children():
		child.queue_free()
	var per_group := {}
	for name in nodes:
		per_group[nodes[name]["group"]] = per_group.get(nodes[name]["group"], 0) + nodes[name]["here"]
	var groups := Lineages.groups(rules)
	if per_group.get("", 0) > 0:
		groups.append("")
	for g in groups:
		_legend.add_child(_chip(g, Lineages.color_for(g, rules), per_group.get(g, 0)))
	_tree.clear()
	var root := _tree.create_item()
	if nodes.is_empty():
		var empty := _tree.create_item(root)
		empty.set_text(0, "No rules yet")
		return
	var top := _tree.create_item(root)
	top.set_text(0, "Root of the human mtDNA tree")
	top.set_custom_color(0, Color(0.7, 0.75, 0.8))
	top.set_selectable(0, false)
	for name in data["roots"]:
		_add(top, name, nodes, rules)


func _add(parent: TreeItem, name: String, nodes: Dictionary, rules: Array) -> void:
	var n: Dictionary = nodes[name]
	var item := _tree.create_item(parent)
	var color := Lineages.color_for(n["group"], rules)
	item.set_text(0, name)
	item.set_custom_color(0, color.lightened(0.25) if n["own_rule"] else Color(0.85, 0.87, 0.9))
	var label: String = n["group"] if n["group"] else "no lineage"
	if n["split"]:
		var from: String = nodes[n["parent"]]["group"]
		label += "  (split from %s)" % (from if from else "no lineage")
	item.set_text(1, label)
	item.set_custom_bg_color(1, Color(color, 0.85 if n["own_rule"] else 0.28))
	item.set_custom_color(1, Color.WHITE if n["own_rule"] else color.lightened(0.35))
	item.set_tooltip_text(1, "%s: %s%s" % [name, ("rule %s = %s" % [name, n["group"]]) if n["own_rule"] else ("inherits %s" % (n["group"] if n["group"] else "no lineage")),
		"\nSplits from its parent's lineage here." if n["split"] else ""])
	if n["total"] > 0:
		item.set_text(2, str(n["total"]) if n["here"] == n["total"] else "%d (%d here)" % [n["total"], n["here"]])
		item.set_tooltip_text(2, "%d sample(s) in %s, %d not in a clade shown below it" % [n["total"], name, n["here"]])
	else:
		item.set_custom_color(0, Color(0.55, 0.58, 0.62) if not n["own_rule"] else color.lightened(0.1))
	item.set_selectable(0, false)
	item.set_selectable(1, false)
	item.set_selectable(2, false)
	for child in n["children"]:
		_add(item, child, nodes, rules)


static func _chip(group: String, color: Color, count: int) -> HBoxContainer:
	var box := HBoxContainer.new()
	var swatch := ColorRect.new()
	swatch.color = color
	swatch.custom_minimum_size = Vector2(14, 14)
	swatch.size_flags_vertical = Control.SIZE_SHRINK_CENTER
	box.add_child(swatch)
	var label := Label.new()
	label.text = "%s %d" % [group if group else "No lineage", count]
	box.add_child(label)
	return box
