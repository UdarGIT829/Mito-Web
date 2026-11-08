class_name GraphEvaluator
extends RefCounted
## Evaluates an analysis graph against a VariantStore.
##
## Results are cached per node and reused while the node's settings, its
## inputs, and the data revision are unchanged, so editing one node only
## recomputes it and what is downstream of it.
##
## Result per node: {"frame": MutationFrame or null, "error": String,
## "summary": String}.

var _cache := {}  # node id -> {"key": String, "result": Dictionary}
var _revision := -1


func evaluate(graph: Dictionary, store: VariantStore, revision: int) -> Dictionary:
	if revision != _revision:
		_cache.clear()
		_revision = revision
	var results := {}
	var keys := {}
	for node_id in GraphModel.topological_order(graph):
		var node := GraphModel.node(graph, node_id)
		var input_keys: Array = []
		for port in GraphModel.TYPES[node["type"]]["inputs"].size():
			# Input names too: renaming a group relabels it downstream (statistics series).
			input_keys.append(GraphModel.inputs_of(graph, node_id, port).map(func(id): return [keys.get(id, "?"), GraphModel.node(graph, id).get("title", "")]))
		var key := str([node["type"], node["config"], input_keys]).sha256_text()
		keys[node_id] = key
		var cached: Dictionary = _cache.get(node_id, {})
		if cached.get("key") == key:
			results[node_id] = cached["result"]
			continue
		var result := _evaluate_node(graph, node, store, results)
		_cache[node_id] = {"key": key, "result": result}
		results[node_id] = result
	for node_id in _cache.keys():
		if not results.has(node_id):
			_cache.erase(node_id)
	return results


func _evaluate_node(graph: Dictionary, node: Dictionary, store: VariantStore, results: Dictionary) -> Dictionary:
	var config: Dictionary = node["config"]
	var inputs: Array = []
	for port in GraphModel.TYPES[node["type"]]["inputs"].size():
		var frames: Array = []
		for source_id in GraphModel.inputs_of(graph, node["id"], port):
			var upstream: Dictionary = results.get(source_id, {})
			if upstream.get("frame") == null:
				return _error("An input has a problem; fix the node feeding this one.")
			frames.append(upstream["frame"])
		inputs.append(frames)
	match node["type"]:
		"data_provider":
			var sample_ids := SampleQuery.select(store, config.get("criteria", {}))
			if sample_ids.is_empty():
				return _error("No samples match. Import files or change the selection.")
			var options := {"filter_values": ["PASS"]} if config.get("pass_only", false) else {}
			return _ok(MutationFrame.from_samples(store, sample_ids, options))
		"tag_filter":
			if inputs[0].is_empty():
				return _error("Connect an input.")
			var tags: Array = config.get("tags", [])
			if tags.is_empty():
				return _error("Choose at least one tag.")
			return _ok(filter_by_tags(inputs[0][0], tags, config.get("match", "any"), config.get("restrict_samples", false)))
		"comparison":
			if inputs[0].is_empty():
				return _error("Connect at least one Include input.")
			var combined: MutationFrame = inputs[0][0]
			for i in range(1, inputs[0].size()):
				combined = combined.union(inputs[0][i]) if config.get("combine") == "union" else combined.intersect(inputs[0][i])
			for excluded in inputs[1]:
				combined = combined.subtract(excluded)
			return _ok(combined)
		"logical_filter":
			if inputs[0].is_empty():
				return _error("Connect an input.")
			var filters: Dictionary = config.get("filters", {})
			var problem := MutationFilters.validate(filters)
			if problem:
				return _error(problem)
			var flag: String = config.get("flag_name", "passes").strip_edges()
			return _ok(MutationFilters.apply(inputs[0][0], filters, config.get("mode", "drop"), flag if flag else "passes"))
		"germline_somatic":
			if inputs[0].is_empty():
				return _error("Connect an input.")
			var thresholds: Dictionary = config.get("thresholds", AfClasses.DEFAULT_THRESHOLDS)
			var problem := AfClasses.validate(thresholds)
			if problem:
				return _error(problem)
			var keep: Array = config.get("keep", [])
			if config.get("mode", "keep") == "keep" and keep.is_empty():
				return _error("Choose which call classes to keep.")
			return _ok(AfClasses.apply(inputs[0][0], thresholds, keep, config.get("mode", "keep")))
		"annotation_filter":
			if inputs[0].is_empty():
				return _error("Connect an input.")
			var load_problem := AnnotationDb.ensure_loaded()
			if load_problem:
				return _error(load_problem)
			var rules: Dictionary = config.get("rules", {})
			var problem := AnnotationFilters.validate(rules)
			if problem:
				return _error(problem)
			var flag: String = config.get("flag_name", "passes").strip_edges()
			return _ok(AnnotationFilters.apply(inputs[0][0], rules, config.get("mode", "drop"), flag if flag else "passes"))
		"enrich":
			if inputs[0].is_empty():
				return _error("Connect an input.")
			var problem := AnnotationDb.ensure_loaded()
			if problem:
				return _error(problem)
			var enriched: MutationFrame = inputs[0][0].take_rows(range(inputs[0][0].row_count()))
			for source in config.get("sources", AnnotationDb.SOURCES):
				if not source in enriched.annotation_sources:
					enriched.annotation_sources.append(source)
			return _ok(enriched)
		"table":
			if inputs[0].is_empty():
				return _error("Connect an input.")
			return _ok(inputs[0][0])
		"statistics":
			if inputs[0].is_empty():
				return _error("Connect one or more groups.")
			var spec: Dictionary = config.get("window", {})
			var problem := WindowStats.validate(spec)
			if problem:
				return _error(problem)
			var window_list := WindowStats.windows(spec)
			var groups: Array = []
			var sources := GraphModel.inputs_of(graph, node["id"], 0)
			for i in inputs[0].size():
				var source := GraphModel.node(graph, sources[i])
				var label := group_label(source)
				var base := label
				var n := 2
				while groups.any(func(g): return g["label"] == label):
					label = "%s (%d)" % [base, n]
					n += 1
				var counted := WindowStats.count(inputs[0][i], window_list)
				counted["label"] = label
				counted["node_id"] = sources[i]
				counted["frame"] = inputs[0][i]
				groups.append(counted)
			return {"frame": null, "error": "", "stats": {"windows": window_list, "groups": groups},
				"summary": "%d windows · %d group%s" % [window_list.size(), groups.size(), "" if groups.size() == 1 else "s"]}
	return _error("Unknown node type.")


## A readable name for a group arriving from `source`.
static func group_label(source: Dictionary) -> String:
	var custom := GraphModel.custom_title(source)
	if custom:
		return custom
	var config: Dictionary = source.get("config", {})
	match source.get("type"):
		"tag_filter":
			if config.get("tags"):
				return "Tags: " + ", ".join(config["tags"])
		"data_provider":
			var criteria: Dictionary = config.get("criteria", {})
			if criteria.has("tags_none"):
				return "Samples without " + ", ".join(criteria["tags_none"])
			var tags: Array = criteria.get("tags_any", criteria.get("tags_all", []))
			return "Samples: " + (", ".join(tags) if tags else "all")
	return source.get("title", "Group")


## Keeps mutations whose carriers include samples with the tags.
## match "any": some carrier has at least one of the tags.
## match "all": every tag is present on at least one carrier.
## restrict_samples narrows the universe and carriers to samples that have
## any of the tags first.
static func filter_by_tags(frame: MutationFrame, tags: Array, match_mode: String, restrict_samples: bool) -> MutationFrame:
	var store := frame.store
	var has_tag := {}  # sample id -> PackedStringArray of matching tags
	for sample_id in frame.samples:
		var sample_tags: PackedStringArray = store.sample_by_id(sample_id).get("tags", PackedStringArray())
		var matching := PackedStringArray()
		for tag in tags:
			if tag in sample_tags:
				matching.append(tag)
		has_tag[sample_id] = matching
	var source := frame
	if restrict_samples:
		source = frame.take_rows(range(frame.row_count()))
		source.samples = PackedInt32Array(Array(frame.samples).filter(func(id): return not has_tag[id].is_empty()))
		source.carriers = frame.carriers.map(func(c: Dictionary) -> Dictionary:
			var kept := {}
			for id in c:
				if not has_tag[id].is_empty():
					kept[id] = c[id]
			return kept)
	var rows: Array = []
	for r in source.row_count():
		var seen := {}
		for sample_id in source.carriers[r]:
			for tag in has_tag[sample_id]:
				seen[tag] = true
		var ok := seen.size() == tags.size() if match_mode == "all" else not seen.is_empty()
		if ok:
			rows.append(r)
	return source.take_rows(rows)


static func _ok(frame: MutationFrame) -> Dictionary:
	var summary := "%d mutation%s · %d sample%s" % [
		frame.row_count(), "" if frame.row_count() == 1 else "s",
		frame.samples.size(), "" if frame.samples.size() == 1 else "s"]
	if frame.af_thresholds:
		var totals := {}
		for r in frame.row_count():
			var counts := AfClasses.row_counts(frame, r)
			for c in counts:
				totals[c] = totals.get(c, 0) + counts[c]
		summary += "\ncalls: " + (AfClasses.describe_counts(totals) if totals else "none")
	if frame.annotation_sources:
		var annotated := 0
		for allele in frame.alleles:
			if not AnnotationDb.lookup(frame.store.allele_pos[allele], frame.store.allele_ref[allele], frame.store.allele_alt[allele], frame.annotation_sources).is_empty():
				annotated += 1
		summary += " · %d annotated" % annotated
	for name in frame.flags:
		var passing := 0
		for v in frame.flags[name]:
			passing += v
		summary += " · %s: %d" % [name, passing]
	return {"frame": frame, "error": "", "summary": summary}


static func _error(message: String) -> Dictionary:
	return {"frame": null, "error": message, "summary": ""}
