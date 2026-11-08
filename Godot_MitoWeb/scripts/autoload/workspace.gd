extends Node
## Autoload "Workspace": owns the open WorkspaceModel, runs imports, and
## persists it.
##
## user:// is backed by the browser's IndexedDB in web exports, so the
## autosave survives reloads. Browser storage can still be cleared or evicted,
## which is why the exported workspace file is the durable save.

signal changed
signal error(message: String)
signal import_progress(fraction: float, message: String)
signal import_finished(summary: String)
## Asks the graph editor to rebuild its node views after an outside edit.
signal graph_rebuild_requested
signal haplogroup_progress(done: int, total: int)

const AUTOSAVE_PATH := "user://workspace.mitoweb"
const LEGACY_AUTOSAVE_PATH := "user://workspace.json"
const AUTOSAVE_DELAY_SECONDS := 1.0

var model: WorkspaceModel
var importing := false
## Increases whenever imported samples or calls change, so cached analysis
## results know to recompute. Not saved.
var data_revision := 0
## Increases when live annotations are added, so tables re-read them.
var annotations_revision := 0
var _evaluator := GraphEvaluator.new()
var _last_results := {}
var _importer: VcfImporter
var _autosave_timer: Timer


func _ready() -> void:
	_autosave_timer = Timer.new()
	_autosave_timer.one_shot = true
	_autosave_timer.wait_time = AUTOSAVE_DELAY_SECONDS
	_autosave_timer.timeout.connect(save_autosave)
	add_child(_autosave_timer)
	if not _load_autosave():
		model = WorkspaceModel.new()
	AnnotationDb.live = model.live_annotations()


func is_storage_persistent() -> bool:
	return OS.is_userfs_persistent()


func new_project() -> void:
	model = WorkspaceModel.new()
	data_revision += 1
	_mark_changed()


func rename(value: String) -> void:
	model.set_name(value)
	_mark_changed()


## Adds uploaded files: VCFs, zips of VCFs, and manifests (see
## UploadBundle). Returns a summary message.
func add_files(files: Array) -> String:
	var expanded := UploadBundle.expand(files)
	var added := 0
	var duplicates := 0
	var unreadable := 0
	for file in expanded["files"]:
		var result := model.add_source(file["name"], file["bytes"])
		if not result[1]:
			duplicates += 1
		elif result[0]["status"] == WorkspaceModel.STATUS_ERROR:
			unreadable += 1
		else:
			added += 1
	var manifest := UploadBundle.apply_manifests(model.sources(), expanded["manifests"])
	if added + unreadable > 0 or manifest["applied"] > 0:
		model.touch()
		_mark_changed()
	var parts: Array = ["Added %d file(s)." % added]
	if duplicates:
		parts.append("%d skipped: already in the project (same content)." % duplicates)
	var flagged := model.duplicate_report().size()
	if flagged:
		parts.append("%d file(s) look like samples already in the project and are marked Duplicate." % flagged)
	if unreadable:
		parts.append("%d could not be read as VCF (see the Files list)." % unreadable)
	if expanded["manifests"]:
		parts.append("Manifest settings applied to %d row(s)." % manifest["applied"])
	if manifest["unmatched"]:
		parts.append("Manifest names %d file(s) not uploaded: %s." % [manifest["unmatched"].size(), ", ".join(manifest["unmatched"].slice(0, 3))])
	for e in expanded["errors"] + manifest["errors"]:
		parts.append(e)
	if expanded["skipped"]:
		parts.append("Skipped %d other file(s) in the zip." % expanded["skipped"].size())
	return " ".join(parts)


## Removes files that repeat data already in the project. Returns how many.
func remove_duplicates() -> int:
	var removed := model.remove_duplicates()
	if removed:
		data_revision += 1
		_mark_changed()
	return removed


func remove_sources(source_ids: Array) -> void:
	for source_id in source_ids:
		model.remove_source(source_id)
	data_revision += 1
	_mark_changed()


func undo_imports(source_ids: Array) -> void:
	for source_id in source_ids:
		model.undo_import(source_id)
	data_revision += 1
	_mark_changed()


## Evaluated analysis graph: {node id: {frame, error, summary}}. Cached per
## node, so calling this often is cheap.
## Adds tags from a parsed table to imported samples, and to the import
## settings of their files so a re-import keeps them.
func import_tags(table: Dictionary, match_column: int, tag_columns: Array, prefix: bool, subject_column: int = -1) -> Dictionary:
	return import_table(table, match_column, tag_columns, prefix, subject_column)


## import_tags plus attribute columns (per-sample values such as clinical
## scores) and a haplogroup column that sets haplogroups by hand.
## The result adds "haplogroups_set" and "haplogroup_errors".
func import_table(table: Dictionary, match_column: int, tag_columns: Array, prefix: bool, subject_column: int = -1,
		attribute_columns: Array = [], haplogroup_column: int = -1) -> Dictionary:
	var result := TagImport.apply(model.store.samples, table, match_column, tag_columns, prefix, subject_column, attribute_columns, haplogroup_column)
	var planned: Array = []
	for source in model.sources():
		planned.append_array(source["plan"].get("samples", []))
	TagImport.apply(planned, table, match_column, tag_columns, prefix, subject_column, attribute_columns)
	var set_count := 0
	var errors: Array = []
	for pair in result["haplogroups"]:
		var manual := HaplogroupClassifier.manual_result(pair[1])
		if manual.has("error"):
			if not manual["error"] in errors:
				errors.append(manual["error"])
			continue
		HaplogroupClassifier.apply_to_sample(pair[0], manual, NAN)
		set_count += 1
	result["haplogroups_set"] = set_count
	result["haplogroup_errors"] = errors
	if set_count:
		_apply_lineages()
	data_revision += 1
	model.touch()
	_mark_changed()
	return result


## Sets the haplogroup of samples by hand (overrides the classifier until
## cleared). An empty name clears the override. Returns a message.
func set_haplogroups(sample_ids: Array, text: String) -> String:
	var samples: Array = sample_ids.map(func(id): return model.store.sample_by_id(id)).filter(func(s): return not s.is_empty())
	if samples.is_empty():
		return "Select one or more samples first."
	if text.strip_edges().is_empty():
		for sample in samples:
			if HaplogroupClassifier.is_manual(sample):
				sample.erase("haplogroup")
				sample["tags"] = PackedStringArray(Array(sample["tags"]).filter(func(t): return not HaplogroupClassifier.TAG_PREFIXES.any(func(p): return t.begins_with(p)) and not t.begins_with(Lineages.PREFIX)))
		data_revision += 1
		model.touch()
		_mark_changed()
		return "Cleared the haplogroup set by hand on %d sample(s). Assign haplogroups to classify them again." % samples.size()
	var manual := HaplogroupClassifier.manual_result(text)
	if manual.has("error"):
		return manual["error"]
	for sample in samples:
		HaplogroupClassifier.apply_to_sample(sample, manual, NAN)
	_apply_lineages()
	data_revision += 1
	model.touch()
	_mark_changed()
	return "Set %d sample(s) to %s%s." % [samples.size(), manual["haplogroup"], " (from %s)" % manual["given"] if manual.has("note") else ""]


## Saves lineage rules and tags every sample with its lineage.
## Returns a summary message.
func set_lineage_rules(rules: Array) -> String:
	var stored := model.lineage_rules()
	stored.clear()
	stored.append_array(rules)
	var counts := _apply_lineages()
	data_revision += 1
	model.touch()
	_mark_changed()
	if rules.is_empty():
		return "Lineage rules cleared."
	var parts: Array = Lineages.groups(rules).map(func(g): return "%s %d" % [g, counts.get(g, 0)])
	if counts.get("", 0):
		parts.append("no rule %d" % counts[""])
	return "Lineages: " + ", ".join(parts) + "."


func _apply_lineages() -> Dictionary:
	return Lineages.apply(model.store.samples, model.lineage_rules())


## Assigns every imported sample a PhyloTree 17 haplogroup from its calls
## with AF >= af_min (or no AF), in chunks so the page stays responsive.
## Returns a summary message.
func assign_haplogroups(af_min: float) -> String:
	var problem := HaplogroupClassifier.ensure_loaded()
	if problem:
		return problem
	var store := model.store
	var ids := SampleQuery.select(store)
	if ids.is_empty():
		return "Import samples first."
	var manual := Array(ids).filter(func(id): return HaplogroupClassifier.is_manual(store.sample_by_id(id))).size()
	ids = PackedInt32Array(Array(ids).filter(func(id): return not HaplogroupClassifier.is_manual(store.sample_by_id(id))))
	var variants := HaplogroupClassifier.all_sample_variants(store, ids, af_min)
	var last := Time.get_ticks_msec()
	var counts := {}
	for k in ids.size():
		var result := HaplogroupClassifier.classify(variants[ids[k]])
		HaplogroupClassifier.apply_to_sample(store.sample_by_id(ids[k]), result, af_min)
		counts[result["macro"]] = counts.get(result["macro"], 0) + 1
		if Time.get_ticks_msec() - last > 25:
			haplogroup_progress.emit(k + 1, ids.size())
			await get_tree().process_frame
			last = Time.get_ticks_msec()
	var all_ids := SampleQuery.select(store)
	var checks := HaplogroupClassifier.mixture_check(store, all_ids)
	var flagged := 0
	for id in checks:
		store.sample_by_id(id)["haplogroup_check"] = checks[id]
		if checks[id]["flag"]:
			flagged += 1
	_apply_lineages()
	data_revision += 1
	model.touch()
	_mark_changed()
	var disagree := Array(all_ids).filter(func(id): return HaplogroupClassifier.manual_disagreement(store.sample_by_id(id)) != "").size()
	var checked := "Mixture check: %s%s." % ["%d possible mixture(s)" % flagged if flagged else "no mixtures",
		"; %d file(s) disagree with the haplogroup set by hand" % disagree if disagree else ""]
	if ids.is_empty():
		return "All %d samples have haplogroups set by hand, so none were reclassified. %s See the Mixture check column." % [manual, checked]
	var macros := counts.keys()
	macros.sort_custom(func(a, b): return counts[a] > counts[b])
	return checked + " Assigned haplogroups to %d samples (%s)%s. Tags: haplogroup=…, macro haplogroup=…, haplogroup clade=…" % [ids.size(),
		", ".join(macros.slice(0, 6).map(func(m): return "%s %d" % [m, counts[m]])),
		"; kept %d set by hand" % manual if manual else ""]


## Stores live Ensembl results in the workspace.
func add_live_annotations(results: Dictionary) -> void:
	var live := model.live_annotations()
	for k in results:
		live[k] = {"ensembl": results[k]}
	annotations_revision += 1
	model.touch()
	_mark_changed()


func graph_results() -> Dictionary:
	if not importing:
		_last_results = _evaluator.evaluate(model.graph(), model.store, data_revision * 1000000 + annotations_revision)
	return _last_results


## Call after editing the analysis graph in place.
func graph_edited() -> void:
	model.touch()
	_mark_changed()


## Call after editing a source's plan in place.
func plan_edited() -> void:
	model.touch()
	_mark_changed()


## Imports the given pending sources one after another.
## Imports the given pending sources. With skip_duplicates, files the
## duplicate report flags are left out ("Import all ready").
func import_sources(source_ids: Array, skip_duplicates: bool = false) -> void:
	if importing:
		return
	var queue: Array = []
	var missing_build: Array = []
	var duplicates: Array = []
	var report := model.duplicate_report() if skip_duplicates else {}
	for source_id in source_ids:
		var source := model.source_by_id(source_id)
		if source.get("status") != WorkspaceModel.STATUS_PENDING:
			continue
		if report.has(source_id):
			duplicates.append(source["name"])
			continue
		if source["plan"].get("build", "").is_empty():
			missing_build.append(source["name"])
		else:
			queue.append(source)
	if queue.is_empty():
		import_finished.emit(
			"Choose a reference build for: %s." % ", ".join(missing_build) if missing_build
			else ("Nothing to import: %d duplicate file(s) skipped (use Remove duplicates, or select them and Import selected)." % duplicates.size()
				if duplicates else "Nothing to import: select pending files."))
		return
	importing = true
	changed.emit()
	var imported := 0
	var calls := 0
	var failures: Array = []
	var cancelled := false
	for k in queue.size():
		var source: Dictionary = queue[k]
		_importer = VcfImporter.new()
		_importer.progress.connect(func(fraction: float, message: String) -> void:
			import_progress.emit((k + fraction) / queue.size(), message))
		var result: Dictionary = await _importer.import_source(source, model.store)
		if result.get("cancelled"):
			cancelled = true
			break
		if result["ok"]:
			model.mark_imported(source["id"], result)
			imported += 1
			calls += result["calls"]
		else:
			failures.append("%s: %s" % [source["name"], result["error"]])
	_importer = null
	importing = false
	data_revision += 1
	_mark_changed()
	var parts: Array = ["Imported %d file(s), %d calls." % [imported, calls]]
	if cancelled:
		parts.append("Cancelled; the file in progress was not imported.")
	if missing_build:
		parts.append("Skipped (no reference build chosen): %s." % ", ".join(missing_build))
	if duplicates:
		parts.append("Skipped %d duplicate file(s)." % duplicates.size())
	if failures:
		parts.append("Failed: %s" % "; ".join(failures))
	import_finished.emit(" ".join(parts))


func cancel_import() -> void:
	if _importer:
		_importer.cancelled = true


func export_bytes() -> PackedByteArray:
	return model.to_bytes()


func export_file_name() -> String:
	var safe := model.get_name().strip_edges().validate_filename().replace(" ", "_")
	if safe.is_empty():
		safe = "workspace"
	return "%s.%s" % [safe, WorkspaceModel.FILE_EXTENSION]


## Replaces the open workspace with an imported one. Returns an error message,
## or "" on success.
func import_bytes(bytes: PackedByteArray) -> String:
	if importing:
		return "Wait for the current import to finish."
	var result := WorkspaceModel.from_bytes(bytes)
	if not result["ok"]:
		return result["error"]
	model = result["model"]
	data_revision += 1
	_mark_changed()
	return ""


func save_autosave() -> void:
	var file := FileAccess.open(AUTOSAVE_PATH, FileAccess.WRITE)
	if file == null:
		error.emit("Autosave failed: %s" % error_string(FileAccess.get_open_error()))
		return
	file.store_buffer(model.to_bytes())
	file.close()


func _load_autosave() -> bool:
	for path in [AUTOSAVE_PATH, LEGACY_AUTOSAVE_PATH]:
		if not FileAccess.file_exists(path):
			continue
		var result := WorkspaceModel.from_bytes(FileAccess.get_file_as_bytes(path))
		if result["ok"]:
			model = result["model"]
			return true
		push_warning("Ignoring unreadable autosave %s: %s" % [path, result["error"]])
	return false


func _mark_changed() -> void:
	AnnotationDb.live = model.live_annotations()
	changed.emit()
	if not importing:
		_autosave_timer.start()
