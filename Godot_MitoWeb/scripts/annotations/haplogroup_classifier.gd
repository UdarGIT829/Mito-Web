class_name HaplogroupClassifier
extends RefCounted
## Assigns each sample the best-matching PhyloTree 17 haplogroup, using
## HaploGrep's Kulczynski score with its per-variant weights:
##   score = 1/2 * ( w(found) / w(expected) + w(found) / w(sample) )
## found = sample variants the haplogroup expects; expected = the
## haplogroup's profile; sample = the sample's variants that occur anywhere
## in the tree. Hotspots (309.1C, 315.1C, 16182C, ...) are ignored and indels
## are aligned with HaploGrep's rules, as HaploGrep does. Bundle built by
## tools/build_phylotree.py.

const BUNDLE_PATH := "res://annotations/phylotree17.json.gz"
const TOP := 3

static var _loaded := false
static var _error := ""
static var tree_name := ""
static var sources: Array = []
static var _poly_index := {}
static var _polys := PackedStringArray()
static var _weights := PackedFloat64Array()
static var _names := PackedStringArray()
static var _parent := PackedInt32Array()
static var _adds: Array = []
static var _removes: Array = []
static var _expected_w := PackedFloat64Array()
static var _hotspots := {}
static var _clusters := {}
static var _cluster_nodes := {}
static var _rules: Array = []


static func ensure_loaded() -> String:
	if _loaded or _error:
		return _error
	var raw := FileAccess.get_file_as_bytes(BUNDLE_PATH)
	var unzipped := GzipReader.decompress(raw)
	var bundle = JSON.parse_string(unzipped["bytes"].get_string_from_utf8()) if unzipped["ok"] and raw.size() > 0 else null
	if typeof(bundle) != TYPE_DICTIONARY or bundle.get("format") != "mito-web-phylotree":
		_error = "The haplogroup tree is missing from this build."
		return _error
	tree_name = bundle["tree"]
	sources = bundle["sources"]
	_polys = PackedStringArray(bundle["polys"])
	for i in _polys.size():
		_poly_index[_polys[i]] = i
	_weights = PackedFloat64Array(bundle["weights"])
	for h in bundle["hotspots"]:
		_hotspots[str(h).to_upper()] = true
	for c in bundle["clusters"]:
		_clusters[c] = true
	_cluster_nodes = bundle.get("cluster_nodes", {})
	_rules = bundle["rules"]
	for node in bundle["nodes"]:
		_names.append(node[0])
		_parent.append(int(node[1]))
		_adds.append(PackedInt32Array(node[2]))
		_removes.append(PackedInt32Array(node[3]))
	# Nodes are stored parents-first, so one forward pass accumulates weights.
	_expected_w.resize(_names.size())
	for n in _names.size():
		var w := _expected_w[_parent[n]] if _parent[n] >= 0 else 0.0
		for p in _adds[n]:
			w += _weights[p]
		for p in _removes[n]:
			w -= _weights[p]
		_expected_w[n] = w
	_loaded = true
	return ""


static func haplogroup_count() -> int:
	ensure_loaded()
	return _names.size()


## PhyloTree-style notation for one rCRS-normalized allele: 73G, 59D,
## 309.1C. Multi-base substitutions become one entry per changed base.
## Imported alleles are left-normalized (the VCF convention), but PhyloTree
## places an indel at the last equivalent position in a repeat (309.1C, not
## 302.1C; 59d 60d, not 57d 58d), so indels are shifted right first.
static func notation(position: int, ref: String, alt: String) -> PackedStringArray:
	var out := PackedStringArray()
	var rcrs := ReferenceSequence.get_build(ReferenceSequence.RCRS).sequence
	var length := rcrs.length()
	if ref.length() == alt.length():
		for i in ref.length():
			if ref[i] != alt[i]:
				out.append("%d%s" % [position + i, alt[i]])
	elif ref.length() > alt.length() and ref.begins_with(alt):
		var k := ref.length() - alt.length()
		var first := position + alt.length()  # first deleted position
		# Deleting first..first+k-1 equals deleting one base later when the
		# base just after the deletion repeats its first base.
		while first + k <= length and rcrs[first + k - 1] == rcrs[first - 1]:
			first += 1
		for i in k:
			out.append("%dD" % (first + i))
	elif alt.length() > ref.length() and alt.begins_with(ref):
		var after := position + ref.length() - 1  # inserted after this position
		var inserted := alt.substr(ref.length())
		while after < length and rcrs[after] == inserted[0]:
			inserted = inserted.substr(1) + inserted[0]
			after += 1
		out.append("%d.1%s" % [after, inserted])
	return out


## A sample's variants in tree notation: calls placed on rCRS whose AF is at
## least `af_min`, or that have no AF (genotype-only data such as 1000
## Genomes). Alignment rules are applied and hotspots removed.
static func sample_variants(store: VariantStore, sample_id: int, af_min: float) -> PackedStringArray:
	return all_sample_variants(store, PackedInt32Array([sample_id]), af_min).get(sample_id, PackedStringArray())


## sample_variants for many samples in one pass over the calls:
## {sample_id: PackedStringArray}.
static func all_sample_variants(store: VariantStore, sample_ids: PackedInt32Array, af_min: float) -> Dictionary:
	var found := {}
	for id in sample_ids:
		found[id] = {}
	var notations := {}  # cache: allele id (mapped calls) or VCF allele text (native calls)
	var c := store.calls
	var sample_col: PackedInt32Array = c["sample"]
	var allele_col: PackedInt32Array = c["allele"]
	var af_col: PackedFloat64Array = c["af"]
	var status_col: PackedByteArray = c["status"]
	for i in sample_col.size():
		var set: Dictionary = found.get(sample_col[i], {})
		var a := allele_col[i]
		if a == VariantStore.NO_ALLELE or not found.has(sample_col[i]):
			continue
		var af := af_col[i]
		if not is_nan(af) and af < af_min:
			continue
		# Native rCRS calls: the VCF's own allele, read like HaploGrep does.
		# Calls mapped from another build: the rCRS allele, shifted to
		# PhyloTree's indel positions (their VCF positions are not rCRS).
		var key
		if status_col[i] == VariantStore.CallStatus.NATIVE:
			key = "%d:%s:%s" % [c["vcf_pos"][i], store.strings[c["vcf_ref"][i]], store.strings[c["vcf_alt"][i]]]
			if not notations.has(key):
				notations[key] = literal_notation(c["vcf_pos"][i], store.strings[c["vcf_ref"][i]], store.strings[c["vcf_alt"][i]])
		else:
			key = a
			if not notations.has(key):
				notations[key] = notation(store.allele_pos[a], store.allele_ref[a], store.allele_alt[a])
		for v in notations[key]:
			set[v] = true
	var out := {}
	for id in sample_ids:
		out[id] = _finish(found[id])
	return out


## Notation for an allele exactly as written in an rCRS VCF, the way
## HaploGrep reads VCFs: equal lengths are per-base substitutions; otherwise
## the tail of the longer allele is the indel (308 CCTCCCC>CTCCCC is 314d,
## 3106 CN>C is 3107d), with positions kept as given.
static func literal_notation(position: int, ref: String, alt: String) -> PackedStringArray:
	ref = ref.to_upper()
	alt = alt.to_upper()
	var out := PackedStringArray()
	if ref.length() == alt.length():
		for i in ref.length():
			if ref[i] != alt[i]:
				out.append("%d%s" % [position + i, alt[i]])
	elif ref.length() > alt.length():
		for p in range(position + alt.length(), position + ref.length()):
			out.append("%dD" % p)
	else:
		out.append("%d.1%s" % [position + ref.length() - 1, alt.substr(ref.length())])
	return out


## Same as sample_variants, from a list of variant strings (e.g. a profile).
static func normalize_variants(variants: Array) -> PackedStringArray:
	var found := {}
	for v in variants:
		found[str(v).strip_edges().to_upper()] = true
	return _finish(found)


static func _finish(found: Dictionary) -> PackedStringArray:
	ensure_loaded()
	for rule in _rules:
		var error: Array = rule[0]
		if error.all(func(e): return found.has(str(e).to_upper())):
			for e in error:
				found.erase(str(e).to_upper())
			for e in rule[1]:
				found[str(e).to_upper()] = true
	var out := PackedStringArray()
	for v in found:
		if not _hotspots.has(v):
			out.append(v)
	return out


## Returns {"haplogroup", "score", "macro", "top": [[name, score], ...],
## "found": [...], "missing": [...], "extra": [...], "not_in_tree": [...]}.
static func classify(variants: PackedStringArray) -> Dictionary:
	var problem := ensure_loaded()
	if problem:
		return {"error": problem}
	var in_sample := PackedByteArray()
	in_sample.resize(_polys.size())
	var sample_w := 0.0
	var not_in_tree := PackedStringArray()
	for v in variants:
		var id = _poly_index.get(v)
		if id == null:
			not_in_tree.append(v)
		else:
			in_sample[id] = 1
			sample_w += _weights[id]
	var found_w := PackedFloat64Array()
	found_w.resize(_names.size())
	var best: Array = []  # [[score, node], ...] best first
	for n in _names.size():
		var f := found_w[_parent[n]] if _parent[n] >= 0 else 0.0
		for p in _adds[n]:
			if in_sample[p]:
				f += _weights[p]
		for p in _removes[n]:
			if in_sample[p]:
				f -= _weights[p]
		found_w[n] = f
		var e := _expected_w[n]
		var score := 0.5 * ((f / e if e > 0 else 1.0) + (f / sample_w if sample_w > 0 else 1.0))
		if best.size() < TOP or score > best[-1][0]:
			best.append([score, n])
			best.sort_custom(func(a, b): return a[0] > b[0])
			if best.size() > TOP:
				best.pop_back()
	var top: int = best[0][1]
	var expected := profile(top)
	var found := PackedStringArray()
	var missing := PackedStringArray()
	for p in expected:
		if p in variants:
			found.append(p)
		else:
			missing.append(p)
	var extra := PackedStringArray()
	for v in variants:
		if not v in expected and _poly_index.has(v):
			extra.append(v)
	return {
		"haplogroup": _names[top], "score": best[0][0], "macro": macro(top),
		"top": best.map(func(b): return [_names[b[1]], b[0]]),
		"found": found, "missing": missing, "extra": extra, "not_in_tree": not_in_tree,
	}


## The expected variant profile of a haplogroup node.
static func profile(node: int) -> PackedStringArray:
	var chain: Array = []
	var n := node
	while n >= 0:
		chain.append(n)
		n = _parent[n]
	var current := {}
	chain.reverse()
	for i in chain:
		for p in _adds[i]:
			current[_polys[p]] = true
		for p in _removes[i]:
			current.erase(_polys[p])
	return PackedStringArray(current.keys())


static func node_index(name: String) -> int:
	ensure_loaded()
	return _names.find(name)


## The nearest enclosing major haplogroup (e.g. J1c2 -> J, HV0a -> HV).
static func macro(node: int) -> String:
	var n := node
	while n >= 0:
		if _cluster_nodes.has(_names[n]):
			return _cluster_nodes[_names[n]]
		if _clusters.has(_names[n]):
			return _names[n]
		n = _parent[n]
	return _names[node].left(1)


const TAG_PREFIXES := ["haplogroup=", "macro haplogroup="]


## Stores a classification on a sample and replaces its haplogroup tags:
## "haplogroup=J1c2" and "macro haplogroup=J" (usable in tag filters).
static func apply_to_sample(sample: Dictionary, result: Dictionary, af_min: float) -> void:
	sample["haplogroup"] = {
		"name": result["haplogroup"], "score": result["score"], "macro": result["macro"],
		"top": result["top"], "found": result["found"].size(),
		"missing": Array(result["missing"]), "extra": Array(result["extra"]),
		"af_min": af_min, "tree": tree_name,
	}
	var tags := PackedStringArray()
	for tag in sample.get("tags", PackedStringArray()):
		if not TAG_PREFIXES.any(func(p): return tag.begins_with(p)):
			tags.append(tag)
	tags.append("haplogroup=" + result["haplogroup"])
	tags.append("macro haplogroup=" + result["macro"])
	sample["tags"] = tags


## Classifies samples (synchronously). Returns the number classified.
static func assign(store: VariantStore, sample_ids: PackedInt32Array, af_min: float) -> int:
	if ensure_loaded():
		return 0
	var variants := all_sample_variants(store, sample_ids, af_min)
	for id in sample_ids:
		apply_to_sample(store.sample_by_id(id), classify(variants[id]), af_min)
	return sample_ids.size()
