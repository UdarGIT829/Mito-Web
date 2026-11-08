extends RefCounted
## Haplogroup classifier (PhyloTree 17, HaploGrep scoring).

## A variant that occurs nowhere in PhyloTree 17 (checked when written).
const PRIVATE := "3000T"

var t


func run(runner) -> void:
	t = runner
	check(HaplogroupClassifier.ensure_loaded() == "", "tree bundle loads")
	check(HaplogroupClassifier.haplogroup_count() == 5435, "5,435 haplogroups")
	test_notation()
	test_profiles()
	test_rules_and_hotspots()
	await test_from_calls()


func check(condition: bool, label: String) -> void:
	t.check(condition, label)


func test_notation() -> void:
	check(Array(HaplogroupClassifier.notation(73, "A", "G")) == ["73G"], "SNV")
	check(Array(HaplogroupClassifier.notation(522, "ACA", "A")) == ["523D", "524D"], "deletion: one entry per deleted base")
	check(Array(HaplogroupClassifier.notation(309, "C", "CC")) == ["309.1C"], "insertion after the anchor")
	check(Array(HaplogroupClassifier.notation(73, "AT", "GC")) == ["73G", "74C"], "multi-base substitution split per base")


func test_profiles() -> void:
	var j1c := HaplogroupClassifier.node_index("J1c")
	var profile := HaplogroupClassifier.profile(j1c)
	for v in ["185A", "228A", "295T", "462T", "489C", "3010A", "12612G", "13708A", "14798C", "16069T", "16126C"]:
		check(v in profile, "J1c expects %s" % v)
	var r := HaplogroupClassifier.classify(profile)
	check(r["haplogroup"] == "J1c" and is_equal_approx(r["score"], 1.0) and r["macro"] == "J", "a haplogroup's own profile classifies as itself (score 1)")
	check(r["missing"].is_empty() and r["extra"].is_empty(), "nothing missing or extra")
	var rcrs := HaplogroupClassifier.classify(PackedStringArray())
	check(rcrs["haplogroup"] == "H2a2a1", "no variants = rCRS = H2a2a1")
	check(not "263G" in HaplogroupClassifier.profile(HaplogroupClassifier.node_index("L0")), "back mutations remove variants (L0 has 263G!)")
	var b4 := HaplogroupClassifier.node_index("B4")
	check(HaplogroupClassifier.macro(b4) == "B" and HaplogroupClassifier.macro(HaplogroupClassifier.node_index("J1c2")) == "J", "major haplogroups (B4 -> B via cluster members)")
	var partial := Array(profile).filter(func(v): return not v in ["185A", "228A", "14798C"])  # J1c's own markers
	check(HaplogroupClassifier.classify(PackedStringArray(partial))["haplogroup"] == "J1", "without J1c's own markers the parent J1 wins")
	var extra := profile.duplicate()
	extra.append(PRIVATE)
	var with_private := HaplogroupClassifier.classify(extra)
	check(with_private["haplogroup"] == "J1c" and PRIVATE in with_private["not_in_tree"] and is_equal_approx(with_private["score"], 1.0), "variants outside the tree are reported, not scored")


func test_rules_and_hotspots() -> void:
	var left := HaplogroupClassifier.normalize_variants(["8271D", "8272D", "8273D", "8274D", "8275D", "8276D", "8277D", "8278D", "8279D", "73G"])
	check("8281D" in left and not "8271D" in left, "alignment rules move the 9-bp deletion to the tree's positions")
	var hot := HaplogroupClassifier.normalize_variants(["309.1C", "315.1C", "16182C", "73G"])
	check(Array(hot) == ["73G"], "hotspots are ignored")


func test_from_calls() -> void:
	# S1 carries the J1c profile homoplasmically; S2 carries it at AF 0.3.
	var j1c := HaplogroupClassifier.profile(HaplogroupClassifier.node_index("J1c"))
	var rcrs := ReferenceSequence.get_build(ReferenceSequence.RCRS).sequence
	var lines := PackedStringArray(["##fileformat=VCFv4.2", "##contig=<ID=chrM,length=16569>", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2"])
	var snvs: Array = Array(j1c).filter(func(v): return v.left(-1).is_valid_int())
	snvs.sort_custom(func(a, b): return a.left(-1).to_int() < b.left(-1).to_int())
	for v in snvs:
		var pos: int = v.left(-1).to_int()
		lines.append("chrM\t%d\t.\t%s\t%s\t.\tPASS\t.\tGT:AF\t1:0.99\t1:0.3" % [pos, rcrs[pos - 1], v[-1]])
	var raw := ("\n".join(lines) + "\n").to_utf8_buffer()
	var store := VariantStore.new()
	var inspected := VcfImporter.inspect("j.vcf", raw)
	await VcfImporter.new().import_source({"id": 1, "name": "j.vcf", "raw": raw, "plan": inspected["plan"]}, store)
	var ids := SampleQuery.select(store)
	var variants := HaplogroupClassifier.all_sample_variants(store, ids, 0.75)
	check(variants[ids[1]].is_empty(), "calls below the AF cut-off are not used")
	HaplogroupClassifier.assign(store, ids, 0.75)
	var s1 := store.sample_by_id(ids[0])
	var s2 := store.sample_by_id(ids[1])
	check(s1["haplogroup"]["name"] == "J1c" and "haplogroup=J1c" in s1["tags"] and "macro haplogroup=J" in s1["tags"], "J1c sample tagged")
	check(s2["haplogroup"]["name"] == "H2a2a1", "heteroplasmic-only sample stays at rCRS")
	HaplogroupClassifier.assign(store, ids, 0.2)
	check(store.sample_by_id(ids[1])["haplogroup"]["name"] == "J1c", "a lower cut-off lets AF 0.3 calls in")
	check(Array(store.sample_by_id(ids[0])["tags"]).filter(func(x): return x.begins_with("haplogroup=")).size() == 1, "re-running replaces haplogroup tags")
	check(SampleQuery.select(store, {"tags_any": ["macro haplogroup=J"]}).size() == 2, "tags are usable in sample selection")
