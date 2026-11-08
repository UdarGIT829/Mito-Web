extends RefCounted
## Patient-level comparisons: statistical tests against reference values
## (scipy, and p-values a published Stata analysis reported for the reference
## table below), lineage rules, haplogroups set by hand, and the Compare
## Patients node.

var t

## A published 25-patient reference table: [subject, haplogroup, sex, age, duration, activity score, damage score]
const TABLE1 := [
	["P01", "L1c2b1b", "F", 15, 6, 8, 2], ["P02", "L1c3", "M", 19, 1, 10, 3], ["P03", "L3", "F", 12, 2, 8, 0],
	["P04", "L3", "F", 21, 15, 2, 0], ["P05", "L3", "F", 8, 1, 8, 0],
	["P06", "A2u", "M", 17, 6, 14, 5], ["P07", "A2j", "F", 18, 2, 16, 2], ["P08", "A2", "M", 14, 6, 18, 1],
	["P09", "A2h1", "F", 12, 1.5, 12, 1], ["P10", "D1c", "F", 14, 2, 10, 2], ["P11", "A2o", "F", 17, 1, 6, 0],
	["P12", "A2w", "F", 23, 17, 8, 16], ["P13", "C1c*4", "F", 12, 4, 6, 3], ["P14", "B2d", "F", 12, 0.5, 12, 0],
	["P15", "A2", "F", 18, 16, 2, 0], ["P16", "A2", "F", 9, 1, 29, 0], ["P17", "C1b", "F", 16, 1, 8, 0],
	["P18", "H1aq1", "F", 17, 0.25, 6, 0], ["P19", "U5a1d1", "F", 16, 6, 0, 0], ["P20", "H3", "F", 12, 8, 8, 0],
	["P21", "T2b21", "F", 16, 4, 5, 0],
	["P22", "B4c1b2a", "M", 20, 12, 4, 0], ["P23", "N", "F", 19, 5, 4, 0], ["P24", "D4", "F", 20, 8, 2, 0], ["P25", "D5c2", "F", 9, 5, 2, 0],
]


func run(runner) -> void:
	t = runner
	test_tests()
	test_manual_and_lineages()
	test_compare()
	test_table_import()
	test_lineage_tree()
	await test_mixture()


func check(condition: bool, label: String) -> void:
	t.check(condition, label)


func near(a: float, b: float, tol: float = 0.0015) -> bool:
	return absf(a - b) <= tol


func test_tests() -> void:
	var af := [8, 10, 8, 2, 8]
	var am := [14, 16, 18, 12, 10, 6, 8, 6, 12, 2, 29, 8]
	var ca := [6, 0, 8, 5]
	var asian := [4, 4, 2, 2]
	# Published pairwise p-values for groups that are complete in the reference table.
	check(near(GroupTests.rank_sum(af, ca)["p"], 0.159), "Wilcoxon African vs Caucasian = 0.159 (published)")
	check(near(GroupTests.rank_sum(af, asian)["p"], 0.075), "Wilcoxon African vs Asian = 0.075 (published 0.074)")
	check(near(GroupTests.rank_sum(ca, asian)["p"], 0.243), "Wilcoxon Caucasian vs Asian = 0.243 (published 0.242)")
	check(near(GroupTests.rank_sum([2, 3, 0, 0, 0], [5, 2, 1, 1, 2, 0, 16, 3, 0, 0, 0, 0, 0])["p"], 0.672), "damage-score Wilcoxon African vs Amerindian = 0.672 (published)")
	check(near(GroupTests.kruskal_wallis([[2, 3, 0, 0, 0], [5, 2, 1, 1, 2, 0, 16, 3, 0, 0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])["p"], 0.121, 0.002), "damage-score Kruskal-Wallis = 0.122 (published 0.121)")
	check(near(GroupTests.kruskal_wallis([af, am, ca, asian])["p"], 0.01574, 0.0001), "Kruskal-Wallis with tie correction (scipy 0.01574)")
	check(near(GroupTests.fisher_exact([[3, 1], [1, 3]])["p"], 0.48571, 0.0001), "Fisher 2x2 (scipy 0.4857)")
	check(near(GroupTests.fisher_exact([[2, 3], [7, 6], [0, 4], [0, 4]])["p"], 0.10397, 0.0001), "Fisher r x c, non-zero damage score with 13 Amerindian = 0.104 (published 0.103)")
	check(near(GroupTests.fisher_exact([[11, 1], [0, 4]])["p"], 0.002747, 0.00001), "Fisher Amerindian vs Asian activity >= 6 (p < 0.05, published)")
	var anova := GroupTests.anova([[15, 19, 12, 21, 8], [17, 18, 14, 12, 14, 17, 23, 12, 12, 18, 9, 16], [17, 16, 12, 16], [20, 19, 20, 9]])
	check(near(anova["f"], 0.22434, 0.0001) and near(anova["p"], 0.87839, 0.0001), "ANOVA (scipy F 0.2243, p 0.8784)")
	var dur := [6, 1, 2, 15, 1, 6, 2, 6, 1.5, 2, 1, 17, 4, 0.5, 16, 1, 1, 0.25, 6, 8, 4, 12, 5, 8, 5]
	var sle := [8, 10, 8, 2, 8, 14, 16, 18, 12, 10, 6, 8, 6, 12, 2, 29, 8, 6, 0, 8, 5, 4, 4, 2, 2]
	var sp := GroupTests.spearman(dur, sle)
	check(near(sp["rho"], -0.43221, 0.0001) and near(sp["p"], 0.03095, 0.0002), "Spearman (scipy rho -0.432, p 0.031)")
	check(GroupTests.median_iqr(ca) == [5.5, 2.5, 7.0] and GroupTests.median_iqr(af) == [8.0, 8.0, 8.0] and GroupTests.median_iqr(asian) == [3.0, 2.0, 4.0], "medians and quartiles as published")


func test_manual_and_lineages() -> void:
	var rules: Array = Lineages.PRESETS["four_lineages"]["rules"]
	check(Lineages.validate(rules).is_empty(), "preset lineage rules name real haplogroups")
	var expected := {"African": 5, "Amerindian": 12, "Caucasian": 4, "Asian": 4}
	var counts := {}
	for row in TABLE1:
		var node := HaplogroupClassifier.resolve_name(row[1])
		var group := Lineages.group_for(node, rules)
		counts[group] = counts.get(group, 0) + 1
	check(counts == expected, "every reference haplogroup falls in its published lineage (got %s)" % counts)
	check(HaplogroupClassifier.node_name(HaplogroupClassifier.resolve_name("C1c*4")) == "C1c", "C1c*4 resolves to C1c")
	check(HaplogroupClassifier.manual_result("xyz").has("error"), "non-haplogroups are rejected")
	var parsed := Lineages.parse_text("# comment\nL3 = African\nA2 = Amerindian\nnot-a-group = X\nbad line")
	check(parsed["rules"].size() == 3 and parsed["errors"].size() == 2, "rule text parses, unknown clades and bad lines are reported (%s)" % [parsed["errors"]])
	var sample := {"tags": PackedStringArray(["cohort"])}
	HaplogroupClassifier.apply_to_sample(sample, HaplogroupClassifier.manual_result("A2o"), NAN)
	check(HaplogroupClassifier.is_manual(sample) and "haplogroup=A2o" in sample["tags"] and "haplogroup clade=A" in sample["tags"], "a haplogroup set by hand tags like a classified one")
	Lineages.apply([sample], rules)
	check("lineage=Amerindian" in sample["tags"], "lineage tag from rules")
	Lineages.apply([sample], [])
	check(not Array(sample["tags"]).any(func(x): return x.begins_with("lineage=")), "clearing rules removes lineage tags")


## A store of the reference patients (one sample each, two for P03) with values.
func _store() -> VariantStore:
	var store := VariantStore.new()
	var samples: Array = []
	var id := 1
	for row in TABLE1:
		var copies := 2 if row[0] == "P03" else 1
		for k in copies:
			var sample := {"id": id, "source_id": 1, "column": "", "name": "%s-%d" % [row[0], k], "subject": row[0],
				"tags": PackedStringArray(["cohort"]), "build": "rcrs",
				"attributes": {"Sex": row[2], "Age": float(row[3]), "Duration": float(row[4]), "Activity": float(row[5]), "Damage": float(row[6])}}
			HaplogroupClassifier.apply_to_sample(sample, HaplogroupClassifier.manual_result(row[1]), NAN)
			samples.append(sample)
			id += 1
	store.append(samples, VariantStore.new_call_columns())
	Lineages.apply(store.samples, Lineages.PRESETS["four_lineages"]["rules"])
	return store


## Groups by lineage and five Compare Patients nodes, like a paper's tables.
func _cohort_graph() -> Dictionary:
	var g := GraphModel.new_graph()
	var groups: Array = []
	for lineage in ["African", "Amerindian", "Caucasian", "Asian"]:
		var n := GraphModel.add_node(g, "data_provider", Vector2.ZERO, {"criteria": {"tags_all": ["cohort", Lineages.PREFIX + lineage]}})
		n["title"] = lineage
		groups.append(n)
	var yes := {"op": "=", "value": "yes"}
	var tables := {
		"Demographics": [{"label": "Female", "attribute": "Sex", "summary": "percent", "condition": {"op": "=", "value": "F"}},
			{"label": "Age", "attribute": "Age", "summary": "mean"}],
		"Activity and damage": [{"label": "Activity", "attribute": "Activity", "summary": "median", "pairwise": true},
			{"label": "Damage", "attribute": "Damage", "summary": "median", "pairwise": true}],
		"Thresholds": [{"label": "Activity >= 6", "attribute": "Activity", "summary": "percent", "condition": {"op": ">=", "value": 6.0}, "pairwise": true}],
		"Outcomes": [{"label": "Outcome 1", "attribute": "Outcome1", "summary": "percent", "condition": yes},
			{"label": "Outcome 2", "attribute": "Outcome2", "summary": "percent", "condition": yes},
			{"label": "Outcome 3", "attribute": "Outcome3", "summary": "percent", "condition": yes},
			{"label": "Outcome 4", "attribute": "Outcome4", "summary": "percent", "condition": yes}],
		"Duration": [{"label": "Duration vs activity", "attribute": "Duration", "summary": "correlation", "with": "Activity"}],
	}
	for name in tables:
		var compare := GraphModel.add_node(g, "patient_compare", Vector2.ZERO, {"name": name, "measures": tables[name]})
		for group in groups:
			GraphModel.connect_nodes(g, group["id"], 0, compare["id"], 0)
	return g


func test_compare() -> void:
	var store := _store()
	var g := _cohort_graph()
	check(GraphModel.validate(g) == "", "cohort comparison graph is valid")
	var results := GraphEvaluator.new().evaluate(g, store, 1)
	var by_name := {}
	for n in g["nodes"]:
		if n["type"] == "patient_compare":
			by_name[n["config"]["name"]] = results[n["id"]]
	var t3: Dictionary = by_name["Activity and damage"]
	check(t3["error"] == "", "activity and damage evaluate (%s)" % t3["error"])
	var p: Dictionary = t3["patients"]
	check(p["groups"].map(func(x): return x["n"]) == [5, 12, 4, 4], "patients counted once per subject (P03 has two samples)")
	var activity: Dictionary = p["rows"][0]
	check(activity["cells"].map(func(c): return c["text"]) == ["8 [8, 8]", "11 [7, 15]", "5.5 [2.5, 7]", "3 [2, 4]"], "activity score medians (got %s)" % [activity["cells"].map(func(c): return c["text"])])
	var pair: Array = activity["pairwise"].filter(func(x): return x["a"] == "African" and x["b"] == "Caucasian")
	check(pair.size() == 1 and near(pair[0]["p"], 0.159), "pairwise Wilcoxon in the results table")
	var t4: Dictionary = by_name["Thresholds"]["patients"]
	check(t4["rows"][0]["cells"].map(func(c): return c["text"]) == ["4 (80%)", "11 (92%)", "2 (50%)", "0 (0%)"], "activity >= 6 percentages match the published 80, 92, 50, 0")
	var t2: Dictionary = by_name["Demographics"]["patients"]
	check(t2["rows"][0]["cells"][0]["text"] == "4 (80%)" and near(t2["rows"][1]["p"], 0.8784, 0.0005), "female share and age ANOVA")
	var corr: Dictionary = by_name["Duration"]["patients"]["rows"][0]
	check(near(corr["p"], 0.031, 0.001) and "rho -0.43" in corr["note"], "overall Spearman duration vs activity (%s)" % corr["note"])
	check(by_name["Outcomes"]["patients"]["rows"][0]["missing"] == 25, "missing values are counted, not guessed")
	var t5: Dictionary = by_name["Outcomes"]["patients"]
	check(t5["missing_attributes"] == ["Outcome1", "Outcome2", "Outcome3", "Outcome4"] and "No patient in these groups has a value for Outcome1" in t5["warnings"][0], "missing value columns are named in a warning")
	var csv := PatientCompare.to_csv(p)
	check(csv.begins_with("Measure,African (n = 5)") and "African vs Caucasian" in csv, "results export as CSV")
	# Disagreeing samples are reported.
	store.sample_by_id(4)["attributes"]["Activity"] = 3.0
	var again: Dictionary = GraphEvaluator.new().evaluate(g, store, 2)
	var warned: bool = again.values().any(func(r): return r.get("patients", {}).get("warnings", []).any(func(w): return "different Activity" in w))
	check(warned, "a patient whose samples disagree is reported")


func test_table_import() -> void:
	var table := TagImport.parse("Subject,Haplogroup,Score,Sex\nP01,L1c2b1b,8,F\nP06,A2u,,M\n")
	var samples := [{"name": "a", "subject": "P01", "tags": PackedStringArray()}, {"name": "b", "subject": "P06", "tags": PackedStringArray()}]
	var r := TagImport.apply(samples, table, 0, [], false, -1, [2, 3], 1)
	check(samples[0]["attributes"] == {"Score": 8.0, "Sex": "F"} and samples[1]["attributes"] == {"Sex": "M"}, "value columns: numbers as numbers, blanks skipped")
	check(r["haplogroups"].size() == 2 and r["haplogroups"][1][1] == "A2u", "haplogroup column returned for setting by hand")
	check(PatientCompare.matches(1.0, {"op": "=", "value": "yes"}) and PatientCompare.matches("Yes", {"op": "=", "value": "yes"}) and not PatientCompare.matches("no", {"op": "=", "value": "yes"}), "yes/1 conditions")
	check(PatientCompare.parse_condition(">= 6") == {"op": ">=", "value": 6.0} and PatientCompare.parse_condition("F") == {"op": "=", "value": "F"}, "conditions parse")


## Synthetic VCF: every SNV of the major haplogroup's profile at AF 1, plus
## the minor haplogroup's own SNVs at AF 0.1.
static func _mix_vcf(major: String, minor: String, minor_af: float) -> String:
	var ref := ReferenceSequence.get_build(ReferenceSequence.RCRS)
	var lines := ["##fileformat=VCFv4.2", "##contig=<ID=chrM,length=16569>", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"]
	var major_set := {}
	for v in HaplogroupClassifier.profile(HaplogroupClassifier.node_index(major)):
		major_set[v] = true
	var records := {}
	for pair in [[major_set.keys(), 1.0], [Array(HaplogroupClassifier.profile(HaplogroupClassifier.node_index(minor))).filter(func(v): return not major_set.has(v)), minor_af]]:
		for v in pair[0]:
			var digits := ""
			for ch in v:
				if ch >= "0" and ch <= "9":
					digits += ch
				else:
					break
			var alt: String = v.substr(digits.length())
			if alt.length() != 1 or not alt in ["A", "C", "G", "T"]:
				continue
			var pos := int(digits)
			var base := ref.supplied_bases(pos, 1)
			if base == alt or records.has(pos):
				continue
			records[pos] = "chrM\t%d\t.\t%s\t%s\t.\tPASS\t.\tGT:AF\t0/1:%s" % [pos, base, alt, pair[1]]
	var keys := records.keys()
	keys.sort()
	for k in keys:
		lines.append(records[k])
	return "\n".join(lines) + "\n"


func test_mixture() -> void:
	check(HaplogroupClassifier.tree_distance(HaplogroupClassifier.node_index("H1"), HaplogroupClassifier.node_index("H1")) == 0
		and HaplogroupClassifier.tree_distance(HaplogroupClassifier.node_index("H1"), HaplogroupClassifier.node_index("H")) == 1, "tree distance")
	var store := VariantStore.new()
	var id := 1
	for case in [["U5b1b1", "J1c5a1", 0.1], ["A2", "A2", 0.1]]:
		var raw: PackedByteArray = _mix_vcf(case[0], case[1], case[2]).to_utf8_buffer()
		var inspected := VcfImporter.inspect("m%d.vcf" % id, raw)
		await VcfImporter.new().import_source({"id": id, "name": "m%d.vcf" % id, "raw": raw, "plan": inspected["plan"]}, store)
		id += 1
	var ids := SampleQuery.select(store)
	var checks := HaplogroupClassifier.mixture_check(store, ids)
	var mixed: Dictionary = checks[ids[0]]
	check(mixed["flag"] and mixed["minor"].begins_with("J1c5") and mixed["major"].begins_with("U5b1b1"), "a low-AF J1c5a1 component in a U5b1b1 sample is flagged (got %s in %s)" % [mixed["minor"], mixed["major"]])
	check(not checks[ids[1]]["flag"], "a clean sample is not flagged")
	var exact_big := GroupTests.fisher_exact([[5, 20], [0, 64], [33, 250], [646, 1510]])
	check(exact_big["method"] == "exact" and exact_big["p"] < 0.001, "a cohort-vs-1000G sized table is still exact")
	var big := GroupTests.fisher_exact([[30, 10, 5], [20, 20, 20], [300, 320, 280], [400, 410, 390]])
	check(big["method"] == "Monte Carlo" and big["p"] < 0.01 and big["p"] > 0.0, "very large tables use a seeded Monte Carlo (p %s)" % big["p"])
	var even := GroupTests.fisher_exact([[100, 100, 100], [100, 100, 100], [100, 100, 100]])
	check(even["method"] == "Monte Carlo" and even["p"] > 0.5, "Monte Carlo p is high when nothing differs (p %s)" % even["p"])
	check(GroupTests.fisher_exact([[3, 1], [1, 3]])["method"] == "exact", "small tables stay exact")


func test_lineage_tree() -> void:
	var store := _store()
	var rules: Array = Lineages.PRESETS["four_lineages"]["rules"]
	var t := Lineages.tree(rules, store.samples)
	var nodes: Dictionary = t["nodes"]
	check(t["roots"].has("L0") and t["roots"].has("L1") and not t["roots"].has("N"), "roots are the deepest clades (got %s)" % [t["roots"]])
	check(nodes["N"]["parent"] == "L3" and nodes["M"]["parent"] == "L3" and nodes["A2"]["parent"] == "A" and nodes["A"]["parent"] == "N", "clades nest as in the tree")
	check(nodes["A"]["group"] == "Asian" and not nodes["A"]["own_rule"], "A outside A2 inherits Asian from N")
	check(nodes["B2"]["parent"] == "B" and nodes["B"]["parent"] == "R", "B (a cluster without its own node) sits between R and B2")
	check(nodes["A2"]["split"] and nodes["A2"]["group"] == "Amerindian" and not nodes["R"]["split"], "A2 splits from Asian; R inherits Asian")
	check(nodes["A2"]["total"] == 8 and nodes["A2"]["own_rule"], "A2 counts its 8 samples")
	check(nodes["L3"]["total"] == 24 and nodes["L3"]["here"] == 4, "L3 counts every sample below it (24 of 26; 2 are L1), 4 directly (got %d, %d)" % [nodes["L3"]["total"], nodes["L3"]["here"]])
	check(nodes["B"]["total"] == 2, "B samples counted once although B matches several nodes (got %d)" % nodes["B"]["total"])
	check(Lineages.color_for("African", rules) != Lineages.color_for("Asian", rules), "lineages get distinct colours")
