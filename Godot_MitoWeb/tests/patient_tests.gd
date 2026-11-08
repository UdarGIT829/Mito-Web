extends RefCounted
## Patient-level comparisons: statistical tests against reference values
## (scipy and the pSLE 2026 paper's Table 3), lineage rules, haplogroups set
## by hand, and the Compare Patients node on the paper's Table 1 data.

var t

## Udar et al. 2026 Table 1 (Amerindian n = 12 as listed): [subject, haplogroup, sex, age, duration, SLEDAI, SDI]
const TABLE1 := [
	["20-03", "L1c2b1b", "F", 15, 6, 8, 2], ["20-16", "L1c3", "M", 19, 1, 10, 3], ["20-19", "L3", "F", 12, 2, 8, 0],
	["20-23", "L3", "F", 21, 15, 2, 0], ["21-55", "L3", "F", 8, 1, 8, 0],
	["20-05", "A2u", "M", 17, 6, 14, 5], ["20-07", "A2j", "F", 18, 2, 16, 2], ["20-09", "A2", "M", 14, 6, 18, 1],
	["20-11", "A2h1", "F", 12, 1.5, 12, 1], ["20-13", "D1c", "F", 14, 2, 10, 2], ["20-15", "A2o", "F", 17, 1, 6, 0],
	["20-21", "A2w", "F", 23, 17, 8, 16], ["20-25", "C1c*4", "F", 12, 4, 6, 3], ["20-34", "B2d", "F", 12, 0.5, 12, 0],
	["21-37", "A2", "F", 18, 16, 2, 0], ["21-39", "A2", "F", 9, 1, 29, 0], ["21-57", "C1b", "F", 16, 1, 8, 0],
	["20-01", "H1aq1", "F", 17, 0.25, 6, 0], ["20-17", "U5a1d1", "F", 16, 6, 0, 0], ["20-32", "H3", "F", 12, 8, 8, 0],
	["20-35", "T2b21", "F", 16, 4, 5, 0],
	["20-29", "B4c1b2a", "M", 20, 12, 4, 0], ["21-44", "N", "F", 19, 5, 4, 0], ["21-53", "D4", "F", 20, 8, 2, 0], ["21-59", "D5c2", "F", 9, 5, 2, 0],
]


func run(runner) -> void:
	t = runner
	test_tests()
	test_manual_and_lineages()
	test_compare()
	test_table_import()


func check(condition: bool, label: String) -> void:
	t.check(condition, label)


func near(a: float, b: float, tol: float = 0.0015) -> bool:
	return absf(a - b) <= tol


func test_tests() -> void:
	var af := [8, 10, 8, 2, 8]
	var am := [14, 16, 18, 12, 10, 6, 8, 6, 12, 2, 29, 8]
	var ca := [6, 0, 8, 5]
	var asian := [4, 4, 2, 2]
	# Table 3 pairs that do not involve the 13th Amerindian patient missing from Table 1.
	check(near(GroupTests.rank_sum(af, ca)["p"], 0.159), "Wilcoxon African vs Caucasian = 0.159 (paper)")
	check(near(GroupTests.rank_sum(af, asian)["p"], 0.075), "Wilcoxon African vs Asian = 0.075 (paper 0.074)")
	check(near(GroupTests.rank_sum(ca, asian)["p"], 0.243), "Wilcoxon Caucasian vs Asian = 0.243 (paper 0.242)")
	check(near(GroupTests.rank_sum([2, 3, 0, 0, 0], [5, 2, 1, 1, 2, 0, 16, 3, 0, 0, 0, 0, 0])["p"], 0.672), "SDI Wilcoxon African vs Amerindian = 0.672 (paper)")
	check(near(GroupTests.kruskal_wallis([[2, 3, 0, 0, 0], [5, 2, 1, 1, 2, 0, 16, 3, 0, 0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])["p"], 0.121, 0.002), "SDI Kruskal-Wallis = 0.122 (paper 0.121)")
	check(near(GroupTests.kruskal_wallis([af, am, ca, asian])["p"], 0.01574, 0.0001), "Kruskal-Wallis with tie correction (scipy 0.01574)")
	check(near(GroupTests.fisher_exact([[3, 1], [1, 3]])["p"], 0.48571, 0.0001), "Fisher 2x2 (scipy 0.4857)")
	check(near(GroupTests.fisher_exact([[2, 3], [7, 6], [0, 4], [0, 4]])["p"], 0.10397, 0.0001), "Fisher r x c, non-zero SDI with 13 Amerindian = 0.104 (paper 0.103)")
	check(near(GroupTests.fisher_exact([[11, 1], [0, 4]])["p"], 0.002747, 0.00001), "Fisher Amerindian vs Asian SLEDAI >= 6 (p < 0.05, paper)")
	var anova := GroupTests.anova([[15, 19, 12, 21, 8], [17, 18, 14, 12, 14, 17, 23, 12, 12, 18, 9, 16], [17, 16, 12, 16], [20, 19, 20, 9]])
	check(near(anova["f"], 0.22434, 0.0001) and near(anova["p"], 0.87839, 0.0001), "ANOVA (scipy F 0.2243, p 0.8784)")
	var dur := [6, 1, 2, 15, 1, 6, 2, 6, 1.5, 2, 1, 17, 4, 0.5, 16, 1, 1, 0.25, 6, 8, 4, 12, 5, 8, 5]
	var sle := [8, 10, 8, 2, 8, 14, 16, 18, 12, 10, 6, 8, 6, 12, 2, 29, 8, 6, 0, 8, 5, 4, 4, 2, 2]
	var sp := GroupTests.spearman(dur, sle)
	check(near(sp["rho"], -0.43221, 0.0001) and near(sp["p"], 0.03095, 0.0002), "Spearman (scipy rho -0.432, p 0.031)")
	check(GroupTests.median_iqr(ca) == [5.5, 2.5, 7.0] and GroupTests.median_iqr(af) == [8.0, 8.0, 8.0] and GroupTests.median_iqr(asian) == [3.0, 2.0, 4.0], "medians and quartiles as in Table 3")


func test_manual_and_lineages() -> void:
	var rules: Array = Lineages.PRESETS["psle_2026"]["rules"]
	check(Lineages.validate(rules).is_empty(), "paper lineage rules name real haplogroups")
	var expected := {"African": 5, "Amerindian": 12, "Caucasian": 4, "Asian": 4}
	var counts := {}
	for row in TABLE1:
		var node := HaplogroupClassifier.resolve_name(row[1])
		var group := Lineages.group_for(node, rules)
		counts[group] = counts.get(group, 0) + 1
	check(counts == expected, "every Table 1 haplogroup falls in the paper's lineage (got %s)" % counts)
	check(HaplogroupClassifier.node_name(HaplogroupClassifier.resolve_name("C1c*4")) == "C1c", "C1c*4 resolves to C1c")
	check(HaplogroupClassifier.manual_result("xyz").has("error"), "non-haplogroups are rejected")
	var parsed := Lineages.parse_text("# comment\nL3 = African\nA2 = Amerindian\nnot-a-group = X\nbad line")
	check(parsed["rules"].size() == 3 and parsed["errors"].size() == 2, "rule text parses, unknown clades and bad lines are reported (%s)" % [parsed["errors"]])
	var sample := {"tags": PackedStringArray(["lupusN"])}
	HaplogroupClassifier.apply_to_sample(sample, HaplogroupClassifier.manual_result("A2o"), NAN)
	check(HaplogroupClassifier.is_manual(sample) and "haplogroup=A2o" in sample["tags"] and "haplogroup clade=A" in sample["tags"], "a haplogroup set by hand tags like a classified one")
	Lineages.apply([sample], rules)
	check("lineage=Amerindian" in sample["tags"], "lineage tag from rules")
	Lineages.apply([sample], [])
	check(not Array(sample["tags"]).any(func(x): return x.begins_with("lineage=")), "clearing rules removes lineage tags")


## A store of Table 1 patients (one sample each, two for 20-19) with values.
func _store() -> VariantStore:
	var store := VariantStore.new()
	var samples: Array = []
	var id := 1
	for row in TABLE1:
		var copies := 2 if row[0] == "20-19" else 1
		for k in copies:
			var sample := {"id": id, "source_id": 1, "column": "", "name": "%s-%d" % [row[0], k], "subject": row[0],
				"tags": PackedStringArray(["lupusN"]), "build": "rcrs",
				"attributes": {"Sex": row[2], "Age": float(row[3]), "Duration": float(row[4]), "SLEDAI": float(row[5]), "SDI": float(row[6])}}
			HaplogroupClassifier.apply_to_sample(sample, HaplogroupClassifier.manual_result(row[1]), NAN)
			samples.append(sample)
			id += 1
	store.append(samples, VariantStore.new_call_columns())
	Lineages.apply(store.samples, Lineages.PRESETS["psle_2026"]["rules"])
	return store


func test_compare() -> void:
	var store := _store()
	var g := Workflows.build(Workflows.find_builtin("lupus_lineages"))
	check(GraphModel.validate(g) == "", "lupus workflow graph is valid")
	var results := GraphEvaluator.new().evaluate(g, store, 1)
	var by_name := {}
	for n in g["nodes"]:
		if n["type"] == "patient_compare":
			by_name[n["config"]["name"]] = results[n["id"]]
	var t3: Dictionary = by_name["Table 3: Activity and damage"]
	check(t3["error"] == "", "Table 3 evaluates (%s)" % t3["error"])
	var p: Dictionary = t3["patients"]
	check(p["groups"].map(func(x): return x["n"]) == [5, 12, 4, 4], "patients counted once per subject (20-19 has two samples)")
	var sledai: Dictionary = p["rows"][0]
	check(sledai["cells"].map(func(c): return c["text"]) == ["8 [8, 8]", "11 [7, 15]", "5.5 [2.5, 7]", "3 [2, 4]"], "SLEDAI medians (got %s)" % [sledai["cells"].map(func(c): return c["text"])])
	var pair: Array = sledai["pairwise"].filter(func(x): return x["a"] == "African" and x["b"] == "Caucasian")
	check(pair.size() == 1 and near(pair[0]["p"], 0.159), "pairwise Wilcoxon in the results table")
	var t4: Dictionary = by_name["Table 4: Activity thresholds"]["patients"]
	check(t4["rows"][0]["cells"].map(func(c): return c["text"]) == ["4 (80%)", "11 (92%)", "2 (50%)", "0 (0%)"], "SLEDAI >= 6 percentages match Table 4 (80, 92, 50, 0)")
	var t2: Dictionary = by_name["Table 2: Demographics"]["patients"]
	check(t2["rows"][0]["cells"][0]["text"] == "4 (80%)" and near(t2["rows"][1]["p"], 0.8784, 0.0005), "Table 2: female share and age ANOVA")
	var corr: Dictionary = by_name["Duration and activity"]["patients"]["rows"][0]
	check(near(corr["p"], 0.031, 0.001) and "rho -0.43" in corr["note"], "overall Spearman duration vs SLEDAI (%s)" % corr["note"])
	check(by_name["Table 5: Lupus nephritis"]["patients"]["rows"][0]["missing"] == 25, "missing values are counted, not guessed")
	var csv := PatientCompare.to_csv(p)
	check(csv.begins_with("Measure,African (n = 5)") and "African vs Caucasian" in csv, "results export as CSV")
	# Disagreeing samples are reported.
	store.sample_by_id(4)["attributes"]["SLEDAI"] = 3.0
	var again: Dictionary = GraphEvaluator.new().evaluate(g, store, 2)
	var warned: bool = again.values().any(func(r): return r.get("patients", {}).get("warnings", []).any(func(w): return "different SLEDAI" in w))
	check(warned, "a patient whose samples disagree is reported")


func test_table_import() -> void:
	var table := TagImport.parse("Subject,Haplogroup,SLEDAI,Sex\n20-03,L1c2b1b,8,F\n20-05,A2u,,M\n")
	var samples := [{"name": "a", "subject": "20-03", "tags": PackedStringArray()}, {"name": "b", "subject": "20-05", "tags": PackedStringArray()}]
	var r := TagImport.apply(samples, table, 0, [], false, -1, [2, 3], 1)
	check(samples[0]["attributes"] == {"SLEDAI": 8.0, "Sex": "F"} and samples[1]["attributes"] == {"Sex": "M"}, "value columns: numbers as numbers, blanks skipped")
	check(r["haplogroups"].size() == 2 and r["haplogroups"][1][1] == "A2u", "haplogroup column returned for setting by hand")
	check(PatientCompare.matches(1.0, {"op": "=", "value": "yes"}) and PatientCompare.matches("Yes", {"op": "=", "value": "yes"}) and not PatientCompare.matches("no", {"op": "=", "value": "yes"}), "yes/1 conditions")
	check(PatientCompare.parse_condition(">= 6") == {"op": ">=", "value": 6.0} and PatientCompare.parse_condition("F") == {"op": "=", "value": "F"}, "conditions parse")
