class_name AnnotationFilters
extends RefCounted
## Mutation filters on public annotations (rule "annotation-rules-v1").
##
## Config (every key optional; all set conditions must hold for a mutation
## to "match"):
##   clinvar_categories: ["pathogenic", "benign", "uncertain", "conflicting",
##                        "other", "absent"]   ClinVar wording is grouped as
##                        documented in clinvar_category(); "absent" = no
##                        ClinVar record (never treated as benign)
##   clinvar_min_stars: 0..4   review status, ClinVar's star scheme; mutations
##                        without a record never satisfy a stars minimum
##   gnomad_af_hom_min / gnomad_af_hom_max: 0..1, gnomad_keep_absent: bool
##   haplogroup_marker: "any" | "yes" | "no"    (PhyloTree flag via gnomAD)
##   consequences: [term, ...]   any of these VEP terms (gnomAD or live Ensembl)
##   trna_predicted_pathogenic: true   MitoTIP or PON-mt-tRNA says pathogenic
##   action: "keep" (keep matching) | "remove" (remove matching)
## Use with a drop/flag mode like MutationFilters.

const RULE := "annotation-rules-v1"
const CATEGORIES := ["pathogenic", "benign", "uncertain", "conflicting", "other", "absent"]
const COMMON_CONSEQUENCES := [
	"missense_variant", "synonymous_variant", "stop_gained", "stop_lost", "start_lost",
	"frameshift_variant", "inframe_insertion", "inframe_deletion",
	"non_coding_transcript_exon_variant", "intergenic_variant", "incomplete_terminal_codon_variant",
]


## Groups ClinVar significance wording into one category.
##   conflicting: "Conflicting ..." or both pathogenic and benign wording
##   pathogenic: Pathogenic / Likely pathogenic (any combination)
##   benign: Benign / Likely benign (any combination)
##   uncertain: Uncertain significance
##   other: Affects, risk factor, association, drug response, not provided ...
static func clinvar_category(annotations: Dictionary) -> String:
	if not annotations.has("clinvar"):
		return "absent"
	var text := str(annotations["clinvar"].get("significance", "")).to_lower()
	var pathogenic := "pathogenic" in text
	var benign := "benign" in text
	if "conflicting" in text or (pathogenic and benign):
		return "conflicting"
	if pathogenic:
		return "pathogenic"
	if benign:
		return "benign"
	if "uncertain" in text:
		return "uncertain"
	return "other"


## ClinVar review status as stars (0-4).
static func clinvar_stars(annotations: Dictionary) -> int:
	var review := str(annotations.get("clinvar", {}).get("review", "")).to_lower()
	# "no assertion criteria provided" contains "criteria provided": check it first.
	if review.is_empty() or "no assertion" in review or "no classification" in review:
		return 0
	if "practice guideline" in review:
		return 4
	if "expert panel" in review:
		return 3
	if "multiple submitters" in review and "no conflicts" in review:
		return 2
	if "criteria provided" in review:
		return 1
	return 0


static func validate(config: Dictionary) -> String:
	for c in config.get("clinvar_categories", []):
		if not c in CATEGORIES:
			return "Unknown ClinVar category: %s." % c
	var stars = config.get("clinvar_min_stars", 0)
	if typeof(stars) not in [TYPE_INT, TYPE_FLOAT] or stars < 0 or stars > 4:
		return "ClinVar review stars must be 0 to 4."
	for k in ["gnomad_af_hom_min", "gnomad_af_hom_max"]:
		var v = config.get(k)
		if v != null and (typeof(v) not in [TYPE_INT, TYPE_FLOAT] or v < 0 or v > 1):
			return "gnomAD frequencies must be between 0 and 1."
	if config.get("haplogroup_marker", "any") not in ["any", "yes", "no"]:
		return "Haplogroup marker must be any, yes or no."
	if config.get("action", "keep") not in ["keep", "remove"]:
		return "Action must be keep or remove."
	if not has_condition(config):
		return "Set at least one annotation condition."
	return ""


static func has_condition(config: Dictionary) -> bool:
	return not config.get("clinvar_categories", []).is_empty() or config.get("clinvar_min_stars", 0) > 0 \
		or config.get("gnomad_af_hom_min") != null or config.get("gnomad_af_hom_max") != null \
		or config.get("haplogroup_marker", "any") != "any" or not config.get("consequences", []).is_empty() \
		or config.get("trna_predicted_pathogenic", false)


## Which conditions a mutation fails (empty = it matches every condition).
static func failed_conditions(annotations: Dictionary, config: Dictionary) -> Array:
	var failed: Array = []
	var categories: Array = config.get("clinvar_categories", [])
	if categories and not clinvar_category(annotations) in categories:
		failed.append("ClinVar category")
	if config.get("clinvar_min_stars", 0) > 0 and (not annotations.has("clinvar") or clinvar_stars(annotations) < config["clinvar_min_stars"]):
		failed.append("ClinVar review")
	var g: Dictionary = annotations.get("gnomad", {})
	var lo = config.get("gnomad_af_hom_min")
	var hi = config.get("gnomad_af_hom_max")
	if lo != null or hi != null:
		if g.is_empty() or g.get("af_hom") == null:
			if not config.get("gnomad_keep_absent", false):
				failed.append("not in gnomAD")
		else:
			var af := float(g["af_hom"])
			if (lo != null and af < lo) or (hi != null and af > hi):
				failed.append("gnomAD frequency")
	match config.get("haplogroup_marker", "any"):
		"yes":
			if not g.get("hap_defining", false):
				failed.append("not a haplogroup marker")
		"no":
			if g.get("hap_defining", false):
				failed.append("haplogroup marker")
	var wanted: Array = config.get("consequences", [])
	if wanted:
		var terms := AnnotationDb.consequences(annotations)
		if not wanted.any(func(w): return w in terms):
			failed.append("consequence")
	if config.get("trna_predicted_pathogenic", false):
		var calls := (str(g.get("mitotip", "")) + " " + str(g.get("pon_mt", ""))).to_lower()
		if not "pathogenic" in calls:
			failed.append("tRNA prediction")
	return failed


## Applies the rules. The node's action decides whether matching mutations
## are kept or removed; mode "flag" keeps all rows and adds a pass column.
static func apply(frame: MutationFrame, config: Dictionary, mode: String = "drop", flag_name: String = "passes") -> MutationFrame:
	var store := frame.store
	var remove: bool = config.get("action", "keep") == "remove"
	var passes := PackedByteArray()
	var reasons := PackedStringArray()
	for allele in frame.alleles:
		var annotations := AnnotationDb.lookup(store.allele_pos[allele], store.allele_ref[allele], store.allele_alt[allele])
		var failed := failed_conditions(annotations, config)
		var matches := failed.is_empty()
		var ok := not matches if remove else matches
		passes.append(1 if ok else 0)
		if ok:
			reasons.append("")
		elif remove:
			reasons.append("matches the removal rule (ClinVar: %s)" % clinvar_category(annotations))
		else:
			reasons.append(", ".join(failed))
	if mode == "flag":
		var out := frame.take_rows(range(frame.row_count()))
		out.set_flag(flag_name, passes, reasons)
		return out
	return frame.take_rows(range(frame.row_count()).filter(func(r: int) -> bool: return passes[r] == 1))


static func describe(config: Dictionary) -> String:
	var parts: Array = []
	if config.get("clinvar_categories"):
		parts.append("ClinVar: " + ", ".join(config["clinvar_categories"]))
	if config.get("clinvar_min_stars", 0) > 0:
		parts.append("review >= %d star(s)" % config["clinvar_min_stars"])
	var lo = config.get("gnomad_af_hom_min")
	var hi = config.get("gnomad_af_hom_max")
	if lo != null or hi != null:
		parts.append("gnomAD AF hom %s–%s%s" % ["0" if lo == null else str(lo), "1" if hi == null else str(hi),
			" (or absent)" if config.get("gnomad_keep_absent", false) else ""])
	if config.get("haplogroup_marker", "any") != "any":
		parts.append("haplogroup marker: " + config["haplogroup_marker"])
	if config.get("consequences"):
		parts.append("consequence: " + ", ".join(config["consequences"]))
	if config.get("trna_predicted_pathogenic", false):
		parts.append("tRNA predicted pathogenic")
	if parts.is_empty():
		return "No conditions set"
	return ("Keep" if config.get("action", "keep") == "keep" else "Remove") + " mutations matching:\n" + "\n".join(parts)
