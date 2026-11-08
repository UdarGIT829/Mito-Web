class_name MitoGenes
extends RefCounted
## Pinned NC_012920.1 gene ranges (the 37 NCBI RefSeq GFF3 "gene" features),
## ported from mito_viewer/reference_features.py. 1-based, inclusive.

const LENGTH := 16569

## [symbol, ncbi_gene_id, start, end, strand, biotype]
const GENES := [
	["MT-TF", "4558", 577, 647, "+", "tRNA"],
	["MT-RNR1", "4549", 648, 1601, "+", "rRNA"],
	["MT-TV", "4577", 1602, 1670, "+", "tRNA"],
	["MT-RNR2", "4550", 1671, 3229, "+", "rRNA"],
	["MT-TL1", "4567", 3230, 3304, "+", "tRNA"],
	["MT-ND1", "4535", 3307, 4262, "+", "protein_coding"],
	["MT-TI", "4565", 4263, 4331, "+", "tRNA"],
	["MT-TQ", "4572", 4329, 4400, "-", "tRNA"],
	["MT-TM", "4569", 4402, 4469, "+", "tRNA"],
	["MT-ND2", "4536", 4470, 5511, "+", "protein_coding"],
	["MT-TW", "4578", 5512, 5579, "+", "tRNA"],
	["MT-TA", "4553", 5587, 5655, "-", "tRNA"],
	["MT-TN", "4570", 5657, 5729, "-", "tRNA"],
	["MT-TC", "4511", 5761, 5826, "-", "tRNA"],
	["MT-TY", "4579", 5826, 5891, "-", "tRNA"],
	["MT-CO1", "4512", 5904, 7445, "+", "protein_coding"],
	["MT-TS1", "4574", 7446, 7514, "-", "tRNA"],
	["MT-TD", "4555", 7518, 7585, "+", "tRNA"],
	["MT-CO2", "4513", 7586, 8269, "+", "protein_coding"],
	["MT-TK", "4566", 8295, 8364, "+", "tRNA"],
	["MT-ATP8", "4509", 8366, 8572, "+", "protein_coding"],
	["MT-ATP6", "4508", 8527, 9207, "+", "protein_coding"],
	["MT-CO3", "4514", 9207, 9990, "+", "protein_coding"],
	["MT-TG", "4563", 9991, 10058, "+", "tRNA"],
	["MT-ND3", "4537", 10059, 10404, "+", "protein_coding"],
	["MT-TR", "4573", 10405, 10469, "+", "tRNA"],
	["MT-ND4L", "4539", 10470, 10766, "+", "protein_coding"],
	["MT-ND4", "4538", 10760, 12137, "+", "protein_coding"],
	["MT-TH", "4564", 12138, 12206, "+", "tRNA"],
	["MT-TS2", "4575", 12207, 12265, "+", "tRNA"],
	["MT-TL2", "4568", 12266, 12336, "+", "tRNA"],
	["MT-ND5", "4540", 12337, 14148, "+", "protein_coding"],
	["MT-ND6", "4541", 14149, 14673, "-", "protein_coding"],
	["MT-TE", "4556", 14674, 14742, "-", "tRNA"],
	["MT-CYB", "4519", 14747, 15887, "+", "protein_coding"],
	["MT-TT", "4576", 15888, 15953, "+", "tRNA"],
	["MT-TP", "4571", 15956, 16023, "-", "tRNA"],
]

static var _cache := {}


static func symbols() -> PackedStringArray:
	var out := PackedStringArray()
	for gene in GENES:
		out.append(gene[0])
	return out


## Genes overlapped by the reference interval of an allele, wrapping around
## the circular origin. Returns {"genes": PackedStringArray, "region":
## "genic" | "control_region" | "intergenic"}.
static func locus(position: int, ref: String = "") -> Dictionary:
	var key := "%d:%d" % [position, maxi(ref.length(), 1)]
	if _cache.has(key):
		return _cache[key]
	var end := position + maxi(ref.length(), 1) - 1
	var intervals := [[position, end]] if end <= LENGTH else [[position, LENGTH], [1, end - LENGTH]]
	var genes := PackedStringArray()
	for gene in GENES:
		for interval in intervals:
			if interval[0] <= gene[3] and interval[1] >= gene[2]:
				genes.append(gene[0])
				break
	var control := false
	for interval in intervals:
		if interval[0] <= 576 or interval[1] >= 16024:
			control = true
	var result := {
		"genes": genes,
		"region": "genic" if genes else ("control_region" if control else "intergenic"),
	}
	_cache[key] = result
	return result
