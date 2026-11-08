class_name SequenceFeatures
extends RefCounted
## Sequencing-slippage repeat profiles and reference context, ported from
## main.py (longest_repeat, repeat_profile, reference_context).

static var _profile_cache := {}


static func is_dna_motif(motif: String) -> bool:
	if motif.is_empty():
		return false
	for i in motif.length():
		if not motif[i] in "ACGTN":
			return false
	return true


## True when a motif is not made from a smaller repeated motif.
static func is_minimal_motif(motif: String) -> bool:
	for size in range(1, motif.length()):
		if motif.length() % size:
			continue
		if motif.substr(0, size).repeat(motif.length() / size) == motif:
			return false
	return true


## Longest repeated run of a fixed-length motif: [motif, count], or ["", 0].
static func longest_repeat(sequence: String, motif_length: int) -> Array:
	sequence = sequence.to_upper()
	var best_motif := ""
	var best_count := 0
	for start in maxi(sequence.length() - motif_length + 1, 0):
		var motif := sequence.substr(start, motif_length)
		if not is_dna_motif(motif) or not is_minimal_motif(motif):
			continue
		var count := 1
		var next_start := start + motif_length
		while sequence.substr(next_start, motif_length) == motif:
			count += 1
			next_start += motif_length
		if count > best_count and count > 1:
			best_motif = motif
			best_count = count
	return [best_motif, best_count]


## {1: [motif, count], 2: [...], 3: [...]} for 1-, 2- and 3-base repeats.
static func repeat_profile(sequence: String) -> Dictionary:
	if not _profile_cache.has(sequence):
		_profile_cache[sequence] = {
			1: longest_repeat(sequence, 1),
			2: longest_repeat(sequence, 2),
			3: longest_repeat(sequence, 3),
		}
	return _profile_cache[sequence]


## Reference bases flanking a 1-based position, wrapping around the circular
## genome: [before, after].
static func reference_context(reference: String, position: int, flank_size: int = 6) -> Array:
	var length := reference.length()
	var before := ""
	var after := ""
	for offset in range(flank_size, 0, -1):
		before += reference[posmod(position - offset - 1, length)]
	for offset in range(1, flank_size + 1):
		after += reference[posmod(position + offset - 1, length)]
	return [before, after]
