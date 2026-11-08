class_name MitoMapper
extends RefCounted
## Conservative forward-strand mitochondrial allele projection, ported from
## mito_viewer/mapping/mitochondrial.py (algorithm
## "mitochondrial-unit-edit-all-optimal-v1").
##
## All optimal unit-edit global alignments are considered. Any nonunique or
## gapped mapping at an edit is rejected; no arbitrary repeat alignment is
## chosen. The same mapper with source == target trims VCF anchors and
## left-normalizes alleles, so native rCRS calls get identical identities.

const VERSION := "mitochondrial-unit-edit-all-optimal-v1"
const BAND := 64

enum Status { MAPPED, TARGET_REFERENCE, AMBIGUOUS, REFERENCE_MISMATCH, UNSUPPORTED, INVALID_INPUT }
const STATUS_NAMES := ["MAPPED", "TARGET_REFERENCE", "AMBIGUOUS", "REFERENCE_MISMATCH", "UNSUPPORTED", "INVALID_INPUT"]

const _UNSET := -3

var source: ReferenceSequence
var target: ReferenceSequence
var score := -1
var is_ready := false
var _map: PackedInt32Array  ## source index -> unique target index, or -1
var _supplied: String
var _cache := {}


func _init(source_ref: ReferenceSequence, target_ref: ReferenceSequence) -> void:
	source = source_ref
	target = target_ref
	_supplied = source.sequence + source.sequence.substr(0, source.circular_extension)


## Computes the alignment. Yields to the main loop every `yield_rows` rows
## when > 0 so a browser tab stays responsive. Returns "" or an error message.
func build(yield_rows: int = 0) -> String:
	if target.circular_extension:
		return "Target must be an unextended reference."
	var n := source.length()
	var m := target.length()
	if absi(n - m) >= BAND:
		return "Reference length difference exceeds the safe alignment band."
	_map = PackedInt32Array()
	_map.resize(n)
	if source.sequence == target.sequence:
		for i in n:
			_map[i] = i
		score = 0
		is_ready = true
		return ""

	var a := source.sequence.to_ascii_buffer()
	var b := target.sequence.to_ascii_buffer()
	var width := 2 * BAND + 1
	var inf := n + m + 1
	# Cell (i, j) lives at i * width + (j - i + BAND).
	var dp := PackedInt32Array()
	dp.resize((n + 1) * width)
	dp.fill(inf)
	for i in n + 1:
		var lo := maxi(0, i - BAND)
		var hi := mini(m, i + BAND)
		var base := i * width - i + BAND
		var pbase := base - width + 1
		var plo := maxi(0, i - 1 - BAND)
		var phi := mini(m, i - 1 + BAND)
		for j in range(lo, hi + 1):
			var value: int
			if i == 0:
				value = j
			elif j == 0:
				value = i
			else:
				var up := dp[pbase + j] + 1 if j >= plo and j <= phi else inf
				var diag := inf
				if j - 1 >= plo and j - 1 <= phi:
					diag = dp[pbase + j - 1] + (0 if a[i - 1] == b[j - 1] else 1)
				var left := dp[base + j - 1] + 1 if j > lo else inf
				value = mini(up, mini(diag, left))
			dp[base + j] = value
		if yield_rows > 0 and i % yield_rows == 0:
			await Engine.get_main_loop().process_frame
	score = dp[n * width - n + BAND + m]
	if score >= BAND:
		return "Cannot prove alignment completeness; the references differ too much."

	# Walk back through every optimal path, recording which target base(s)
	# each source base aligns to. Any gap or second option makes it ambiguous.
	var first := PackedInt32Array()
	first.resize(n)
	first.fill(_UNSET)
	var ambiguous := PackedByteArray()
	ambiguous.resize(n)
	var cur := PackedByteArray()
	cur.resize(m + 1)
	var nxt := PackedByteArray()
	nxt.resize(m + 1)
	cur[m] = 1
	for i in range(n, -1, -1):
		var lo := maxi(0, i - BAND)
		var hi := mini(m, i + BAND)
		var base := i * width - i + BAND
		var pbase := base - width + 1
		var plo := maxi(0, i - 1 - BAND)
		var phi := mini(m, i - 1 + BAND)
		for j in range(hi, lo - 1, -1):
			if not cur[j]:
				continue
			var value := dp[base + j]
			if j > lo and dp[base + j - 1] + 1 == value:
				cur[j - 1] = 1
			if i > 0:
				if j >= plo and j <= phi and dp[pbase + j] + 1 == value:
					ambiguous[i - 1] = 1  # source base deleted on this path
					nxt[j] = 1
				if j > 0 and j - 1 >= plo and j - 1 <= phi \
						and dp[pbase + j - 1] + (0 if a[i - 1] == b[j - 1] else 1) == value:
					if first[i - 1] == _UNSET:
						first[i - 1] = j - 1
					elif first[i - 1] != j - 1:
						ambiguous[i - 1] = 1
					nxt[j - 1] = 1
		for j in range(lo, hi + 1):
			cur[j] = 0
		var swap := cur
		cur = nxt
		nxt = swap
		if yield_rows > 0 and i % yield_rows == 0:
			await Engine.get_main_loop().process_frame
	for i in n:
		_map[i] = -1 if ambiguous[i] or first[i] < 0 else first[i]
	is_ready = true
	return ""


## Projects one allele. Returns {"status": Status, "detail": String,
## "position": int, "ref": String, "alt": String} (converted fields are set
## only for MAPPED). Results are memoized.
func map_variant(position: int, ref: String, alt: String) -> Dictionary:
	var key := "%d:%s:%s" % [position, ref, alt]
	if not _cache.has(key):
		_cache[key] = _map_variant(position, ref.to_upper(), alt.to_upper())
	return _cache[key]


func _map_variant(position: int, ref: String, alt: String) -> Dictionary:
	if position < 1:
		return _result(Status.INVALID_INPUT, "Expected a positive one-based position.")
	if ref.is_empty() or alt.is_empty() or ref == alt or not _is_acgt(ref) or not _is_acgt(alt):
		return _result(Status.UNSUPPORTED, "Expected distinct nonempty literal A/C/G/T alleles, one ALT at a time.")
	var seq := source.sequence
	var tgt := target.sequence
	var start := position - 1
	if start + ref.length() > _supplied.length() or _supplied.substr(start, ref.length()) != ref:
		return _result(Status.REFERENCE_MISMATCH, "Original REF does not match the declared source reference.")
	start = start % seq.length()
	if start + ref.length() > seq.length():
		return _result(Status.UNSUPPORTED, "Allele spans the circular origin.")
	# Remove unchanged VCF anchoring before projecting the biological edit.
	var deleted := ref
	var inserted := alt
	while deleted and inserted and deleted[0] == inserted[0]:
		start += 1
		deleted = deleted.substr(1)
		inserted = inserted.substr(1)
	while deleted and inserted and deleted[-1] == inserted[-1]:
		deleted = deleted.left(-1)
		inserted = inserted.left(-1)
	var end := start + deleted.length()
	var indices: Array = range(start, end) if deleted else [start - 1, start]
	for idx in indices:
		if idx < 0 or idx >= seq.length():
			return _result(Status.UNSUPPORTED, "Insertion at the circular origin is not supported.")
	var mapped: Array = []
	for idx in indices:
		if _map[idx] < 0:
			return _result(Status.AMBIGUOUS, "Edit touches a reference gap or nonunique optimal alignment.")
		mapped.append(_map[idx])
	for k in range(1, mapped.size()):
		if mapped[k] != mapped[k - 1] + 1:
			return _result(Status.AMBIGUOUS, "Edit spans a target-reference gap.")
	var left: int = mapped[0] if deleted else mapped[1]
	var right: int = mapped[-1] + 1 if deleted else left
	if "N" in tgt.substr(left, right - left):
		return _result(Status.UNSUPPORTED, "Target edit includes an unknown reference base.")
	var haplotype := tgt.substr(0, left) + inserted + tgt.substr(right)
	if haplotype == tgt:
		return _result(Status.TARGET_REFERENCE, "Projected allele equals the target reference; this is not a benign classification.")
	# Maximal common suffix first gives a deterministic leftmost edit in repeats.
	# Bases after the edit are identical by construction, so start past them.
	var tl := tgt.length()
	var hl := haplotype.length()
	var suffix := mini(tl - right, mini(tl, hl))
	while suffix < mini(tl, hl) and tgt.unicode_at(tl - 1 - suffix) == haplotype.unicode_at(hl - 1 - suffix):
		suffix += 1
	var te := tl - suffix
	var he := hl - suffix
	var prefix := mini(left, mini(te, he))
	while prefix < mini(te, he) and tgt.unicode_at(prefix) == haplotype.unicode_at(prefix):
		prefix += 1
	var r := tgt.substr(prefix, te - prefix)
	var a := haplotype.substr(prefix, he - prefix)
	if r.is_empty() or a.is_empty():
		if prefix:
			prefix -= 1
			r = tgt[prefix] + r
			a = tgt[prefix] + a
		elif suffix:
			r += tgt[te]
			a += haplotype[he]
		else:
			return _result(Status.UNSUPPORTED, "Cannot represent a whole-reference deletion as an anchored allele.")
	if tgt.substr(prefix, r.length()) != r or tgt.substr(0, prefix) + a + tgt.substr(prefix + r.length()) != haplotype:
		push_error("Converted allele failed haplotype reconstruction at %d %s>%s" % [position, ref, alt])
		return _result(Status.UNSUPPORTED, "Converted allele failed haplotype reconstruction.")
	var result := _result(Status.MAPPED, "Unique edit projection.")
	result["position"] = prefix + 1
	result["ref"] = r
	result["alt"] = a
	return result


static func _result(status: Status, detail: String) -> Dictionary:
	return {"status": status, "detail": detail}


static func _is_acgt(text: String) -> bool:
	for i in text.length():
		if not text[i] in "ACGT":
			return false
	return true
