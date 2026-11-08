class_name EnsemblClient
extends Node
## Looks up consequences for alleles the bundle does not cover, using the
## Ensembl VEP REST API (which allows requests from web pages). Requests are
## batched (200 alleles each) and sent one at a time.

signal progress(done: int, total: int)
signal finished(results: Dictionary, error: String)

## distance=0: no upstream/downstream calls, which on the compact
## mitochondrial genome would list nearly every neighboring gene.
const URL := "https://rest.ensembl.org/vep/human/region?distance=0"
const BATCH := 200
const PAUSE_SECONDS := 1.0

var busy := false
var _http: HTTPRequest


func _ready() -> void:
	_http = HTTPRequest.new()
	_http.timeout = 60
	add_child(_http)


## `alleles`: Array of [position, ref, alt] on rCRS. Emits `finished` with
## {key: {"consequence": [...], "genes": [...], "protein": [...], ...}}.
func lookup(alleles: Array) -> void:
	if busy:
		return
	busy = true
	var results := {}
	var error := ""
	for start in range(0, alleles.size(), BATCH):
		var batch := alleles.slice(start, start + BATCH)
		var body := JSON.stringify({"variants": batch.map(func(a): return to_vep_input(a[0], a[1], a[2]))})
		var sent := _http.request(URL, ["Content-Type: application/json", "Accept: application/json"], HTTPClient.METHOD_POST, body)
		if sent != OK:
			error = "Could not start the request (%s)." % error_string(sent)
			break
		var response: Array = await _http.request_completed
		var code: int = response[1]
		if response[0] != HTTPRequest.RESULT_SUCCESS or code != 200:
			error = "Ensembl did not answer (HTTP %d)%s." % [code, "; it may be rate limiting, try again shortly" if code == 429 else ""]
			break
		var parsed = JSON.parse_string(response[3].get_string_from_utf8())
		if typeof(parsed) != TYPE_ARRAY:
			error = "Ensembl returned an unexpected response."
			break
		results.merge(parse_response(parsed, batch), true)
		progress.emit(mini(start + BATCH, alleles.size()), alleles.size())
		if start + BATCH < alleles.size():
			await get_tree().create_timer(PAUSE_SECONDS).timeout
	busy = false
	finished.emit(results, error)


## VCF-style region input; the "MT" contig is rCRS in Ensembl (GRCh38).
static func to_vep_input(position: int, ref: String, alt: String) -> String:
	return "MT %d . %s %s . . ." % [position, ref, alt]


## Maps VEP output back to requested alleles by their input string.
static func parse_response(entries: Array, batch: Array) -> Dictionary:
	var by_input := {}
	for a in batch:
		by_input[to_vep_input(a[0], a[1], a[2])] = AnnotationDb.key(a[0], a[1], a[2])
	var date := Time.get_date_string_from_system(true)
	var out := {}
	for entry in entries:
		if typeof(entry) != TYPE_DICTIONARY or not by_input.has(entry.get("input")):
			continue
		var consequence: Array = []
		var genes: Array = []
		var protein: Array = []
		var most_severe: String = entry.get("most_severe_consequence", "")
		if most_severe:
			consequence.append(most_severe)
		for tc in entry.get("transcript_consequences", []):
			for term in tc.get("consequence_terms", []):
				if not term in consequence:
					consequence.append(term)
			var symbol: String = tc.get("gene_symbol", "")
			if symbol and not symbol in genes:
				genes.append(symbol)
			if tc.get("amino_acids") and tc.get("protein_start"):
				var change := "%s p.%s" % [symbol, protein_change(tc["amino_acids"], int(tc["protein_start"]))]
				if not change in protein:
					protein.append(change)
		out[by_input[entry["input"]]] = {"consequence": consequence, "genes": genes, "protein": protein,
			"source": "Ensembl VEP REST", "fetched": date}
	for a in batch:
		var k := AnnotationDb.key(a[0], a[1], a[2])
		if not out.has(k):
			out[k] = {"consequence": [], "genes": [], "protein": [], "source": "Ensembl VEP REST", "fetched": date, "note": "no result"}
	return out


## "L/R" at 156 -> "L156R"; synonymous "L" at 243 -> "L243=".
static func protein_change(amino_acids: String, position: int) -> String:
	if "/" in amino_acids:
		var parts := amino_acids.split("/", true, 1)
		return "%s%d%s" % [parts[0], position, parts[1]]
	return "%s%d=" % [amino_acids, position]
