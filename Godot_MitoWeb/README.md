# Mito Web

Mitochondrial variant analysis that runs entirely in the browser. It is built
with Godot 4.7 (GDScript, single-threaded web export) and hosted on GitHub
Pages. There is no server: uploaded files and projects stay in the user's
browser (`user://`, which is backed by IndexedDB). The exported
`*.mitoweb.json` workspace file is the durable save.

The Python code elsewhere in this repository (`mito_viewer/`, `main.py`,
`web_viewer.py`) is the reference implementation being ported.

## Layout

| Path | Purpose |
|---|---|
| `scripts/workspace_model.gd` | Workspace data and the versioned `.mitoweb` file format |
| `scripts/autoload/workspace.gd` | Open workspace, imports, autosave to `user://` |
| `scripts/autoload/web_files.gd` | Browser upload/download bridge; desktop file-dialog fallback |
| `scripts/io/` | gzip/BGZF decompression, VCF reading |
| `scripts/genome/` | Bundled references, rCRS mapper (port of `mito_viewer/mapping`), repeat profiles |
| `scripts/data/` | Columnar `VariantStore` and the VCF importer |
| `scripts/engine/` | Analysis engine: sample selection, `MutationFrame`, filters, comparison, table rows, gene loci |
| `scripts/graph/` | Analysis graph data (`GraphModel`) and cached evaluation (`GraphEvaluator`) |
| `scripts/ui/` | Files, Samples and Analysis (GraphEdit) tabs; filter and table dialogs |
| `reference/` | rCRS (NC_012920.1) and UCSC hg19 chrM FASTA |
| `tests/` | Headless tests; `fixtures/mapper_cases.json` is generated from the Python mapper |

## Import model

Every ALT allele of every sample becomes one call. Calls keep their original
VCF position/REF/ALT, AF (`AF`, else `VF`, else `HF`), FILTER and raw sample
value. Each call is projected onto rCRS with the same conservative mapper as
the Python code. Native rCRS calls are trimmed and left-normalized the same
way, so identities agree across builds. Calls that cannot be placed
(ambiguous, reference mismatch, unsupported such as the `N` at 3107) are kept
with their status and excluded from analysis.

The reference build comes from the `##contig` length (16569 rCRS, 16571
hg19, 16629 hg19 + 58). When a file does not state it, the user must choose;
REF agreement is shown as evidence only. In multi-sample VCFs a sample has a
call only when its genotype contains that ALT.

## Analysis engine

A `MutationFrame` is the table that flows between analysis steps: one row per
rCRS mutation carried by a sample in its universe, with the call each carrier
contributes. It is built with `MutationFrame.from_samples()` from samples
chosen by `SampleQuery` (tags all/any/none, `<NONE>`, subjects, files).

- `MutationFilters`: position, genes, rCRS reference base, ALT, slippage
  metadata, presence percentage, and AF ranges with a carrier percentage.
  Ported from the Python Table pipeline; percentages compare counts, not
  rounded values. `apply(..., "flag")` keeps rows and adds a pass/fail column.
- `Comparison`: Python's common/partial/unique statuses and per-sample
  present / not_in / unique rules.
- Set operations: `intersect`, `subtract`, `union`.
- `MutationTable.build()`: display rows with genes, region, slippage
  footprint, per-sample AF, sorting and paging.

One intentional difference from the Python Table: the `polymorphism`
metadata filter is applied (the Python Table accepted it but ignored it).

## Analysis graph

The Analysis tab is a node graph saved in the workspace. Every connection
carries a `MutationFrame`.

| Node | Does |
|---|---|
| Data Provider | Samples (all, or by tags) and their mutations; optional PASS-only calls |
| Filter on Tag | Mutations with a carrier having any of the tags, or every tag across carriers; can also narrow the samples |
| Logical Comparison | Include inputs (in every / in any), minus every Exclude input |
| Logical Filter | `MutationFilters` settings; removes failing rows or adds a pass/fail column |
| Enrich with Annotations | Pass-through until annotations are added |
| Table | Views, copies or downloads the mutations at that point |

Results are cached per node; an edit recomputes only that node and what is
downstream of it. "Insert example" builds the treated-vs-untreated design.

## Commands

From this folder, with `godot` pointing at Godot 4.7.1:

```bash
godot --headless --script res://tests/run_tests.gd               # tests
godot --headless --export-release "Web" ../build/web/index.html  # web build
python3 -m http.server -d ../build/web 8000                      # try it locally
```

Pushing changes under `Godot_MitoWeb/` to `main` runs
`.github/workflows/pages.yml`, which tests, exports and deploys to Pages.
