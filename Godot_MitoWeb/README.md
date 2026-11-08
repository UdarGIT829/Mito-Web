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
| `scripts/data/` | Columnar `VariantStore`, VCF importer, zip/manifest uploads, tag and window imports |
| `scripts/engine/` | Analysis engine: sample selection, `MutationFrame`, filters, comparison, table rows, gene loci, patient tests (`GroupTests`, `PatientCompare`) |
| `scripts/graph/` | Analysis graph data (`GraphModel`), cached evaluation (`GraphEvaluator`), workflows (`Workflows`) |
| `scripts/table/` | `TableViewModel`: columns, column filters, sorting, paging, export |
| `scripts/ui/` | Files, Samples, Analysis (GraphEdit), Tables, Statistics and Patients tabs; dialogs |
| `scripts/annotations/` | Annotation bundle lookup (`AnnotationDb`) and live Ensembl VEP client |
| `annotations/` | Public annotation bundle (`tools/build_annotations.py`) and haplogroup tree (`tools/build_phylotree.py`) |
| `reference/` | rCRS (NC_012920.1) and UCSC hg19 chrM FASTA |
| `tests/` | Headless tests; `fixtures/mapper_cases.json` is generated from the Python mapper |

## Uploading: zips and manifests

Upload VCFs, or a `.zip` of them. A CSV/TSV whose name contains `manifest`
(inside the zip, or uploaded next to the VCFs) sets each file's import
settings, so files arrive ready to import:

```csv
# lines starting with # are comments
file,sample,subject,build,tags
AP4_EV_PCR.vcf,,AP4,rcrs,EV;PCR
PL20-19-A_S16.vcf,,PL-20-19,hg19_plus58,lupusN;A
ALL.chrMT.phase3.vcf.gz,HG00096,HG00096,rcrs,1000G;EUR;GBR
```

`sample` selects a column of a multi-sample VCF (blank = all samples);
`build` accepts rcrs/hg38/GRCh37, hg19, hg19_plus58 (hg19+58); `tags` are
separated by `;`.

## Duplicates

- Same content: a file identical to one in the project once decompressed
  (e.g. `.vcf` vs `.vcf.gz`, CRLF line endings) is not added again.
- Same sample: a single-sample file whose VCF sample column matches another
  file's (e.g. a strict and a filtered call set of one sample) is marked
  Duplicate, with the share of calls in common. Generic column names (e.g.
  `16`) also need half of the calls shared. Multi-sample files are flagged
  when most sample columns match.
- Imported files count as the originals. "Import all ready" skips
  duplicates ("Import selected" imports them deliberately), and "Remove
  duplicates" removes the later copies with their imported samples.

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
| Germline / Somatic | Classifies every call by AF (see below); keeps only calls of chosen classes, or labels them |
| Annotation Filter | Keeps or removes mutations by ClinVar category/review, gnomAD frequency, haplogroup marker, consequence, tRNA prediction |
| Enrich with Annotations | Adds gnomAD, ClinVar and live Ensembl annotation columns to the tables after it |
| Table | Names a point in the graph to view on the Tables tab |
| Statistics | Mutation density along the genome per connected group; view on the Statistics tab |
| Compare Patients | Per-patient values (clinical scores) compared across groups; view on the Patients tab |

Results are cached per node; an edit recomputes only that node and what is
downstream of it.

Every node can be given a name and a description: double-click its title bar
or press "Edit". The description is shown on the node and in its tooltip. A
named group labels its series on the Statistics tab, and named Table and
Statistics nodes are listed by that name when they have no table name.
Names and descriptions are saved with the workspace and in workflows.

## Workflows

Preset analyses, from "Workflows…" on the Analysis tab. A workflow asks for
its settings, then inserts its graph (or replaces the current one).

- **Treated vs untreated** (group tags A and B): exclusive and shared
  mutations, annotated, with density statistics including all somatic calls.
- **Recreate haplogroup markers** (haplogroup, whole clade or major
  haplogroup only, common %, rare %;
  assigns haplogroups first if no sample has them):
  *gained* = carried by >= common% of the haplogroup and < rare% of everyone
  else; *absent* = the reverse (how H, the rCRS haplogroup, is defined: it
  lacks 2706G and 7028T). Both sets are flagged against PhyloTree markers.
  On 1000 Genomes (defaults: whole clade, 90% / 5%) all 27 clades
  present find markers, 175 in all, and every one is a PhyloTree marker:
  e.g. J (295T, 12612G, 13708A, 16069T + JT markers), U (11467G, 12308G,
  12372A), M (489C, 10400T, 14783C, 15043A), B (the 8281-8289 9-bp
  deletion), and as absent markers H (2706G, 7028T), HV (11719A, 14766T),
  R (12705T, 16223T), N (8701G, 9540C, 10873C) and L3 (769A, 1018A, ...);
  ancestral clades show up as "absent" because rCRS itself is H, inside R,
  N and L3. With "major haplogroup only" ancestral groups find little,
  because their descendants count as everyone else. Recurrent markers (e.g.
  10398G) are excluded.

"Save current graph as workflow…" keeps the graph in the project;
saved workflows can be exported as `.mitoflow.json` and imported elsewhere.

## Tables tab

Shows a Table node's mutations, 100–1000 rows per page. Click a column title
to sort or filter it (number ranges, text contains, genes, flag yes/no, or
per sample: carries / does not carry / AF range). Column filters only change
the view, not the graph. Up to 60 samples get one AF column each (optional up
to 300); larger sets use a compact samples column. Flag columns from a
Logical Filter in flag mode come with a "why not" column. Double-click a row
for its calls: AF, FILTER, the allele as written in the VCF, placement on
rCRS, repeats and reference context. Copy (TSV) and Download CSV export
every filtered row, not just the current page.

## Germline / somatic classification

Rule `af-quarters-v1`, applied to each sample's call (never to a mutation as
a whole): AF >= 0.75 germline, 0.50 <= AF < 0.75 indeterminate, AF < 0.50
somatic, no AF unknown. These are fixed quarters of the AF scale, not
quantiles of the data: mitochondrial heteroplasmy is bimodal (one hump near
1, one near 0), so the middle band is nearly empty and calls there are
flagged rather than labelled. Cut-offs are editable per node. Keeping only
somatic calls removes germline carriers from a mutation; mutations left with
no carriers drop out, while the sample universe (presence denominator) stays.

## Annotation rules

Rule `annotation-rules-v1`. ClinVar wording is grouped into pathogenic,
benign, uncertain, conflicting (conflicting or mixed), other (risk factor,
affects, association, ...) and absent (no record; never treated as benign).
Review status maps to ClinVar's stars (0-4).

## Tags, values and haplogroups from a table

Samples tab, "Import from a table": a CSV/TSV with a header row. One column
names the sample (matched to sample name, then subject, so one row can cover
all of a patient's samples). Other columns can become:

- **values**: per-sample attributes such as clinical scores (SLEDAI 8, SDI 2,
  Sex F). Numbers stay numbers. Used by Compare Patients nodes.
- **tags**: labels for grouping (optionally `column=value`).
- **a haplogroup set by hand** (one column): overrides the classifier, e.g.
  after resolving mixed-up samples. Names resolve to the nearest PhyloTree 17
  haplogroup (`C1c*4` -> C1c).

Tables with numeric columns default to values, text-only tables to tags.
Values and tags are also written into the files' import settings, so they
survive Undo import and re-import. "Set haplogroup…" sets or clears a
haplogroup by hand for the selected samples (or a patient's row);
"Assign haplogroups" leaves hand-set ones alone.

## Lineages

Samples tab, "Lineages…": rules `clade = lineage`, one per line. The most
specific rule on a sample's path to the root wins, so `N = Asian` with
`A2 = Amerindian` puts A2 in Amerindian and other N in Asian. Samples get
`lineage=` tags (re-applied whenever haplogroups change). The built-in preset
is the pSLE 2026 paper's: African = L; Amerindian = A2, B2, C1b, C1c, C1d,
C4c, D1, D4h3a, X2a; Caucasian = H, HV, J, T, V, U, K, W, X, I; Asian = the
rest of M and N.

## Compare Patients and the Patients tab

The Compare Patients node compares connected groups patient by patient (a
subject's samples count once) on imported values. Each measure is one row:

| Summary | Cell | Omnibus test | Pairwise |
|---|---|---|---|
| Median | median [25th, 75th] | Kruskal-Wallis (tie-corrected) | Wilcoxon rank-sum |
| Mean | mean (SD) | one-way ANOVA | Wilcoxon rank-sum |
| Percent where (e.g. `>= 6`, `= F`, `= yes`) | n (%) | Fisher's exact (r x c) | Fisher's exact |
| Correlation with another value | Spearman rho per group | overall Spearman | |

Tests follow Stata 18's defaults (Wilcoxon uses the normal approximation with
tie correction, as `ranksum`; percentiles as `summarize, detail`), and are
checked against the paper and scipy in `tests/patient_tests.gd`. The Patients
tab lists every patient's haplogroup, lineage and values (disagreeing samples
highlighted) and shows each node's results table, exportable as CSV.

The **Lupus lineages (pSLE 2026)** workflow recreates Udar et al., J Clin Med
2026;15:86, Tables 2-5 plus duration-vs-activity correlations, from
haplogroups and a clinical table with columns Sex, Age, Age_at_Dx, Duration,
SLEDAI, SDI, LN_II, LN_III_IV, LN_V, CKD_ESRD. On the 25 patients of the
paper's Table 1 it reproduces every value that does not depend on the 13th
Amerindian patient missing from that table.

### Lupus expansions and quality checks

- **Computed values** for Compare Patients: per-patient mutation counts from
  the group's own mutations (all; somatic AF < 0.50; somatic with a 5% floor
  that drops amplicon noise; indeterminate; germline), and the patient's
  haplogroup, major haplogroup and lineage, all usable like imported values.
- **Lupus: mutations and outcomes** workflow: mutation counts by lineage,
  their Spearman correlation with SLEDAI-2K and SDI, and the somatic
  mutations found only in Amerindian or only in African patients, annotated.
- **Cohort vs 1000 Genomes** workflow: lineage shares in the cohort against
  1000 Genomes MXL (Mexican Ancestry, Los Angeles), the other admixed
  American populations and the rest, in non-overlapping groups.
- **Mixture check** (with Assign haplogroups, for every sample): a major
  profile (homoplasmic + AF >= 0.5) and a minor profile (homoplasmic + AF <
  0.5) are classified separately, in the spirit of haplocheck. A minor
  haplogroup at least 2 nodes away, supported by 3+ low-AF variants it
  expects and the major one does not, is flagged as possible contamination
  or a mix-up. A haplogroup set by hand whose major haplogroup differs from
  the file's own calls is flagged too ("file looks like ..."), which catches
  swapped files. Shown on the Samples and Patients tabs.
- Fisher's exact test switches to a seeded Monte Carlo (20,000 tables) when
  exact enumeration would exceed 3 million tables.

## Statistics (density windows)

A Statistics node takes any number of groups (each connected input) and a
window set:

- **Sliding windows**: size and step in bp (step < size overlaps), limited to
  a range (e.g. one gene), optionally wrapping past 16,569 to 1.
- **Genes and regions**: the 37 genes, the control region (16024–576, across
  the origin) and the intergenic stretches.

- **Imported list**: a CSV/TSV with start and end columns (1-based; start >
  end wraps the origin) or a BED file, or windows kept from a scan. Per window and group the Statistics tab shows distinct
mutations, carrier calls, or AF-weighted calls (heteroplasmy load), optionally
per kb and per sample (divided by the group's own sample count; use "keep
only these samples" on tag filters so each group's count is its own). Chart
and table, with CSV export. Short regions (tRNAs ~70 bp) give noisy per-kb
values.

## Feature matrix and permutation tests

"Download matrix" writes one row per sample (all groups) with subject, group,
every tag as a 0/1 column (labels for models) and the chosen count per
window: the input for models built outside this tool.

"Compare within" tests one group's samples split by two tags (A vs B), per
window: difference in mean (A minus B), permutation p, and family-wise
adjusted p (single-step max-T). Comparing sets that were built from the same
tags (e.g. Treated-exclusive vs Untreated-exclusive) would be circular, so the
test always splits one group. Units are subjects (samples of a subject
averaged; set subjects with "Set subject from column" when importing tags) or
samples. When all labellings fit in the permutation budget they are
enumerated exactly, and the smallest possible p is reported: with 3 vs 1
cell lines it is 0.25. The seed is recorded for reproducibility.

"Scan for windows" tests 100/200/500/1000/2000 bp windows at quarter-size
steps around the genome, adjusting for the whole scan, and can keep the best
non-overlapping windows as the node's window list.

## Haplogroups

Samples tab, "Assign haplogroups": every sample is matched against PhyloTree
Build 17 (5,435 haplogroups, rCRS orientation) with HaploGrep's weighted
Kulczynski score, using its variant weights, hotspots (309.1C, 315.1C,
16182C, ...) and indel alignment rules. Calls count when AF >= the chosen
cut-off (default 0.75, the germline cut-off) or when there is no AF. Results
are stored on the sample (score, runner-up, missing and extra variants) and as
tags `haplogroup=J1c2`, `macro haplogroup=J`, and `haplogroup clade=` for
every major haplogroup it belongs to (J, R, N, L3), usable in tag filters.
`macro` is the nearest major haplogroup only; `clade` includes descendants,
so `haplogroup clade=R` selects H, J, T, U, K... samples too.

Tree data: github.com/genepi/phylotree-rcrs-17 (MIT), from van Oven M.
(2015) PhyloTree Build 17, Forensic Sci Int Genet Suppl 5; method from
Schönherr S et al. (2023) Haplogrep 3, Nucleic Acids Res 51:W263-W268.

Checked against HaploGrep 3.3.2 itself (same tree, phylotree-rcrs@17.3) on
all 2,534 1000 Genomes samples: identical haplogroup for every sample, score
within 0.0005 for 99.4%. Native rCRS calls are read from the VCF the way
HaploGrep reads them (positions as given; a shorter ALT deletes the tail of
REF); calls mapped from hg19 are shifted to PhyloTree's indel positions.
Also: 99.8% of HaploGrep's 5,426 canonical profiles classify as themselves.

## Annotations

`annotations/mito_annotations.json.gz` (~0.35 MB) is built only from public
bulk releases, never from user data:

| Source | Fields | License |
|---|---|---|
| gnomAD v3.1 mtDNA, 56,434 genomes (Laricchia et al. 2022) | AF hom/het, max heteroplasmy, PhyloTree haplogroup-marker flag, per-haplogroup AF hom, VEP consequence and protein change, MitoTIP and PON-mt-tRNA | CC0 1.0 |
| ClinVar GRCh38 VCF, MT records | significance, review status, condition, variation ID (release date stored) | Public domain (NCBI) |

Alleles are normalized exactly like imported calls, which a test checks for
every bundle key. Refresh with `python3 tools/build_annotations.py`.

For mutations gnomAD does not cover, the Tables tab can ask the Ensembl VEP
REST API (`distance=0`, batches of 200); results are saved in the workspace.
MITOMAP is not included: it has no bulk download and its service does not
accept requests from web pages. Absence from a source is never treated as a
classification.

## Commands

From this folder, with `godot` pointing at Godot 4.7.1:

```bash
godot --headless --script res://tests/run_tests.gd               # tests
godot --headless --export-release "Web" ../build/web/index.html  # web build
python3 -m http.server -d ../build/web 8000                      # try it locally
```

Pushing changes under `Godot_MitoWeb/` to `main` runs
`.github/workflows/pages.yml`, which tests, exports and deploys to Pages.
