# Mito Web

**Open it: https://udargit829.github.io/Mito-Web/**

Mitochondrial DNA variant analysis that runs entirely in your browser. Upload
VCF files (or a zip of them) and Mito Web places every call on the rCRS
reference, assigns PhyloTree 17 haplogroups (matching HaploGrep 3), and lets
you build analyses as a node graph: filters, comparisons between groups,
germline/somatic classes by allele fraction, gnomAD and ClinVar annotation,
mutation density along the genome with permutation tests, and patient-level
comparisons of clinical values (Kruskal-Wallis, Wilcoxon, Fisher's exact,
ANOVA, Spearman).

No account and nothing to install. To try it without your own data, open the
Files tab and press **Load example data**: the 1000 Genomes Project's
mitochondrial genotypes for 2,534 people, with population labels. Then try
**Analysis > Workflows…**, for example "Recreate haplogroup markers".

## Your data stays on your computer

There is no server. Files you upload are read in the browser and kept only in
that browser's storage. **Export workspace** saves a project to a file;
**New project** clears it. The only network requests are loading the site and,
if you press "Look up missing consequences" on the Tables tab, the positions of
the listed variants, sent to the Ensembl VEP service.

## Credits

- Example data: The 1000 Genomes Project Consortium, *A global reference for
  human genetic variation*, Nature 2015;526:68-74 (IGSR, release 20130502).
- Haplogroups: PhyloTree Build 17 (van Oven 2015) via
  [genepi/phylotree-rcrs-17](https://github.com/genepi/phylotree-rcrs-17) (MIT),
  scored with the HaploGrep 3 method (Schönherr et al., NAR 2023).
- Annotations: gnomAD v3.1 mtDNA (CC0; Laricchia et al., Genome Research 2022),
  ClinVar (public domain, NCBI), Ensembl VEP.
- Built with the [Godot Engine](https://godotengine.org) (MIT).

## For developers

The app is the Godot project in [`Godot_MitoWeb/`](Godot_MitoWeb/README.md),
which documents its features, file formats and layout. Each push to `main`
runs the tests and deploys the site with GitHub Actions.

This repository is public: never commit sequencing data, databases, clinical
tables or upload bundles. `Godot_MitoWeb/tools/check_public.py` enforces this
in CI.

## License

[MIT](LICENSE)
