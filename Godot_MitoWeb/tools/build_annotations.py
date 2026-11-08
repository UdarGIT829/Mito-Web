#!/usr/bin/env python3
"""Build the public annotation bundle shipped with Mito Web.

Sources (downloaded unless local paths are given):
  gnomAD v3.1 mtDNA sites VCF  (56,434 genomes; CC0)
  ClinVar GRCh38 VCF, MT records (NCBI; public domain)

Every allele is normalized onto rCRS exactly as the app's MitoMapper does
for rCRS input (trim VCF anchors, then the leftmost anchored representation),
so bundle keys match imported calls. Nothing from any user data set is used.

Usage:
  python3 tools/build_annotations.py                  # download and build
  python3 tools/build_annotations.py --gnomad g.vcf.bgz --clinvar clinvar.vcf.gz
Output: annotations/mito_annotations.json.gz
"""
from __future__ import annotations

import argparse
import datetime
import gzip
import io
import json
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GNOMAD_URL = "https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz"
CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
HAP_ORDER = ["A", "B", "C", "D", "E", "F", "G", "H", "HV", "I", "J", "K", "L0", "L1", "L2", "L3", "L4",
             "L5", "M", "N", "P", "R", "T", "U", "V", "W", "X", "Y", "Z"]
HAP_AF_MIN = 0.01  # per-haplogroup homoplasmic AF kept when at least this


def read_rcrs() -> str:
    lines = (ROOT / "reference" / "rcrs_chrM.fa").read_text().splitlines()
    seq = "".join(l.strip() for l in lines if not l.startswith(">")).upper()
    assert len(seq) == 16569
    return seq


def normalize(rcrs: str, position: int, ref: str, alt: str):
    """Port of MitoMapper._map_variant for source == target == rCRS.
    Returns (position, ref, alt) or None when the allele cannot be placed."""
    ref, alt = ref.upper(), alt.upper()
    if position < 1 or not ref or not alt or ref == alt or set(ref + alt) - set("ACGT"):
        return None
    start = position - 1
    if start + len(ref) > len(rcrs) or rcrs[start:start + len(ref)] != ref:
        return None
    deleted, inserted = ref, alt
    while deleted and inserted and deleted[0] == inserted[0]:
        start += 1
        deleted, inserted = deleted[1:], inserted[1:]
    while deleted and inserted and deleted[-1] == inserted[-1]:
        deleted, inserted = deleted[:-1], inserted[:-1]
    if not deleted and (start - 1 < 0 or start >= len(rcrs)):
        return None
    left = start
    right = start + len(deleted)
    if "N" in rcrs[left:right]:
        return None
    hap = rcrs[:left] + inserted + rcrs[right:]
    if hap == rcrs:
        return None
    tl, hl = len(rcrs), len(hap)
    suffix = min(tl - right, tl, hl)
    while suffix < min(tl, hl) and rcrs[tl - 1 - suffix] == hap[hl - 1 - suffix]:
        suffix += 1
    te, he = tl - suffix, hl - suffix
    prefix = min(left, te, he)
    while prefix < min(te, he) and rcrs[prefix] == hap[prefix]:
        prefix += 1
    r, a = rcrs[prefix:te], hap[prefix:he]
    if not r or not a:
        if prefix:
            prefix -= 1
            r, a = rcrs[prefix] + r, rcrs[prefix] + a
        elif suffix:
            r, a = r + rcrs[te], a + hap[he]
        else:
            return None
    return prefix + 1, r, a


def open_text(path_or_url: str) -> io.TextIOBase:
    if path_or_url.startswith("http"):
        print(f"downloading {path_or_url}")
        data = urllib.request.urlopen(path_or_url).read()
    else:
        data = Path(path_or_url).read_bytes()
    if data[:2] == b"\x1f\x8b":
        # gzip.decompress handles the multiple members of BGZF files.
        data = gzip.decompress(data)
    return io.StringIO(data.decode("utf-8"))


def info_dict(text: str) -> dict:
    out = {}
    for part in text.split(";"):
        key, _, value = part.partition("=")
        out[key] = value if _ else True
    return out


def protein_change(amino_acids: str, position: str) -> str:
    """'L/R' at 156 -> 'L156R'; synonymous 'L' at 243 -> 'L243='."""
    if "/" in amino_acids:
        before, after = amino_acids.split("/", 1)
        return f"{before}{position}{after}"
    return f"{amino_acids}{position}="


def num(value, digits=6):
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return float(f"{x:.{digits}g}")


def read_gnomad(source: str, rcrs: str, records: dict) -> dict:
    meta = {"name": "gnomAD v3.1 mtDNA", "samples": 56434, "url": GNOMAD_URL,
            "license": "CC0 1.0", "citation": "Laricchia et al. 2022, Genome Research 32:569-582"}
    vep_fields = []
    kept = skipped = 0
    for line in open_text(source):
        if line.startswith("##INFO=<ID=vep"):
            vep_fields = line.split("Format: ")[1].rstrip('">\n').split("|")
        if line.startswith("#"):
            continue
        chrom, pos, _id, ref, alt, _q, filt, info = line.rstrip("\n").split("\t")[:8]
        key = normalize(rcrs, int(pos), ref, alt)
        if key is None:
            skipped += 1
            continue
        i = info_dict(info)
        consequences, genes, protein = [], [], []
        for entry in str(i.get("vep", "")).split(","):
            parts = dict(zip(vep_fields, entry.split("|")))
            for term in parts.get("Consequence", "").split("&"):
                if term and term not in consequences:
                    consequences.append(term)
            if parts.get("SYMBOL") and parts["SYMBOL"] not in genes:
                genes.append(parts["SYMBOL"])
            if parts.get("Amino_acids") and parts.get("Protein_position"):
                change = f"{parts['SYMBOL']} p.{protein_change(parts['Amino_acids'], parts['Protein_position'])}"
                if change not in protein:
                    protein.append(change)
        haps = []
        if "hap_AF_hom" in i:
            for name, af in zip(HAP_ORDER, str(i["hap_AF_hom"]).split("|")):
                value = num(af, 4)
                if value is not None and value >= HAP_AF_MIN:
                    haps.append([name, value])
            haps.sort(key=lambda h: -h[1])
        g = {
            "af_hom": num(i.get("AF_hom")), "af_het": num(i.get("AF_het")),
            "ac_hom": int(i.get("AC_hom", 0)), "ac_het": int(i.get("AC_het", 0)), "an": int(i.get("AN", 0)),
            "max_hl": num(i.get("max_hl"), 3), "filters": "" if filt == "PASS" else filt,
            "hap_defining": "hap_defining_variant" in i, "haplogroups": haps,
            "consequence": consequences, "genes": genes, "protein": protein,
        }
        for k in ("mitotip_trna_prediction", "pon_mt_trna_prediction"):
            if i.get(k):
                g[k.replace("_trna_prediction", "")] = str(i[k]).replace("_", " ")
        records.setdefault("%d:%s:%s" % key, {})["gnomad"] = g
        kept += 1
    meta.update({"records": kept, "not_placed": skipped})
    return meta


def read_clinvar(source: str, rcrs: str, records: dict) -> dict:
    meta = {"name": "ClinVar", "url": CLINVAR_URL, "license": "Public domain (NCBI)",
            "citation": "Landrum et al., ClinVar, Nucleic Acids Research"}
    kept = skipped = 0
    for line in open_text(source):
        if line.startswith("##fileDate="):
            meta["release"] = line.split("=", 1)[1].strip()
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if fields[0] not in ("MT", "chrM"):
            continue
        key = normalize(rcrs, int(fields[1]), fields[3], fields[4])
        if key is None:
            skipped += 1
            continue
        i = info_dict(fields[7])
        clean = lambda v: str(v or "").replace("_", " ").replace("|", "; ")
        entry = {
            "variation_id": fields[2],
            "significance": clean(i.get("CLNSIG")),
            "review": clean(i.get("CLNREVSTAT")),
            "condition": "; ".join(c.replace("_", " ") for c in str(i.get("CLNDN") or "").split("|")
                                   if c and c not in ("not_provided", ".")),
        }
        bucket = records.setdefault("%d:%s:%s" % key, {})
        if "clinvar" in bucket:  # rare: two records normalize to one allele
            bucket["clinvar"]["significance"] += " / " + entry["significance"]
        else:
            bucket["clinvar"] = entry
        kept += 1
    meta.update({"records": kept, "not_placed": skipped})
    return meta


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gnomad", default=GNOMAD_URL)
    parser.add_argument("--clinvar", default=CLINVAR_URL)
    parser.add_argument("--out", default=str(ROOT / "annotations" / "mito_annotations.json.gz"))
    args = parser.parse_args()
    rcrs = read_rcrs()
    records: dict = {}
    sources = {"gnomad": read_gnomad(args.gnomad, rcrs, records), "clinvar": read_clinvar(args.clinvar, rcrs, records)}
    bundle = {
        "format": "mito-web-annotations", "version": 1,
        "built": datetime.date.today().isoformat(),
        "reference": "rCRS NC_012920.1", "normalization": "mitochondrial-unit-edit-all-optimal-v1 (rCRS identity)",
        "haplogroup_order": HAP_ORDER, "sources": sources, "records": records,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    data = json.dumps(bundle, separators=(",", ":"), sort_keys=True).encode()
    out.write_bytes(gzip.compress(data, compresslevel=9, mtime=0))
    print(f"{len(records)} alleles; gnomAD {sources['gnomad']['records']}, ClinVar {sources['clinvar']['records']}; "
          f"{len(data) / 1e6:.1f} MB JSON, {out.stat().st_size / 1e6:.2f} MB gzipped -> {out}")


if __name__ == "__main__":
    main()
