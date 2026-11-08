#!/usr/bin/env python3
"""Build the haplogroup tree bundle shipped with Mito Web.

Source: genepi/phylotree-rcrs-17 (MIT license): PhyloTree Build 17 in rCRS
orientation, with HaploGrep's variant weights, hotspots and alignment rules.
  van Oven M. (2015) PhyloTree Build 17. Forensic Sci Int Genet Suppl 5.
  Schönherr S et al. (2023) Haplogrep 3. Nucleic Acids Res 51:W263-W268.

The XML tree is rooted at H2a2a1 (the rCRS sequence). Each edge lists the
variants, relative to rCRS, gained on the way to a haplogroup; a trailing "!"
is a back mutation (an even number of "!" re-applies the variant). Here every
haplogroup's full expected profile is resolved and stored as additions and
removals relative to its parent, so the app can score a sample against all
haplogroups in one pass down the tree.

Usage: python3 tools/build_phylotree.py [--src DIR]   (default: download)
Output: annotations/phylotree17.json.gz
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPO = "https://raw.githubusercontent.com/genepi/phylotree-rcrs-17/HEAD/src/"
FILES = ["tree.xml", "weights.txt", "tree.yaml", "rules.csv"]


def fetch(src: Path | None, name: str) -> str:
    if src:
        return (src / name).read_text()
    print("downloading", REPO + name)
    return urllib.request.urlopen(REPO + name).read().decode()


def position_key(poly: str) -> str:
    """Site of a variant: '263G' -> '263', '309.1C' -> '309.1', '523d' -> '523'."""
    m = re.match(r"^(\d+(?:\.\d+)?)", poly)
    return m.group(1) if m else poly


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--src", type=Path, help="folder holding tree.xml, weights.txt, tree.yaml, rules.csv")
    parser.add_argument("--out", type=Path, default=ROOT / "annotations" / "phylotree17.json.gz")
    args = parser.parse_args()
    texts = {name: fetch(args.src, name) for name in FILES}

    yaml = texts["tree.yaml"]
    hotspots = re.findall(r'^\s*-\s*"([^"]+)"', yaml.split("hotspots:")[1].split("genes:")[0], re.M)
    clusters = re.findall(r'^\s*-\s*label:\s*"([^"]+)"', yaml.split("clusters:")[1], re.M)
    # Clusters may name extra member nodes (e.g. B -> B4'5, B4, B5, B6) when
    # the tree has no node with the cluster's own name.
    cluster_nodes = {}
    for block in re.split(r'^\s*-\s*label:', yaml.split("clusters:")[1].split("template:")[0], flags=re.M)[1:]:
        label = re.match(r'\s*"([^"]+)"', block).group(1)
        if "nodes:" in block:
            members = block.split("nodes:")[1].split("frequencies:")[0]
            for node in re.findall(r'-\s*"([^"]+)"', members):
                cluster_nodes[node] = label
    version = re.search(r"^version:\s*(\S+)", yaml, re.M).group(1)

    weights = {}
    for line in texts["weights.txt"].splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            weights[parts[0].strip().upper()] = float(parts[1])

    rules = []
    for line in texts["rules.csv"].splitlines()[1:]:
        if "," in line:
            error, expected = line.split(",", 1)
            rules.append([error.split(), expected.split()])

    root = ET.fromstring(texts["tree.xml"])
    hot = {h.upper() for h in hotspots}
    names, parents, profiles = [], [], []

    def visit(node, parent_index: int, profile: dict) -> None:
        for h in node.findall("haplogroup"):
            current = dict(profile)
            for p in h.findall("details/poly"):
                raw = (p.text or "").strip().upper()
                if not raw:
                    continue
                bangs = len(raw) - len(raw.rstrip("!"))
                poly = raw.rstrip("!")
                key = position_key(poly)
                if bangs % 2 == 1:
                    if current.get(key) == poly:
                        del current[key]
                else:
                    current[key] = poly
            index = len(names)
            names.append(h.get("name"))
            parents.append(parent_index)
            profiles.append({k: v for k, v in current.items() if v not in hot})
            visit(h, index, current)

    visit(root, -1, {})

    poly_ids: dict[str, int] = {}
    def pid(poly: str) -> int:
        if poly not in poly_ids:
            poly_ids[poly] = len(poly_ids)
        return poly_ids[poly]

    nodes = []
    for i, name in enumerate(names):
        own = set(profiles[i].values())
        parent = set(profiles[parents[i]].values()) if parents[i] >= 0 else set()
        nodes.append([name, parents[i], sorted(pid(p) for p in own - parent), sorted(pid(p) for p in parent - own)])
    polys = sorted(poly_ids, key=poly_ids.get)
    bundle = {
        "format": "mito-web-phylotree", "version": 1, "tree": "PhyloTree Build 17 (Haplogrep tree %s)" % version,
        "reference": "rCRS NC_012920.1",
        "sources": [
            "van Oven M. (2015). PhyloTree Build 17: Growing the human mitochondrial DNA tree. Forensic Sci Int Genet Suppl Ser 5.",
            "Schönherr S, Weissensteiner H, Kronenberg F, Forer L. Haplogrep 3. Nucleic Acids Res 51:W263-W268 (2023).",
            "Tree, weights, hotspots and alignment rules: github.com/genepi/phylotree-rcrs-17 (MIT license).",
        ],
        "hotspots": hotspots, "clusters": clusters, "cluster_nodes": cluster_nodes, "rules": rules,
        "polys": polys, "weights": [weights.get(p, 1.0) for p in polys],
        "nodes": nodes,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    data = json.dumps(bundle, separators=(",", ":"), ensure_ascii=False).encode()
    args.out.write_bytes(gzip.compress(data, compresslevel=9, mtime=0))
    missing = sum(1 for p in polys if p not in weights)
    print(f"{len(nodes)} haplogroups, {len(polys)} variants ({missing} without a weight -> 1.0), "
          f"{len(hotspots)} hotspots, {len(rules)} rules; {args.out.stat().st_size / 1e3:.0f} kB -> {args.out}")


if __name__ == "__main__":
    main()
