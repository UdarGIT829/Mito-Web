#!/usr/bin/env python3
"""Build examples/1000_genomes_chrMT.zip, the example data the app loads.

Public data from the 1000 Genomes Project (IGSR), phase 3 release 20130502:
  ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz
  integrated_call_samples_v3.20130502.ALL.panel
Each sample is tagged 1000G, its continental group and its population.

Usage: python3 tools/build_1000g_example.py [--vcf FILE] [--panel FILE]
Without --vcf/--panel the files are downloaded from the IGSR FTP site.
"""
import argparse
import io
import pathlib
import urllib.request
import zipfile

BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
VCF = "ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz"
PANEL = "integrated_call_samples_v3.20130502.ALL.panel"
OUT = pathlib.Path(__file__).resolve().parent.parent / "examples" / "1000_genomes_chrMT.zip"
FIXED_TIME = (2026, 1, 1, 0, 0, 0)  # reproducible zip


def fetch(path, name):
    if path:
        return pathlib.Path(path).read_bytes()
    with urllib.request.urlopen(BASE + name) as response:
        return response.read()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf")
    parser.add_argument("--panel")
    args = parser.parse_args()
    vcf = fetch(args.vcf, VCF)
    panel = fetch(args.panel, PANEL).decode()
    lines = [
        "# 1000 Genomes phase 3 chrMT genotypes (public; ftp.1000genomes.ebi.ac.uk, release 20130502).",
        "# Tags: 1000G, continental group (AFR, AMR, EAS, EUR, SAS) and population code, from %s." % PANEL,
        "file,sample,subject,build,tags",
    ]
    for row in panel.splitlines()[1:]:
        parts = row.split("\t")
        if len(parts) < 3 or not parts[0]:
            continue
        sample, pop, super_pop = parts[0], parts[1], parts[2]
        lines.append(",".join([VCF, sample, sample, "rcrs", "1000G;%s;%s" % (super_pop, pop)]))
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as z:
        for name, data in [(VCF, vcf), ("manifest.csv", ("\n".join(lines) + "\n").encode())]:
            info = zipfile.ZipInfo("1000_genomes_chrMT/" + name, FIXED_TIME)
            info.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(info, data)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_bytes(buffer.getvalue())
    print("%d samples tagged -> %s (%d kB)" % (len(lines) - 3, OUT, len(buffer.getvalue()) // 1000))


if __name__ == "__main__":
    main()
