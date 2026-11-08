#!/usr/bin/env python3
"""Fails if the repository would publish data it must not.

This repository is public and the site is built from it. Lab data
(sequencing files, databases, clinical tables) and lab-identifying names must
never be committed. Run from the repository root; CI runs it on every push.

Checks every tracked file:
  1. Data-type files are allowed only where listed in ALLOWED.
  2. The example zip holds only the 1000 Genomes chrMT VCF (pinned by hash)
     and a manifest naming 1000 Genomes samples.
  3. No VCF records pasted into other files.
  4. No lab-identifying names: words are hashed and compared to DENY_HASHES
     (hashes, so this file does not list the names), plus a generic pattern
     for patient-style IDs.
"""
import hashlib
import io
import re
import subprocess
import sys
import zipfile

DATA_EXTENSIONS = (".vcf", ".vcf.gz", ".bgz", ".bam", ".bai", ".cram", ".fastq", ".fastq.gz", ".fq.gz",
                   ".sqlite", ".db", ".xlsx", ".xls", ".csv", ".tsv", ".pdf", ".zip", ".gz", ".mitoweb")
ALLOWED = [
    re.compile(r"^Godot_MitoWeb/examples/1000_genomes_chrMT\.zip$"),
    re.compile(r"^Godot_MitoWeb/annotations/[\w.-]+\.json\.gz$"),
]
EXAMPLE_ZIP = "Godot_MitoWeb/examples/1000_genomes_chrMT.zip"
EXAMPLE_VCF = "1000_genomes_chrMT/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz"
EXAMPLE_VCF_SHA256 = "ab1a8264c49c15728905d60dd11c881697ffcfe88f0437cfe01ad07a4e74c726"
EXAMPLE_MANIFEST = "1000_genomes_chrMT/manifest.csv"
DENY_HASHES = {'a0ba39df9d55a7479252d169cc5c3ce23e769bfcc2b1180b57ef02a66766d21b', '5a17fad1d7ee4c46676572afe4e5275b829d9e1133bcd9555f41e5ce9da3d256', '3d934c1eae590079406d8b3fbb39e267abff73a804e3be3e6f52a74dc0fc936c', '71fce443333775896c1d0c48d181ad4fe99221433c0a2210dea227c48315b33b', 'bf849cf0fa5a8074eed0166fbb850a151de91d373058b74d1faec6804fa70998', '844fcaa6f06c50a312a0e2c91c44fe9787a6f536ec646f759cc313e0fa5eaa35', 'fe93ba24aa758ff8068366b215327b1374f72af89f1926d2d280a23e57741a19', '1f311134efe1f98d5a9e049b42d7c73c3f9d825abe3319e16e21cfa65fa368b7', '63a98318c415bf5b00c9d867e4109a9d9caf8238e450fa4b810457b3f9412248', '4c5652c8553711737553a48da8207d6d3c67434f5739d442b66378c6884369f1', '47085116d81e750b96901bf33c5d1139395243d13eb8dee085ee71e48987f584', }
PATIENT_ID = re.compile(r"\bPL-?2\d-?\d\d\b", re.IGNORECASE)
VCF_RECORD = re.compile(r"^(chrM|chrMT|MT|M)\t\d+\t", re.MULTILINE)
CODE_SUFFIXES = (".gd", ".py", ".json")  # tests build synthetic VCF text in code


def main():
    files = subprocess.run(["git", "ls-files", "-z"], capture_output=True, check=True).stdout.decode().split("\0")
    problems = []
    for path in filter(None, files):
        lower = path.lower()
        if lower.endswith(DATA_EXTENSIONS) and not any(p.match(path) for p in ALLOWED):
            problems.append("%s: data files may not be committed (only %s)" % (path, ", ".join(p.pattern for p in ALLOWED)))
            continue
        raw = open(path, "rb").read()
        if path == EXAMPLE_ZIP:
            problems.extend(check_example(raw))
            continue
        if b"\0" in raw[:8192]:
            continue  # binary (icons, compressed annotation bundles)
        text = raw.decode("utf-8", "replace")
        if not lower.endswith(CODE_SUFFIXES) and VCF_RECORD.search(text):
            problems.append("%s: contains VCF records" % path)
        for word in set(re.findall(r"[A-Za-z0-9]+", text)):
            if hashlib.sha256(word.lower().encode()).hexdigest() in DENY_HASHES:
                problems.append("%s: contains a lab-identifying name" % path)
                break
        if PATIENT_ID.search(text):
            problems.append("%s: contains a patient-style ID" % path)
    for p in problems:
        print("NOT PUBLIC: " + p)
    print("check_public: %d tracked files, %d problem(s)" % (len([f for f in files if f]), len(problems)))
    return 1 if problems else 0


def check_example(raw):
    try:
        z = zipfile.ZipFile(io.BytesIO(raw))
    except zipfile.BadZipFile:
        return ["%s: not a zip" % EXAMPLE_ZIP]
    names = sorted(n for n in z.namelist() if not n.endswith("/"))
    if names != sorted([EXAMPLE_VCF, EXAMPLE_MANIFEST]):
        return ["%s: may only contain the 1000 Genomes VCF and its manifest (has %s)" % (EXAMPLE_ZIP, names)]
    out = []
    if hashlib.sha256(z.read(EXAMPLE_VCF)).hexdigest() != EXAMPLE_VCF_SHA256:
        out.append("%s: the VCF is not the 1000 Genomes release file" % EXAMPLE_ZIP)
    for line in z.read(EXAMPLE_MANIFEST).decode().splitlines():
        if line.startswith("#") or line.startswith("file,"):
            continue
        sample = line.split(",")[1] if line.count(",") >= 4 else ""
        if not re.fullmatch(r"(HG|NA)\d{5}", sample):
            out.append("%s: manifest row is not a 1000 Genomes sample: %s" % (EXAMPLE_ZIP, line[:60]))
            break
    return out


if __name__ == "__main__":
    sys.exit(main())
