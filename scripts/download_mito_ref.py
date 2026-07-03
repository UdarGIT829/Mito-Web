#!/usr/bin/env python3
"""Download the hg38 mitochondrial reference FASTA from UCSC."""

import argparse
import gzip
from pathlib import Path
from urllib.request import urlretrieve


DEFAULT_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/"
    "chromosomes/chrM.fa.gz"
)
DEFAULT_OUTPUT = Path(__file__).resolve().parent.parent / "reference" / "hg38_chrM.fa"
EXPECTED_LENGTH = 16569


def read_fasta_sequence(path):
    sequence_parts = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            sequence_parts.append(line.upper())
    return "".join(sequence_parts)


def download_mito_reference(url=DEFAULT_URL, output_path=DEFAULT_OUTPUT):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    compressed_path = output_path.with_suffix(output_path.suffix + ".gz")
    urlretrieve(url, compressed_path)

    with gzip.open(compressed_path, "rt") as source:
        with output_path.open("w") as destination:
            for line in source:
                destination.write(line)

    sequence = read_fasta_sequence(output_path)
    if len(sequence) != EXPECTED_LENGTH:
        raise ValueError(
            f"Expected chrM length {EXPECTED_LENGTH}, got {len(sequence)} "
            f"from {output_path}"
        )

    return output_path


def main():
    parser = argparse.ArgumentParser(
        description="Download UCSC hg38 chrM FASTA for master CSV ref bases."
    )
    parser.add_argument("--url", default=DEFAULT_URL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    output_path = download_mito_reference(args.url, args.output)
    print(output_path)


if __name__ == "__main__":
    main()
