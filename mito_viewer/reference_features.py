"""Pinned mitochondrial gene ranges for fast, reproducible locus lookup.

The coordinates below are the 37 ``gene`` features from the NCBI RefSeq
NC_012920.1 GFF3 record. Coordinates are 1-based and inclusive. Official HGNC
symbols are used for display while NCBI Gene IDs provide stable identifiers.

Source:
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
    ?db=nuccore&id=NC_012920.1&rettype=gff3&retmode=text
"""

from __future__ import annotations

from dataclasses import asdict, dataclass


MITO_REFERENCE_ACCESSION = "NC_012920.1"
MITO_REFERENCE_LENGTH = 16569
MITO_REFERENCE_SOURCE = "NCBI RefSeq"
MITO_REFERENCE_URL = (
    "https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1"
)


@dataclass(frozen=True)
class MitoGene:
    symbol: str
    ncbi_gene_id: str
    hgnc_id: str
    start: int
    end: int
    strand: str
    biotype: str

    def to_dict(self) -> dict[str, str | int]:
        return asdict(self)


MITO_GENES = (
    MitoGene("MT-TF", "4558", "HGNC:7481", 577, 647, "+", "tRNA"),
    MitoGene("MT-RNR1", "4549", "HGNC:7470", 648, 1601, "+", "rRNA"),
    MitoGene("MT-TV", "4577", "HGNC:7500", 1602, 1670, "+", "tRNA"),
    MitoGene("MT-RNR2", "4550", "HGNC:7471", 1671, 3229, "+", "rRNA"),
    MitoGene("MT-TL1", "4567", "HGNC:7490", 3230, 3304, "+", "tRNA"),
    MitoGene("MT-ND1", "4535", "HGNC:7455", 3307, 4262, "+", "protein_coding"),
    MitoGene("MT-TI", "4565", "HGNC:7488", 4263, 4331, "+", "tRNA"),
    MitoGene("MT-TQ", "4572", "HGNC:7495", 4329, 4400, "-", "tRNA"),
    MitoGene("MT-TM", "4569", "HGNC:7492", 4402, 4469, "+", "tRNA"),
    MitoGene("MT-ND2", "4536", "HGNC:7456", 4470, 5511, "+", "protein_coding"),
    MitoGene("MT-TW", "4578", "HGNC:7501", 5512, 5579, "+", "tRNA"),
    MitoGene("MT-TA", "4553", "HGNC:7475", 5587, 5655, "-", "tRNA"),
    MitoGene("MT-TN", "4570", "HGNC:7493", 5657, 5729, "-", "tRNA"),
    MitoGene("MT-TC", "4511", "HGNC:7477", 5761, 5826, "-", "tRNA"),
    MitoGene("MT-TY", "4579", "HGNC:7502", 5826, 5891, "-", "tRNA"),
    MitoGene("MT-CO1", "4512", "HGNC:7419", 5904, 7445, "+", "protein_coding"),
    MitoGene("MT-TS1", "4574", "HGNC:7497", 7446, 7514, "-", "tRNA"),
    MitoGene("MT-TD", "4555", "HGNC:7478", 7518, 7585, "+", "tRNA"),
    MitoGene("MT-CO2", "4513", "HGNC:7421", 7586, 8269, "+", "protein_coding"),
    MitoGene("MT-TK", "4566", "HGNC:7489", 8295, 8364, "+", "tRNA"),
    MitoGene("MT-ATP8", "4509", "HGNC:7415", 8366, 8572, "+", "protein_coding"),
    MitoGene("MT-ATP6", "4508", "HGNC:7414", 8527, 9207, "+", "protein_coding"),
    MitoGene("MT-CO3", "4514", "HGNC:7422", 9207, 9990, "+", "protein_coding"),
    MitoGene("MT-TG", "4563", "HGNC:7486", 9991, 10058, "+", "tRNA"),
    MitoGene("MT-ND3", "4537", "HGNC:7458", 10059, 10404, "+", "protein_coding"),
    MitoGene("MT-TR", "4573", "HGNC:7496", 10405, 10469, "+", "tRNA"),
    MitoGene("MT-ND4L", "4539", "HGNC:7460", 10470, 10766, "+", "protein_coding"),
    MitoGene("MT-ND4", "4538", "HGNC:7459", 10760, 12137, "+", "protein_coding"),
    MitoGene("MT-TH", "4564", "HGNC:7487", 12138, 12206, "+", "tRNA"),
    MitoGene("MT-TS2", "4575", "HGNC:7498", 12207, 12265, "+", "tRNA"),
    MitoGene("MT-TL2", "4568", "HGNC:7491", 12266, 12336, "+", "tRNA"),
    MitoGene("MT-ND5", "4540", "HGNC:7461", 12337, 14148, "+", "protein_coding"),
    MitoGene("MT-ND6", "4541", "HGNC:7462", 14149, 14673, "-", "protein_coding"),
    MitoGene("MT-TE", "4556", "HGNC:7479", 14674, 14742, "-", "tRNA"),
    MitoGene("MT-CYB", "4519", "HGNC:7427", 14747, 15887, "+", "protein_coding"),
    MitoGene("MT-TT", "4576", "HGNC:7499", 15888, 15953, "+", "tRNA"),
    MitoGene("MT-TP", "4571", "HGNC:7494", 15956, 16023, "-", "tRNA"),
)


def _circular_reference_intervals(
    position: int,
    length: int,
) -> tuple[tuple[int, int], ...]:
    if length > MITO_REFERENCE_LENGTH:
        raise ValueError("Reference allele cannot exceed the mitochondrial genome.")
    end = position + length - 1
    if end <= MITO_REFERENCE_LENGTH:
        return ((position, end),)
    return (
        (position, MITO_REFERENCE_LENGTH),
        (1, end - MITO_REFERENCE_LENGTH),
    )


def mitochondrial_gene_locus(
    position: int,
    ref: str = "",
) -> dict:
    """Return genes overlapped by one NC_012920.1 reference interval."""
    try:
        position = int(position)
    except (TypeError, ValueError) as exc:
        raise ValueError("Mitochondrial position must be an integer.") from exc
    if not 1 <= position <= MITO_REFERENCE_LENGTH:
        raise ValueError(
            f"Mitochondrial position must be between 1 and "
            f"{MITO_REFERENCE_LENGTH}."
        )

    ref = str(ref or "").strip().upper()
    intervals = _circular_reference_intervals(position, max(len(ref), 1))
    genes = [
        gene
        for gene in MITO_GENES
        if any(
            interval_start <= gene.end and interval_end >= gene.start
            for interval_start, interval_end in intervals
        )
    ]
    in_control_region = any(
        interval_start <= 576 or interval_end >= 16024
        for interval_start, interval_end in intervals
    )
    region = (
        "genic"
        if genes
        else "control_region"
        if in_control_region
        else "intergenic"
    )
    region_label = {
        "genic": "Gene",
        "control_region": "Control region / D-loop",
        "intergenic": "Intergenic",
    }[region]

    return {
        "reference": {
            "accession": MITO_REFERENCE_ACCESSION,
            "length": MITO_REFERENCE_LENGTH,
            "source": MITO_REFERENCE_SOURCE,
            "url": MITO_REFERENCE_URL,
        },
        "query": {
            "position": position,
            "ref": ref,
            "intervals": [
                {"start": start, "end": end}
                for start, end in intervals
            ],
        },
        "region": region,
        "region_label": region_label,
        "genes": [gene.to_dict() for gene in genes],
    }
