"""Canonical domain models shared by the importer and web viewer."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import vcf_parser


@dataclass(frozen=True, order=True)
class AlleleKey:
    """Identity of one mutation allele.

    Position, reference, and alternate allele all participate in equality and
    hashing. Display values and sample-specific call data belong on other
    models rather than on the allele identity.
    """

    position: int
    ref: str
    alt: str

    def __post_init__(self) -> None:
        try:
            position = int(self.position)
        except (TypeError, ValueError) as exc:
            raise ValueError("Allele position must be an integer.") from exc

        ref = str(self.ref).strip().upper()
        alt = str(self.alt).strip().upper()
        if position < 1:
            raise ValueError("Allele position must be positive.")
        if not ref:
            raise ValueError("Allele reference cannot be empty.")
        if not alt:
            raise ValueError("Allele alternate cannot be empty.")

        object.__setattr__(self, "position", position)
        object.__setattr__(self, "ref", ref)
        object.__setattr__(self, "alt", alt)


@dataclass
class Sample:
    """One imported VCF sample and its population tags."""

    sample_id: str
    population: list[str]
    source_path: Path
    mutations: list[vcf_parser.VCFMutation] = field(default_factory=list)

    def __contains__(self, tags: str | list[str] | tuple[str, ...]) -> bool:
        """Return True when this sample has the given population tag(s)."""
        if isinstance(tags, str):
            return tags in self.population
        return all(tag in self.population for tag in tags)

    def __iter__(self):
        """Iterate over population tags."""
        return iter(self.population)

    @property
    def population_key(self) -> str:
        """Return a stable serialized population tag key."""
        return "|".join(self.population)

    @property
    def label(self) -> str:
        """Return a compact subject/population label."""
        if not self.population:
            return self.sample_id
        return f"{self.sample_id}_{'_'.join(self.population)}"

    def has_any(self, tags: list[str] | tuple[str, ...]) -> bool:
        """Return True when this sample has at least one of the given tags."""
        return any(tag in self for tag in tags)

    def has_all(self, tags: list[str] | tuple[str, ...]) -> bool:
        """Return True when this sample has every given tag."""
        return tags in self

    def is_subject(self, sample_id: str) -> bool:
        """Return True when this sample belongs to the given subject."""
        return self.sample_id == sample_id

    @property
    def mutation_alleles(self) -> set[AlleleKey]:
        """Return comparable mutation alleles for this sample."""
        return {
            AlleleKey(
                position=mutation.position,
                ref=mutation.ref,
                alt=alt,
            )
            for mutation in self.mutations
            for alt in mutation.alts
        }

    def __and__(self, other: "Sample") -> set[AlleleKey]:
        """Return mutation alleles common to both samples."""
        return self.mutation_alleles & other.mutation_alleles

    def __or__(self, other: "Sample") -> set[AlleleKey]:
        """Return mutation alleles present in either sample."""
        return self.mutation_alleles | other.mutation_alleles

    def __sub__(self, other: "Sample") -> set[AlleleKey]:
        """Return mutation alleles present in this sample but not the other."""
        return self.mutation_alleles - other.mutation_alleles

    def __xor__(self, other: "Sample") -> set[AlleleKey]:
        """Return mutation alleles that differ between two samples."""
        return self.mutation_alleles ^ other.mutation_alleles

    def common_mutations(self, *others: "Sample") -> set[AlleleKey]:
        """Return mutation alleles common to this sample and all others."""
        common = self.mutation_alleles
        for other in others:
            common &= other.mutation_alleles
        return common

    def different_mutations(self, *others: "Sample") -> set[AlleleKey]:
        """Return mutation alleles not shared by every provided sample."""
        samples = (self, *others)
        union: set[AlleleKey] = set()
        common = samples[0].mutation_alleles
        for sample in samples:
            union |= sample.mutation_alleles
            common &= sample.mutation_alleles
        return union - common


@dataclass
class SampleAlleleCall:
    """One sample's observed data for a canonical allele."""

    allele: AlleleKey
    sample_id: int | str
    label: str
    af: float | None
    af_text: str
    filter: str
    vcf_ref: str
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_json(self) -> dict[str, Any]:
        """Return the existing browser-facing representation."""
        return {
            "sample_id": self.sample_id,
            "label": self.label,
            "af": self.af,
            "af_text": self.af_text,
            "filter": self.filter,
            "vcf_ref": self.vcf_ref,
            "metadata": self.metadata,
        }


@dataclass
class DerivedSample:
    """In-memory sample materialized from a comparison result."""

    id: str
    label: str
    calls: list[SampleAlleleCall]
    mutations: list[vcf_parser.VCFMutation]
    source_description: str

    @property
    def subject_id(self) -> str:
        return "Derived"

    @property
    def population_key(self) -> str:
        return self.label

    @property
    def source_file(self) -> str:
        return self.source_description

    @property
    def mutation_count(self) -> int:
        return len({call.allele for call in self.calls})

    @property
    def vcf_iterator(self) -> vcf_parser.VCFIterator:
        return vcf_parser.VCFIterator.from_mutations(
            self.mutations,
            path=self.id,
        )

    def sample_row(self) -> dict[str, Any]:
        """Return the existing browser-facing sample representation."""
        return {
            "id": self.id,
            "subject_id": self.subject_id,
            "population_key": self.population_key,
            "source_file": self.source_file,
            "mutation_count": self.mutation_count,
            "is_derived": True,
        }


@dataclass
class SampleAlleleSet:
    """Allele set plus per-allele call details for one sample group."""

    calls_by_allele: dict[AlleleKey, list[SampleAlleleCall]] = field(
        default_factory=dict
    )

    def add(self, call: SampleAlleleCall) -> None:
        self.calls_by_allele.setdefault(call.allele, []).append(call)

    def __contains__(self, allele: AlleleKey) -> bool:
        return allele in self.calls_by_allele

    def __iter__(self):
        return iter(self.calls_by_allele)

    def __and__(self, other: "SampleAlleleSet") -> set[AlleleKey]:
        return set(self) & set(other)

    def __sub__(self, other: "SampleAlleleSet") -> set[AlleleKey]:
        return set(self) - set(other)

    def __or__(self, other: "SampleAlleleSet") -> set[AlleleKey]:
        return set(self) | set(other)

    def calls(self, allele: AlleleKey) -> list[SampleAlleleCall]:
        return self.calls_by_allele.get(allele, [])
