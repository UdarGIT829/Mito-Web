"""Domain values for durable derived samples and immutable runs."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

from mito_viewer.domain import AFRule, AlleleKey, MetadataFilter, MutationFilters

from .models import Visibility


DERIVED_COMPARISON_STATUSES = frozenset(
    {"common", "partial", "unique", "__none__"}
)
DERIVED_SAMPLE_STATUSES = frozenset({"present", "unique", "not_in"})


class DerivedInputKind(str, Enum):
    GROUP = "group"
    SAMPLE = "sample"
    RUN = "run"


class PresenceRequirement(str, Enum):
    ANY = "any"
    ALL = "all"
    NONE = "none"
    EXACTLY_ONE = "exactly_one"


class DerivedRunStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


class EvidenceStatus(str, Enum):
    QUALIFYING_PRESENT = "qualifying_present"
    FILTERED_OUT = "filtered_out"
    NOT_OBSERVED = "not_observed"


@dataclass(frozen=True)
class DerivedInput:
    """One group, observed sample, or exact parent run clause."""

    kind: DerivedInputKind
    input_id: str
    role: str
    requirement: PresenceRequirement

    def __post_init__(self) -> None:
        input_id = str(self.input_id or "").strip()
        role = str(self.role or "").strip()
        if not input_id:
            raise ValueError("Derived input ID cannot be empty.")
        if not role:
            raise ValueError("Derived input role cannot be empty.")
        object.__setattr__(self, "kind", DerivedInputKind(self.kind))
        object.__setattr__(self, "input_id", input_id)
        object.__setattr__(self, "role", role)
        object.__setattr__(
            self,
            "requirement",
            PresenceRequirement(self.requirement),
        )

    def to_dict(self) -> dict[str, str]:
        return {
            "kind": self.kind.value,
            "input_id": self.input_id,
            "role": self.role,
            "requirement": self.requirement.value,
        }


@dataclass(frozen=True)
class DerivedComparison:
    """The exact viewer comparison predicate attached to a durable run."""

    statuses: frozenset[str]
    input_statuses: tuple[frozenset[str], ...]

    def __post_init__(self) -> None:
        statuses = frozenset(str(status).strip() for status in self.statuses)
        unknown_statuses = statuses - DERIVED_COMPARISON_STATUSES
        if not statuses:
            raise ValueError(
                "Durable comparison requires at least one global status."
            )
        if unknown_statuses:
            raise ValueError(
                "Unknown durable comparison status: "
                f"{sorted(unknown_statuses)[0]}."
            )

        input_statuses = tuple(
            frozenset(str(status).strip() for status in allowed)
            for allowed in self.input_statuses
        )
        for allowed in input_statuses:
            unknown = allowed - DERIVED_SAMPLE_STATUSES
            if unknown:
                raise ValueError(
                    "Unknown durable sample comparison status: "
                    f"{sorted(unknown)[0]}."
                )

        object.__setattr__(self, "statuses", statuses)
        object.__setattr__(self, "input_statuses", input_statuses)

    def to_dict(self) -> dict[str, Any]:
        return {
            "statuses": sorted(self.statuses),
            "input_statuses": [
                sorted(allowed) for allowed in self.input_statuses
            ],
        }

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> "DerivedComparison":
        if not isinstance(value, dict):
            raise ValueError("Durable comparison definition must be an object.")
        return cls(
            statuses=frozenset(value.get("statuses") or ()),
            input_statuses=tuple(
                frozenset(allowed)
                for allowed in value.get("input_statuses") or ()
            ),
        )


@dataclass(frozen=True)
class DerivedDefinition:
    """Normalized logical definition for one derived calculation."""

    inputs: tuple[DerivedInput, ...]
    filters: MutationFilters = field(default_factory=MutationFilters)
    version: int = 1
    comparison: DerivedComparison | None = None

    def __post_init__(self) -> None:
        inputs = tuple(
            item if isinstance(item, DerivedInput) else DerivedInput(**item)
            for item in self.inputs
        )
        if not inputs:
            raise ValueError("Derived definition requires at least one input.")
        identities = [(item.kind, item.input_id) for item in inputs]
        if len(set(identities)) != len(identities):
            raise ValueError("Derived definition inputs must be unique.")
        if all(item.requirement is PresenceRequirement.NONE for item in inputs):
            raise ValueError("Derived definition requires at least one positive input.")
        filters = (
            self.filters
            if isinstance(self.filters, MutationFilters)
            else MutationFilters(**self.filters)
        )
        version = int(self.version)
        comparison = self.comparison
        if comparison is not None and not isinstance(
            comparison,
            DerivedComparison,
        ):
            comparison = DerivedComparison.from_dict(comparison)
        if version not in {1, 2}:
            raise ValueError(f"Unsupported derived definition version: {version}")
        if version == 1 and comparison is not None:
            raise ValueError(
                "Version 1 derived definitions cannot contain a comparison."
            )
        if version == 2 and comparison is None:
            raise ValueError(
                "Version 2 derived definitions require a comparison."
            )
        if (
            comparison is not None
            and len(comparison.input_statuses) != len(inputs)
        ):
            raise ValueError(
                "Durable comparison requires one status set per input."
            )
        if comparison is not None and len(inputs) < 2:
            raise ValueError(
                "Durable comparison requires at least two inputs."
            )
        object.__setattr__(self, "inputs", inputs)
        object.__setattr__(self, "filters", filters)
        object.__setattr__(self, "version", version)
        object.__setattr__(self, "comparison", comparison)

    def to_dict(self) -> dict[str, Any]:
        value = {
            "version": self.version,
            "inputs": [item.to_dict() for item in self.inputs],
            "filters": {
                "position": self.filters.position,
                "alt": self.filters.alt,
                "af_rules": [
                    {
                        "operator": rule.operator,
                        "threshold": rule.threshold,
                    }
                    for rule in self.filters.af_rules
                ],
                "metadata_filters": [
                    {
                        "field": item.field,
                        "value": item.value,
                    }
                    for item in self.filters.metadata_filters
                ],
            },
        }
        if self.comparison is not None:
            value["comparison"] = self.comparison.to_dict()
        return value

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> "DerivedDefinition":
        filters = value.get("filters") or {}
        comparison = value.get("comparison")
        return cls(
            version=value.get("version", 2 if comparison is not None else 1),
            inputs=tuple(DerivedInput(**item) for item in value.get("inputs") or ()),
            filters=MutationFilters(
                position=filters.get("position"),
                alt=filters.get("alt", ""),
                af_rules=tuple(
                    AFRule(
                        item["operator"],
                        item["threshold"],
                    )
                    for item in filters.get("af_rules") or ()
                ),
                metadata_filters=tuple(
                    MetadataFilter(
                        item["field"],
                        item["value"],
                    )
                    for item in filters.get("metadata_filters") or ()
                ),
            ),
            comparison=(
                DerivedComparison.from_dict(comparison)
                if comparison is not None
                else None
            ),
        )

    @property
    def normalized_json(self) -> str:
        return _canonical_json(self.to_dict())

    @property
    def fingerprint(self) -> str:
        return _sha256_text(self.normalized_json)


@dataclass(frozen=True)
class DerivedSampleRecord:
    id: str
    catalog_sample_id: str
    perspective_id: str
    dataset_id: str
    name: str
    description: str
    visibility: Visibility
    current_run_id: str | None
    created_at: str
    updated_at: str


@dataclass(frozen=True)
class DerivedRunRecord:
    id: str
    derived_sample_id: str
    definition: DerivedDefinition
    input_snapshot_fingerprint: str
    reference_id: str
    normalization_version: str
    status: DerivedRunStatus
    output_count: int
    output_fingerprint: str
    error_text: str
    created_at: str
    completed_at: str | None


@dataclass(frozen=True)
class MaterializedAllele:
    id: str
    run_id: str
    allele: AlleleKey
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class AlleleEvidence:
    id: str
    derived_allele_id: str
    source_catalog_sample_id: str | None
    parent_derived_allele_id: str | None
    original_catalog_sample_id: str | None
    source_mutation_id: str | None
    source_alt_index: int | None
    input_role: str
    evaluation_status: EvidenceStatus
    af: float | None
    af_text: str
    filter_text: str
    metadata: dict[str, Any]
    lineage: dict[str, Any]


@dataclass(frozen=True)
class DerivedCalculationResult:
    derived_sample: DerivedSampleRecord
    run: DerivedRunRecord
    alleles: tuple[MaterializedAllele, ...]


@dataclass(frozen=True)
class DerivedStaleness:
    run_id: str
    stale_group_ids: tuple[str, ...] = ()
    missing_group_ids: tuple[str, ...] = ()
    changed_sample_ids: tuple[str, ...] = ()
    changed_source_cohort_ids: tuple[str, ...] = ()
    parent_runs_with_updates: tuple[str, ...] = ()

    @property
    def stale(self) -> bool:
        return any(
            (
                self.stale_group_ids,
                self.missing_group_ids,
                self.changed_sample_ids,
                self.changed_source_cohort_ids,
                self.parent_runs_with_updates,
            )
        )


@dataclass(frozen=True)
class ResolvedSampleInput:
    clause_index: int
    catalog_sample_id: str
    source_group_id: str | None
    group_membership_fingerprint: str
    sample_fingerprint: str
    source_database_fingerprint: str
    input_role: str
    display_order: int


@dataclass(frozen=True)
class ResolvedParentInput:
    clause_index: int
    parent_run_id: str
    source_group_id: str | None
    group_membership_fingerprint: str
    input_role: str
    display_order: int


@dataclass(frozen=True)
class PendingAlleleEvidence:
    source_catalog_sample_id: str | None
    parent_derived_allele_id: str | None
    original_catalog_sample_id: str | None
    source_mutation_id: str | None
    source_alt_index: int | None
    input_role: str
    evaluation_status: EvidenceStatus
    af: float | None = None
    af_text: str = ""
    filter_text: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)
    lineage: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class PendingMaterializedAllele:
    allele: AlleleKey
    metadata: dict[str, Any] = field(default_factory=dict)
    evidence: tuple[PendingAlleleEvidence, ...] = ()


def allele_output_fingerprint(alleles: list[AlleleKey] | tuple[AlleleKey, ...]) -> str:
    payload = [
        [allele.position, allele.ref, allele.alt] for allele in sorted(set(alleles))
    ]
    return _sha256_text(_canonical_json(payload))


def snapshot_fingerprint(value: Any) -> str:
    return _sha256_text(_canonical_json(value))


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    )


def _sha256_text(value: str) -> str:
    return "sha256:" + hashlib.sha256(value.encode("utf-8")).hexdigest()
