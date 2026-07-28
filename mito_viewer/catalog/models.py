"""Domain models for the writable multi-perspective catalog."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from enum import Enum
from typing import Any


class Visibility(str, Enum):
    PRIVATE = "private"
    PERSPECTIVE = "perspective"
    SHARED = "shared"
    PUBLIC = "public"


class CohortType(str, Enum):
    SOURCE = "source"
    DERIVED = "derived"


class AccessLevel(str, Enum):
    VIEW = "view"
    USE = "use"
    MANAGE = "manage"


class SampleType(str, Enum):
    OBSERVED = "observed"
    DERIVED = "derived"


def _required_text(value: object, label: str) -> str:
    text = str(value or "").strip()
    if not text:
        raise ValueError(f"{label} cannot be empty.")
    return text


def _optional_text(value: object) -> str:
    return str(value or "").strip()


@dataclass(frozen=True)
class StudyPerspective:
    """Editable workspace that owns datasets, groups, and derived results."""

    id: str
    name: str
    description: str = ""
    visibility: Visibility = Visibility.PRIVATE
    created_at: str = ""
    updated_at: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "id", _required_text(self.id, "Perspective ID"))
        object.__setattr__(
            self,
            "name",
            _required_text(self.name, "Perspective name"),
        )
        object.__setattr__(
            self,
            "description",
            _optional_text(self.description),
        )
        object.__setattr__(self, "visibility", Visibility(self.visibility))


@dataclass(frozen=True)
class Cohort:
    """Reusable collection of observed samples or derived results."""

    id: str
    name: str
    cohort_type: CohortType
    owner_perspective_id: str | None = None
    description: str = ""
    source_database_identifier: str = ""
    source_database_fingerprint: str = ""
    reference_id: str = ""
    normalization_version: str = ""
    provenance: dict[str, Any] = field(default_factory=dict)
    visibility: Visibility = Visibility.PRIVATE
    created_at: str = ""
    updated_at: str = ""

    def __post_init__(self) -> None:
        cohort_type = CohortType(self.cohort_type)
        owner_id = _optional_text(self.owner_perspective_id) or None
        source_identifier = _optional_text(self.source_database_identifier)
        source_fingerprint = _optional_text(self.source_database_fingerprint)

        if cohort_type is CohortType.SOURCE:
            if not source_identifier:
                raise ValueError(
                    "Source cohort requires a database identifier."
                )
            if not source_fingerprint:
                raise ValueError(
                    "Source cohort requires a database fingerprint."
                )
        if cohort_type is CohortType.DERIVED and owner_id is None:
            raise ValueError(
                "Derived cohort requires an owning perspective."
            )

        object.__setattr__(self, "id", _required_text(self.id, "Cohort ID"))
        object.__setattr__(self, "name", _required_text(self.name, "Cohort name"))
        object.__setattr__(self, "cohort_type", cohort_type)
        object.__setattr__(self, "owner_perspective_id", owner_id)
        object.__setattr__(
            self,
            "description",
            _optional_text(self.description),
        )
        object.__setattr__(
            self,
            "source_database_identifier",
            source_identifier,
        )
        object.__setattr__(
            self,
            "source_database_fingerprint",
            source_fingerprint,
        )
        object.__setattr__(self, "reference_id", _optional_text(self.reference_id))
        object.__setattr__(
            self,
            "normalization_version",
            _optional_text(self.normalization_version),
        )
        object.__setattr__(self, "provenance", dict(self.provenance or {}))
        object.__setattr__(self, "visibility", Visibility(self.visibility))


@dataclass(frozen=True)
class Dataset:
    """Perspective-owned analysis assembly of accessible cohorts."""

    id: str
    perspective_id: str
    name: str
    description: str = ""
    derived_results_cohort_id: str | None = None
    visibility: Visibility = Visibility.PRIVATE
    created_at: str = ""
    updated_at: str = ""

    def __post_init__(self) -> None:
        derived_id = _optional_text(self.derived_results_cohort_id) or None
        object.__setattr__(self, "id", _required_text(self.id, "Dataset ID"))
        object.__setattr__(
            self,
            "perspective_id",
            _required_text(self.perspective_id, "Perspective ID"),
        )
        object.__setattr__(self, "name", _required_text(self.name, "Dataset name"))
        object.__setattr__(
            self,
            "description",
            _optional_text(self.description),
        )
        object.__setattr__(
            self,
            "derived_results_cohort_id",
            derived_id,
        )
        object.__setattr__(self, "visibility", Visibility(self.visibility))


@dataclass(frozen=True)
class CatalogSample:
    """Stable catalog identity for an observed or derived sample."""

    id: str
    cohort_id: str
    sample_type: SampleType
    display_label: str
    source_sample_id: str | None = None
    source_fingerprint: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)
    created_at: str = ""

    def __post_init__(self) -> None:
        sample_type = SampleType(self.sample_type)
        source_sample_id = _optional_text(self.source_sample_id) or None
        if sample_type is SampleType.OBSERVED and source_sample_id is None:
            raise ValueError(
                "Observed catalog sample requires a source sample ID."
            )
        if sample_type is SampleType.DERIVED and source_sample_id is not None:
            raise ValueError(
                "Derived catalog sample cannot have a source sample ID."
            )

        object.__setattr__(self, "id", _required_text(self.id, "Sample ID"))
        object.__setattr__(
            self,
            "cohort_id",
            _required_text(self.cohort_id, "Cohort ID"),
        )
        object.__setattr__(self, "sample_type", sample_type)
        object.__setattr__(
            self,
            "display_label",
            _required_text(self.display_label, "Sample display label"),
        )
        object.__setattr__(self, "source_sample_id", source_sample_id)
        object.__setattr__(
            self,
            "source_fingerprint",
            _optional_text(self.source_fingerprint),
        )
        object.__setattr__(self, "metadata", dict(self.metadata or {}))


def sample_membership_fingerprint(sample_ids: list[str] | tuple[str, ...]) -> str:
    """Return an order-independent fingerprint for exact group membership."""
    normalized = sorted({str(sample_id).strip() for sample_id in sample_ids})
    serialized = json.dumps(
        normalized,
        separators=(",", ":"),
        ensure_ascii=False,
    )
    return "sha256:" + hashlib.sha256(serialized.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class SampleGroup:
    """Perspective-owned editable selection within one dataset."""

    id: str
    perspective_id: str
    dataset_id: str
    name: str
    membership_fingerprint: str
    description: str = ""
    visibility: Visibility = Visibility.PRIVATE
    created_at: str = ""
    updated_at: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "id", _required_text(self.id, "Group ID"))
        object.__setattr__(
            self,
            "perspective_id",
            _required_text(self.perspective_id, "Perspective ID"),
        )
        object.__setattr__(
            self,
            "dataset_id",
            _required_text(self.dataset_id, "Dataset ID"),
        )
        object.__setattr__(self, "name", _required_text(self.name, "Group name"))
        object.__setattr__(
            self,
            "membership_fingerprint",
            _required_text(
                self.membership_fingerprint,
                "Membership fingerprint",
            ),
        )
        object.__setattr__(
            self,
            "description",
            _optional_text(self.description),
        )
        object.__setattr__(self, "visibility", Visibility(self.visibility))


@dataclass(frozen=True)
class SampleGroupSnapshot:
    """Immutable membership value ready to be attached to a calculation."""

    group_id: str
    perspective_id: str
    dataset_id: str
    group_name: str
    sample_ids: tuple[str, ...]
    membership_fingerprint: str
    captured_at: str

    def __post_init__(self) -> None:
        sample_ids = tuple(
            _required_text(sample_id, "Catalog sample ID")
            for sample_id in self.sample_ids
        )
        if len(set(sample_ids)) != len(sample_ids):
            raise ValueError("Group snapshot cannot contain duplicate samples.")
        fingerprint = _required_text(
            self.membership_fingerprint,
            "Membership fingerprint",
        )
        if fingerprint != sample_membership_fingerprint(sample_ids):
            raise ValueError(
                "Group snapshot fingerprint does not match its samples."
            )
        object.__setattr__(
            self,
            "group_id",
            _required_text(self.group_id, "Group ID"),
        )
        object.__setattr__(
            self,
            "perspective_id",
            _required_text(self.perspective_id, "Perspective ID"),
        )
        object.__setattr__(
            self,
            "dataset_id",
            _required_text(self.dataset_id, "Dataset ID"),
        )
        object.__setattr__(
            self,
            "group_name",
            _required_text(self.group_name, "Group name"),
        )
        object.__setattr__(self, "sample_ids", sample_ids)
        object.__setattr__(self, "membership_fingerprint", fingerprint)
        object.__setattr__(
            self,
            "captured_at",
            _required_text(self.captured_at, "Snapshot timestamp"),
        )
