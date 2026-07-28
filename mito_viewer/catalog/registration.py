"""Read-only study database registration for the writable catalog."""

from __future__ import annotations

import hashlib
import json
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from mito_viewer.repositories import StudyRepository, inspect_study_database

from .models import CatalogSample, Cohort, CohortType, SampleType, Visibility
from .repository import CatalogRegistrationConflictError, CatalogRepository


REGISTRATION_VERSION = 1
DEFAULT_REFERENCE_ID = "NC_012920.1"
DEFAULT_GENOME_BUILD = "GRCh38/hg38"
DEFAULT_CONTIG = "chrM"
DEFAULT_COORDINATE_SYSTEM = "1-based"
DEFAULT_NORMALIZATION_VERSION = "vcf-as-imported-v1"


@dataclass(frozen=True)
class SourceRegistrationSpec:
    """Configuration required to register one source study database."""

    path: Path
    identifier: str
    display_name: str
    description: str = ""
    reference_id: str = DEFAULT_REFERENCE_ID
    genome_build: str = DEFAULT_GENOME_BUILD
    contig: str = DEFAULT_CONTIG
    coordinate_system: str = DEFAULT_COORDINATE_SYSTEM
    normalization_version: str = DEFAULT_NORMALIZATION_VERSION
    reference_sequence_fingerprint: str = ""
    visibility: Visibility = Visibility.PUBLIC
    provenance: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        path = Path(self.path).resolve()
        identifier = str(self.identifier or "").strip()
        display_name = str(self.display_name or "").strip()
        if not identifier:
            raise ValueError("Source database identifier cannot be empty.")
        if not display_name:
            raise ValueError("Source cohort display name cannot be empty.")
        object.__setattr__(self, "path", path)
        object.__setattr__(self, "identifier", identifier)
        object.__setattr__(self, "display_name", display_name)
        object.__setattr__(self, "visibility", Visibility(self.visibility))
        object.__setattr__(self, "provenance", dict(self.provenance or {}))


@dataclass(frozen=True)
class SourceRegistrationResult:
    """Catalog entities and status produced by source registration."""

    cohort: Cohort
    samples: tuple[CatalogSample, ...]
    created: bool


def sha256_file(path: str | Path) -> str:
    """Return a prefixed SHA-256 fingerprint without modifying the file."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        while chunk := source.read(1024 * 1024):
            digest.update(chunk)
    return f"sha256:{digest.hexdigest()}"


def fasta_sequence_fingerprint(path: str | Path) -> str:
    """Hash normalized FASTA sequence content, excluding its header/wrapping."""
    sequence_parts = []
    with Path(path).open(encoding="utf-8") as source:
        for line in source:
            line = line.strip()
            if line and not line.startswith(">"):
                sequence_parts.append(line.upper())
    sequence = "".join(sequence_parts)
    if not sequence:
        raise ValueError(f"Reference FASTA has no sequence: {path}")
    return "sha256:" + hashlib.sha256(sequence.encode("ascii")).hexdigest()


def stable_source_cohort_id(identifier: str) -> str:
    """Derive a portable catalog ID from the logical source identifier."""
    return str(
        uuid.uuid5(
            uuid.NAMESPACE_URL,
            f"urn:mito-viewer:source-cohort:{identifier}",
        )
    )


def stable_catalog_sample_id(cohort_id: str, source_sample_id: str) -> str:
    """Derive a catalog-wide ID without coupling it to source content."""
    return str(
        uuid.uuid5(
            uuid.NAMESPACE_URL,
            "urn:mito-viewer:observed-sample:"
            f"{cohort_id}:{source_sample_id}",
        )
    )


def source_sample_fingerprint(
    repository: StudyRepository,
    source_sample_id: str,
) -> str:
    """Hash the complete stored snapshot for one source sample."""
    digest = hashlib.sha256()
    for record in repository.iter_sample_snapshot_records(source_sample_id):
        serialized = json.dumps(
            record,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        )
        digest.update(serialized.encode("utf-8"))
        digest.update(b"\n")
    return f"sha256:{digest.hexdigest()}"


def register_study_database(
    catalog: CatalogRepository,
    spec: SourceRegistrationSpec,
) -> SourceRegistrationResult:
    """Register one schema-valid source database and its observed samples."""
    if not spec.path.is_file():
        raise FileNotFoundError(f"Study database not found: {spec.path}")

    database_fingerprint = sha256_file(spec.path)
    schema_report = inspect_study_database(spec.path)
    schema_report.require_valid()
    cohort_id = stable_source_cohort_id(spec.identifier)
    registered_at = _utc_now()

    with StudyRepository.open(spec.path) as study:
        sample_rows = study.registration_samples()
        samples = tuple(
            _catalog_sample(
                study,
                row,
                cohort_id=cohort_id,
                database_fingerprint=database_fingerprint,
                registered_at=registered_at,
            )
            for row in sample_rows
        )

    if sha256_file(spec.path) != database_fingerprint:
        raise CatalogRegistrationConflictError(
            f"Source database changed during registration: {spec.path}"
        )

    source_stat = spec.path.stat()
    cohort = Cohort(
        id=cohort_id,
        name=spec.display_name,
        description=spec.description,
        cohort_type=CohortType.SOURCE,
        source_database_identifier=spec.identifier,
        source_database_fingerprint=database_fingerprint,
        reference_id=spec.reference_id,
        normalization_version=spec.normalization_version,
        provenance={
            "registration_version": REGISTRATION_VERSION,
            "source_kind": "legacy-study-sqlite",
            "source_path": str(spec.path),
            "source_size_bytes": source_stat.st_size,
            "study_schema_version": schema_report.user_version,
            "study_schema_label": schema_report.version_label,
            "reference": {
                "accession": spec.reference_id,
                "genome_build": spec.genome_build,
                "contig": spec.contig,
                "coordinate_system": spec.coordinate_system,
                "sequence_fingerprint": spec.reference_sequence_fingerprint,
            },
            "mutation_storage": "read-only-source-database",
            **spec.provenance,
        },
        visibility=spec.visibility,
        created_at=registered_at,
        updated_at=registered_at,
    )
    created = catalog.register_source_snapshot(cohort, samples)
    stored_cohort = catalog.get_cohort(cohort.id)
    if stored_cohort is None:
        raise RuntimeError("Source registration did not persist its cohort.")
    return SourceRegistrationResult(
        cohort=stored_cohort,
        samples=tuple(catalog.list_catalog_samples(cohort.id)),
        created=created,
    )


def current_source_specs(
    project_root: str | Path,
) -> tuple[SourceRegistrationSpec, ...]:
    """Return the two approved source registrations for this project."""
    project_root = Path(project_root).resolve()
    reference_fingerprint = fasta_sequence_fingerprint(
        project_root / "reference" / "hg38_chrM.fa"
    )
    common = {
        "reference_sequence_fingerprint": reference_fingerprint,
        "provenance": {
            "registration_scope": "initial-local-study-cohorts",
        },
    }
    return (
        SourceRegistrationSpec(
            path=project_root / "EV_Study.sqlite",
            identifier="EV_Study.sqlite",
            display_name="EV Study",
            description="Existing EV mitochondrial study cohort.",
            **common,
        ),
        SourceRegistrationSpec(
            path=project_root / "lupusN.sqlite",
            identifier="lupusN.sqlite",
            display_name="Lupus N",
            description="Existing lupus mitochondrial study cohort.",
            **common,
        ),
    )


def register_current_source_cohorts(
    catalog: CatalogRepository,
    project_root: str | Path,
) -> tuple[SourceRegistrationResult, ...]:
    """Register EV and lupusN using the approved reference metadata."""
    with catalog.transaction():
        return tuple(
            register_study_database(catalog, spec)
            for spec in current_source_specs(project_root)
        )


def _catalog_sample(
    study: StudyRepository,
    row: dict,
    *,
    cohort_id: str,
    database_fingerprint: str,
    registered_at: str,
) -> CatalogSample:
    source_sample_id = str(row["id"])
    population_label = str(row["population_key"]).replace("|", "_").strip()
    display_label = " ".join(
        part
        for part in (str(row["subject_id"]).strip(), population_label)
        if part
    )
    return CatalogSample(
        id=stable_catalog_sample_id(cohort_id, source_sample_id),
        cohort_id=cohort_id,
        sample_type=SampleType.OBSERVED,
        display_label=display_label,
        source_sample_id=source_sample_id,
        source_fingerprint=source_sample_fingerprint(
            study,
            source_sample_id,
        ),
        metadata={
            "subject_id": row["subject_id"],
            "population_key": row["population_key"],
            "population_tags": list(row["population_tags"]),
            "source_file": row["source_file"],
            "source_database_fingerprint": database_fingerprint,
        },
        created_at=registered_at,
    )


def _utc_now() -> str:
    return (
        datetime.now(timezone.utc)
        .isoformat(timespec="microseconds")
        .replace("+00:00", "Z")
    )
