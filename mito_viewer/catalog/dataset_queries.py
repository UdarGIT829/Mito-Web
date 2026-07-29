"""Request-scoped Dataset identity and source resolution for viewer queries."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Mapping

from mito_viewer.domain import SampleAlleleCall
from mito_viewer.domain.comparison import comparison_rows
from mito_viewer.repositories import (
    NO_TAGS_FILTER,
    StudyRepository,
    inspect_study_database,
)

from .derived_models import (
    DerivedRunRecord,
    DerivedRunStatus,
    DerivedSampleRecord,
    EvidenceStatus,
)
from .derived_repository import DerivedCatalogRepository
from .derived_service import ReferenceCompatibilityError
from .models import CatalogSample, Cohort, CohortType, Dataset, SampleType
from .repository import CatalogAccessError, CatalogRepository
from .workspaces import CatalogWorkspaceService


class DatasetSourceResolutionError(ValueError):
    """Raised when an attached source cohort cannot be queried safely."""


@dataclass(frozen=True)
class DatasetReferenceIdentity:
    """Reference properties required for safe cross-cohort allele merging."""

    reference_id: str
    normalization_version: str
    genome_build: str = ""
    contig: str = ""
    coordinate_system: str = ""
    sequence_fingerprint: str = ""

    @classmethod
    def from_cohort(cls, cohort: Cohort) -> "DatasetReferenceIdentity":
        reference = cohort.provenance.get("reference") or {}
        return cls(
            reference_id=str(cohort.reference_id or "").strip(),
            normalization_version=str(
                cohort.normalization_version or ""
            ).strip(),
            genome_build=str(reference.get("genome_build") or "").strip(),
            contig=str(reference.get("contig") or "").strip(),
            coordinate_system=str(
                reference.get("coordinate_system") or ""
            ).strip(),
            sequence_fingerprint=str(
                reference.get("sequence_fingerprint") or ""
            ).strip(),
        )

    @classmethod
    def from_run(cls, run: DerivedRunRecord) -> "DatasetReferenceIdentity":
        return cls(
            reference_id=str(run.reference_id or "").strip(),
            normalization_version=str(
                run.normalization_version or ""
            ).strip(),
        )

    @property
    def core(self) -> tuple[str, str]:
        return self.reference_id, self.normalization_version

    @property
    def label(self) -> str:
        return f"{self.reference_id}/{self.normalization_version}"


@dataclass(frozen=True)
class DatasetSourceRef:
    """One Dataset source cohort resolved to an allowed local database."""

    cohort: Cohort
    database_path: Path
    reference: DatasetReferenceIdentity

    def __post_init__(self) -> None:
        if self.cohort.cohort_type is not CohortType.SOURCE:
            raise ValueError("DatasetSourceRef requires a source cohort.")
        object.__setattr__(self, "database_path", self.database_path.resolve())

    @property
    def cohort_id(self) -> str:
        return self.cohort.id

    @property
    def database_id(self) -> str:
        return self.cohort.source_database_identifier


@dataclass(frozen=True)
class DatasetSampleRef:
    """Catalog-wide sample identity plus its Dataset query provenance."""

    catalog_sample: CatalogSample
    cohort: Cohort
    source: DatasetSourceRef | None = None
    derived_sample: DerivedSampleRecord | None = None
    current_run: DerivedRunRecord | None = None

    def __post_init__(self) -> None:
        if self.catalog_sample.cohort_id != self.cohort.id:
            raise ValueError("Dataset sample and cohort identities do not match.")
        if self.catalog_sample.sample_type is SampleType.OBSERVED:
            if self.source is None or self.derived_sample is not None:
                raise ValueError(
                    "Observed Dataset samples require source provenance."
                )
            if self.source.cohort_id != self.cohort.id:
                raise ValueError(
                    "Observed Dataset sample source cohort does not match."
                )
        elif self.source is not None or self.derived_sample is None:
            raise ValueError(
                "Derived Dataset samples require a derived catalog record."
            )
        if (
            self.derived_sample is not None
            and self.derived_sample.catalog_sample_id != self.catalog_sample.id
        ):
            raise ValueError(
                "Derived Dataset and catalog sample identities do not match."
            )

    @property
    def id(self) -> str:
        """Return the opaque catalog-wide ID exposed to the browser."""
        return self.catalog_sample.id

    @property
    def sample_type(self) -> SampleType:
        return self.catalog_sample.sample_type

    @property
    def source_sample_id(self) -> str | None:
        return self.catalog_sample.source_sample_id

    def to_browser_payload(self) -> dict:
        """Return the stable identity portion of the viewer sample payload."""
        metadata = self.catalog_sample.metadata
        is_derived = self.sample_type is SampleType.DERIVED
        payload = {
            "id": self.id,
            "catalog_sample_id": self.catalog_sample.id,
            "sample_type": self.sample_type.value,
            "cohort_id": self.cohort.id,
            "cohort_name": self.cohort.name,
            "display_label": self.catalog_sample.display_label,
            "subject_id": (
                "Derived"
                if is_derived
                else str(metadata.get("subject_id") or "")
            ),
            "population_key": (
                self.catalog_sample.display_label
                if is_derived
                else str(metadata.get("population_key") or "")
            ),
            "population_tags": (
                []
                if is_derived
                else list(metadata.get("population_tags") or [])
            ),
            "source_file": (
                ""
                if is_derived
                else str(metadata.get("source_file") or "")
            ),
            "source_database_id": (
                self.source.database_id if self.source is not None else ""
            ),
            "source_sample_id": self.source_sample_id or "",
            "is_derived": is_derived,
            "mutation_count": (
                self.current_run.output_count
                if self.current_run is not None
                else None
            ),
            "current_run_id": (
                self.current_run.id if self.current_run is not None else None
            ),
        }
        if self.derived_sample is not None:
            payload["catalog_derived_sample_id"] = self.derived_sample.id
        return payload


@dataclass(frozen=True)
class DatasetQueryScope:
    """Complete request-local identity boundary for one selected Dataset."""

    perspective_id: str
    dataset: Dataset
    cohorts: tuple[Cohort, ...]
    sources: tuple[DatasetSourceRef, ...]
    samples: tuple[DatasetSampleRef, ...]
    reference: DatasetReferenceIdentity | None
    _samples_by_id: Mapping[str, DatasetSampleRef] = field(
        init=False,
        repr=False,
        compare=False,
    )
    _sources_by_cohort_id: Mapping[str, DatasetSourceRef] = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if self.dataset.perspective_id != self.perspective_id:
            raise ValueError("Dataset query scope perspective does not match.")
        sample_map = {sample.id: sample for sample in self.samples}
        if len(sample_map) != len(self.samples):
            raise ValueError("Dataset query scope contains duplicate sample IDs.")
        source_map = {source.cohort_id: source for source in self.sources}
        if len(source_map) != len(self.sources):
            raise ValueError("Dataset query scope contains duplicate sources.")
        object.__setattr__(
            self,
            "_samples_by_id",
            MappingProxyType(sample_map),
        )
        object.__setattr__(
            self,
            "_sources_by_cohort_id",
            MappingProxyType(source_map),
        )

    @property
    def dataset_id(self) -> str:
        return self.dataset.id

    @property
    def observed_samples(self) -> tuple[DatasetSampleRef, ...]:
        return tuple(
            sample
            for sample in self.samples
            if sample.sample_type is SampleType.OBSERVED
        )

    @property
    def derived_samples(self) -> tuple[DatasetSampleRef, ...]:
        return tuple(
            sample
            for sample in self.samples
            if sample.sample_type is SampleType.DERIVED
        )

    def require_sample(self, sample_id: str) -> DatasetSampleRef:
        """Resolve one opaque ID and reject anything outside this Dataset."""
        sample_id = str(sample_id or "").strip()
        sample = self._samples_by_id.get(sample_id)
        if sample is None:
            raise CatalogAccessError(
                f"Sample {sample_id!r} is not available in "
                f"dataset {self.dataset.id!r}."
            )
        return sample

    def resolve_samples(
        self,
        sample_ids: list[str] | tuple[str, ...],
    ) -> tuple[DatasetSampleRef, ...]:
        """Resolve IDs in caller order without accepting duplicates."""
        normalized = tuple(str(sample_id or "").strip() for sample_id in sample_ids)
        if any(not sample_id for sample_id in normalized):
            raise ValueError("Dataset sample IDs cannot be empty.")
        if len(set(normalized)) != len(normalized):
            raise ValueError("Dataset sample IDs must be unique.")
        return tuple(self.require_sample(sample_id) for sample_id in normalized)

    def samples_for_source(
        self,
        cohort_id: str,
    ) -> tuple[DatasetSampleRef, ...]:
        if cohort_id not in self._sources_by_cohort_id:
            raise CatalogAccessError(
                f"Source cohort {cohort_id!r} is not available in "
                f"dataset {self.dataset.id!r}."
            )
        return tuple(
            sample
            for sample in self.observed_samples
            if sample.cohort.id == cohort_id
        )

    def browser_sample_payloads(self) -> list[dict]:
        return [sample.to_browser_payload() for sample in self.samples]


class DatasetQueryService:
    """Build Dataset query scopes from catalog state and allowed sources."""

    def __init__(
        self,
        catalog: CatalogRepository,
        *,
        source_database_paths: Mapping[str, str | Path],
    ) -> None:
        self.catalog = catalog
        self.workspaces = CatalogWorkspaceService(catalog)
        self.derived = DerivedCatalogRepository(catalog)
        self.source_database_paths = {
            str(identifier): Path(path).resolve()
            for identifier, path in source_database_paths.items()
        }

    def open_scope(
        self,
        perspective_id: str,
        dataset_id: str,
    ) -> DatasetQueryScope:
        """Resolve one Dataset without mutating process-global viewer state."""
        perspective_id = str(perspective_id or "").strip()
        dataset_id = str(dataset_id or "").strip()
        if not perspective_id or not dataset_id:
            raise ValueError(
                "Dataset queries require a Study Perspective and Dataset."
            )

        workspace = self.workspaces.open_dataset(
            dataset_id,
            acting_perspective_id=perspective_id,
        )
        source_refs = tuple(
            self._resolve_source(cohort)
            for cohort in workspace.cohorts
            if cohort.cohort_type is CohortType.SOURCE
        )
        source_by_cohort = {
            source.cohort_id: source
            for source in source_refs
        }
        derived_records = {
            record.catalog_sample_id: record
            for record in self.derived.list_derived_samples(
                dataset_id,
                acting_perspective_id=perspective_id,
            )
        }
        sample_refs = []
        derived_run_references = []
        for sample in self.workspaces.dataset_samples(
            dataset_id,
            acting_perspective_id=perspective_id,
        ):
            cohort = next(
                (
                    item
                    for item in workspace.cohorts
                    if item.id == sample.cohort_id
                ),
                None,
            )
            if cohort is None:
                raise RuntimeError(
                    f"Dataset sample {sample.id!r} has no attached cohort."
                )
            if sample.sample_type is SampleType.OBSERVED:
                source = source_by_cohort.get(cohort.id)
                if source is None:
                    raise RuntimeError(
                        f"Observed sample {sample.id!r} has no source."
                    )
                sample_refs.append(
                    DatasetSampleRef(
                        catalog_sample=sample,
                        cohort=cohort,
                        source=source,
                    )
                )
                continue

            derived_sample = derived_records.get(sample.id)
            if derived_sample is None:
                raise RuntimeError(
                    f"Derived catalog sample {sample.id!r} has no "
                    "derived record."
                )
            current_run = None
            if derived_sample.current_run_id:
                current_run = self.derived.get_run(
                    derived_sample.current_run_id
                )
                if current_run is None:
                    raise RuntimeError(
                        f"Derived sample {derived_sample.id!r} has a "
                        "missing current run."
                    )
                if current_run.status is not DerivedRunStatus.COMPLETED:
                    raise ValueError(
                        f"Derived sample {derived_sample.id!r} current "
                        "run is not completed."
                    )
                derived_run_references.append(
                    DatasetReferenceIdentity.from_run(current_run)
                )
            sample_refs.append(
                DatasetSampleRef(
                    catalog_sample=sample,
                    cohort=cohort,
                    derived_sample=derived_sample,
                    current_run=current_run,
                )
            )

        reference = _compatible_reference(
            tuple(source.reference for source in source_refs),
            tuple(derived_run_references),
        )
        return DatasetQueryScope(
            perspective_id=perspective_id,
            dataset=workspace.dataset,
            cohorts=workspace.cohorts,
            sources=source_refs,
            samples=tuple(sample_refs),
            reference=reference,
        )

    def sample_payloads(
        self,
        scope: DatasetQueryScope,
        *,
        subject_id: str = "",
        tags: list[str] | tuple[str, ...] = (),
    ) -> list[dict]:
        """Return Dataset samples after source subject/tag filtering."""
        subject_id = str(subject_id or "").strip()
        tags = tuple(str(tag).strip() for tag in tags if str(tag).strip())
        observed_rows: dict[str, dict] = {}

        for source in scope.sources:
            allowed_by_source_id = {
                sample.source_sample_id: sample
                for sample in scope.samples_for_source(source.cohort_id)
            }
            with StudyRepository.open(source.database_path) as study:
                rows = study.samples(
                    subject_id=subject_id or None,
                    tags=list(tags),
                )
            for row in rows:
                sample = allowed_by_source_id.get(str(row["id"]))
                if sample is not None:
                    observed_rows[sample.id] = row

        payloads = []
        for sample in scope.samples:
            if sample.sample_type is SampleType.OBSERVED:
                source_row = observed_rows.get(sample.id)
                if source_row is None:
                    continue
                payload = sample.to_browser_payload()
                payload.update(
                    {
                        "subject_id": source_row["subject_id"],
                        "population_key": source_row["population_key"],
                        "source_file": source_row["source_file"],
                        "mutation_count": source_row["mutation_count"],
                    }
                )
                payloads.append(payload)
            else:
                # Preserve existing viewer behavior: source filters narrow
                # observed samples without hiding durable results.
                payloads.append(sample.to_browser_payload())
        return payloads

    @staticmethod
    def subjects(scope: DatasetQueryScope) -> list[dict]:
        """Return subjects counted only across observed Dataset samples."""
        sample_ids_by_subject: dict[str, set[str]] = {}
        for sample in scope.observed_samples:
            subject_id = str(
                sample.catalog_sample.metadata.get("subject_id") or ""
            ).strip()
            if subject_id:
                sample_ids_by_subject.setdefault(subject_id, set()).add(
                    sample.id
                )
        return [
            {
                "id": subject_id,
                "subject_id": subject_id,
                "sample_count": len(sample_ids),
            }
            for subject_id, sample_ids in sorted(
                sample_ids_by_subject.items()
            )
        ]

    @staticmethod
    def population_tags(scope: DatasetQueryScope) -> list[dict]:
        """Return source tags counted only across observed Dataset samples."""
        sample_ids_by_tag: dict[str, set[str]] = {}
        untagged_sample_ids = set()
        for sample in scope.observed_samples:
            tags = [
                str(tag).strip()
                for tag in (
                    sample.catalog_sample.metadata.get("population_tags")
                    or []
                )
                if str(tag).strip()
            ]
            if not tags:
                untagged_sample_ids.add(sample.id)
            for tag in tags:
                sample_ids_by_tag.setdefault(tag, set()).add(sample.id)

        rows = [
            {
                "tag": tag,
                "sample_count": len(sample_ids),
            }
            for tag, sample_ids in sorted(sample_ids_by_tag.items())
        ]
        if untagged_sample_ids:
            rows.append(
                {
                    "tag": NO_TAGS_FILTER,
                    "label": "<NONE>",
                    "sample_count": len(untagged_sample_ids),
                }
            )
        return rows

    def counts(self, scope: DatasetQueryScope) -> dict[str, int]:
        """Return Dataset-wide subject, sample, and mutation counts."""
        payloads = self.sample_payloads(scope)
        subject_ids = {
            str(
                sample.catalog_sample.metadata.get("subject_id") or ""
            ).strip()
            for sample in scope.observed_samples
        }
        subject_ids.discard("")
        return {
            "subjects": len(subject_ids),
            "samples": len(payloads),
            "mutations": sum(
                int(payload.get("mutation_count") or 0)
                for payload in payloads
            ),
        }

    def mutation_rows(
        self,
        scope: DatasetQueryScope,
        *,
        sample_id: str = "",
        subject_id: str = "",
        tags: list[str] | tuple[str, ...] = (),
        position: int | str | None = None,
        alt: str = "",
        af_rules=(),
        metadata_filters=(),
        limit: int = 500,
        offset: int = 0,
        include_total: bool = False,
    ) -> list[dict] | dict:
        """Query, normalize, and merge mutation rows within a Dataset."""
        try:
            limit = int(limit)
        except (TypeError, ValueError) as exc:
            raise ValueError("Mutation limit must be an integer.") from exc
        if limit < 1:
            raise ValueError("Mutation limit must be positive.")
        try:
            offset = int(offset)
        except (TypeError, ValueError) as exc:
            raise ValueError("Mutation offset must be an integer.") from exc
        if offset < 0:
            raise ValueError("Mutation offset cannot be negative.")
        subject_id = str(subject_id or "").strip()
        tags = tuple(str(tag).strip() for tag in tags if str(tag).strip())

        selected = scope.require_sample(sample_id) if sample_id else None
        if selected is not None and selected.sample_type is SampleType.DERIVED:
            derived_rows = self._derived_mutation_rows(
                selected,
                position=position,
                alt=alt,
                af_rules=af_rules,
                metadata_filters=metadata_filters,
                limit=None if include_total else offset + limit,
            )
            if include_total:
                return _mutation_result_payload(
                    derived_rows[offset : offset + limit],
                    matched_count=len(derived_rows),
                    limit=limit,
                    offset=offset,
                )
            return derived_rows[offset : offset + limit]

        selected_by_cohort: dict[str, tuple[DatasetSampleRef, ...]] = {}
        if selected is not None:
            if _observed_sample_matches(selected, subject_id, tags):
                selected_by_cohort[selected.cohort.id] = (selected,)
        else:
            for source in scope.sources:
                selected_by_cohort[source.cohort_id] = tuple(
                    sample
                    for sample in scope.samples_for_source(
                        source.cohort_id
                    )
                    if _observed_sample_matches(
                        sample,
                        subject_id,
                        tags,
                    )
                )

        rows = []
        matched_count = 0
        for source in scope.sources:
            samples = selected_by_cohort.get(source.cohort_id, ())
            if not samples:
                continue
            by_source_id = {
                sample.source_sample_id: sample
                for sample in samples
            }
            with StudyRepository.open(source.database_path) as study:
                source_rows = study.mutation_rows(
                    sample_ids=tuple(by_source_id),
                    position=position,
                    alt=alt or None,
                    af_rules=af_rules,
                    metadata_filters=metadata_filters,
                    limit=offset + limit,
                )
                if include_total:
                    matched_count += study.mutation_count(
                        sample_ids=tuple(by_source_id),
                        position=position,
                        alt=alt or None,
                        af_rules=af_rules,
                        metadata_filters=metadata_filters,
                    )
            for source_row in source_rows:
                sample = by_source_id[str(source_row["sample_id"])]
                rows.append(
                    self._observed_mutation_payload(
                        source_row,
                        sample,
                    )
                )

        cohort_order = {
            cohort.id: index
            for index, cohort in enumerate(scope.cohorts)
        }
        sample_order = {
            sample.id: index
            for index, sample in enumerate(scope.samples)
        }
        rows.sort(
            key=lambda row: (
                int(row["pos"]),
                str(row["ref"]),
                str(row["alt"]),
                cohort_order[row["cohort_id"]],
                sample_order[row["catalog_sample_id"]],
                str(row["source_mutation_id"]),
            )
        )
        rows = rows[offset : offset + limit]
        if include_total:
            return _mutation_result_payload(
                rows,
                matched_count=matched_count,
                limit=limit,
                offset=offset,
            )
        return rows

    def compare_rows(
        self,
        scope: DatasetQueryScope,
        *,
        sample_ids: list[str] | tuple[str, ...],
        position: int | str | None = None,
        alt: str = "",
        af_rules=(),
        metadata_filters=(),
        statuses=(),
        sample_statuses: Mapping[str, set[str] | frozenset[str]] | None = None,
        limit: int = 2000,
    ) -> list[dict]:
        """Compare Dataset samples across observed sources and derived runs."""
        try:
            limit = int(limit)
        except (TypeError, ValueError) as exc:
            raise ValueError("Comparison limit must be an integer.") from exc
        if limit < 1:
            raise ValueError("Comparison limit must be positive.")

        selected = scope.resolve_samples(tuple(sample_ids))
        if len(selected) < 2:
            return []
        selected_ids = {sample.id for sample in selected}
        normalized_sample_statuses = {
            str(sample_id): set(allowed_statuses)
            for sample_id, allowed_statuses in (sample_statuses or {}).items()
        }
        unknown_constraints = (
            set(normalized_sample_statuses) - selected_ids
        )
        if unknown_constraints:
            unknown = sorted(unknown_constraints)[0]
            scope.require_sample(unknown)
            raise ValueError(
                f"Comparison constraint sample {unknown!r} is not selected."
            )

        calls = []
        for source in scope.sources:
            source_samples = tuple(
                sample
                for sample in selected
                if (
                    sample.sample_type is SampleType.OBSERVED
                    and sample.cohort.id == source.cohort_id
                )
            )
            if not source_samples:
                continue
            by_source_id = {
                sample.source_sample_id: sample
                for sample in source_samples
            }
            with StudyRepository.open(source.database_path) as study:
                source_calls = study.allele_calls(
                    list(by_source_id),
                    position=position,
                    alt=alt or None,
                    af_rules=af_rules,
                    metadata_filters=metadata_filters,
                )
            for call in source_calls:
                sample = by_source_id[str(call.sample_id)]
                call.sample_id = sample.id
                call.label = sample.catalog_sample.display_label
                calls.append(call)

        for sample in selected:
            if sample.sample_type is not SampleType.DERIVED:
                continue
            calls.extend(
                self._derived_allele_calls(
                    sample,
                    position=position,
                    alt=alt,
                    af_rules=af_rules,
                    metadata_filters=metadata_filters,
                )
            )

        labels = {
            sample.id: sample.catalog_sample.display_label
            for sample in selected
        }
        provenance = {
            sample.id: {
                "catalog_sample_id": sample.id,
                "sample_type": sample.sample_type.value,
                "cohort_id": sample.cohort.id,
                "cohort_name": sample.cohort.name,
                "source_database_id": (
                    sample.source.database_id
                    if sample.source is not None
                    else ""
                ),
                "source_sample_id": sample.source_sample_id or "",
            }
            for sample in selected
        }
        return comparison_rows(
            [sample.id for sample in selected],
            calls,
            sample_labels=labels,
            statuses=statuses,
            sample_statuses=normalized_sample_statuses,
            sample_provenance=provenance,
            limit=limit,
        )

    @staticmethod
    def _observed_mutation_payload(
        source_row: dict,
        sample: DatasetSampleRef,
    ) -> dict:
        source_mutation_id = str(source_row["id"])
        payload = {
            **source_row,
            "id": f"{sample.id}|source-mutation:{source_mutation_id}",
            "sample_id": sample.id,
            "catalog_sample_id": sample.id,
            "sample_type": SampleType.OBSERVED.value,
            "cohort_id": sample.cohort.id,
            "cohort_name": sample.cohort.name,
            "source_database_id": sample.source.database_id,
            "source_sample_id": sample.source_sample_id,
            "source_mutation_id": source_mutation_id,
        }
        return payload

    def _derived_mutation_rows(
        self,
        sample: DatasetSampleRef,
        *,
        position: int | str | None,
        alt: str,
        af_rules,
        metadata_filters,
        limit: int | None,
    ) -> list[dict]:
        if sample.current_run is None or sample.derived_sample is None:
            return []
        run = sample.current_run
        record = sample.derived_sample
        expression = _definition_expression(run)
        rows = []
        for materialized in self.derived.list_run_alleles(run.id):
            allele = materialized.allele
            if position and str(allele.position) != str(position):
                continue
            if alt and allele.alt != str(alt).strip().upper():
                continue
            evidence = self.derived.list_allele_evidence(materialized.id)
            qualifying = [
                item
                for item in evidence
                if item.evaluation_status
                is EvidenceStatus.QUALIFYING_PRESENT
            ]
            af_texts = _unique(
                item.af_text for item in qualifying if item.af_text
            )
            filter_texts = _unique(
                item.filter_text for item in qualifying if item.filter_text
            )
            if not _af_rules_match_text(",".join(af_texts), af_rules):
                continue
            original_labels = _unique(
                self._catalog_sample_label(
                    item.original_catalog_sample_id
                )
                for item in qualifying
                if item.original_catalog_sample_id
            )
            metadata = {
                "DERIVED_EXPRESSION": expression,
                "DERIVED_LABEL": record.name,
                "DERIVED_SAMPLE_ID": record.id,
                "DERIVED_RUN_ID": run.id,
                "DERIVED_SOURCE_SAMPLES": ";".join(original_labels),
                "DERIVED_PRESENT_SAMPLES": ";".join(original_labels),
                "DERIVED_SET_STATUSES": "durable",
                "CATALOG_DEFINITION_JSON": run.definition.normalized_json,
            }
            if not _metadata_filters_match(
                metadata,
                metadata_filters,
                allele.alt,
            ):
                continue
            rows.append(
                {
                    "id": (
                        f"{sample.id}|derived-allele:{materialized.id}"
                    ),
                    "sample_id": sample.id,
                    "catalog_sample_id": sample.id,
                    "sample_type": SampleType.DERIVED.value,
                    "cohort_id": sample.cohort.id,
                    "cohort_name": sample.cohort.name,
                    "source_database_id": "",
                    "source_sample_id": "",
                    "source_mutation_id": "",
                    "subject_id": "Derived",
                    "population_key": record.name,
                    "source_file": f"Catalog run {run.id}",
                    "pos": allele.position,
                    "ref": allele.ref,
                    "vcf_ref": allele.ref,
                    "alt": allele.alt,
                    "af": ",".join(af_texts),
                    "filter": ",".join(filter_texts),
                    "metadata_json": json.dumps(
                        metadata,
                        sort_keys=True,
                    ),
                }
            )
            if limit is not None and len(rows) >= limit:
                break
        return rows

    def _derived_allele_calls(
        self,
        sample: DatasetSampleRef,
        *,
        position: int | str | None,
        alt: str,
        af_rules,
        metadata_filters,
    ) -> list[SampleAlleleCall]:
        if sample.current_run is None or sample.derived_sample is None:
            raise ValueError(
                f"Derived sample {sample.id!r} has no usable current run."
            )
        run = sample.current_run
        record = sample.derived_sample
        expression = _definition_expression(run)
        calls = []
        for materialized in self.derived.list_run_alleles(run.id):
            allele = materialized.allele
            if position and str(allele.position) != str(position):
                continue
            if alt and allele.alt != str(alt).strip().upper():
                continue
            evidence = self.derived.list_allele_evidence(materialized.id)
            qualifying = [
                item
                for item in evidence
                if item.evaluation_status
                is EvidenceStatus.QUALIFYING_PRESENT
            ]
            af_texts = _unique(
                item.af_text for item in qualifying if item.af_text
            )
            filter_texts = _unique(
                item.filter_text for item in qualifying if item.filter_text
            )
            if not _af_rules_match_text(",".join(af_texts), af_rules):
                continue
            original_labels = _unique(
                self._catalog_sample_label(
                    item.original_catalog_sample_id
                )
                for item in qualifying
                if item.original_catalog_sample_id
            )
            metadata = {
                "DERIVED_EXPRESSION": expression,
                "DERIVED_LABEL": record.name,
                "DERIVED_SAMPLE_ID": record.id,
                "DERIVED_RUN_ID": run.id,
                "DERIVED_SOURCE_SAMPLES": ";".join(original_labels),
                "DERIVED_PRESENT_SAMPLES": ";".join(original_labels),
                "DERIVED_SET_STATUSES": "durable",
                "CATALOG_DEFINITION_JSON": run.definition.normalized_json,
            }
            if not _metadata_filters_match(
                metadata,
                metadata_filters,
                allele.alt,
            ):
                continue
            calls.append(
                SampleAlleleCall(
                    allele=allele,
                    sample_id=sample.id,
                    label=sample.catalog_sample.display_label,
                    af=(
                        qualifying[0].af
                        if len(qualifying) == 1
                        else None
                    ),
                    af_text=",".join(af_texts),
                    filter=",".join(filter_texts),
                    vcf_ref=allele.ref,
                    metadata=metadata,
                )
            )
        return calls

    def _catalog_sample_label(self, catalog_sample_id: str) -> str:
        sample = self.catalog.get_catalog_sample(catalog_sample_id)
        return sample.display_label if sample is not None else catalog_sample_id

    def _resolve_source(self, cohort: Cohort) -> DatasetSourceRef:
        path = self.source_database_paths.get(
            cohort.source_database_identifier
        )
        if path is None:
            registered_path_text = str(
                cohort.provenance.get("source_path") or ""
            ).strip()
            registered_path = (
                Path(registered_path_text).resolve()
                if registered_path_text
                else None
            )
            matches = {
                candidate
                for candidate in self.source_database_paths.values()
                if registered_path is not None
                and candidate == registered_path
            }
            if len(matches) == 1:
                path = matches.pop()
            elif len(matches) > 1:
                raise DatasetSourceResolutionError(
                    "Registered source database resolves ambiguously for "
                    f"{cohort.source_database_identifier!r}."
                )

        if path is None or not path.is_file():
            raise DatasetSourceResolutionError(
                "Registered source database is not available for "
                f"cohort {cohort.name!r} "
                f"({cohort.source_database_identifier!r})."
            )

        report = inspect_study_database(path)
        if not report.valid:
            raise DatasetSourceResolutionError(
                "Registered source database is incompatible for "
                f"cohort {cohort.name!r}: {path}."
            )
        return DatasetSourceRef(
            cohort=cohort,
            database_path=path,
            reference=DatasetReferenceIdentity.from_cohort(cohort),
        )


def _compatible_reference(
    source_references: tuple[DatasetReferenceIdentity, ...],
    derived_references: tuple[DatasetReferenceIdentity, ...],
) -> DatasetReferenceIdentity | None:
    """Return one compatible identity or reject unsafe allele merging."""
    for reference in (*source_references, *derived_references):
        if not all(reference.core):
            raise ReferenceCompatibilityError(
                "Every Dataset source must declare a reference and "
                "normalization identity."
            )

    source_identity_set = set(source_references)
    if len(source_identity_set) > 1:
        labels = ", ".join(
            sorted(reference.label for reference in source_identity_set)
        )
        raise ReferenceCompatibilityError(
            f"Dataset source cohorts use incompatible references: {labels}"
        )

    core_set = {
        reference.core
        for reference in (*source_references, *derived_references)
    }
    if len(core_set) > 1:
        labels = ", ".join(
            f"{reference}/{normalization}"
            for reference, normalization in sorted(core_set)
        )
        raise ReferenceCompatibilityError(
            f"Dataset samples use incompatible references: {labels}"
        )

    if source_references:
        return source_references[0]
    if derived_references:
        return derived_references[0]
    return None


def _definition_expression(run: DerivedRunRecord) -> str:
    clauses = [
        f"{item.requirement.value.upper()}({item.role}:{item.input_id})"
        for item in run.definition.inputs
    ]
    expression = " AND ".join(clauses)
    filters = run.definition.filters
    filter_values = []
    if filters.position:
        filter_values.append(f"position={filters.position}")
    if filters.alt:
        filter_values.append(f"ALT={filters.alt}")
    filter_values.extend(
        f"AF {rule.operator} {rule.threshold}"
        for rule in filters.af_rules
    )
    filter_values.extend(
        f"{item.field}={item.value}"
        for item in filters.metadata_filters
    )
    if filter_values:
        expression += " FILTER " + " AND ".join(filter_values)
    return expression


def _unique(values) -> list[str]:
    return list(dict.fromkeys(str(value) for value in values if value))


def _mutation_result_payload(
    rows: list[dict],
    *,
    matched_count: int,
    limit: int,
    offset: int = 0,
) -> dict:
    next_offset = int(offset) + len(rows)
    return {
        "rows": rows,
        "matched_count": int(matched_count),
        "shown_count": len(rows),
        "limit": int(limit),
        "offset": int(offset),
        "next_offset": next_offset,
        "has_more": next_offset < int(matched_count),
        "truncated": next_offset < int(matched_count),
    }


def _af_rules_match_text(af_text: str, af_rules) -> bool:
    rules = tuple(af_rules or ())
    if not rules:
        return True
    values = []
    for raw_value in str(af_text or "").split(","):
        try:
            values.append(float(raw_value))
        except ValueError:
            continue
    for rule in rules:
        operator, threshold = rule
        threshold = float(threshold)
        if not any(
            _af_value_matches(value, operator, threshold)
            for value in values
        ):
            return False
    return True


def _af_value_matches(value: float, operator: str, threshold: float) -> bool:
    return {
        "gt": value > threshold,
        "gte": value >= threshold,
        "lt": value < threshold,
        "lte": value <= threshold,
        "eq": value == threshold,
        "neq": value != threshold,
    }[str(operator)]


def _metadata_filters_match(
    metadata: dict,
    metadata_filters,
    alt: str,
) -> bool:
    """Preserve the legacy derived-row behavior for active table filters."""
    for item in metadata_filters or ():
        field, raw_value = item
        if field == "polymorphism":
            if str(metadata.get("POLYMORPHISM", "")) != raw_value:
                return False
        elif field == "reference_contains_alt":
            alt_text = str(alt).upper()
            before = str(
                metadata.get("REFERENCE_6_BEFORE", "")
            ).upper()
            after = str(
                metadata.get("REFERENCE_6_AFTER", "")
            ).upper()
            contains = bool(alt_text) and (
                alt_text in before or alt_text in after
            )
            if raw_value == "contains" and not contains:
                return False
            if raw_value == "not_contains" and contains:
                return False
        elif field == "reference_context":
            context = (
                str(metadata.get("REFERENCE_6_BEFORE", ""))
                + str(metadata.get("REFERENCE_6_AFTER", ""))
            )
            if raw_value not in context:
                return False
        elif field == "reference_repeat":
            before_seen = _single_base_repeat_seen(
                metadata.get("REFERENCE_6_BEFORE", "")
            )
            after_seen = _single_base_repeat_seen(
                metadata.get("REFERENCE_6_AFTER", "")
            )
            if raw_value == "before" and not (
                before_seen and not after_seen
            ):
                return False
            if raw_value == "after" and not (
                after_seen and not before_seen
            ):
                return False
            if raw_value == "one" and before_seen == after_seen:
                return False
            if raw_value == "both" and not (
                before_seen and after_seen
            ):
                return False
            if raw_value == "none" and (
                before_seen or after_seen
            ):
                return False
            if raw_value == "either" and not (
                before_seen or after_seen
            ):
                return False
        else:
            key = {
                "repeat_base": "REPEAT_1_BASE",
                "repeat_count": "REPEAT_1_BASE_COUNT",
                "repeat_2_bases": "REPEAT_2_BASES",
                "repeat_2_count": "REPEAT_2_BASES_COUNT",
                "repeat_3_bases": "REPEAT_3_BASES",
                "repeat_3_count": "REPEAT_3_BASES_COUNT",
            }.get(field)
            if key is None:
                continue
            if field.endswith("count"):
                operator, separator, threshold = raw_value.partition("|")
                if not separator:
                    continue
                if key not in metadata or metadata.get(key, "") == "":
                    return False
                try:
                    value = int(metadata.get(key))
                    threshold_value = int(threshold)
                except (TypeError, ValueError):
                    return False
                if not _af_value_matches(
                    value,
                    operator,
                    threshold_value,
                ):
                    return False
            elif str(metadata.get(key, "")) != raw_value:
                return False
    return True


def _single_base_repeat_seen(sequence) -> bool:
    sequence = str(sequence or "").upper()
    return any(base * 2 in sequence for base in ("A", "C", "G", "T", "N"))


def _observed_sample_matches(
    sample: DatasetSampleRef,
    subject_id: str,
    tags: tuple[str, ...],
) -> bool:
    metadata = sample.catalog_sample.metadata
    if (
        subject_id
        and str(metadata.get("subject_id") or "") != subject_id
    ):
        return False
    sample_tags = {
        str(tag)
        for tag in metadata.get("population_tags") or ()
    }
    for tag in tags:
        if tag == NO_TAGS_FILTER:
            if sample_tags:
                return False
        elif tag not in sample_tags:
            return False
    return True
