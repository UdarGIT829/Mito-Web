"""Minimal adapter between the legacy viewer and the durable catalog."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import vcf_parser
from mito_viewer.domain import (
    DerivedSample,
    MutationFilters,
    SampleAlleleCall,
)

from .derived_models import (
    DerivedDefinition,
    DerivedInput,
    DerivedInputKind,
    EvidenceStatus,
    PresenceRequirement,
)
from .derived_repository import DerivedCatalogRepository
from .derived_service import DerivedAnalysisService
from .models import CatalogSample, Cohort, CohortType, SampleType
from .repository import CatalogNotFoundError, CatalogRepository
from .workspaces import CatalogWorkspaceService


DURABLE_DERIVED_PREFIX = "derived:"


class CatalogViewerService:
    """Expose catalog workspaces through legacy-viewer-shaped values."""

    def __init__(self, catalog: CatalogRepository) -> None:
        self.catalog = catalog
        self.workspaces = CatalogWorkspaceService(catalog)
        self.derived = DerivedCatalogRepository(catalog)
        self.analysis = DerivedAnalysisService(catalog)

    def list_perspectives(self) -> list[dict]:
        return [
            {
                "id": item.id,
                "name": item.name,
                "description": item.description,
                "visibility": item.visibility.value,
            }
            for item in self.catalog.list_perspectives()
        ]

    def create_perspective(self, name: str) -> dict:
        workspace = self.workspaces.create_perspective(name)
        return {
            "id": workspace.perspective.id,
            "name": workspace.perspective.name,
            "description": workspace.perspective.description,
            "visibility": workspace.perspective.visibility.value,
        }

    def list_datasets(self, perspective_id: str) -> list[dict]:
        workspace = self.workspaces.select_perspective(perspective_id)
        return [
            {
                "id": item.id,
                "perspective_id": item.perspective_id,
                "name": item.name,
                "description": item.description,
                "derived_results_cohort_id": item.derived_results_cohort_id,
            }
            for item in workspace.datasets
        ]

    def create_dataset(
        self,
        perspective_id: str,
        name: str,
        *,
        database_id: str,
        database_path: str | Path,
    ) -> dict:
        cohort = self.source_cohort_for_database(
            database_id,
            database_path,
        )
        if cohort is None:
            raise ValueError(f"{database_id!r} is not registered as a source cohort.")
        workspace = self.workspaces.create_dataset(
            perspective_id,
            name,
            cohort_ids=(cohort.id,),
        )
        return {
            "id": workspace.dataset.id,
            "perspective_id": workspace.dataset.perspective_id,
            "name": workspace.dataset.name,
            "description": workspace.dataset.description,
            "derived_results_cohort_id": (workspace.dataset.derived_results_cohort_id),
        }

    def workspace(
        self,
        perspective_id: str,
        dataset_id: str,
        *,
        database_id: str,
        database_path: str | Path,
    ) -> dict:
        workspace = self.workspaces.open_dataset(
            dataset_id,
            acting_perspective_id=perspective_id,
        )
        source_cohort = self.source_cohort_for_database(
            database_id,
            database_path,
            cohorts=workspace.cohorts,
        )
        groups = self.catalog.list_sample_groups(
            dataset_id,
            acting_perspective_id=perspective_id,
        )
        return {
            "perspective_id": perspective_id,
            "dataset": {
                "id": workspace.dataset.id,
                "name": workspace.dataset.name,
                "derived_results_cohort_id": (
                    workspace.dataset.derived_results_cohort_id
                ),
            },
            "cohorts": [
                {
                    "id": cohort.id,
                    "name": cohort.name,
                    "cohort_type": cohort.cohort_type.value,
                    "database_id": (
                        cohort.source_database_identifier
                        if cohort.cohort_type is CohortType.SOURCE
                        else ""
                    ),
                    "current_database": (
                        source_cohort is not None and cohort.id == source_cohort.id
                    ),
                }
                for cohort in workspace.cohorts
            ],
            "groups": [
                self._group_payload(
                    group.id,
                    perspective_id=perspective_id,
                    source_cohort=source_cohort,
                )
                for group in groups
            ],
            "derived_samples": [
                self._derived_payload(item)
                for item in self.derived.list_derived_samples(
                    dataset_id,
                    acting_perspective_id=perspective_id,
                )
            ],
        }

    def create_group(
        self,
        perspective_id: str,
        dataset_id: str,
        name: str,
        viewer_sample_ids: Iterable[str],
        *,
        database_id: str,
        database_path: str | Path,
    ) -> dict:
        source_cohort = self.require_dataset_source_cohort(
            perspective_id,
            dataset_id,
            database_id=database_id,
            database_path=database_path,
        )
        sample_ids = [
            self._catalog_sample_for_viewer_id(
                source_cohort,
                viewer_sample_id,
                allow_derived=False,
            ).id
            for viewer_sample_id in viewer_sample_ids
        ]
        group = self.workspaces.create_sample_group(
            perspective_id,
            dataset_id,
            name,
            sample_ids,
        )
        return self._group_payload(
            group.id,
            perspective_id=perspective_id,
            source_cohort=source_cohort,
        )

    def attach_database(
        self,
        perspective_id: str,
        dataset_id: str,
        *,
        database_id: str,
        database_path: str | Path,
    ) -> dict:
        workspace = self.workspaces.open_dataset(
            dataset_id,
            acting_perspective_id=perspective_id,
        )
        cohort = self.source_cohort_for_database(
            database_id,
            database_path,
        )
        if cohort is None:
            raise ValueError(
                f"{database_id!r} is not registered as a source cohort."
            )
        self.catalog.add_cohort_to_dataset(
            dataset_id,
            cohort.id,
            acting_perspective_id=perspective_id,
            display_order=len(workspace.cohorts),
        )
        return {
            "id": cohort.id,
            "name": cohort.name,
            "database_id": cohort.source_database_identifier,
            "cohort_type": cohort.cohort_type.value,
        }

    def durable_derived_samples(
        self,
        perspective_id: str,
        dataset_id: str,
    ) -> dict[str, DerivedSample]:
        records = self.derived.list_derived_samples(
            dataset_id,
            acting_perspective_id=perspective_id,
        )
        return {
            durable_viewer_sample_id(record.id): self._legacy_derived(record)
            for record in records
            if record.current_run_id
        }

    def save_comparison(
        self,
        perspective_id: str,
        dataset_id: str,
        name: str,
        viewer_sample_ids: Iterable[str],
        *,
        database_id: str,
        database_path: str | Path,
        sample_statuses: dict[str, set[str]],
        global_statuses: Iterable[str],
        filters: MutationFilters,
    ) -> dict:
        source_cohort = self.require_dataset_source_cohort(
            perspective_id,
            dataset_id,
            database_id=database_id,
            database_path=database_path,
        )
        inputs = []
        positive_count = 0
        viewer_sample_ids = tuple(str(item) for item in viewer_sample_ids)
        for viewer_sample_id in viewer_sample_ids:
            statuses = set(sample_statuses.get(viewer_sample_id, {"present"}))
            if statuses == {"present"}:
                requirement = PresenceRequirement.ANY
                positive_count += 1
            elif statuses == {"not_in"}:
                requirement = PresenceRequirement.NONE
            else:
                raise ValueError(
                    "Durable saves currently require each selected sample "
                    "to be exactly Present or Not In. Mixed Unique/None "
                    "rules can still be viewed but cannot yet be persisted."
                )
            kind, input_id = self._derived_input_for_viewer_sample(
                source_cohort,
                viewer_sample_id,
                perspective_id=perspective_id,
            )
            inputs.append(
                DerivedInput(
                    kind=kind,
                    input_id=input_id,
                    role=(
                        "required-present"
                        if requirement is PresenceRequirement.ANY
                        else "required-absent"
                    ),
                    requirement=requirement,
                )
            )
        if positive_count < 1:
            raise ValueError(
                "A durable derived sample requires at least one Present input."
            )
        expected_status = _fixed_output_status(
            positive_count,
            len(viewer_sample_ids),
        )
        selected_global_statuses = set(global_statuses)
        if selected_global_statuses and expected_status not in selected_global_statuses:
            raise ValueError(
                "The selected global comparison statuses exclude the "
                f"{expected_status!r} result implied by the sample rules."
            )

        result = self.analysis.create_derived_sample(
            perspective_id=perspective_id,
            dataset_id=dataset_id,
            name=name,
            definition=DerivedDefinition(
                inputs=tuple(inputs),
                filters=filters,
            ),
        )
        payload = self._derived_payload(result.derived_sample)
        payload.update(
            {
                "id": durable_viewer_sample_id(result.derived_sample.id),
                "catalog_derived_sample_id": result.derived_sample.id,
                "run_id": result.run.id,
                "mutation_count": result.run.output_count,
                "is_derived": True,
                "durable": True,
                "subject_id": "Derived",
                "population_key": result.derived_sample.name,
                "source_file": f"Catalog run {result.run.id}",
            }
        )
        return payload

    def source_cohort_for_database(
        self,
        database_id: str,
        database_path: str | Path,
        *,
        cohorts: Iterable[Cohort] | None = None,
    ) -> Cohort | None:
        candidates = (
            tuple(cohorts)
            if cohorts is not None
            else tuple(self.catalog.list_cohorts())
        )
        resolved_path = Path(database_path).resolve()
        for cohort in candidates:
            if cohort.cohort_type is not CohortType.SOURCE:
                continue
            provenance_path = cohort.provenance.get("source_path")
            if cohort.source_database_identifier == database_id:
                return cohort
            if provenance_path and Path(provenance_path).resolve() == resolved_path:
                return cohort
        return None

    def require_dataset_source_cohort(
        self,
        perspective_id: str,
        dataset_id: str,
        *,
        database_id: str,
        database_path: str | Path,
    ) -> Cohort:
        workspace = self.workspaces.open_dataset(
            dataset_id,
            acting_perspective_id=perspective_id,
        )
        cohort = self.source_cohort_for_database(
            database_id,
            database_path,
            cohorts=workspace.cohorts,
        )
        if cohort is None:
            raise ValueError(
                f"Database {database_id!r} is not part of dataset "
                f"{workspace.dataset.name!r}."
            )
        return cohort

    def _derived_input_for_viewer_sample(
        self,
        source_cohort: Cohort,
        viewer_sample_id: str,
        *,
        perspective_id: str,
    ) -> tuple[DerivedInputKind, str]:
        if is_durable_viewer_sample_id(viewer_sample_id):
            derived_id = parse_durable_viewer_sample_id(viewer_sample_id)
            record = self.derived.get_derived_sample(
                derived_id,
                acting_perspective_id=perspective_id,
            )
            if record is None or record.current_run_id is None:
                raise CatalogNotFoundError(
                    f"Durable derived sample not found: {viewer_sample_id}"
                )
            return DerivedInputKind.RUN, record.current_run_id
        sample = self._catalog_sample_for_viewer_id(
            source_cohort,
            viewer_sample_id,
            allow_derived=False,
        )
        return DerivedInputKind.SAMPLE, sample.id

    def _catalog_sample_for_viewer_id(
        self,
        source_cohort: Cohort,
        viewer_sample_id: str,
        *,
        allow_derived: bool,
    ) -> CatalogSample:
        if is_durable_viewer_sample_id(viewer_sample_id):
            if not allow_derived:
                raise ValueError(
                    "Sample groups can contain observed samples only in "
                    "the current viewer integration."
                )
        else:
            for sample in self.catalog.list_catalog_samples(source_cohort.id):
                if sample.source_sample_id == str(viewer_sample_id):
                    return sample
        raise CatalogNotFoundError(
            f"Viewer sample not found in the current source cohort: {viewer_sample_id}"
        )

    def _legacy_derived(self, record) -> DerivedSample:
        run = self.derived.get_run(record.current_run_id)
        if run is None:
            raise CatalogNotFoundError(
                f"Current derived run not found: {record.current_run_id}"
            )
        viewer_id = durable_viewer_sample_id(record.id)
        calls = []
        mutations = []
        expression = _definition_expression(run.definition)
        for materialized in self.derived.list_run_alleles(run.id):
            evidence = self.derived.list_allele_evidence(materialized.id)
            qualifying = [
                item
                for item in evidence
                if item.evaluation_status is EvidenceStatus.QUALIFYING_PRESENT
            ]
            af_texts = _unique(item.af_text for item in qualifying if item.af_text)
            filter_texts = _unique(
                item.filter_text for item in qualifying if item.filter_text
            )
            original_labels = _unique(
                self._sample_label(item.original_catalog_sample_id)
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
            calls.append(
                SampleAlleleCall(
                    allele=materialized.allele,
                    sample_id=viewer_id,
                    label=f"Derived {record.name}",
                    af=(qualifying[0].af if len(qualifying) == 1 else None),
                    af_text=",".join(af_texts),
                    filter=",".join(filter_texts),
                    vcf_ref=materialized.allele.ref,
                    metadata=metadata,
                )
            )
            mutations.append(
                vcf_parser.VCFMutation(
                    position=materialized.allele.position,
                    alt=materialized.allele.alt,
                    metadata=metadata,
                    ref=materialized.allele.ref,
                    filter=",".join(filter_texts),
                )
            )
        return DerivedSample(
            id=viewer_id,
            label=record.name,
            calls=calls,
            mutations=mutations,
            source_description=f"Catalog run {run.id}: {expression}",
        )

    def _group_payload(
        self,
        group_id: str,
        *,
        perspective_id: str,
        source_cohort: Cohort | None,
    ) -> dict:
        group = self.catalog.get_sample_group(
            group_id,
            acting_perspective_id=perspective_id,
        )
        members = self.catalog.sample_group_members(
            group_id,
            acting_perspective_id=perspective_id,
        )
        visible_ids = [
            sample.source_sample_id
            for sample in members
            if (
                source_cohort is not None
                and sample.cohort_id == source_cohort.id
                and sample.sample_type is SampleType.OBSERVED
            )
        ]
        return {
            "id": group.id,
            "name": group.name,
            "description": group.description,
            "sample_ids": visible_ids,
            "visible_sample_count": len(visible_ids),
            "total_sample_count": len(members),
            "membership_fingerprint": group.membership_fingerprint,
        }

    @staticmethod
    def _derived_payload(record) -> dict:
        return {
            "id": durable_viewer_sample_id(record.id),
            "catalog_derived_sample_id": record.id,
            "name": record.name,
            "description": record.description,
            "current_run_id": record.current_run_id,
            "visibility": record.visibility.value,
            "durable": True,
        }

    def _sample_label(self, catalog_sample_id: str) -> str:
        sample = self.catalog.get_catalog_sample(catalog_sample_id)
        return sample.display_label if sample else catalog_sample_id


def durable_viewer_sample_id(derived_sample_id: str) -> str:
    return DURABLE_DERIVED_PREFIX + str(derived_sample_id)


def is_durable_viewer_sample_id(viewer_sample_id: str) -> bool:
    return str(viewer_sample_id).startswith(DURABLE_DERIVED_PREFIX)


def parse_durable_viewer_sample_id(viewer_sample_id: str) -> str:
    value = str(viewer_sample_id)
    if not is_durable_viewer_sample_id(value):
        raise ValueError(f"Not a durable derived sample ID: {value}")
    derived_id = value[len(DURABLE_DERIVED_PREFIX) :]
    if not derived_id:
        raise ValueError("Durable derived sample ID cannot be empty.")
    return derived_id


def _fixed_output_status(positive_count: int, total_count: int) -> str:
    if positive_count == total_count:
        return "common"
    if positive_count == 1:
        return "unique"
    return "partial"


def _definition_expression(definition: DerivedDefinition) -> str:
    clauses = [
        (f"{item.requirement.value.upper()}({item.role}:{item.input_id})")
        for item in definition.inputs
    ]
    filters = definition.filters
    filter_values = []
    if filters.position is not None:
        filter_values.append(f"POS={filters.position}")
    if filters.alt:
        filter_values.append(f"ALT={filters.alt}")
    filter_values.extend(
        f"AF {item.operator} {item.threshold}" for item in filters.af_rules
    )
    filter_values.extend(
        f"{item.field}={item.value}" for item in filters.metadata_filters
    )
    expression = " AND ".join(clauses)
    if filter_values:
        expression += " FILTER " + " AND ".join(filter_values)
    return expression


def _unique(values: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(str(value) for value in values if value))
