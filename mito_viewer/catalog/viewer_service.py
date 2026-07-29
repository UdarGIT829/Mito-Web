"""Minimal adapter between the legacy viewer and the durable catalog."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Mapping

import vcf_parser
from mito_viewer.domain import (
    DerivedSample,
    MutationFilters,
    SampleAlleleCall,
)

from .derived_models import (
    DERIVED_COMPARISON_STATUSES,
    DerivedComparison,
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

    def __init__(
        self,
        catalog: CatalogRepository,
        *,
        source_database_paths: Mapping[str, str | Path] | None = None,
    ) -> None:
        self.catalog = catalog
        self.workspaces = CatalogWorkspaceService(catalog)
        self.derived = DerivedCatalogRepository(catalog)
        self.analysis = DerivedAnalysisService(
            catalog,
            source_database_paths=source_database_paths,
        )

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
    ) -> dict:
        workspace = self.workspaces.open_dataset(
            dataset_id,
            acting_perspective_id=perspective_id,
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
                }
                for cohort in workspace.cohorts
            ],
            "groups": [
                self._group_payload(
                    group.id,
                    perspective_id=perspective_id,
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
    ) -> dict:
        dataset_samples = {
            sample.id: sample
            for sample in self.workspaces.dataset_samples(
                dataset_id,
                acting_perspective_id=perspective_id,
            )
        }
        sample_ids = []
        for viewer_sample_id in viewer_sample_ids:
            viewer_sample_id = str(viewer_sample_id)
            sample = dataset_samples.get(viewer_sample_id)
            if sample is None:
                raise CatalogNotFoundError(
                    "Catalog sample is not available in the Dataset: "
                    f"{viewer_sample_id}"
                )
            if sample.sample_type is not SampleType.OBSERVED:
                raise ValueError(
                    "Sample groups can contain observed samples only."
                )
            sample_ids.append(sample.id)
        group = self.workspaces.create_sample_group(
            perspective_id,
            dataset_id,
            name,
            sample_ids,
        )
        return self._group_payload(
            group.id,
            perspective_id=perspective_id,
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
        sample_statuses: dict[str, set[str]],
        global_statuses: Iterable[str],
        filters: MutationFilters,
    ) -> dict:
        inputs = []
        dataset_samples = {
            sample.id: sample
            for sample in self.workspaces.dataset_samples(
                dataset_id,
                acting_perspective_id=perspective_id,
            )
        }
        viewer_sample_ids = tuple(str(item) for item in viewer_sample_ids)
        normalized_sample_statuses = tuple(
            frozenset(
                sample_statuses.get(viewer_sample_id, {"present"})
            )
            for viewer_sample_id in viewer_sample_ids
        )
        for viewer_sample_id in viewer_sample_ids:
            kind, input_id = self._derived_input_for_viewer_sample(
                viewer_sample_id,
                dataset_samples=dataset_samples,
            )
            inputs.append(
                DerivedInput(
                    kind=kind,
                    input_id=input_id,
                    role="comparison-input",
                    requirement=PresenceRequirement.ANY,
                )
            )
        selected_global_statuses = frozenset(global_statuses)
        if not selected_global_statuses:
            selected_global_statuses = (
                DERIVED_COMPARISON_STATUSES - {"__none__"}
            )

        result = self.analysis.create_derived_sample(
            perspective_id=perspective_id,
            dataset_id=dataset_id,
            name=name,
            definition=DerivedDefinition(
                inputs=tuple(inputs),
                filters=filters,
                version=2,
                comparison=DerivedComparison(
                    statuses=selected_global_statuses,
                    input_statuses=normalized_sample_statuses,
                ),
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
    ) -> Cohort | None:
        resolved_path = Path(database_path).resolve()
        for cohort in self.catalog.list_cohorts():
            if cohort.cohort_type is not CohortType.SOURCE:
                continue
            provenance_path = cohort.provenance.get("source_path")
            if cohort.source_database_identifier == database_id:
                return cohort
            if (
                provenance_path
                and Path(provenance_path).resolve() == resolved_path
            ):
                return cohort
        return None

    def _derived_input_for_viewer_sample(
        self,
        viewer_sample_id: str,
        *,
        dataset_samples: dict[str, CatalogSample] | None = None,
    ) -> tuple[DerivedInputKind, str]:
        catalog_sample = (dataset_samples or {}).get(viewer_sample_id)
        if catalog_sample is not None:
            if catalog_sample.sample_type is SampleType.OBSERVED:
                return DerivedInputKind.SAMPLE, catalog_sample.id
            record = self.derived.find_by_catalog_sample(
                catalog_sample.id
            )
            if record is None or record.current_run_id is None:
                raise CatalogNotFoundError(
                    "Durable derived sample not found for catalog sample: "
                    f"{viewer_sample_id}"
                )
            return DerivedInputKind.RUN, record.current_run_id
        raise CatalogNotFoundError(
            "Catalog sample is not available in the Dataset: "
            f"{viewer_sample_id}"
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
    ) -> dict:
        group = self.catalog.get_sample_group(
            group_id,
            acting_perspective_id=perspective_id,
        )
        members = self.catalog.sample_group_members(
            group_id,
            acting_perspective_id=perspective_id,
        )
        visible_ids = [sample.id for sample in members]
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
            "catalog_sample_id": record.catalog_sample_id,
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


def _definition_expression(definition: DerivedDefinition) -> str:
    if definition.comparison is None:
        clauses = [
            (f"{item.requirement.value.upper()}({item.role}:{item.input_id})")
            for item in definition.inputs
        ]
        expression = " AND ".join(clauses)
    else:
        expression = _comparison_expression(definition)
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
    if filter_values:
        expression += " FILTER " + " AND ".join(filter_values)
    return expression


def _comparison_expression(definition: DerivedDefinition) -> str:
    comparison = definition.comparison
    if comparison is None:
        return ""
    present = []
    unique = []
    absent = []
    for item, statuses in zip(
        definition.inputs,
        comparison.input_statuses,
    ):
        if "present" in statuses:
            present.append(item.input_id)
        if "unique" in statuses:
            unique.append(item.input_id)
        if "not_in" in statuses:
            absent.append(item.input_id)

    branches = []
    if present:
        branches.append(f"PRESENT_ALL({', '.join(present)})")
    if unique:
        branches.append(f"UNIQUE_ANY({', '.join(unique)})")
    expression = " OR ".join(branches) if branches else "ANY_SELECTED_INPUT"
    if len(branches) > 1:
        expression = f"({expression})"
    if absent:
        expression += f" AND NOT_IN_ALL({', '.join(absent)})"
    statuses = ", ".join(sorted(comparison.statuses))
    return f"({expression}) AND STATUS_IN({statuses})"


def _unique(values: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(str(value) for value in values if value))
