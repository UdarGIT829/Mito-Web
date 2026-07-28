"""Perspective-scoped workflows built on the catalog repository."""

from __future__ import annotations

from dataclasses import dataclass

from .models import (
    CatalogSample,
    Cohort,
    CohortType,
    Dataset,
    SampleGroup,
    SampleGroupSnapshot,
    StudyPerspective,
    Visibility,
)
from .repository import (
    CatalogAccessError,
    CatalogNotFoundError,
    CatalogRepository,
)


@dataclass(frozen=True)
class PerspectiveWorkspace:
    """Request-scoped view of one perspective's available catalog state."""

    perspective: StudyPerspective
    datasets: tuple[Dataset, ...]
    accessible_cohorts: tuple[Cohort, ...]


@dataclass(frozen=True)
class DatasetWorkspace:
    """Dataset together with its attached and system-managed cohorts."""

    dataset: Dataset
    cohorts: tuple[Cohort, ...]
    derived_results_cohort: Cohort


class CatalogWorkspaceService:
    """Compose repository calls into complete perspective workflows."""

    def __init__(self, repository: CatalogRepository) -> None:
        self.repository = repository

    def create_perspective(
        self,
        name: str,
        *,
        perspective_id: str | None = None,
        description: str = "",
        visibility: Visibility = Visibility.PRIVATE,
    ) -> PerspectiveWorkspace:
        perspective = self.repository.create_perspective(
            name,
            perspective_id=perspective_id,
            description=description,
            visibility=visibility,
        )
        return self.select_perspective(perspective.id)

    def select_perspective(
        self,
        perspective_id: str,
    ) -> PerspectiveWorkspace:
        """Resolve perspective state without setting process-global context."""
        perspective = self.repository.get_perspective(perspective_id)
        if perspective is None:
            raise CatalogNotFoundError(
                f"Perspective not found: {perspective_id}"
            )
        return PerspectiveWorkspace(
            perspective=perspective,
            datasets=tuple(
                self.repository.list_datasets(perspective.id)
            ),
            accessible_cohorts=tuple(
                self.repository.list_cohorts(
                    perspective_id=perspective.id,
                    usable_only=True,
                )
            ),
        )

    def create_dataset(
        self,
        perspective_id: str,
        name: str,
        *,
        cohort_ids: list[str] | tuple[str, ...] = (),
        dataset_id: str | None = None,
        derived_results_cohort_id: str | None = None,
        description: str = "",
        visibility: Visibility = Visibility.PRIVATE,
    ) -> DatasetWorkspace:
        """Create a complete dataset with exactly one results cohort."""
        cohort_ids = _unique_ids(cohort_ids, "Cohort")
        with self.repository.transaction():
            dataset = self.repository.create_dataset(
                perspective_id,
                name,
                dataset_id=dataset_id,
                description=description,
                visibility=visibility,
            )
            results_cohort = self.repository.create_cohort(
                f"{dataset.name} — Derived Results",
                CohortType.DERIVED,
                cohort_id=derived_results_cohort_id,
                owner_perspective_id=perspective_id,
                description=(
                    "System-managed derived results for dataset "
                    f"{dataset.name}."
                ),
                provenance={
                    "system_managed": True,
                    "purpose": "dataset-derived-results",
                    "dataset_id": dataset.id,
                },
                visibility=visibility,
            )
            dataset = self.repository.set_derived_results_cohort(
                dataset.id,
                results_cohort.id,
                acting_perspective_id=perspective_id,
            )
            for display_order, cohort_id in enumerate(cohort_ids):
                self.repository.add_cohort_to_dataset(
                    dataset.id,
                    cohort_id,
                    acting_perspective_id=perspective_id,
                    display_order=display_order,
                )
            self.repository.add_cohort_to_dataset(
                dataset.id,
                results_cohort.id,
                acting_perspective_id=perspective_id,
                display_order=len(cohort_ids),
            )

        return DatasetWorkspace(
            dataset=dataset,
            cohorts=tuple(
                self.repository.list_dataset_cohorts(dataset.id)
            ),
            derived_results_cohort=results_cohort,
        )

    def open_dataset(
        self,
        dataset_id: str,
        *,
        acting_perspective_id: str,
    ) -> DatasetWorkspace:
        dataset = self.repository.get_dataset(dataset_id)
        if dataset is None:
            raise CatalogNotFoundError(f"Dataset not found: {dataset_id}")
        if dataset.perspective_id != acting_perspective_id:
            raise CatalogAccessError(
                f"Perspective {acting_perspective_id!r} does not own "
                f"dataset {dataset_id!r}."
            )
        if dataset.derived_results_cohort_id is None:
            raise RuntimeError(
                f"Dataset {dataset_id!r} has no Derived Results cohort."
            )
        results_cohort = self.repository.get_cohort(
            dataset.derived_results_cohort_id
        )
        if results_cohort is None:
            raise RuntimeError(
                f"Dataset {dataset_id!r} has a missing results cohort."
            )
        return DatasetWorkspace(
            dataset=dataset,
            cohorts=tuple(
                self.repository.list_dataset_cohorts(dataset.id)
            ),
            derived_results_cohort=results_cohort,
        )

    def dataset_samples(
        self,
        dataset_id: str,
        *,
        acting_perspective_id: str,
    ) -> tuple[CatalogSample, ...]:
        return tuple(
            self.repository.list_dataset_samples(
                dataset_id,
                acting_perspective_id=acting_perspective_id,
            )
        )

    def create_sample_group(
        self,
        perspective_id: str,
        dataset_id: str,
        name: str,
        sample_ids: list[str] | tuple[str, ...],
        *,
        group_id: str | None = None,
        description: str = "",
        visibility: Visibility = Visibility.PRIVATE,
    ) -> SampleGroup:
        return self.repository.create_sample_group(
            perspective_id,
            dataset_id,
            name,
            sample_ids,
            group_id=group_id,
            description=description,
            visibility=visibility,
        )

    def update_sample_group(
        self,
        group_id: str,
        *,
        acting_perspective_id: str,
        name: str | None = None,
        description: str | None = None,
        visibility: Visibility | None = None,
        sample_ids: list[str] | tuple[str, ...] | None = None,
        expected_membership_fingerprint: str | None = None,
    ) -> SampleGroup:
        return self.repository.update_sample_group(
            group_id,
            acting_perspective_id=acting_perspective_id,
            name=name,
            description=description,
            visibility=visibility,
            sample_ids=sample_ids,
            expected_membership_fingerprint=(
                expected_membership_fingerprint
            ),
        )

    def snapshot_sample_group(
        self,
        group_id: str,
        *,
        acting_perspective_id: str,
    ) -> SampleGroupSnapshot:
        return self.repository.snapshot_sample_group(
            group_id,
            acting_perspective_id=acting_perspective_id,
        )


def _unique_ids(
    values: list[str] | tuple[str, ...],
    label: str,
) -> tuple[str, ...]:
    normalized = tuple(str(value or "").strip() for value in values)
    if any(not value for value in normalized):
        raise ValueError(f"{label} IDs cannot be empty.")
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"{label} IDs cannot contain duplicates.")
    return normalized
