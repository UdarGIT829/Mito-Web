"""Persistence operations for durable derived samples and immutable runs."""

from __future__ import annotations

import json
import sqlite3
import uuid
from datetime import datetime, timezone

from mito_viewer.domain import AlleleKey

from .derived_models import (
    AlleleEvidence,
    DerivedCalculationResult,
    DerivedDefinition,
    DerivedRunRecord,
    DerivedRunStatus,
    DerivedSampleRecord,
    DerivedStaleness,
    EvidenceStatus,
    MaterializedAllele,
    PendingMaterializedAllele,
    ResolvedParentInput,
    ResolvedSampleInput,
    allele_output_fingerprint,
)
from .models import SampleType, Visibility
from .repository import (
    CatalogAccessError,
    CatalogNotFoundError,
    CatalogRepository,
)


class DerivedCycleError(ValueError):
    """Raised when a run dependency would introduce a cycle."""


class DerivedCatalogRepository:
    """Store and read derived state using one catalog connection."""

    def __init__(self, catalog: CatalogRepository) -> None:
        self.catalog = catalog
        self.connection = catalog.connection

    def persist_completed_run(
        self,
        *,
        perspective_id: str,
        dataset_id: str,
        name: str,
        definition: DerivedDefinition,
        input_snapshot_fingerprint: str,
        reference_id: str,
        normalization_version: str,
        sample_inputs: tuple[ResolvedSampleInput, ...],
        parent_inputs: tuple[ResolvedParentInput, ...],
        materialized_alleles: tuple[PendingMaterializedAllele, ...],
        description: str = "",
        visibility: Visibility = Visibility.PRIVATE,
        derived_sample_id: str | None = None,
        catalog_sample_id: str | None = None,
        run_id: str | None = None,
    ) -> DerivedCalculationResult:
        name = str(name or "").strip()
        if not name:
            raise ValueError("Derived sample name cannot be empty.")
        definition = (
            definition
            if isinstance(definition, DerivedDefinition)
            else DerivedDefinition.from_dict(definition)
        )
        visibility = Visibility(visibility)
        run_id = run_id or _new_id()
        derived_sample_id = derived_sample_id or _new_id()
        output_fingerprint = allele_output_fingerprint(
            tuple(item.allele for item in materialized_alleles)
        )
        now = _utc_now()

        with self.catalog.transaction():
            dataset = self.catalog.get_dataset(dataset_id)
            if dataset is None:
                raise CatalogNotFoundError(f"Dataset not found: {dataset_id}")
            if dataset.perspective_id != perspective_id:
                raise CatalogAccessError(
                    f"Perspective {perspective_id!r} does not own "
                    f"dataset {dataset_id!r}."
                )
            if dataset.derived_results_cohort_id is None:
                raise ValueError(
                    f"Dataset {dataset_id!r} has no Derived Results cohort."
                )

            existing = self._find_derived_sample(derived_sample_id)
            if existing is None:
                catalog_sample_id = catalog_sample_id or _new_id()
                self.catalog.create_catalog_sample(
                    dataset.derived_results_cohort_id,
                    SampleType.DERIVED,
                    name,
                    acting_perspective_id=perspective_id,
                    catalog_sample_id=catalog_sample_id,
                    metadata={
                        "derived_sample_id": derived_sample_id,
                        "dataset_id": dataset_id,
                    },
                )
                self.connection.execute(
                    """
                    INSERT INTO derived_samples(
                        id,
                        catalog_sample_id,
                        perspective_id,
                        dataset_id,
                        name,
                        description,
                        visibility,
                        current_run_id,
                        created_at,
                        updated_at
                    )
                    VALUES (?, ?, ?, ?, ?, ?, ?, NULL, ?, ?)
                    """,
                    (
                        derived_sample_id,
                        catalog_sample_id,
                        perspective_id,
                        dataset_id,
                        name,
                        str(description or "").strip(),
                        visibility.value,
                        now,
                        now,
                    ),
                )
            else:
                self._require_derived_owner(existing, perspective_id)
                if existing.dataset_id != dataset_id:
                    raise ValueError("A derived sample cannot move between datasets.")
                name = existing.name

            self.assert_acyclic_dependencies(
                run_id,
                tuple(item.parent_run_id for item in parent_inputs),
            )
            self._revalidate_inputs(sample_inputs, parent_inputs)
            self.connection.execute(
                """
                INSERT INTO derived_runs(
                    id,
                    derived_sample_id,
                    definition_json,
                    input_snapshot_fingerprint,
                    reference_id,
                    normalization_version,
                    status,
                    output_count,
                    output_fingerprint,
                    error_text,
                    created_at,
                    completed_at
                )
                VALUES (?, ?, ?, ?, ?, ?, 'completed', ?, ?, '', ?, ?)
                """,
                (
                    run_id,
                    derived_sample_id,
                    definition.normalized_json,
                    input_snapshot_fingerprint,
                    reference_id,
                    normalization_version,
                    len(materialized_alleles),
                    output_fingerprint,
                    now,
                    now,
                ),
            )
            self._insert_sample_inputs(run_id, sample_inputs)
            self._insert_parent_inputs(run_id, parent_inputs)
            self._insert_materialized_alleles(
                run_id,
                materialized_alleles,
            )
            self.connection.execute(
                """
                UPDATE derived_samples
                SET current_run_id = ?, updated_at = ?
                WHERE id = ?
                """,
                (run_id, now, derived_sample_id),
            )

        derived_sample = self.get_derived_sample(
            derived_sample_id,
            acting_perspective_id=perspective_id,
        )
        run = self.get_run(run_id)
        if derived_sample is None or run is None:
            raise RuntimeError("Completed derived run could not be reloaded.")
        return DerivedCalculationResult(
            derived_sample=derived_sample,
            run=run,
            alleles=tuple(self.list_run_alleles(run_id)),
        )

    def get_derived_sample(
        self,
        derived_sample_id: str,
        *,
        acting_perspective_id: str,
    ) -> DerivedSampleRecord | None:
        record = self._find_derived_sample(derived_sample_id)
        if record is None:
            return None
        self._require_derived_owner(record, acting_perspective_id)
        return record

    def find_by_catalog_sample(
        self,
        catalog_sample_id: str,
    ) -> DerivedSampleRecord | None:
        row = self.connection.execute(
            """
            SELECT *
            FROM derived_samples
            WHERE catalog_sample_id = ?
            """,
            (catalog_sample_id,),
        ).fetchone()
        return self._derived_sample_from_row(row) if row else None

    def list_derived_samples(
        self,
        dataset_id: str,
        *,
        acting_perspective_id: str,
    ) -> list[DerivedSampleRecord]:
        dataset = self.catalog.get_dataset(dataset_id)
        if dataset is None:
            raise CatalogNotFoundError(f"Dataset not found: {dataset_id}")
        if dataset.perspective_id != acting_perspective_id:
            raise CatalogAccessError(
                f"Perspective {acting_perspective_id!r} does not own "
                f"dataset {dataset_id!r}."
            )
        rows = self.connection.execute(
            """
            SELECT *
            FROM derived_samples
            WHERE dataset_id = ?
            ORDER BY name COLLATE NOCASE, id
            """,
            (dataset_id,),
        ).fetchall()
        return [self._derived_sample_from_row(row) for row in rows]

    def get_run(self, run_id: str) -> DerivedRunRecord | None:
        row = self.connection.execute(
            """
            SELECT *
            FROM derived_runs
            WHERE id = ?
            """,
            (run_id,),
        ).fetchone()
        return self._run_from_row(row) if row else None

    def derived_sample_for_run(
        self,
        run_id: str,
    ) -> DerivedSampleRecord | None:
        row = self.connection.execute(
            """
            SELECT derived_samples.*
            FROM derived_runs
            JOIN derived_samples
                ON derived_samples.id = derived_runs.derived_sample_id
            WHERE derived_runs.id = ?
            """,
            (run_id,),
        ).fetchone()
        return self._derived_sample_from_row(row) if row else None

    def list_runs(
        self,
        derived_sample_id: str,
        *,
        acting_perspective_id: str,
    ) -> list[DerivedRunRecord]:
        sample = self._find_derived_sample(derived_sample_id)
        if sample is None:
            raise CatalogNotFoundError(f"Derived sample not found: {derived_sample_id}")
        self._require_derived_owner(sample, acting_perspective_id)
        rows = self.connection.execute(
            """
            SELECT *
            FROM derived_runs
            WHERE derived_sample_id = ?
            ORDER BY created_at, id
            """,
            (derived_sample_id,),
        ).fetchall()
        return [self._run_from_row(row) for row in rows]

    def list_run_alleles(self, run_id: str) -> list[MaterializedAllele]:
        if self.get_run(run_id) is None:
            raise CatalogNotFoundError(f"Derived run not found: {run_id}")
        rows = self.connection.execute(
            """
            SELECT *
            FROM derived_run_alleles
            WHERE run_id = ?
            ORDER BY position, ref, alt
            """,
            (run_id,),
        ).fetchall()
        return [self._allele_from_row(row) for row in rows]

    def list_run_sample_inputs(
        self,
        run_id: str,
    ) -> list[ResolvedSampleInput]:
        """Return the immutable observed-sample inputs captured for a run."""
        if self.get_run(run_id) is None:
            raise CatalogNotFoundError(f"Derived run not found: {run_id}")
        rows = self.connection.execute(
            """
            SELECT *
            FROM derived_run_input_samples
            WHERE run_id = ?
            ORDER BY input_clause_index, display_order, catalog_sample_id
            """,
            (run_id,),
        ).fetchall()
        return [
            ResolvedSampleInput(
                clause_index=row["input_clause_index"],
                catalog_sample_id=row["catalog_sample_id"],
                source_group_id=row["source_group_id"],
                group_membership_fingerprint=(
                    row["group_membership_fingerprint"]
                ),
                sample_fingerprint=row["sample_fingerprint"],
                source_database_fingerprint=(
                    row["source_database_fingerprint"]
                ),
                input_role=row["input_role"],
                display_order=row["display_order"],
            )
            for row in rows
        ]

    def list_run_parent_inputs(
        self,
        run_id: str,
    ) -> list[ResolvedParentInput]:
        """Return the immutable parent-run inputs captured for a run."""
        if self.get_run(run_id) is None:
            raise CatalogNotFoundError(f"Derived run not found: {run_id}")
        rows = self.connection.execute(
            """
            SELECT *
            FROM derived_run_parent_runs
            WHERE run_id = ?
            ORDER BY input_clause_index, display_order, parent_run_id
            """,
            (run_id,),
        ).fetchall()
        return [
            ResolvedParentInput(
                clause_index=row["input_clause_index"],
                parent_run_id=row["parent_run_id"],
                source_group_id=row["source_group_id"],
                group_membership_fingerprint=(
                    row["group_membership_fingerprint"]
                ),
                input_role=row["input_role"],
                display_order=row["display_order"],
            )
            for row in rows
        ]

    def list_allele_evidence(
        self,
        derived_allele_id: str,
    ) -> list[AlleleEvidence]:
        rows = self.connection.execute(
            """
            SELECT *
            FROM derived_allele_sources
            WHERE derived_allele_id = ?
            ORDER BY input_role, evaluation_status, id
            """,
            (derived_allele_id,),
        ).fetchall()
        return [self._evidence_from_row(row) for row in rows]

    def staleness(self, run_id: str) -> DerivedStaleness:
        if self.get_run(run_id) is None:
            raise CatalogNotFoundError(f"Derived run not found: {run_id}")
        stale_groups: set[str] = set()
        missing_groups: set[str] = set()
        changed_samples: set[str] = set()
        changed_cohorts: set[str] = set()

        sample_rows = self.connection.execute(
            """
            SELECT
                inputs.*,
                catalog_samples.source_fingerprint AS current_sample_fingerprint,
                catalog_samples.cohort_id,
                cohorts.source_database_fingerprint AS current_database_fingerprint
            FROM derived_run_input_samples inputs
            JOIN catalog_samples
                ON catalog_samples.id = inputs.catalog_sample_id
            JOIN cohorts ON cohorts.id = catalog_samples.cohort_id
            WHERE inputs.run_id = ?
            """,
            (run_id,),
        ).fetchall()
        for row in sample_rows:
            if row["sample_fingerprint"] != row["current_sample_fingerprint"]:
                changed_samples.add(row["catalog_sample_id"])
            if (
                row["source_database_fingerprint"]
                != row["current_database_fingerprint"]
            ):
                changed_cohorts.add(row["cohort_id"])
            self._collect_group_staleness(
                row,
                stale_groups,
                missing_groups,
            )

        parent_rows = self.connection.execute(
            """
            SELECT
                inputs.*,
                derived_samples.current_run_id
            FROM derived_run_parent_runs inputs
            JOIN derived_runs parent
                ON parent.id = inputs.parent_run_id
            JOIN derived_samples
                ON derived_samples.id = parent.derived_sample_id
            WHERE inputs.run_id = ?
            """,
            (run_id,),
        ).fetchall()
        parent_updates = {
            row["parent_run_id"]
            for row in parent_rows
            if row["current_run_id"] != row["parent_run_id"]
        }
        for row in parent_rows:
            self._collect_group_staleness(
                row,
                stale_groups,
                missing_groups,
            )

        return DerivedStaleness(
            run_id=run_id,
            stale_group_ids=tuple(sorted(stale_groups)),
            missing_group_ids=tuple(sorted(missing_groups)),
            changed_sample_ids=tuple(sorted(changed_samples)),
            changed_source_cohort_ids=tuple(sorted(changed_cohorts)),
            parent_runs_with_updates=tuple(sorted(parent_updates)),
        )

    def assert_acyclic_dependencies(
        self,
        run_id: str,
        parent_run_ids: tuple[str, ...],
    ) -> None:
        if len(set(parent_run_ids)) != len(parent_run_ids):
            raise ValueError("Parent run inputs cannot contain duplicates.")
        if run_id in parent_run_ids:
            raise DerivedCycleError("A derived run cannot depend on itself.")
        for parent_run_id in parent_run_ids:
            parent = self.get_run(parent_run_id)
            if parent is None:
                raise CatalogNotFoundError(
                    f"Parent derived run not found: {parent_run_id}"
                )
            if parent.status is not DerivedRunStatus.COMPLETED:
                raise ValueError(f"Parent run {parent_run_id!r} is not completed.")

        if not parent_run_ids:
            return
        placeholders = ",".join("?" for _ in parent_run_ids)
        row = self.connection.execute(
            f"""
            WITH RECURSIVE ancestors(run_id) AS (
                SELECT parent_run_id
                FROM derived_run_parent_runs
                WHERE run_id IN ({placeholders})
                UNION
                SELECT links.parent_run_id
                FROM derived_run_parent_runs links
                JOIN ancestors ON ancestors.run_id = links.run_id
            )
            SELECT 1
            FROM ancestors
            WHERE run_id = ?
            LIMIT 1
            """,
            (*parent_run_ids, run_id),
        ).fetchone()
        if row:
            raise DerivedCycleError("Derived run dependency would introduce a cycle.")

    def _revalidate_inputs(
        self,
        sample_inputs: tuple[ResolvedSampleInput, ...],
        parent_inputs: tuple[ResolvedParentInput, ...],
    ) -> None:
        for item in sample_inputs:
            sample = self.catalog.get_catalog_sample(item.catalog_sample_id)
            if sample is None:
                raise CatalogNotFoundError(
                    f"Catalog sample not found: {item.catalog_sample_id}"
                )
            cohort = self.catalog.get_cohort(sample.cohort_id)
            if cohort is None:
                raise CatalogNotFoundError(f"Cohort not found: {sample.cohort_id}")
            if sample.source_fingerprint != item.sample_fingerprint:
                raise ValueError(
                    f"Catalog sample {sample.id!r} changed during calculation."
                )
            if cohort.source_database_fingerprint != item.source_database_fingerprint:
                raise ValueError(
                    f"Source cohort {cohort.id!r} changed during calculation."
                )
            self._require_group_fingerprint(
                item.source_group_id,
                item.group_membership_fingerprint,
            )
        for item in parent_inputs:
            self._require_group_fingerprint(
                item.source_group_id,
                item.group_membership_fingerprint,
            )

    def _require_group_fingerprint(
        self,
        group_id: str | None,
        expected_fingerprint: str,
    ) -> None:
        if group_id is None:
            return
        row = self.connection.execute(
            """
            SELECT membership_fingerprint
            FROM sample_groups
            WHERE id = ?
            """,
            (group_id,),
        ).fetchone()
        if row is None:
            raise ValueError(f"Input sample group {group_id!r} was deleted.")
        if row["membership_fingerprint"] != expected_fingerprint:
            raise ValueError(
                f"Input sample group {group_id!r} changed during calculation."
            )

    def _insert_sample_inputs(
        self,
        run_id: str,
        values: tuple[ResolvedSampleInput, ...],
    ) -> None:
        self.connection.executemany(
            """
            INSERT INTO derived_run_input_samples(
                run_id,
                input_clause_index,
                catalog_sample_id,
                source_group_id,
                group_membership_fingerprint,
                sample_fingerprint,
                source_database_fingerprint,
                input_role,
                display_order
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            [
                (
                    run_id,
                    item.clause_index,
                    item.catalog_sample_id,
                    item.source_group_id,
                    item.group_membership_fingerprint,
                    item.sample_fingerprint,
                    item.source_database_fingerprint,
                    item.input_role,
                    item.display_order,
                )
                for item in values
            ],
        )

    def _insert_parent_inputs(
        self,
        run_id: str,
        values: tuple[ResolvedParentInput, ...],
    ) -> None:
        self.connection.executemany(
            """
            INSERT INTO derived_run_parent_runs(
                run_id,
                input_clause_index,
                parent_run_id,
                source_group_id,
                group_membership_fingerprint,
                input_role,
                display_order
            )
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            [
                (
                    run_id,
                    item.clause_index,
                    item.parent_run_id,
                    item.source_group_id,
                    item.group_membership_fingerprint,
                    item.input_role,
                    item.display_order,
                )
                for item in values
            ],
        )

    def _insert_materialized_alleles(
        self,
        run_id: str,
        values: tuple[PendingMaterializedAllele, ...],
    ) -> None:
        for value in values:
            allele_id = _allele_id(run_id, value.allele)
            self.connection.execute(
                """
                INSERT INTO derived_run_alleles(
                    id,
                    run_id,
                    position,
                    ref,
                    alt,
                    result_metadata_json
                )
                VALUES (?, ?, ?, ?, ?, ?)
                """,
                (
                    allele_id,
                    run_id,
                    value.allele.position,
                    value.allele.ref,
                    value.allele.alt,
                    _json_dump(value.metadata),
                ),
            )
            self.connection.executemany(
                """
                INSERT INTO derived_allele_sources(
                    id,
                    derived_allele_id,
                    source_catalog_sample_id,
                    parent_derived_allele_id,
                    original_catalog_sample_id,
                    source_mutation_id,
                    source_alt_index,
                    input_role,
                    evaluation_status,
                    af,
                    af_text,
                    filter_text,
                    metadata_json,
                    lineage_json
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                [
                    (
                        _new_id(),
                        allele_id,
                        evidence.source_catalog_sample_id,
                        evidence.parent_derived_allele_id,
                        evidence.original_catalog_sample_id,
                        evidence.source_mutation_id,
                        evidence.source_alt_index,
                        evidence.input_role,
                        evidence.evaluation_status.value,
                        evidence.af,
                        evidence.af_text,
                        evidence.filter_text,
                        _json_dump(evidence.metadata),
                        _json_dump(evidence.lineage),
                    )
                    for evidence in value.evidence
                ],
            )

    def _find_derived_sample(
        self,
        derived_sample_id: str,
    ) -> DerivedSampleRecord | None:
        row = self.connection.execute(
            """
            SELECT *
            FROM derived_samples
            WHERE id = ?
            """,
            (derived_sample_id,),
        ).fetchone()
        return self._derived_sample_from_row(row) if row else None

    def _collect_group_staleness(
        self,
        row: sqlite3.Row,
        stale_groups: set[str],
        missing_groups: set[str],
    ) -> None:
        group_id = row["source_group_id"]
        if group_id is None:
            return
        group = self.connection.execute(
            """
            SELECT membership_fingerprint
            FROM sample_groups
            WHERE id = ?
            """,
            (group_id,),
        ).fetchone()
        if group is None:
            missing_groups.add(group_id)
        elif group["membership_fingerprint"] != row["group_membership_fingerprint"]:
            stale_groups.add(group_id)

    @staticmethod
    def _require_derived_owner(
        sample: DerivedSampleRecord,
        perspective_id: str,
    ) -> None:
        if sample.perspective_id != perspective_id:
            raise CatalogAccessError(
                f"Perspective {perspective_id!r} does not own "
                f"derived sample {sample.id!r}."
            )

    @staticmethod
    def _derived_sample_from_row(row: sqlite3.Row) -> DerivedSampleRecord:
        return DerivedSampleRecord(
            id=row["id"],
            catalog_sample_id=row["catalog_sample_id"],
            perspective_id=row["perspective_id"],
            dataset_id=row["dataset_id"],
            name=row["name"],
            description=row["description"],
            visibility=Visibility(row["visibility"]),
            current_run_id=row["current_run_id"],
            created_at=row["created_at"],
            updated_at=row["updated_at"],
        )

    @staticmethod
    def _run_from_row(row: sqlite3.Row) -> DerivedRunRecord:
        return DerivedRunRecord(
            id=row["id"],
            derived_sample_id=row["derived_sample_id"],
            definition=DerivedDefinition.from_dict(json.loads(row["definition_json"])),
            input_snapshot_fingerprint=row["input_snapshot_fingerprint"],
            reference_id=row["reference_id"],
            normalization_version=row["normalization_version"],
            status=DerivedRunStatus(row["status"]),
            output_count=row["output_count"],
            output_fingerprint=row["output_fingerprint"],
            error_text=row["error_text"],
            created_at=row["created_at"],
            completed_at=row["completed_at"],
        )

    @staticmethod
    def _allele_from_row(row: sqlite3.Row) -> MaterializedAllele:
        return MaterializedAllele(
            id=row["id"],
            run_id=row["run_id"],
            allele=AlleleKey(
                position=row["position"],
                ref=row["ref"],
                alt=row["alt"],
            ),
            metadata=_json_load(row["result_metadata_json"]),
        )

    @staticmethod
    def _evidence_from_row(row: sqlite3.Row) -> AlleleEvidence:
        return AlleleEvidence(
            id=row["id"],
            derived_allele_id=row["derived_allele_id"],
            source_catalog_sample_id=row["source_catalog_sample_id"],
            parent_derived_allele_id=row["parent_derived_allele_id"],
            original_catalog_sample_id=row["original_catalog_sample_id"],
            source_mutation_id=row["source_mutation_id"],
            source_alt_index=row["source_alt_index"],
            input_role=row["input_role"],
            evaluation_status=EvidenceStatus(row["evaluation_status"]),
            af=row["af"],
            af_text=row["af_text"],
            filter_text=row["filter_text"],
            metadata=_json_load(row["metadata_json"]),
            lineage=_json_load(row["lineage_json"]),
        )


def _allele_id(run_id: str, allele: AlleleKey) -> str:
    return str(
        uuid.uuid5(
            uuid.NAMESPACE_URL,
            "urn:mito-viewer:derived-allele:"
            f"{run_id}:{allele.position}:{allele.ref}:{allele.alt}",
        )
    )


def _new_id() -> str:
    return str(uuid.uuid4())


def _utc_now() -> str:
    return (
        datetime.now(timezone.utc)
        .isoformat(timespec="microseconds")
        .replace("+00:00", "Z")
    )


def _json_dump(value: dict) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    )


def _json_load(value: str) -> dict:
    loaded = json.loads(value or "{}")
    if not isinstance(loaded, dict):
        raise ValueError("Derived JSON fields must contain objects.")
    return loaded
