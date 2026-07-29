"""Calculation service for durable, provenance-preserving derived runs."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping

from mito_viewer.domain import AlleleKey
from mito_viewer.domain.comparison import (
    comparison_status,
    sample_filters_match,
)
from mito_viewer.repositories import StudyRepository

from .derived_models import (
    DerivedCalculationResult,
    DerivedDefinition,
    DerivedInput,
    DerivedInputKind,
    DerivedRunStatus,
    EvidenceStatus,
    PendingAlleleEvidence,
    PendingMaterializedAllele,
    PresenceRequirement,
    ResolvedParentInput,
    ResolvedSampleInput,
    snapshot_fingerprint,
)
from .derived_repository import DerivedCatalogRepository
from .models import CatalogSample, Cohort, SampleType, Visibility
from .registration import sha256_file
from .repository import (
    CatalogAccessError,
    CatalogNotFoundError,
    CatalogRepository,
)


class ReferenceCompatibilityError(ValueError):
    """Raised when inputs cannot share canonical allele coordinates."""


@dataclass
class _ResolvedClause:
    definition: DerivedInput
    member_sets: list[set[AlleleKey]]


class DerivedAnalysisService:
    """Resolve inputs, calculate allele sets, and persist immutable runs."""

    def __init__(
        self,
        catalog: CatalogRepository,
        *,
        source_database_paths: Mapping[str, str | Path] | None = None,
        maximum_output_count: int = 100_000,
    ) -> None:
        self.catalog = catalog
        self.derived_repository = DerivedCatalogRepository(catalog)
        self.source_database_paths = {
            str(identifier): Path(path).resolve()
            for identifier, path in (source_database_paths or {}).items()
        }
        self.maximum_output_count = int(maximum_output_count)
        if self.maximum_output_count < 1:
            raise ValueError("Maximum output count must be positive.")

    def create_derived_sample(
        self,
        *,
        perspective_id: str,
        dataset_id: str,
        name: str,
        definition: DerivedDefinition,
        description: str = "",
        visibility: Visibility = Visibility.PRIVATE,
        derived_sample_id: str | None = None,
        catalog_sample_id: str | None = None,
        run_id: str | None = None,
    ) -> DerivedCalculationResult:
        if derived_sample_id is not None:
            existing = self.derived_repository.get_derived_sample(
                derived_sample_id,
                acting_perspective_id=perspective_id,
            )
            if existing is not None:
                raise ValueError(f"Derived sample already exists: {derived_sample_id}")
        return self._calculate_and_persist(
            perspective_id=perspective_id,
            dataset_id=dataset_id,
            name=name,
            definition=definition,
            description=description,
            visibility=visibility,
            derived_sample_id=derived_sample_id,
            catalog_sample_id=catalog_sample_id,
            run_id=run_id,
        )

    def recalculate(
        self,
        derived_sample_id: str,
        *,
        perspective_id: str,
        definition: DerivedDefinition | None = None,
        run_id: str | None = None,
    ) -> DerivedCalculationResult:
        sample = self.derived_repository.get_derived_sample(
            derived_sample_id,
            acting_perspective_id=perspective_id,
        )
        if sample is None:
            raise CatalogNotFoundError(f"Derived sample not found: {derived_sample_id}")
        if definition is None:
            if sample.current_run_id is None:
                raise ValueError(
                    f"Derived sample {derived_sample_id!r} has no prior run."
                )
            current_run = self.derived_repository.get_run(sample.current_run_id)
            if current_run is None:
                raise RuntimeError(f"Current run is missing: {sample.current_run_id}")
            definition = current_run.definition
        return self._calculate_and_persist(
            perspective_id=perspective_id,
            dataset_id=sample.dataset_id,
            name=sample.name,
            definition=definition,
            description=sample.description,
            visibility=sample.visibility,
            derived_sample_id=sample.id,
            catalog_sample_id=sample.catalog_sample_id,
            run_id=run_id,
        )

    def _calculate_and_persist(
        self,
        *,
        perspective_id: str,
        dataset_id: str,
        name: str,
        definition: DerivedDefinition,
        description: str,
        visibility: Visibility,
        derived_sample_id: str | None,
        catalog_sample_id: str | None,
        run_id: str | None,
    ) -> DerivedCalculationResult:
        definition = (
            definition
            if isinstance(definition, DerivedDefinition)
            else DerivedDefinition.from_dict(definition)
        )
        available_samples = {
            sample.id: sample
            for sample in self.catalog.list_dataset_samples(
                dataset_id,
                acting_perspective_id=perspective_id,
            )
        }
        (
            sample_inputs,
            parent_inputs,
            clauses,
            references,
            group_snapshots,
        ) = self._resolve_inputs(
            definition,
            perspective_id=perspective_id,
            dataset_id=dataset_id,
            available_samples=available_samples,
        )
        reference_id, normalization_version = _compatible_reference(references)
        if parent_inputs and (
            definition.filters.af_rules or definition.filters.metadata_filters
        ):
            raise ValueError(
                "AF and metadata filters cannot be re-applied to an exact "
                "parent run; filter its observed inputs before chaining."
            )

        observed_evidence, observed_sets = self._load_observed_evidence(
            sample_inputs,
            available_samples,
            definition,
        )
        parent_alleles, parent_sets = self._load_parent_alleles(
            parent_inputs,
            definition,
        )
        clause_lookup = self._populate_clause_sets(
            clauses,
            sample_inputs,
            parent_inputs,
            observed_sets,
            parent_sets,
        )
        result_alleles = self._evaluate_clauses(
            definition,
            clause_lookup,
        )
        if len(result_alleles) > self.maximum_output_count:
            raise ValueError(
                "Derived result exceeds the configured materialization "
                f"limit of {self.maximum_output_count} alleles."
            )
        materialized = tuple(
            self._materialize_allele(
                allele,
                definition=definition,
                clauses=clause_lookup,
                sample_inputs=sample_inputs,
                parent_inputs=parent_inputs,
                observed_evidence=observed_evidence,
                parent_alleles=parent_alleles,
                available_samples=available_samples,
            )
            for allele in sorted(result_alleles)
        )

        parent_output_fingerprints = {
            item.parent_run_id: self._required_run(
                item.parent_run_id
            ).output_fingerprint
            for item in parent_inputs
        }
        input_snapshot = {
            "definition_fingerprint": definition.fingerprint,
            "sample_inputs": [asdict(item) for item in sample_inputs],
            "parent_inputs": [
                {
                    **asdict(item),
                    "parent_output_fingerprint": (
                        parent_output_fingerprints[item.parent_run_id]
                    ),
                }
                for item in parent_inputs
            ],
            "group_snapshots": group_snapshots,
            "reference_id": reference_id,
            "normalization_version": normalization_version,
        }
        return self.derived_repository.persist_completed_run(
            perspective_id=perspective_id,
            dataset_id=dataset_id,
            name=name,
            definition=definition,
            input_snapshot_fingerprint=snapshot_fingerprint(input_snapshot),
            reference_id=reference_id,
            normalization_version=normalization_version,
            sample_inputs=sample_inputs,
            parent_inputs=parent_inputs,
            materialized_alleles=materialized,
            description=description,
            visibility=visibility,
            derived_sample_id=derived_sample_id,
            catalog_sample_id=catalog_sample_id,
            run_id=run_id,
        )

    def _resolve_inputs(
        self,
        definition: DerivedDefinition,
        *,
        perspective_id: str,
        dataset_id: str,
        available_samples: dict[str, CatalogSample],
    ) -> tuple[
        tuple[ResolvedSampleInput, ...],
        tuple[ResolvedParentInput, ...],
        list[_ResolvedClause],
        list[tuple[str, str]],
        list[dict],
    ]:
        sample_inputs = []
        parent_inputs = []
        clauses = []
        references = []
        group_snapshots = []

        for clause_index, input_definition in enumerate(definition.inputs):
            clause = _ResolvedClause(input_definition, [])
            if input_definition.kind is DerivedInputKind.GROUP:
                snapshot = self.catalog.snapshot_sample_group(
                    input_definition.input_id,
                    acting_perspective_id=perspective_id,
                )
                if snapshot.dataset_id != dataset_id:
                    raise ValueError(
                        f"Group {snapshot.group_id!r} belongs to another dataset."
                    )
                if not snapshot.sample_ids:
                    raise ValueError(f"Input group {snapshot.group_id!r} is empty.")
                group_snapshots.append(
                    {
                        "group_id": snapshot.group_id,
                        "membership_fingerprint": (snapshot.membership_fingerprint),
                        "sample_ids": list(snapshot.sample_ids),
                    }
                )
                for display_order, sample_id in enumerate(snapshot.sample_ids):
                    sample = available_samples.get(sample_id)
                    if sample is None:
                        raise ValueError(
                            f"Group sample {sample_id!r} is no longer "
                            "available in the dataset."
                        )
                    if sample.sample_type is not SampleType.OBSERVED:
                        raise ValueError(
                            "Derived samples in groups must be supplied as "
                            "exact parent-run inputs."
                        )
                    resolved, reference = self._resolve_observed_sample(
                        sample,
                        clause_index=clause_index,
                        source_group_id=snapshot.group_id,
                        group_fingerprint=snapshot.membership_fingerprint,
                        role=input_definition.role,
                        display_order=display_order,
                    )
                    sample_inputs.append(resolved)
                    references.append(reference)
            elif input_definition.kind is DerivedInputKind.SAMPLE:
                sample = available_samples.get(input_definition.input_id)
                if sample is None:
                    raise ValueError(
                        f"Input sample {input_definition.input_id!r} is "
                        "not available in the dataset."
                    )
                if sample.sample_type is not SampleType.OBSERVED:
                    raise ValueError(
                        "Derived samples must be supplied as exact parent-run inputs."
                    )
                resolved, reference = self._resolve_observed_sample(
                    sample,
                    clause_index=clause_index,
                    source_group_id=None,
                    group_fingerprint="",
                    role=input_definition.role,
                    display_order=0,
                )
                sample_inputs.append(resolved)
                references.append(reference)
            else:
                run = self._required_run(input_definition.input_id)
                if run.status is not DerivedRunStatus.COMPLETED:
                    raise ValueError(f"Parent run {run.id!r} is not completed.")
                parent_sample = self.derived_repository.derived_sample_for_run(run.id)
                if parent_sample is None:
                    raise RuntimeError(f"Parent run {run.id!r} has no derived sample.")
                if parent_sample.catalog_sample_id not in available_samples:
                    raise CatalogAccessError(
                        f"Parent run {run.id!r} is not available in "
                        f"dataset {dataset_id!r}."
                    )
                parent_inputs.append(
                    ResolvedParentInput(
                        clause_index=clause_index,
                        parent_run_id=run.id,
                        source_group_id=None,
                        group_membership_fingerprint="",
                        input_role=input_definition.role,
                        display_order=0,
                    )
                )
                references.append((run.reference_id, run.normalization_version))
            clauses.append(clause)

        return (
            tuple(sample_inputs),
            tuple(parent_inputs),
            clauses,
            references,
            group_snapshots,
        )

    def _resolve_observed_sample(
        self,
        sample: CatalogSample,
        *,
        clause_index: int,
        source_group_id: str | None,
        group_fingerprint: str,
        role: str,
        display_order: int,
    ) -> tuple[ResolvedSampleInput, tuple[str, str]]:
        cohort = self.catalog.get_cohort(sample.cohort_id)
        if cohort is None:
            raise CatalogNotFoundError(f"Cohort not found: {sample.cohort_id}")
        return (
            ResolvedSampleInput(
                clause_index=clause_index,
                catalog_sample_id=sample.id,
                source_group_id=source_group_id,
                group_membership_fingerprint=group_fingerprint,
                sample_fingerprint=sample.source_fingerprint,
                source_database_fingerprint=(cohort.source_database_fingerprint),
                input_role=role,
                display_order=display_order,
            ),
            (cohort.reference_id, cohort.normalization_version),
        )

    def _load_observed_evidence(
        self,
        sample_inputs: tuple[ResolvedSampleInput, ...],
        available_samples: dict[str, CatalogSample],
        definition: DerivedDefinition,
    ) -> tuple[
        dict[str, dict[AlleleKey, list[dict]]],
        dict[tuple[int, str], set[AlleleKey]],
    ]:
        inputs_by_cohort: dict[str, list[ResolvedSampleInput]] = {}
        for item in sample_inputs:
            sample = available_samples[item.catalog_sample_id]
            inputs_by_cohort.setdefault(sample.cohort_id, []).append(item)

        evidence_by_sample: dict[str, dict[AlleleKey, list[dict]]] = {}
        qualifying_sets: dict[tuple[int, str], set[AlleleKey]] = {}
        for cohort_id, inputs in inputs_by_cohort.items():
            cohort = self.catalog.get_cohort(cohort_id)
            if cohort is None:
                raise CatalogNotFoundError(f"Cohort not found: {cohort_id}")
            database_path = self._source_database_path(cohort)
            if sha256_file(database_path) != cohort.source_database_fingerprint:
                raise ValueError(
                    f"Source database changed since registration: "
                    f"{cohort.source_database_identifier}"
                )
            local_to_catalog = {
                available_samples[item.catalog_sample_id].source_sample_id: (
                    item.catalog_sample_id
                )
                for item in inputs
            }
            with StudyRepository.open(database_path) as study:
                evidence_rows = study.allele_evidence(
                    list(local_to_catalog),
                    position=definition.filters.position,
                    alt=definition.filters.alt,
                    af_rules=definition.filters.af_rules,
                    metadata_filters=definition.filters.metadata_filters,
                )
            for row in evidence_rows:
                catalog_sample_id = local_to_catalog[row["sample_id"]]
                allele = AlleleKey(
                    row["position"],
                    row["ref"],
                    row["alt"],
                )
                evidence_by_sample.setdefault(
                    catalog_sample_id,
                    {},
                ).setdefault(allele, []).append(row)
                if row["qualifies"]:
                    for item in inputs:
                        if item.catalog_sample_id == catalog_sample_id:
                            qualifying_sets.setdefault(
                                (item.clause_index, catalog_sample_id),
                                set(),
                            ).add(allele)
            for item in inputs:
                qualifying_sets.setdefault(
                    (item.clause_index, item.catalog_sample_id),
                    set(),
                )
                evidence_by_sample.setdefault(item.catalog_sample_id, {})
        return evidence_by_sample, qualifying_sets

    def _load_parent_alleles(
        self,
        parent_inputs: tuple[ResolvedParentInput, ...],
        definition: DerivedDefinition,
    ) -> tuple[
        dict[str, dict[AlleleKey, object]],
        dict[tuple[int, str], set[AlleleKey]],
    ]:
        alleles_by_run = {}
        parent_sets = {}
        for item in parent_inputs:
            materialized = self.derived_repository.list_run_alleles(item.parent_run_id)
            selected = {
                entry.allele: entry
                for entry in materialized
                if _identity_filter_matches(
                    entry.allele,
                    definition,
                )
            }
            alleles_by_run[item.parent_run_id] = selected
            parent_sets[(item.clause_index, item.parent_run_id)] = set(selected)
        return alleles_by_run, parent_sets

    @staticmethod
    def _populate_clause_sets(
        clauses: list[_ResolvedClause],
        sample_inputs: tuple[ResolvedSampleInput, ...],
        parent_inputs: tuple[ResolvedParentInput, ...],
        observed_sets: dict[tuple[int, str], set[AlleleKey]],
        parent_sets: dict[tuple[int, str], set[AlleleKey]],
    ) -> list[_ResolvedClause]:
        for item in sample_inputs:
            clauses[item.clause_index].member_sets.append(
                observed_sets[(item.clause_index, item.catalog_sample_id)]
            )
        for item in parent_inputs:
            clauses[item.clause_index].member_sets.append(
                parent_sets[(item.clause_index, item.parent_run_id)]
            )
        for clause in clauses:
            if not clause.member_sets:
                raise ValueError(
                    f"Derived input {clause.definition.input_id!r} "
                    "resolved to no members."
                )
        return clauses

    @staticmethod
    def _evaluate_clauses(
        definition: DerivedDefinition,
        clauses: list[_ResolvedClause],
    ) -> set[AlleleKey]:
        if definition.comparison is not None:
            return _evaluate_comparison(definition, clauses)

        candidates: set[AlleleKey] = set()
        for clause in clauses:
            if clause.definition.requirement is not PresenceRequirement.NONE:
                for member_set in clause.member_sets:
                    candidates.update(member_set)
        return {
            allele
            for allele in candidates
            if all(_clause_matches(clause, allele) for clause in clauses)
        }

    def _materialize_allele(
        self,
        allele: AlleleKey,
        *,
        definition: DerivedDefinition,
        clauses: list[_ResolvedClause],
        sample_inputs: tuple[ResolvedSampleInput, ...],
        parent_inputs: tuple[ResolvedParentInput, ...],
        observed_evidence: dict[str, dict[AlleleKey, list[dict]]],
        parent_alleles: dict[str, dict[AlleleKey, object]],
        available_samples: dict[str, CatalogSample],
    ) -> PendingMaterializedAllele:
        evidence = []
        qualifying_input_count = 0
        for item in sample_inputs:
            sample = available_samples[item.catalog_sample_id]
            cohort = self.catalog.get_cohort(sample.cohort_id)
            rows = observed_evidence[item.catalog_sample_id].get(allele, [])
            if not rows:
                evidence.append(
                    PendingAlleleEvidence(
                        source_catalog_sample_id=item.catalog_sample_id,
                        parent_derived_allele_id=None,
                        original_catalog_sample_id=item.catalog_sample_id,
                        source_mutation_id=None,
                        source_alt_index=None,
                        input_role=item.input_role,
                        evaluation_status=EvidenceStatus.NOT_OBSERVED,
                        metadata={
                            "source_cohort_id": sample.cohort_id,
                            "source_database_identifier": (
                                cohort.source_database_identifier
                            ),
                            "source_database_fingerprint": (
                                item.source_database_fingerprint
                            ),
                        },
                        lineage={"kind": "observed"},
                    )
                )
                continue
            if any(row["qualifies"] for row in rows):
                qualifying_input_count += 1
            for row in rows:
                evidence.append(
                    PendingAlleleEvidence(
                        source_catalog_sample_id=item.catalog_sample_id,
                        parent_derived_allele_id=None,
                        original_catalog_sample_id=item.catalog_sample_id,
                        source_mutation_id=row["mutation_id"],
                        source_alt_index=row["alt_index"],
                        input_role=item.input_role,
                        evaluation_status=(
                            EvidenceStatus.QUALIFYING_PRESENT
                            if row["qualifies"]
                            else EvidenceStatus.FILTERED_OUT
                        ),
                        af=row["af"],
                        af_text=row["af_text"],
                        filter_text=row["filter"],
                        metadata={
                            **row["metadata"],
                            "vcf_ref": row["vcf_ref"],
                            "source_cohort_id": sample.cohort_id,
                            "source_database_identifier": (
                                cohort.source_database_identifier
                            ),
                            "source_database_fingerprint": (
                                item.source_database_fingerprint
                            ),
                        },
                        lineage={"kind": "observed"},
                    )
                )

        for item in parent_inputs:
            parent_allele = parent_alleles[item.parent_run_id].get(allele)
            if parent_allele is None:
                continue
            qualifying_input_count += 1
            ancestor_evidence = self.derived_repository.list_allele_evidence(
                parent_allele.id
            )
            qualifying_ancestors = [
                source
                for source in ancestor_evidence
                if source.evaluation_status is EvidenceStatus.QUALIFYING_PRESENT
            ]
            original_ids = sorted(
                {
                    source.original_catalog_sample_id or source.source_catalog_sample_id
                    for source in qualifying_ancestors
                    if (
                        source.original_catalog_sample_id
                        or source.source_catalog_sample_id
                    )
                }
            )
            if not original_ids:
                original_ids = [None]
            for original_id in original_ids:
                evidence.append(
                    PendingAlleleEvidence(
                        source_catalog_sample_id=None,
                        parent_derived_allele_id=parent_allele.id,
                        original_catalog_sample_id=original_id,
                        source_mutation_id=None,
                        source_alt_index=None,
                        input_role=item.input_role,
                        evaluation_status=(EvidenceStatus.QUALIFYING_PRESENT),
                        metadata={
                            "parent_run_id": item.parent_run_id,
                        },
                        lineage={
                            "kind": "parent_run",
                            "parent_run_id": item.parent_run_id,
                            "parent_derived_allele_id": parent_allele.id,
                            "original_catalog_sample_id": original_id,
                        },
                    )
                )

        return PendingMaterializedAllele(
            allele=allele,
            metadata={
                "definition_fingerprint": definition.fingerprint,
                "qualifying_input_count": qualifying_input_count,
                "clause_evaluations": [
                    {
                        "input_id": clause.definition.input_id,
                        "role": clause.definition.role,
                        "requirement": (clause.definition.requirement.value),
                        "present_count": sum(
                            allele in member_set for member_set in clause.member_sets
                        ),
                        "member_count": len(clause.member_sets),
                    }
                    for clause in clauses
                ],
            },
            evidence=tuple(evidence),
        )

    def _source_database_path(self, cohort: Cohort) -> Path:
        configured = self.source_database_paths.get(cohort.source_database_identifier)
        if configured is not None:
            path = configured
        else:
            path = Path(cohort.provenance.get("source_path") or "").resolve()
        if not path.is_file():
            raise FileNotFoundError(
                "Registered source database not found for "
                f"{cohort.source_database_identifier!r}: {path}"
            )
        return path

    def _required_run(self, run_id: str):
        run = self.derived_repository.get_run(run_id)
        if run is None:
            raise CatalogNotFoundError(f"Derived run not found: {run_id}")
        return run


def _compatible_reference(
    values: list[tuple[str, str]],
) -> tuple[str, str]:
    normalized = {
        (str(reference or "").strip(), str(normalization or "").strip())
        for reference, normalization in values
    }
    if not normalized or any(not all(item) for item in normalized):
        raise ReferenceCompatibilityError(
            "Every derived input must declare a reference and normalization identity."
        )
    if len(normalized) != 1:
        labels = ", ".join(
            f"{reference}/{normalization}"
            for reference, normalization in sorted(normalized)
        )
        raise ReferenceCompatibilityError(
            f"Derived inputs use incompatible references: {labels}"
        )
    return next(iter(normalized))


def _identity_filter_matches(
    allele: AlleleKey,
    definition: DerivedDefinition,
) -> bool:
    filters = definition.filters
    if filters.position is not None and allele.position != filters.position:
        return False
    if filters.alt and allele.alt != filters.alt:
        return False
    return True


def _clause_matches(
    clause: _ResolvedClause,
    allele: AlleleKey,
) -> bool:
    present_count = sum(allele in member_set for member_set in clause.member_sets)
    requirement = clause.definition.requirement
    if requirement is PresenceRequirement.ANY:
        return present_count >= 1
    if requirement is PresenceRequirement.ALL:
        return present_count == len(clause.member_sets)
    if requirement is PresenceRequirement.NONE:
        return present_count == 0
    if requirement is PresenceRequirement.EXACTLY_ONE:
        return present_count == 1
    return False


def _evaluate_comparison(
    definition: DerivedDefinition,
    clauses: list[_ResolvedClause],
) -> set[AlleleKey]:
    """Evaluate a durable viewer comparison with live-table semantics."""
    comparison = definition.comparison
    if comparison is None:
        raise ValueError("Durable comparison definition is missing.")

    candidates = {
        allele
        for clause in clauses
        for member_set in clause.member_sets
        for allele in member_set
    }
    input_ids = tuple(str(index) for index in range(len(clauses)))
    sample_statuses = {
        input_id: set(comparison.input_statuses[index])
        for index, input_id in enumerate(input_ids)
    }
    results = set()
    for allele in candidates:
        present_input_ids = {
            str(index)
            for index, clause in enumerate(clauses)
            if any(allele in member_set for member_set in clause.member_sets)
        }
        present_count = len(present_input_ids)
        status = comparison_status(present_count, len(clauses))
        if status not in comparison.statuses:
            continue
        if not sample_filters_match(
            input_ids,
            sample_statuses,
            present_input_ids,
            present_count,
        ):
            continue
        results.add(allele)
    return results
