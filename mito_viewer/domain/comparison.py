"""Pure comparison-set construction shared by legacy and Dataset queries."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import Any

from .models import SampleAlleleCall


DEFAULT_COMPARE_STATUSES = frozenset({"common", "partial", "unique"})
DEFAULT_SAMPLE_COMPARE_STATUSES = frozenset({"present"})


def comparison_status(present_count: int, total_count: int) -> str:
    """Classify one allele by how many selected samples contain it."""
    if present_count == total_count:
        return "common"
    if present_count == 1:
        return "unique"
    return "partial"


def sample_filters_match(
    compare_sample_ids: Sequence[str],
    sample_statuses: Mapping[str, set[str] | frozenset[str]],
    present_sample_ids: set[str],
    present_count: int,
) -> bool:
    """Return whether direct per-sample constraints match an allele."""
    present_required = set()
    absent_required = set()
    unique_allowed = set()

    for sample_id in compare_sample_ids:
        allowed_statuses = sample_statuses.get(
            sample_id,
            DEFAULT_SAMPLE_COMPARE_STATUSES,
        )
        if not allowed_statuses:
            continue
        if "not_in" in allowed_statuses:
            absent_required.add(sample_id)
        if "present" in allowed_statuses:
            present_required.add(sample_id)
        if "unique" in allowed_statuses:
            unique_allowed.add(sample_id)

    if absent_required & present_sample_ids:
        return False

    present_branch = (
        bool(present_required)
        and present_required.issubset(present_sample_ids)
    )
    unique_branch = (
        bool(unique_allowed)
        and present_count == 1
        and bool(unique_allowed & present_sample_ids)
    )

    if present_required or unique_allowed:
        return present_branch or unique_branch
    return True


def comparison_rows(
    compare_sample_ids: Sequence[str],
    calls: Iterable[SampleAlleleCall],
    *,
    sample_labels: Mapping[str, str],
    statuses: Iterable[str] | None = None,
    sample_statuses: Mapping[str, set[str] | frozenset[str]] | None = None,
    sample_provenance: Mapping[str, Mapping[str, Any]] | None = None,
    limit: int = 2000,
) -> list[dict[str, Any]]:
    """Build browser comparison rows without making storage assumptions."""
    compare_sample_ids = [str(sample_id) for sample_id in compare_sample_ids]
    if len(compare_sample_ids) < 2:
        return []

    global_statuses = set(statuses or DEFAULT_COMPARE_STATUSES)
    normalized_statuses = {
        str(sample_id): set(allowed)
        for sample_id, allowed in (sample_statuses or {}).items()
    }
    provenance = {
        str(sample_id): dict(values)
        for sample_id, values in (sample_provenance or {}).items()
    }
    calls_by_allele: dict[Any, list[SampleAlleleCall]] = {}
    for call in calls:
        calls_by_allele.setdefault(call.allele, []).append(call)

    results = []
    for allele, allele_calls in calls_by_allele.items():
        present_sample_ids = {
            str(call.sample_id)
            for call in allele_calls
        }
        present_count = len(present_sample_ids)
        status = comparison_status(present_count, len(compare_sample_ids))

        if status not in global_statuses:
            continue
        if not sample_filters_match(
            compare_sample_ids,
            normalized_statuses,
            present_sample_ids,
            present_count,
        ):
            continue

        present = []
        for call in allele_calls:
            item = call.to_json()
            item.update(provenance.get(str(call.sample_id), {}))
            present.append(item)
        missing = []
        for sample_id in compare_sample_ids:
            if sample_id in present_sample_ids:
                continue
            item = {
                "sample_id": sample_id,
                "label": sample_labels.get(
                    sample_id,
                    f"Sample {sample_id}",
                ),
            }
            item.update(provenance.get(sample_id, {}))
            missing.append(item)

        results.append(
            {
                "pos": allele.position,
                "ref": allele.ref,
                "alt": allele.alt,
                "group_key": (
                    f"{status}|{allele.position}|{allele.ref}|{allele.alt}"
                ),
                "present": present,
                "missing": missing,
                "status": status,
            }
        )

    results.sort(
        key=lambda row: (
            row["pos"],
            row["alt"],
            row["status"],
        )
    )
    results = results[:limit]

    group_sizes = {}
    for row in results:
        group_sizes[row["group_key"]] = (
            group_sizes.get(row["group_key"], 0) + 1
        )
    for index, row in enumerate(results):
        previous_key = results[index - 1]["group_key"] if index else None
        next_key = (
            results[index + 1]["group_key"]
            if index + 1 < len(results)
            else None
        )
        row["group_size"] = group_sizes[row["group_key"]]
        row["group_start"] = row["group_key"] != previous_key
        row["group_end"] = row["group_key"] != next_key

    return results
