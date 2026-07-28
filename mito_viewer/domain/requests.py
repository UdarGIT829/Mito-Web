"""Validated request objects independent of HTTP and SQLite."""

from __future__ import annotations

from dataclasses import dataclass, field

from .filters import AFRule, MetadataFilter


MITOCHONDRIAL_LENGTH = 16569
COMPARE_STATUSES = frozenset({"common", "partial", "unique"})
SAMPLE_STATUSES = frozenset({"present", "unique", "not_in"})


@dataclass(frozen=True)
class MutationFilters:
    """Validated mutation filters shared by list and comparison requests."""

    position: int | None = None
    alt: str = ""
    af_rules: tuple[AFRule, ...] = ()
    metadata_filters: tuple[MetadataFilter, ...] = ()

    def __post_init__(self) -> None:
        position = self.position
        if position not in (None, ""):
            try:
                position = int(position)
            except (TypeError, ValueError) as exc:
                raise ValueError("Position must be an integer.") from exc
            if not 1 <= position <= MITOCHONDRIAL_LENGTH:
                raise ValueError(
                    f"Position must be between 1 and {MITOCHONDRIAL_LENGTH}."
                )
        else:
            position = None

        alt = str(self.alt or "").strip().upper()
        af_rules = tuple(
            rule if isinstance(rule, AFRule) else AFRule(*rule)
            for rule in self.af_rules
        )
        metadata_filters = tuple(
            item
            if isinstance(item, MetadataFilter)
            else MetadataFilter(*item)
            for item in self.metadata_filters
        )

        object.__setattr__(self, "position", position)
        object.__setattr__(self, "alt", alt)
        object.__setattr__(self, "af_rules", af_rules)
        object.__setattr__(self, "metadata_filters", metadata_filters)


@dataclass(frozen=True)
class SampleConstraint:
    """Allowed presence states for one selected comparison sample."""

    sample_id: str
    statuses: frozenset[str] = field(
        default_factory=lambda: frozenset({"present"})
    )

    def __post_init__(self) -> None:
        sample_id = str(self.sample_id).strip()
        statuses = frozenset(str(status).strip() for status in self.statuses)
        unknown = statuses - SAMPLE_STATUSES
        if not sample_id:
            raise ValueError("Sample constraint requires a sample ID.")
        if unknown:
            raise ValueError(
                f"Unknown sample comparison status: {', '.join(sorted(unknown))}."
            )
        object.__setattr__(self, "sample_id", sample_id)
        object.__setattr__(self, "statuses", statuses)


@dataclass(frozen=True)
class ComparisonRequest:
    """Validated, transport-independent comparison request."""

    sample_ids: tuple[str, ...]
    filters: MutationFilters = field(default_factory=MutationFilters)
    statuses: frozenset[str] = field(default_factory=lambda: COMPARE_STATUSES)
    sample_constraints: tuple[SampleConstraint, ...] = ()
    limit: int = 2000

    def __post_init__(self) -> None:
        sample_ids = tuple(
            str(sample_id).strip()
            for sample_id in self.sample_ids
            if str(sample_id).strip()
        )
        if len(sample_ids) < 2:
            raise ValueError("Comparison requires at least two samples.")
        if len(set(sample_ids)) != len(sample_ids):
            raise ValueError("Comparison sample IDs must be unique.")

        filters = (
            self.filters
            if isinstance(self.filters, MutationFilters)
            else MutationFilters(**self.filters)
        )
        statuses = frozenset(str(status).strip() for status in self.statuses)
        unknown_statuses = statuses - COMPARE_STATUSES
        if unknown_statuses:
            raise ValueError(
                f"Unknown comparison status: {', '.join(sorted(unknown_statuses))}."
            )
        if not statuses:
            raise ValueError("Select at least one comparison status.")

        sample_constraints = tuple(
            constraint
            if isinstance(constraint, SampleConstraint)
            else SampleConstraint(*constraint)
            for constraint in self.sample_constraints
        )
        unknown_sample_ids = {
            constraint.sample_id
            for constraint in sample_constraints
        } - set(sample_ids)
        if unknown_sample_ids:
            raise ValueError(
                "Sample constraints reference unselected samples: "
                + ", ".join(sorted(unknown_sample_ids))
                + "."
            )
        if len({item.sample_id for item in sample_constraints}) != len(
            sample_constraints
        ):
            raise ValueError("Each comparison sample may have one constraint.")

        try:
            limit = int(self.limit)
        except (TypeError, ValueError) as exc:
            raise ValueError("Comparison limit must be an integer.") from exc
        if not 1 <= limit <= 10000:
            raise ValueError("Comparison limit must be between 1 and 10000.")

        object.__setattr__(self, "sample_ids", sample_ids)
        object.__setattr__(self, "filters", filters)
        object.__setattr__(self, "statuses", statuses)
        object.__setattr__(self, "sample_constraints", sample_constraints)
        object.__setattr__(self, "limit", limit)

    @property
    def sample_statuses(self) -> dict[str, set[str]]:
        """Return the mapping expected by the current comparison function."""
        return {
            constraint.sample_id: set(constraint.statuses)
            for constraint in self.sample_constraints
        }
