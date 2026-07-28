"""Validated filter values used by viewer requests and services."""

from __future__ import annotations

import math
from dataclasses import dataclass


AF_OPERATORS = {
    "gt": ">",
    "gte": ">=",
    "lt": "<",
    "lte": "<=",
    "eq": "=",
    "neq": "!=",
}
METADATA_FILTER_FIELDS = {
    "polymorphism",
    "reference_contains_alt",
    "reference_context",
    "reference_repeat",
    "repeat_base",
    "repeat_count",
    "repeat_2_bases",
    "repeat_2_count",
    "repeat_3_bases",
    "repeat_3_count",
}
REFERENCE_CONTAINS_ALT_VALUES = {"contains", "not_contains"}
REFERENCE_REPEAT_VALUES = {"before", "after", "one", "both", "none", "either"}
COUNT_FILTER_FIELDS = {
    "repeat_count",
    "repeat_2_count",
    "repeat_3_count",
}


@dataclass(frozen=True)
class AFRule:
    """One validated allele-fraction comparison."""

    operator: str
    threshold: float

    def __post_init__(self) -> None:
        operator = str(self.operator).strip().lower()
        if operator not in AF_OPERATORS:
            choices = ", ".join(AF_OPERATORS)
            raise ValueError(f"Unknown AF operator {operator!r}; expected {choices}.")
        try:
            threshold = float(self.threshold)
        except (TypeError, ValueError) as exc:
            raise ValueError("AF threshold must be numeric.") from exc
        if not math.isfinite(threshold):
            raise ValueError("AF threshold must be finite.")
        object.__setattr__(self, "operator", operator)
        object.__setattr__(self, "threshold", threshold)

    @classmethod
    def parse(cls, value: str) -> "AFRule":
        """Parse the existing query representation, such as ``gt:0.8``."""
        operator, separator, threshold = str(value).partition(":")
        if not separator:
            raise ValueError("AF rule must use the form operator:threshold.")
        return cls(operator=operator, threshold=threshold)

    def __iter__(self):
        """Allow compatibility with existing ``operator, threshold`` loops."""
        yield self.operator
        yield self.threshold


@dataclass(frozen=True)
class MetadataFilter:
    """One validated mutation-metadata filter."""

    field: str
    value: str

    def __post_init__(self) -> None:
        field = str(self.field).strip()
        value = str(self.value).strip()
        if field not in METADATA_FILTER_FIELDS:
            raise ValueError(f"Unknown metadata filter field: {field!r}.")
        if not value:
            raise ValueError("Metadata filter value cannot be empty.")
        if field == "polymorphism" and value not in {"0", "1"}:
            raise ValueError("Polymorphism filter must be 0 or 1.")
        if (
            field == "reference_contains_alt"
            and value not in REFERENCE_CONTAINS_ALT_VALUES
        ):
            raise ValueError(
                "Reference ALT filter must be contains or not_contains."
            )
        if field == "reference_repeat" and value not in REFERENCE_REPEAT_VALUES:
            raise ValueError("Unknown reference-repeat filter value.")
        if field in COUNT_FILTER_FIELDS:
            operator, separator, threshold = value.partition("|")
            if not separator or operator not in AF_OPERATORS:
                raise ValueError(
                    "Repeat-count filter must use operator|integer."
                )
            try:
                int(threshold)
            except ValueError as exc:
                raise ValueError(
                    "Repeat-count threshold must be an integer."
                ) from exc

        object.__setattr__(self, "field", field)
        object.__setattr__(self, "value", value)

    @classmethod
    def parse(cls, value: str) -> "MetadataFilter":
        """Parse the existing query representation, such as ``polymorphism:1``."""
        field, separator, raw_value = str(value).partition(":")
        if not separator:
            raise ValueError("Metadata filter must use the form field:value.")
        return cls(field=field, value=raw_value)

    def __iter__(self):
        """Allow compatibility with existing ``field, value`` loops."""
        yield self.field
        yield self.value
