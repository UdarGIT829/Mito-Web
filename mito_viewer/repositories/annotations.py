"""Read-only access to the persistent annotation cache."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

from .schema import (
    DatabaseSchemaReport,
    inspect_database_schema,
    read_only_connection,
)


ANNOTATION_REQUIRED_SCHEMA = {
    "annotation_variants": frozenset(
        {
            "id",
            "pos",
            "ref",
            "alt",
            "first_seen_at",
            "last_seen_at",
        }
    ),
    "provider_annotations": frozenset(
        {
            "variant_id",
            "provider",
            "annotation_json",
            "error",
            "retrieved_at",
            "last_attempted_at",
        }
    ),
}


class AnnotationRepository:
    """Query one annotation-cache SQLite database."""

    def __init__(self, path: str | Path, *, validate: bool = True) -> None:
        self.path = Path(path).resolve()
        if not self.path.is_file():
            raise FileNotFoundError(
                f"Annotation database not found: {self.path}"
            )
        if validate:
            self.schema_report().require_valid()

    def schema_report(self) -> DatabaseSchemaReport:
        return inspect_database_schema(
            self.path,
            database_kind="annotation cache",
            required_schema=ANNOTATION_REQUIRED_SCHEMA,
        )

    def fetch_cached(self, position: int, ref: str, alt: str) -> dict | None:
        """Return cached provider envelopes for one mitochondrial allele."""
        position = int(position)
        ref = str(ref).strip().upper()
        alt = str(alt).strip().upper()
        if not ref or not alt:
            raise ValueError(
                "Annotation lookup requires position, ref, and alt."
            )

        connection = read_only_connection(self.path)
        try:
            variant = connection.execute(
                """
                SELECT id, pos, ref, alt, first_seen_at, last_seen_at
                FROM annotation_variants
                WHERE pos = ? AND ref = ? AND alt = ?
                """,
                (position, ref, alt),
            ).fetchone()
            if variant is None:
                return None

            payload: dict[str, Any] = {
                "variant": dict(variant),
                "cache": {"database": self.path.name},
            }
            rows = connection.execute(
                """
                SELECT
                    provider,
                    annotation_json,
                    error,
                    retrieved_at,
                    last_attempted_at
                FROM provider_annotations
                WHERE variant_id = ?
                ORDER BY provider
                """,
                (variant["id"],),
            )
            for row in rows:
                if row["annotation_json"] is not None:
                    try:
                        provider_payload = json.loads(row["annotation_json"])
                    except json.JSONDecodeError as exc:
                        provider_payload = {
                            "source": row["provider"],
                            "error": (
                                "Cached annotation JSON is invalid: "
                                f"{exc}"
                            ),
                        }
                else:
                    provider_payload = {
                        "source": row["provider"],
                        "error": (
                            row["error"] or "No cached annotation response"
                        ),
                    }
                payload[row["provider"]] = provider_payload
            return payload
        finally:
            connection.close()

    def vocabulary(self) -> dict[str, list[dict[str, Any]]]:
        """Return observed ClinVar classifications and conditions with counts."""
        classifications = Counter()
        conditions = Counter()
        connection = read_only_connection(self.path)
        try:
            rows = connection.execute(
                """
                SELECT annotation_json
                FROM provider_annotations
                WHERE provider = 'clinvar' AND annotation_json IS NOT NULL
                """
            )
            for (annotation_json,) in rows:
                try:
                    envelope = json.loads(annotation_json)
                except json.JSONDecodeError:
                    continue
                result = (
                    envelope.get("data", {})
                    .get("summaries", {})
                    .get("result", {})
                )
                for uid in result.get("uids", []):
                    clinical = (
                        result.get(uid, {})
                        .get("germline_classification", {})
                    )
                    classification = clinical.get("description")
                    if classification:
                        classifications[str(classification)] += 1
                    for trait in clinical.get("trait_set", []):
                        condition = trait.get("trait_name")
                        if condition:
                            conditions[str(condition)] += 1
        finally:
            connection.close()

        def ranked(counter: Counter) -> list[dict[str, Any]]:
            return [
                {"value": value, "count": count}
                for value, count in sorted(
                    counter.items(),
                    key=lambda item: (-item[1], item[0].lower()),
                )
            ]

        return {
            "classifications": ranked(classifications),
            "conditions": ranked(conditions),
        }
