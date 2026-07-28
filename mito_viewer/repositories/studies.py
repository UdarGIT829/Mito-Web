"""Read-only access to imported mitochondrial study databases."""

from __future__ import annotations

import json
import sqlite3
from collections.abc import Iterator
from pathlib import Path

from mito_viewer.domain import AlleleKey, SampleAlleleCall
from mito_viewer.domain.filters import AF_OPERATORS

from .schema import (
    DatabaseSchemaReport,
    inspect_database_schema,
    read_only_connection,
)


DATABASE_EXTENSIONS = {".db", ".sqlite", ".sqlite3"}
NO_TAGS_FILTER = "__NO_TAGS__"
REFERENCE_REPEAT_BASES = ("A", "C", "G", "T", "N")
STUDY_TABLES = (
    "subjects",
    "samples",
    "sample_population_tags",
    "mutations",
    "mutation_alts",
)
STUDY_REQUIRED_SCHEMA = {
    "subjects": frozenset({"id", "subject_id"}),
    "samples": frozenset(
        {"id", "subject_id", "population_key", "source_file"}
    ),
    "sample_population_tags": frozenset({"sample_id", "tag", "tag_order"}),
    "mutations": frozenset(
        {
            "id",
            "sample_id",
            "pos",
            "ref",
            "vcf_ref",
            "alt",
            "af",
            "filter",
            "metadata_json",
        }
    ),
    "mutation_alts": frozenset(
        {"mutation_id", "alt_index", "alt", "af", "af_text"}
    ),
}
ALT_FILTER_COLUMNS = {
    "polymorphism": ("mutation_alts.polymorphism", "integer"),
    "repeat_base": ("mutation_alts.repeat_base", "text"),
    "repeat_count": ("mutation_alts.repeat_count", "integer"),
    "repeat_2_bases": ("mutation_alts.repeat_2_bases", "text"),
    "repeat_2_count": ("mutation_alts.repeat_2_count", "integer"),
    "repeat_3_bases": ("mutation_alts.repeat_3_bases", "text"),
    "repeat_3_count": ("mutation_alts.repeat_3_count", "integer"),
}
ALT_METADATA_FIELDS = (
    ("polymorphism", "alt_polymorphism", "POLYMORPHISM"),
    ("repeat_base", "alt_repeat_base", "REPEAT_1_BASE"),
    ("repeat_count", "alt_repeat_count", "REPEAT_1_BASE_COUNT"),
    ("repeat_2_bases", "alt_repeat_2_bases", "REPEAT_2_BASES"),
    ("repeat_2_count", "alt_repeat_2_count", "REPEAT_2_BASES_COUNT"),
    ("repeat_3_bases", "alt_repeat_3_bases", "REPEAT_3_BASES"),
    ("repeat_3_count", "alt_repeat_3_count", "REPEAT_3_BASES_COUNT"),
)


def is_sqlite_database_path(path: str | Path) -> bool:
    database_path = Path(path)
    return (
        database_path.is_file()
        and database_path.suffix.lower() in DATABASE_EXTENSIONS
    )


def inspect_study_database(path: str | Path) -> DatabaseSchemaReport:
    return inspect_database_schema(
        path,
        database_kind="study",
        required_schema=STUDY_REQUIRED_SCHEMA,
    )


def discover_study_databases(
    database_dir: str | Path,
) -> dict[str, Path]:
    """Return only schema-compatible study databases."""
    directory = Path(database_dir)
    databases = {}
    if not directory.exists():
        return databases

    for path in sorted(directory.iterdir(), key=lambda item: item.name.lower()):
        if not is_sqlite_database_path(path):
            continue
        if inspect_study_database(path).valid:
            databases[path.name] = path.resolve()
    return databases


class StudyRepository:
    """Query one schema-validated study database."""

    def __init__(
        self,
        connection: sqlite3.Connection,
        *,
        path: str | Path | None = None,
        owns_connection: bool = False,
    ) -> None:
        self.connection = connection
        self.path = Path(path).resolve() if path is not None else None
        self.owns_connection = owns_connection

    @classmethod
    def open(cls, path: str | Path, *, validate: bool = True) -> "StudyRepository":
        database_path = Path(path).resolve()
        if not database_path.is_file():
            raise FileNotFoundError(f"Study database not found: {database_path}")
        if validate:
            inspect_study_database(database_path).require_valid()
        return cls(
            read_only_connection(database_path),
            path=database_path,
            owns_connection=True,
        )

    def __enter__(self) -> "StudyRepository":
        return self

    def __exit__(self, exc_type, exc, traceback) -> None:
        self.close()

    def close(self) -> None:
        if self.owns_connection:
            self.connection.close()
            self.owns_connection = False

    def schema_report(self) -> DatabaseSchemaReport:
        if self.path is None:
            raise ValueError("Schema reports require a repository path.")
        return inspect_study_database(self.path)

    def table_columns(self, table: str) -> set[str]:
        if table not in STUDY_TABLES:
            raise ValueError(f"Unknown study table: {table}")
        return {
            row["name"]
            for row in self.connection.execute(
                f"PRAGMA table_info({table})"
            )
        }

    def counts(self) -> dict[str, int]:
        return {
            table: self.connection.execute(
                f"SELECT COUNT(*) FROM {table}"
            ).fetchone()[0]
            for table in STUDY_TABLES
        }

    def subjects(self) -> list[dict]:
        rows = self.connection.execute(
            """
            SELECT
                subjects.id,
                subjects.subject_id,
                COUNT(samples.id) AS sample_count
            FROM subjects
            LEFT JOIN samples ON samples.subject_id = subjects.id
            GROUP BY subjects.id
            ORDER BY subjects.subject_id
            """
        ).fetchall()
        return [dict(row) for row in rows]

    def population_tags(self) -> list[dict]:
        rows = self.connection.execute(
            """
            SELECT tag, COUNT(DISTINCT sample_id) AS sample_count
            FROM sample_population_tags
            GROUP BY tag
            ORDER BY tag
            """
        ).fetchall()
        tags = [dict(row) for row in rows]
        if tags:
            untagged_count = self.connection.execute(
                """
                SELECT COUNT(*)
                FROM samples
                WHERE NOT EXISTS (
                    SELECT 1
                    FROM sample_population_tags
                    WHERE sample_population_tags.sample_id = samples.id
                )
                """
            ).fetchone()[0]
            if untagged_count:
                tags.append(
                    {
                        "tag": NO_TAGS_FILTER,
                        "label": "<NONE>",
                        "sample_count": untagged_count,
                    }
                )
        return tags

    def samples(
        self,
        *,
        subject_id: str | None = None,
        tags: list[str] | None = None,
    ) -> list[dict]:
        join_params = []
        where_params = []
        where_clauses = []
        tag_joins = []

        if subject_id:
            where_clauses.append("subjects.subject_id = ?")
            where_params.append(subject_id)

        for index, tag in enumerate(tags or []):
            if tag == NO_TAGS_FILTER:
                where_clauses.append(
                    """
                    NOT EXISTS (
                        SELECT 1
                        FROM sample_population_tags no_tag
                        WHERE no_tag.sample_id = samples.id
                    )
                    """
                )
                continue
            alias = f"tag_{index}"
            tag_joins.append(
                f"""
                JOIN sample_population_tags {alias}
                    ON {alias}.sample_id = samples.id
                    AND {alias}.tag = ?
                """
            )
            join_params.append(tag)

        where_sql = ""
        if where_clauses:
            where_sql = "WHERE " + " AND ".join(where_clauses)

        rows = self.connection.execute(
            f"""
            SELECT
                samples.id,
                subjects.subject_id,
                samples.population_key,
                samples.source_file,
                COUNT(mutations.id) AS mutation_count
            FROM samples
            JOIN subjects ON subjects.id = samples.subject_id
            {" ".join(tag_joins)}
            LEFT JOIN mutations ON mutations.sample_id = samples.id
            {where_sql}
            GROUP BY samples.id
            ORDER BY subjects.subject_id, samples.population_key
            """,
            join_params + where_params,
        ).fetchall()
        return [dict(row) for row in rows]

    def registration_samples(self) -> list[dict]:
        """Return stable source identities and metadata for catalog indexing."""
        rows = self.connection.execute(
            """
            SELECT
                samples.id,
                subjects.subject_id,
                samples.population_key,
                samples.source_file
            FROM samples
            JOIN subjects ON subjects.id = samples.subject_id
            ORDER BY samples.id
            """
        ).fetchall()
        tags_by_sample: dict[str, list[str]] = {}
        for row in self.connection.execute(
            """
            SELECT sample_id, tag
            FROM sample_population_tags
            ORDER BY sample_id, tag_order, tag
            """
        ):
            tags_by_sample.setdefault(str(row["sample_id"]), []).append(
                row["tag"]
            )

        return [
            {
                **dict(row),
                "id": str(row["id"]),
                "population_tags": tags_by_sample.get(str(row["id"]), []),
            }
            for row in rows
        ]

    def iter_sample_snapshot_records(
        self,
        sample_id: str | int,
    ) -> Iterator[dict]:
        """Yield every stored source record contributing to one sample."""
        sample_row = self.connection.execute(
            """
            SELECT
                samples.*,
                subjects.subject_id AS subject_identifier
            FROM samples
            JOIN subjects ON subjects.id = samples.subject_id
            WHERE samples.id = ?
            """,
            (sample_id,),
        ).fetchone()
        if sample_row is None:
            raise KeyError(f"Source sample not found: {sample_id}")
        yield {"record_type": "sample", **dict(sample_row)}

        for row in self.connection.execute(
            """
            SELECT *
            FROM sample_population_tags
            WHERE sample_id = ?
            ORDER BY tag_order, tag
            """,
            (sample_id,),
        ):
            yield {"record_type": "population_tag", **dict(row)}

        for row in self.connection.execute(
            """
            SELECT *
            FROM mutations
            WHERE sample_id = ?
            ORDER BY id
            """,
            (sample_id,),
        ):
            yield {"record_type": "mutation", **dict(row)}

        for row in self.connection.execute(
            """
            SELECT mutation_alts.*
            FROM mutation_alts
            JOIN mutations ON mutations.id = mutation_alts.mutation_id
            WHERE mutations.sample_id = ?
            ORDER BY mutation_alts.mutation_id, mutation_alts.alt_index
            """,
            (sample_id,),
        ):
            yield {"record_type": "mutation_alt", **dict(row)}

    def mutation_samples(
        self,
        position: int,
        ref: str,
        alt: str,
    ) -> list[dict]:
        position = int(position)
        ref = str(ref).strip().upper()
        alt = str(alt).strip().upper()
        if not ref or not alt:
            raise ValueError("Sample lookup requires position, ref, and alt.")

        rows = self.connection.execute(
            """
            SELECT DISTINCT
                samples.id,
                subjects.subject_id,
                samples.population_key,
                samples.source_file
            FROM mutation_alts
            JOIN mutations ON mutations.id = mutation_alts.mutation_id
            JOIN samples ON samples.id = mutations.sample_id
            JOIN subjects ON subjects.id = samples.subject_id
            WHERE mutations.pos = ?
              AND UPPER(mutations.vcf_ref) = ?
              AND UPPER(mutation_alts.alt) = ?
            ORDER BY subjects.subject_id, samples.population_key
            """,
            (position, ref, alt),
        ).fetchall()
        return [
            {
                **dict(row),
                "label": " ".join(
                    filter(
                        None,
                        (
                            row["subject_id"],
                            row["population_key"].replace("|", "_"),
                        ),
                    )
                ),
            }
            for row in rows
        ]

    def mutation_rows(
        self,
        *,
        sample_id: str | None = None,
        position: int | str | None = None,
        alt: str | None = None,
        af_rules=(),
        metadata_filters=(),
        limit: int = 500,
    ) -> list[dict]:
        params = []
        where_clauses = []
        alt_join = ""
        use_alt_rows = False

        if sample_id:
            where_clauses.append("samples.id = ?")
            params.append(sample_id)
        if position:
            where_clauses.append("mutations.pos = ?")
            params.append(position)
        if alt:
            alt_join = (
                "JOIN mutation_alts "
                "ON mutation_alts.mutation_id = mutations.id"
            )
            use_alt_rows = True
            where_clauses.append("mutation_alts.alt = ?")
            params.append(alt)
        if af_rules:
            alt_join = (
                "JOIN mutation_alts "
                "ON mutation_alts.mutation_id = mutations.id"
            )
            use_alt_rows = True
            for operator, threshold in af_rules:
                where_clauses.append(
                    f"mutation_alts.af {AF_OPERATORS[operator]} ?"
                )
                params.append(threshold)
        if metadata_filters and self._add_metadata_filter_sql(
            where_clauses,
            params,
            metadata_filters,
        ):
            alt_join = (
                "JOIN mutation_alts "
                "ON mutation_alts.mutation_id = mutations.id"
            )
            use_alt_rows = True

        where_sql = ""
        if where_clauses:
            where_sql = "WHERE " + " AND ".join(where_clauses)
        alt_select = "mutation_alts.alt" if use_alt_rows else "mutations.alt"
        af_select = (
            "mutation_alts.af_text" if use_alt_rows else "mutations.af"
        )

        params.append(limit)
        rows = self.connection.execute(
            f"""
            SELECT DISTINCT
                mutations.id,
                subjects.subject_id,
                samples.population_key,
                samples.source_file,
                mutations.pos,
                mutations.ref,
                mutations.vcf_ref,
                {alt_select} AS alt,
                {af_select} AS af,
                mutations.filter,
                mutations.metadata_json
            FROM mutations
            JOIN samples ON samples.id = mutations.sample_id
            JOIN subjects ON subjects.id = samples.subject_id
            {alt_join}
            {where_sql}
            ORDER BY
                mutations.pos,
                subjects.subject_id,
                samples.population_key
            LIMIT ?
            """,
            params,
        ).fetchall()
        return [dict(row) for row in rows]

    def allele_calls(
        self,
        sample_ids: list[str],
        *,
        position: int | str | None = None,
        alt: str | None = None,
        af_rules=(),
        metadata_filters=(),
    ) -> list[SampleAlleleCall]:
        if not sample_ids:
            return []

        placeholders = ",".join("?" for _ in sample_ids)
        params = list(sample_ids)
        where_clauses = [f"samples.id IN ({placeholders})"]
        if position:
            where_clauses.append("mutations.pos = ?")
            params.append(position)
        if alt:
            where_clauses.append("mutation_alts.alt = ?")
            params.append(alt)
        if af_rules:
            for operator, threshold in af_rules:
                where_clauses.append(
                    f"mutation_alts.af {AF_OPERATORS[operator]} ?"
                )
                params.append(threshold)
        if metadata_filters:
            self._add_metadata_filter_sql(
                where_clauses,
                params,
                metadata_filters,
            )

        alt_columns = self.table_columns("mutation_alts")
        alt_metadata_selects = [
            f"mutation_alts.{column} AS {alias}"
            for column, alias, _key in ALT_METADATA_FIELDS
            if column in alt_columns
        ]
        metadata_select_sql = ""
        if alt_metadata_selects:
            metadata_select_sql = ",\n                " + ",\n                ".join(
                alt_metadata_selects
            )

        rows = self.connection.execute(
            f"""
            SELECT
                samples.id AS sample_id,
                subjects.subject_id,
                samples.population_key,
                mutations.pos,
                mutations.ref,
                mutations.vcf_ref,
                mutations.filter,
                mutation_alts.alt,
                mutation_alts.af,
                mutation_alts.af_text,
                mutations.af AS mutation_af,
                mutations.metadata_json
                {metadata_select_sql}
            FROM mutation_alts
            JOIN mutations ON mutations.id = mutation_alts.mutation_id
            JOIN samples ON samples.id = mutations.sample_id
            JOIN subjects ON subjects.id = samples.subject_id
            WHERE {" AND ".join(where_clauses)}
            ORDER BY
                mutations.pos,
                mutation_alts.alt,
                subjects.subject_id,
                samples.population_key
            """,
            params,
        ).fetchall()

        calls = []
        for row in rows:
            try:
                source_metadata = json.loads(row["metadata_json"] or "{}")
            except json.JSONDecodeError:
                source_metadata = {}
            call_metadata = {
                key: source_metadata.get(key, "")
                for key in ("REFERENCE_6_BEFORE", "REFERENCE_6_AFTER")
                if key in source_metadata
            }
            for _column, alias, key in ALT_METADATA_FIELDS:
                if alias in row.keys() and row[alias] not in (None, ""):
                    call_metadata[key] = str(row[alias])
            af_text = (
                row["af_text"]
                or source_metadata.get("AF", "")
                or source_metadata.get("VF", "")
                or row["mutation_af"]
                or ""
            )
            calls.append(
                SampleAlleleCall(
                    allele=AlleleKey(
                        position=row["pos"],
                        ref=row["ref"],
                        alt=row["alt"],
                    ),
                    sample_id=row["sample_id"],
                    label=(
                        f"{row['subject_id']} "
                        f"{row['population_key'].replace('|', '_')}"
                    ),
                    af=row["af"],
                    af_text=af_text,
                    filter=row["filter"],
                    vcf_ref=row["vcf_ref"],
                    metadata=call_metadata,
                )
            )
        return calls

    def allele_evidence(
        self,
        sample_ids: list[str],
        *,
        position: int | str | None = None,
        alt: str | None = None,
        af_rules=(),
        metadata_filters=(),
    ) -> list[dict]:
        """Return all source calls and whether each passes selected filters."""
        if not sample_ids:
            return []

        placeholders = ",".join("?" for _ in sample_ids)
        filter_params = list(sample_ids)
        filter_clauses = [f"samples.id IN ({placeholders})"]
        if position:
            filter_clauses.append("mutations.pos = ?")
            filter_params.append(position)
        if alt:
            filter_clauses.append("mutation_alts.alt = ?")
            filter_params.append(alt)
        for operator, threshold in af_rules or ():
            filter_clauses.append(
                f"mutation_alts.af {AF_OPERATORS[operator]} ?"
            )
            filter_params.append(threshold)
        self._add_metadata_filter_sql(
            filter_clauses,
            filter_params,
            metadata_filters,
        )
        qualifying_keys = {
            (
                str(row["sample_id"]),
                str(row["mutation_id"]),
                int(row["alt_index"]),
            )
            for row in self.connection.execute(
                f"""
                SELECT
                    samples.id AS sample_id,
                    mutations.id AS mutation_id,
                    mutation_alts.alt_index
                FROM mutation_alts
                JOIN mutations
                    ON mutations.id = mutation_alts.mutation_id
                JOIN samples ON samples.id = mutations.sample_id
                WHERE {" AND ".join(filter_clauses)}
                """,
                filter_params,
            )
        }

        alt_columns = self.table_columns("mutation_alts")
        alt_metadata_selects = [
            f"mutation_alts.{column} AS {alias}"
            for column, alias, _key in ALT_METADATA_FIELDS
            if column in alt_columns
        ]
        metadata_select_sql = ""
        if alt_metadata_selects:
            metadata_select_sql = ",\n                " + ",\n                ".join(
                alt_metadata_selects
            )
        rows = self.connection.execute(
            f"""
            SELECT
                samples.id AS sample_id,
                mutations.id AS mutation_id,
                mutation_alts.alt_index,
                mutations.pos,
                mutations.ref,
                mutations.vcf_ref,
                mutation_alts.alt,
                mutation_alts.af,
                mutation_alts.af_text,
                mutations.af AS mutation_af,
                mutations.filter,
                mutations.metadata_json
                {metadata_select_sql}
            FROM mutation_alts
            JOIN mutations ON mutations.id = mutation_alts.mutation_id
            JOIN samples ON samples.id = mutations.sample_id
            WHERE samples.id IN ({placeholders})
            ORDER BY
                samples.id,
                mutations.pos,
                mutations.id,
                mutation_alts.alt_index
            """,
            sample_ids,
        ).fetchall()

        evidence = []
        for row in rows:
            raw_metadata = row["metadata_json"] or "{}"
            try:
                metadata = json.loads(raw_metadata)
            except json.JSONDecodeError:
                metadata = {"_invalid_source_json": raw_metadata}
            for _column, alias, key in ALT_METADATA_FIELDS:
                if alias in row.keys() and row[alias] not in (None, ""):
                    metadata[key] = str(row[alias])
            af_text = (
                row["af_text"]
                or metadata.get("AF", "")
                or metadata.get("VF", "")
                or row["mutation_af"]
                or ""
            )
            key = (
                str(row["sample_id"]),
                str(row["mutation_id"]),
                int(row["alt_index"]),
            )
            evidence.append(
                {
                    "sample_id": key[0],
                    "mutation_id": key[1],
                    "alt_index": key[2],
                    "position": row["pos"],
                    "ref": row["ref"],
                    "vcf_ref": row["vcf_ref"],
                    "alt": row["alt"],
                    "af": row["af"],
                    "af_text": af_text,
                    "filter": row["filter"],
                    "metadata": metadata,
                    "qualifies": key in qualifying_keys,
                }
            )
        return evidence

    def sample_labels(self, sample_ids: list[str]) -> dict[str, str]:
        if not sample_ids:
            return {}
        placeholders = ",".join("?" for _ in sample_ids)
        rows = self.connection.execute(
            f"""
            SELECT samples.id, subjects.subject_id, samples.population_key
            FROM samples
            JOIN subjects ON subjects.id = samples.subject_id
            WHERE samples.id IN ({placeholders})
            """,
            sample_ids,
        ).fetchall()
        return {
            str(row["id"]): (
                f"{row['subject_id']} "
                f"{row['population_key'].replace('|', '_')}"
            )
            for row in rows
        }

    @staticmethod
    def _reference_repeat_sql_expression(json_key: str):
        checks = [
            f"json_extract(mutations.metadata_json, '$.{json_key}') LIKE ?"
            for _base in REFERENCE_REPEAT_BASES
        ]
        params = [f"%{base * 2}%" for base in REFERENCE_REPEAT_BASES]
        return "(" + " OR ".join(checks) + ")", params

    @classmethod
    def _reference_repeat_filter_sql(cls, raw_value: str):
        before_sql, before_params = cls._reference_repeat_sql_expression(
            "REFERENCE_6_BEFORE"
        )
        after_sql, after_params = cls._reference_repeat_sql_expression(
            "REFERENCE_6_AFTER"
        )
        if raw_value == "before":
            return (
                f"({before_sql} AND NOT {after_sql})",
                before_params + after_params,
            )
        if raw_value == "after":
            return (
                f"({after_sql} AND NOT {before_sql})",
                after_params + before_params,
            )
        if raw_value == "one":
            return (
                f"(({before_sql} AND NOT {after_sql}) OR "
                f"({after_sql} AND NOT {before_sql}))",
                before_params
                + after_params
                + after_params
                + before_params,
            )
        if raw_value == "both":
            return (
                f"({before_sql} AND {after_sql})",
                before_params + after_params,
            )
        if raw_value == "none":
            return (
                f"(NOT {before_sql} AND NOT {after_sql})",
                before_params + after_params,
            )
        if raw_value == "either":
            return (
                f"({before_sql} OR {after_sql})",
                before_params + after_params,
            )
        return "", []

    def _add_metadata_filter_sql(
        self,
        where_clauses: list[str],
        params: list,
        metadata_filters,
    ) -> bool:
        if not metadata_filters:
            return False

        alt_columns = self.table_columns("mutation_alts")
        joined_alts = False
        for field, raw_value in metadata_filters:
            if field == "reference_contains_alt":
                before_sql = (
                    "UPPER(COALESCE(json_extract(mutations.metadata_json, "
                    "'$.REFERENCE_6_BEFORE'), ''))"
                )
                after_sql = (
                    "UPPER(COALESCE(json_extract(mutations.metadata_json, "
                    "'$.REFERENCE_6_AFTER'), ''))"
                )
                contains_sql = (
                    f"(INSTR({before_sql}, UPPER(mutation_alts.alt)) > 0 OR "
                    f"INSTR({after_sql}, UPPER(mutation_alts.alt)) > 0)"
                )
                if raw_value == "contains":
                    where_clauses.append(f"({contains_sql})")
                elif raw_value == "not_contains":
                    where_clauses.append(f"(NOT ({contains_sql}))")
                else:
                    continue
                joined_alts = True
                continue

            if field == "reference_context":
                where_clauses.append(
                    "(json_extract(mutations.metadata_json, "
                    "'$.REFERENCE_6_BEFORE') LIKE ? OR "
                    "json_extract(mutations.metadata_json, "
                    "'$.REFERENCE_6_AFTER') LIKE ?)"
                )
                params.extend([f"%{raw_value}%", f"%{raw_value}%"])
                continue

            if field == "reference_repeat":
                clause, clause_params = self._reference_repeat_filter_sql(
                    raw_value
                )
                if clause:
                    where_clauses.append(clause)
                    params.extend(clause_params)
                continue

            if field not in ALT_FILTER_COLUMNS:
                continue
            column, field_type = ALT_FILTER_COLUMNS[field]
            if column.split(".")[-1] not in alt_columns:
                continue
            joined_alts = True

            if field_type == "integer":
                if field == "polymorphism":
                    where_clauses.append(f"{column} = ?")
                    params.append(1 if raw_value == "1" else 0)
                else:
                    operator, separator, threshold = raw_value.partition("|")
                    if not separator or operator not in AF_OPERATORS:
                        continue
                    where_clauses.append(
                        f"{column} {AF_OPERATORS[operator]} ?"
                    )
                    params.append(int(threshold))
            else:
                where_clauses.append(f"{column} = ?")
                params.append(raw_value)

        return joined_alts
