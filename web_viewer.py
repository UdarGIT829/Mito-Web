#!/usr/bin/env python3
"""Small web viewer for the SQLite mutation database."""

import argparse
import html
import json
import sqlite3
from dataclasses import dataclass, field
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, unquote, urlparse

import vcf_parser


DEFAULT_DB_PATH = Path(__file__).resolve().parent / "master.sqlite"
VIEWER_TEMPLATE_PATH = Path(__file__).resolve().parent / "templates" / "viewer.html"
DEFAULT_COMPARE_STATUSES = {"common", "partial", "unique"}
DEFAULT_SAMPLE_COMPARE_STATUSES = {"present"}
DERIVED_SAMPLE_PREFIX = "derived:"
AF_OPERATORS = {"gt": ">", "gte": ">=", "lt": "<", "lte": "<=", "eq": "=", "neq": "!="}
METADATA_FILTER_FIELDS = {
    "polymorphism": ("mutation_alts.polymorphism", "integer"),
    "repeat_base": ("mutation_alts.repeat_base", "text"),
    "repeat_count": ("mutation_alts.repeat_count", "integer"),
    "repeat_2_bases": ("mutation_alts.repeat_2_bases", "text"),
    "repeat_2_count": ("mutation_alts.repeat_2_count", "integer"),
    "repeat_3_bases": ("mutation_alts.repeat_3_bases", "text"),
    "repeat_3_count": ("mutation_alts.repeat_3_count", "integer"),
}
REFERENCE_REPEAT_BASES = ("A", "C", "G", "T", "N")


@dataclass(frozen=True)
class MutationAllele:
    """Comparable allele identity for viewer-side set operations."""

    position: int
    ref: str
    alt: str


@dataclass
class SampleAlleleCall:
    """One sample's call for a comparable allele."""

    allele: MutationAllele
    sample_id: int
    label: str
    af: float | None
    af_text: str
    filter: str
    vcf_ref: str
    metadata: dict = field(default_factory=dict)

    def to_json(self):
        return {
            "sample_id": self.sample_id,
            "label": self.label,
            "af": self.af,
            "af_text": self.af_text,
            "filter": self.filter,
            "vcf_ref": self.vcf_ref,
            "metadata": self.metadata,
        }


@dataclass
class DerivedSample:
    """In-memory sample built from a comparison result."""

    id: str
    label: str
    calls: list[SampleAlleleCall]
    mutations: list[vcf_parser.VCFMutation]
    source_description: str

    @property
    def subject_id(self):
        return "Derived"

    @property
    def population_key(self):
        return self.label

    @property
    def source_file(self):
        return self.source_description

    @property
    def mutation_count(self):
        return len({call.allele for call in self.calls})

    @property
    def vcf_iterator(self):
        return vcf_parser.VCFIterator.from_mutations(
            self.mutations,
            path=self.id,
        )

    def sample_row(self):
        return {
            "id": self.id,
            "subject_id": self.subject_id,
            "population_key": self.population_key,
            "source_file": self.source_file,
            "mutation_count": self.mutation_count,
            "is_derived": True,
        }


@dataclass
class SampleAlleleSet:
    """Allele set plus per-allele call details for one sample group."""

    calls_by_allele: dict[MutationAllele, list[SampleAlleleCall]] = field(
        default_factory=dict
    )

    def add(self, call):
        self.calls_by_allele.setdefault(call.allele, []).append(call)

    def __contains__(self, allele):
        return allele in self.calls_by_allele

    def __iter__(self):
        return iter(self.calls_by_allele)

    def __and__(self, other):
        return set(self) & set(other)

    def __sub__(self, other):
        return set(self) - set(other)

    def __or__(self, other):
        return set(self) | set(other)

    def calls(self, allele):
        return self.calls_by_allele.get(allele, [])


def connect_database(db_path):
    connection = sqlite3.connect(db_path)
    connection.row_factory = sqlite3.Row
    return connection


def is_derived_sample_id(sample_id):
    return str(sample_id).startswith(DERIVED_SAMPLE_PREFIX)


def parse_af_rules(values):
    """Parse AF rule query values like gt:0.8 into normalized tuples."""
    rules = []
    for value in values or []:
        operator, separator, threshold = str(value).partition(":")
        if not separator or operator not in AF_OPERATORS:
            continue
        try:
            numeric_threshold = float(threshold)
        except ValueError:
            continue
        rules.append((operator, numeric_threshold))
    return rules


def af_rule_matches(value, operator, threshold):
    if value is None:
        return False
    if operator == "gt":
        return value > threshold
    if operator == "gte":
        return value >= threshold
    if operator == "lt":
        return value < threshold
    if operator == "lte":
        return value <= threshold
    if operator == "eq":
        return value == threshold
    if operator == "neq":
        return value != threshold
    return False


def af_rules_match_values(values, af_rules):
    """Return True when every AF rule matches at least one AF value."""
    if not af_rules:
        return True
    numeric_values = [
        value
        for value in values
        if isinstance(value, (int, float))
    ]
    return all(
        any(af_rule_matches(value, operator, threshold) for value in numeric_values)
        for operator, threshold in af_rules
    )


def af_rules_match_text(af_text, af_rules):
    values = []
    for value in str(af_text or "").split(","):
        try:
            values.append(float(value))
        except ValueError:
            continue
    return af_rules_match_values(values, af_rules)


def parse_metadata_filters(values):
    filters = []
    for value in values or []:
        field, separator, raw_value = str(value).partition(":")
        if not separator:
            continue
        filters.append((field, raw_value))
    return filters


def table_columns(connection, table):
    return {
        row["name"]
        for row in connection.execute(f"PRAGMA table_info({table})").fetchall()
    }


def single_base_repeat_seen(sequence):
    """Return True when a reference context contains an adjacent single-base run."""
    sequence = str(sequence or "").upper()
    return any(base * 2 in sequence for base in REFERENCE_REPEAT_BASES)


def reference_repeat_sql_expression(json_key):
    checks = [
        f"json_extract(mutations.metadata_json, '$.{json_key}') LIKE ?"
        for _base in REFERENCE_REPEAT_BASES
    ]
    params = [f"%{base * 2}%" for base in REFERENCE_REPEAT_BASES]
    return "(" + " OR ".join(checks) + ")", params


def reference_repeat_filter_sql(raw_value):
    before_sql, before_params = reference_repeat_sql_expression("REFERENCE_6_BEFORE")
    after_sql, after_params = reference_repeat_sql_expression("REFERENCE_6_AFTER")

    if raw_value == "before":
        return f"({before_sql} AND NOT {after_sql})", before_params + after_params
    if raw_value == "after":
        return f"({after_sql} AND NOT {before_sql})", after_params + before_params
    if raw_value == "one":
        return (
            f"(({before_sql} AND NOT {after_sql}) OR ({after_sql} AND NOT {before_sql}))",
            before_params + after_params + after_params + before_params,
        )
    if raw_value == "both":
        return f"({before_sql} AND {after_sql})", before_params + after_params
    if raw_value == "none":
        return f"(NOT {before_sql} AND NOT {after_sql})", before_params + after_params
    if raw_value == "either":
        return f"({before_sql} OR {after_sql})", before_params + after_params

    return "", []


def add_metadata_filter_sql(connection, where_clauses, params, metadata_filters):
    if not metadata_filters:
        return False

    alt_columns = table_columns(connection, "mutation_alts")
    joined_alts = False
    for field, raw_value in metadata_filters:
        if field == "reference_context":
            where_clauses.append(
                "(json_extract(mutations.metadata_json, '$.REFERENCE_6_BEFORE') LIKE ? "
                "OR json_extract(mutations.metadata_json, '$.REFERENCE_6_AFTER') LIKE ?)"
            )
            params.extend([f"%{raw_value}%", f"%{raw_value}%"])
            continue

        if field == "reference_repeat":
            clause, clause_params = reference_repeat_filter_sql(raw_value)
            if clause:
                where_clauses.append(clause)
                params.extend(clause_params)
            continue

        if field not in METADATA_FILTER_FIELDS:
            continue
        column, field_type = METADATA_FILTER_FIELDS[field]
        column_name = column.split(".")[-1]
        if column_name not in alt_columns:
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
                where_clauses.append(f"{column} {AF_OPERATORS[operator]} ?")
                params.append(int(threshold))
        else:
            where_clauses.append(f"{column} = ?")
            params.append(raw_value)

    return joined_alts


def metadata_filters_match(metadata, metadata_filters):
    if not metadata_filters:
        return True
    for field, raw_value in metadata_filters:
        if field == "polymorphism":
            if str(metadata.get("POLYMORPHISM", "")) != raw_value:
                return False
        elif field == "reference_context":
            context = (
                str(metadata.get("REFERENCE_6_BEFORE", ""))
                + str(metadata.get("REFERENCE_6_AFTER", ""))
            )
            if raw_value not in context:
                return False
        elif field == "reference_repeat":
            before_seen = single_base_repeat_seen(metadata.get("REFERENCE_6_BEFORE", ""))
            after_seen = single_base_repeat_seen(metadata.get("REFERENCE_6_AFTER", ""))
            if raw_value == "before" and not (before_seen and not after_seen):
                return False
            if raw_value == "after" and not (after_seen and not before_seen):
                return False
            if raw_value == "one" and (before_seen == after_seen):
                return False
            if raw_value == "both" and not (before_seen and after_seen):
                return False
            if raw_value == "none" and (before_seen or after_seen):
                return False
            if raw_value == "either" and not (before_seen or after_seen):
                return False
        else:
            key = {
                "repeat_base": "REPEAT_1_BASE",
                "repeat_count": "REPEAT_1_BASE_COUNT",
                "repeat_2_bases": "REPEAT_2_BASES",
                "repeat_2_count": "REPEAT_2_BASES_COUNT",
                "repeat_3_bases": "REPEAT_3_BASES",
                "repeat_3_count": "REPEAT_3_BASES_COUNT",
            }.get(field)
            if key is None:
                continue
            if field.endswith("count"):
                operator, separator, threshold = raw_value.partition("|")
                if not separator or operator not in AF_OPERATORS:
                    continue
                if key not in metadata or metadata.get(key, "") == "":
                    return False
                try:
                    value = int(metadata.get(key))
                    threshold_value = int(threshold)
                except ValueError:
                    return False
                if not af_rule_matches(value, operator, threshold_value):
                    return False
            elif str(metadata.get(key, "")) != raw_value:
                return False
    return True


def fetch_subjects(connection):
    rows = connection.execute("""
        SELECT
            subjects.id,
            subjects.subject_id,
            COUNT(samples.id) AS sample_count
        FROM subjects
        LEFT JOIN samples ON samples.subject_id = subjects.id
        GROUP BY subjects.id
        ORDER BY subjects.subject_id
    """).fetchall()
    return [dict(row) for row in rows]


def fetch_samples(connection, subject_id=None, tags=None, derived_samples=None):
    join_params = []
    where_params = []
    where_clauses = []
    tag_joins = []

    if subject_id:
        where_clauses.append("subjects.subject_id = ?")
        where_params.append(subject_id)

    for index, tag in enumerate(tags or []):
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

    rows = connection.execute(f"""
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
    """, join_params + where_params).fetchall()
    samples = [dict(row) for row in rows]

    derived_samples = derived_samples or {}
    for sample in derived_samples.values():
        samples.append(sample.sample_row())

    return samples


def derived_mutation_rows(sample, position=None, alt=None, af_rules=None, metadata_filters=None, limit=500):
    rows = []
    seen = set()
    metadata_by_allele = {
        MutationAllele(
            position=mutation.position,
            ref=mutation.ref,
            alt=mutation.alt,
        ): mutation.metadata
        for mutation in sample.mutations
    }
    for call in sorted(sample.calls, key=lambda item: (item.allele.position, item.allele.alt)):
        allele = call.allele
        key = (allele.position, allele.ref, allele.alt)
        if key in seen:
            continue
        if position and str(allele.position) != str(position):
            continue
        if alt and allele.alt != alt:
            continue
        if not af_rules_match_text(call.af_text, af_rules or []):
            continue
        metadata = metadata_by_allele.get(allele, {
            "DERIVED_EXPRESSION": sample.source_description,
            "DERIVED_LABEL": sample.label,
        })
        if not metadata_filters_match(metadata, metadata_filters or []):
            continue

        seen.add(key)
        rows.append({
            "id": f"{sample.id}|{allele.position}|{allele.ref}|{allele.alt}",
            "subject_id": sample.subject_id,
            "population_key": sample.population_key,
            "source_file": sample.source_file,
            "pos": allele.position,
            "ref": allele.ref,
            "vcf_ref": call.vcf_ref,
            "alt": allele.alt,
            "af": call.af_text,
            "filter": call.filter,
            "metadata_json": json.dumps(metadata, sort_keys=True),
        })
        if len(rows) >= limit:
            break

    return rows


def fetch_mutations(
    connection,
    sample_id=None,
    position=None,
    alt=None,
    af_rules=None,
    metadata_filters=None,
    limit=500,
    derived_samples=None,
):
    if sample_id and is_derived_sample_id(sample_id):
        sample = (derived_samples or {}).get(str(sample_id))
        if sample is None:
            return []
        return derived_mutation_rows(
            sample,
            position=position,
            alt=alt,
            af_rules=af_rules,
            metadata_filters=metadata_filters,
            limit=limit,
        )

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
        alt_join = "JOIN mutation_alts ON mutation_alts.mutation_id = mutations.id"
        use_alt_rows = True
        where_clauses.append("mutation_alts.alt = ?")
        params.append(alt)

    if af_rules:
        alt_join = "JOIN mutation_alts ON mutation_alts.mutation_id = mutations.id"
        use_alt_rows = True
        for operator, threshold in af_rules:
            sql_operator = AF_OPERATORS[operator]
            where_clauses.append(f"mutation_alts.af {sql_operator} ?")
            params.append(threshold)

    if metadata_filters:
        if add_metadata_filter_sql(connection, where_clauses, params, metadata_filters):
            alt_join = "JOIN mutation_alts ON mutation_alts.mutation_id = mutations.id"
            use_alt_rows = True

    where_sql = ""
    if where_clauses:
        where_sql = "WHERE " + " AND ".join(where_clauses)

    alt_select = "mutation_alts.alt" if use_alt_rows else "mutations.alt"
    af_select = "mutation_alts.af_text" if use_alt_rows else "mutations.af"

    params.append(limit)
    rows = connection.execute(f"""
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
        ORDER BY mutations.pos, subjects.subject_id, samples.population_key
        LIMIT ?
    """, params).fetchall()
    return [dict(row) for row in rows]


def filter_derived_calls(sample, position=None, alt=None, af_rules=None, metadata_filters=None):
    metadata_by_allele = {
        MutationAllele(
            position=mutation.position,
            ref=mutation.ref,
            alt=mutation.alt,
        ): mutation.metadata
        for mutation in sample.mutations
    }
    calls = []
    for call in sample.calls:
        if position and str(call.allele.position) != str(position):
            continue
        if alt and call.allele.alt != alt:
            continue
        if not af_rules_match_text(call.af_text, af_rules or []):
            continue
        if not metadata_filters_match(
            metadata_by_allele.get(call.allele, {}),
            metadata_filters or [],
        ):
            continue
        calls.append(call)
    return calls


def fetch_allele_calls(
    connection,
    sample_ids,
    position=None,
    alt=None,
    af_rules=None,
    metadata_filters=None,
    derived_samples=None,
):
    """Fetch individual allele calls for the given sample ids."""
    derived_samples = derived_samples or {}
    real_sample_ids = [
        str(sample_id)
        for sample_id in sample_ids
        if not is_derived_sample_id(sample_id)
    ]
    derived_sample_ids = [
        str(sample_id)
        for sample_id in sample_ids
        if is_derived_sample_id(sample_id)
    ]

    calls = []
    for sample_id in derived_sample_ids:
        sample = derived_samples.get(sample_id)
        if sample is not None:
            calls.extend(
                filter_derived_calls(
                    sample,
                    position=position,
                    alt=alt,
                    af_rules=af_rules,
                    metadata_filters=metadata_filters,
                )
            )

    if not real_sample_ids:
        return calls

    placeholders = ",".join("?" for _ in real_sample_ids)
    params = list(real_sample_ids)
    where_clauses = [f"samples.id IN ({placeholders})"]

    if position:
        where_clauses.append("mutations.pos = ?")
        params.append(position)

    if alt:
        where_clauses.append("mutation_alts.alt = ?")
        params.append(alt)

    if af_rules:
        for operator, threshold in af_rules:
            sql_operator = AF_OPERATORS[operator]
            where_clauses.append(f"mutation_alts.af {sql_operator} ?")
            params.append(threshold)

    if metadata_filters:
        add_metadata_filter_sql(connection, where_clauses, params, metadata_filters)

    alt_columns = table_columns(connection, "mutation_alts")
    alt_metadata_fields = [
        ("polymorphism", "alt_polymorphism", "POLYMORPHISM"),
        ("repeat_base", "alt_repeat_base", "REPEAT_1_BASE"),
        ("repeat_count", "alt_repeat_count", "REPEAT_1_BASE_COUNT"),
        ("repeat_2_bases", "alt_repeat_2_bases", "REPEAT_2_BASES"),
        ("repeat_2_count", "alt_repeat_2_count", "REPEAT_2_BASES_COUNT"),
        ("repeat_3_bases", "alt_repeat_3_bases", "REPEAT_3_BASES"),
        ("repeat_3_count", "alt_repeat_3_count", "REPEAT_3_BASES_COUNT"),
    ]
    alt_metadata_selects = [
        f"mutation_alts.{column} AS {alias}"
        for column, alias, _key in alt_metadata_fields
        if column in alt_columns
    ]
    metadata_select_sql = ""
    if alt_metadata_selects:
        metadata_select_sql = ",\n            " + ",\n            ".join(alt_metadata_selects)

    rows = connection.execute(f"""
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
            mutations.metadata_json
            {metadata_select_sql}
        FROM mutation_alts
        JOIN mutations ON mutations.id = mutation_alts.mutation_id
        JOIN samples ON samples.id = mutations.sample_id
        JOIN subjects ON subjects.id = samples.subject_id
        WHERE {" AND ".join(where_clauses)}
        ORDER BY mutations.pos, mutation_alts.alt, subjects.subject_id, samples.population_key
    """, params).fetchall()

    for row in rows:
        allele = MutationAllele(
            position=row["pos"],
            ref=row["ref"],
            alt=row["alt"],
        )
        row_keys = row.keys()
        try:
            source_metadata = json.loads(row["metadata_json"] or "{}")
        except json.JSONDecodeError:
            source_metadata = {}
        call_metadata = {
            key: source_metadata.get(key, "")
            for key in ("REFERENCE_6_BEFORE", "REFERENCE_6_AFTER")
            if key in source_metadata
        }
        for _column, alias, key in alt_metadata_fields:
            if alias in row_keys and row[alias] not in (None, ""):
                call_metadata[key] = str(row[alias])
        calls.append(SampleAlleleCall(
            allele=allele,
            sample_id=row["sample_id"],
            label=f"{row['subject_id']} {row['population_key'].replace('|', '_')}",
            af=row["af"],
            af_text=row["af_text"],
            filter=row["filter"],
            vcf_ref=row["vcf_ref"],
            metadata=call_metadata,
        ))
    return calls


def build_allele_set(calls):
    """Build a comparable allele set from allele calls."""
    allele_set = SampleAlleleSet()
    for call in calls:
        allele_set.add(call)
    return allele_set


def fetch_sample_labels(connection, sample_ids, derived_samples=None):
    """Return display labels for sample ids."""
    if not sample_ids:
        return {}

    derived_samples = derived_samples or {}
    labels = {
        sample_id: f"{sample.subject_id} {sample.population_key.replace('|', '_')}"
        for sample_id, sample in derived_samples.items()
        if sample_id in set(str(item) for item in sample_ids)
    }
    real_sample_ids = [
        str(sample_id)
        for sample_id in sample_ids
        if not is_derived_sample_id(sample_id)
    ]
    if not real_sample_ids:
        return labels

    placeholders = ",".join("?" for _ in real_sample_ids)
    rows = connection.execute(f"""
        SELECT
            samples.id,
            subjects.subject_id,
            samples.population_key
        FROM samples
        JOIN subjects ON subjects.id = samples.subject_id
        WHERE samples.id IN ({placeholders})
    """, real_sample_ids).fetchall()

    labels.update({
        str(row["id"]): f"{row['subject_id']} {row['population_key'].replace('|', '_')}"
        for row in rows
    })
    return labels


def parse_sample_statuses(values):
    """Parse sample-specific status filters from id:status query values."""
    sample_statuses = {}
    for value in values:
        sample_id, separator, status = value.partition(":")
        if not separator or not sample_id:
            continue

        statuses = sample_statuses.setdefault(sample_id, set())
        if status != "__none__":
            statuses.add(status)

    return sample_statuses


def sample_constraint_matches(allowed_statuses, is_present, present_count):
    """Return True when one sample's direct set constraint is satisfied."""
    if not allowed_statuses:
        return True

    return (
        ("present" in allowed_statuses and is_present)
        or ("unique" in allowed_statuses and is_present and present_count == 1)
        or ("not_in" in allowed_statuses and not is_present)
    )


def sample_filters_match(compare_sample_ids, sample_statuses, present_sample_ids, present_count):
    """Return True when direct per-sample constraints match an allele."""
    present_required = set()
    absent_required = set()
    unique_allowed = set()

    for sample_id in compare_sample_ids:
        allowed_statuses = sample_statuses.get(sample_id, DEFAULT_SAMPLE_COMPARE_STATUSES)
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


def compare_row(
    allele,
    status,
    present_calls,
    missing_samples,
):
    """Return one JSON-ready comparison row."""
    return {
        "pos": allele.position,
        "ref": allele.ref,
        "alt": allele.alt,
        "group_key": f"{status}|{allele.position}|{allele.ref}|{allele.alt}",
        "present": [call.to_json() for call in present_calls],
        "missing": missing_samples,
        "status": status,
    }


def fetch_compare(
    connection,
    compare_sample_ids,
    position=None,
    alt=None,
    af_rules=None,
    metadata_filters=None,
    statuses=None,
    sample_statuses=None,
    limit=2000,
    derived_samples=None,
):
    """Fetch a colorizable peer comparison for selected samples."""
    compare_sample_ids = [
        str(sample_id)
        for sample_id in compare_sample_ids
    ]
    if len(compare_sample_ids) < 2:
        return []

    sample_labels = fetch_sample_labels(
        connection,
        compare_sample_ids,
        derived_samples=derived_samples,
    )
    calls = fetch_allele_calls(
        connection,
        compare_sample_ids,
        position=position,
        alt=alt,
        af_rules=af_rules,
        metadata_filters=metadata_filters,
        derived_samples=derived_samples,
    )

    global_statuses = set(statuses or DEFAULT_COMPARE_STATUSES)
    sample_statuses = sample_statuses or {}
    calls_by_allele = {}
    for call in calls:
        calls_by_allele.setdefault(call.allele, []).append(call)

    results = []
    for allele, allele_calls in calls_by_allele.items():
        present_sample_ids = {
            str(call.sample_id)
            for call in allele_calls
        }
        present_count = len(present_sample_ids)
        if present_count == len(compare_sample_ids):
            status = "common"
        elif present_count == 1:
            status = "unique"
        else:
            status = "partial"

        if status not in global_statuses:
            continue

        if not sample_filters_match(
            compare_sample_ids,
            sample_statuses,
            present_sample_ids,
            present_count,
        ):
            continue

        missing_samples = [
            {
                "sample_id": sample_id,
                "label": sample_labels.get(sample_id, f"Sample {sample_id}"),
            }
            for sample_id in compare_sample_ids
            if sample_id not in present_sample_ids
        ]
        results.append(compare_row(
            allele,
            status,
            allele_calls,
            missing_samples,
        ))

    results.sort(key=lambda row: (
        row["pos"],
        row["alt"],
        row["status"],
    ))
    results = results[:limit]

    group_sizes = {}
    for row in results:
        group_sizes[row["group_key"]] = group_sizes.get(row["group_key"], 0) + 1
    for index, row in enumerate(results):
        previous_key = results[index - 1]["group_key"] if index else None
        next_key = results[index + 1]["group_key"] if index + 1 < len(results) else None
        row["group_size"] = group_sizes[row["group_key"]]
        row["group_start"] = row["group_key"] != previous_key
        row["group_end"] = row["group_key"] != next_key

    return results


def database_counts(connection):
    counts = {}
    for table in [
        "subjects",
        "samples",
        "sample_population_tags",
        "mutations",
        "mutation_alts",
    ]:
        counts[table] = connection.execute(
            f"SELECT COUNT(*) FROM {table}"
        ).fetchone()[0]
    return counts


def json_response(handler, payload, status=200):
    body = json.dumps(payload, indent=2).encode("utf-8")
    handler.send_response(status)
    handler.send_header("Content-Type", "application/json; charset=utf-8")
    handler.send_header("Content-Length", str(len(body)))
    handler.end_headers()
    handler.wfile.write(body)


def html_response(handler, body, status=200):
    encoded = body.encode("utf-8")
    handler.send_response(status)
    handler.send_header("Content-Type", "text/html; charset=utf-8")
    handler.send_header("Content-Length", str(len(encoded)))
    handler.end_headers()
    handler.wfile.write(encoded)


def parse_tags(query):
    tags = []
    for value in query.get("tag", []):
        tags.extend(tag for tag in value.split(",") if tag)
    return tags


def unique_nonempty(values):
    unique_values = []
    for value in values:
        if value in (None, ""):
            continue
        text = str(value)
        if text not in unique_values:
            unique_values.append(text)
    return unique_values


def ordered_statuses(statuses):
    """Return compare statuses in a stable, readable order."""
    order = ["present", "unique", "not_in", "common", "partial", "__none__"]
    status_set = set(statuses or [])
    return [
        status for status in order if status in status_set
    ] + sorted(status_set - set(order))


def status_expression(status, labels):
    joined_labels = ", ".join(labels)
    if status == "common":
        return f"AND({joined_labels})"
    if status == "partial":
        return f"SOME_NOT_ALL({joined_labels})"
    if status == "unique":
        operator = "XOR" if len(labels) == 2 else "EXACTLY_ONE"
        return f"{operator}({joined_labels})"
    if status == "__none__":
        return "EMPTY_SET"
    return f"{status.upper()}({joined_labels})"


def sample_presence_sets(compare_sample_ids, source_labels, sample_statuses):
    """Return explicit present/absent sample filters when they are simple."""
    present = []
    absent = []
    unique = []
    complex_filters = []

    for sample_id in compare_sample_ids:
        allowed_statuses = sample_statuses.get(str(sample_id)) if sample_statuses else None
        if allowed_statuses is None:
            continue

        label = source_labels.get(str(sample_id), f"Sample {sample_id}")
        status_set = set(allowed_statuses)
        if "present" in status_set:
            present.append(label)
        if "unique" in status_set:
            unique.append(label)
        if "not_in" in status_set:
            absent.append(label)
        if status_set - {"present", "unique", "not_in"}:
            complex_filters.append(
                f"{label} IS {'/'.join(status.upper() for status in ordered_statuses(status_set)) or 'NONE'}"
            )

    return present, unique_nonempty(absent), unique, complex_filters


def comparison_metadata(
    compare_sample_ids,
    source_labels,
    statuses=None,
    sample_statuses=None,
    position=None,
    alt=None,
):
    """Return readable metadata describing a materialized comparison set."""
    labels = [
        source_labels.get(str(sample_id), f"Sample {sample_id}")
        for sample_id in compare_sample_ids
    ]
    selected_statuses = ordered_statuses(statuses or DEFAULT_COMPARE_STATUSES)
    default_statuses = ordered_statuses(DEFAULT_COMPARE_STATUSES)
    if selected_statuses == default_statuses:
        expression = f"OR({', '.join(labels)})"
    elif len(selected_statuses) == 1:
        expression = status_expression(selected_statuses[0], labels)
    else:
        expression = " OR ".join(
            status_expression(status, labels)
            for status in selected_statuses
        )

    present_samples, absent_samples, unique_samples, sample_filters = sample_presence_sets(
        compare_sample_ids,
        source_labels,
        sample_statuses or {},
    )

    if unique_samples:
        expression_parts = []
        if present_samples:
            if len(present_samples) == 1:
                expression_parts.append(present_samples[0])
            else:
                expression_parts.append(f"AND({', '.join(present_samples)})")
        expression_parts.extend(f"UNIQUE({sample})" for sample in unique_samples)
        expression = " OR ".join(expression_parts)
        if absent_samples:
            expression = f"({expression}) AND NOT({', '.join(absent_samples)})"
    elif present_samples or absent_samples:
        included_samples = present_samples or [
            label for label in labels if label not in absent_samples
        ]
        if len(included_samples) == 1:
            expression = included_samples[0]
        else:
            expression = f"AND({', '.join(included_samples)})"
        if absent_samples:
            expression = f"{expression} AND NOT({', '.join(absent_samples)})"

    if sample_filters:
        expression = f"{expression} WHERE {' AND '.join(sample_filters)}"

    filters = []
    if position:
        filters.append(f"POS={position}")
    if alt:
        filters.append(f"ALT={alt}")
    if filters:
        expression = f"{expression} FILTER {' AND '.join(filters)}"

    return {
        "expression": expression,
        "source_samples": labels,
        "statuses": selected_statuses,
        "sample_filters": sample_filters,
        "present_samples": present_samples,
        "absent_samples": absent_samples,
        "unique_samples": unique_samples,
        "position_filter": str(position or ""),
        "alt_filter": str(alt or ""),
    }


def create_derived_sample(
    connection,
    derived_samples,
    derived_id,
    label,
    compare_sample_ids,
    position=None,
    alt=None,
    af_rules=None,
    metadata_filters=None,
    statuses=None,
    sample_statuses=None,
    limit=2000,
):
    rows = fetch_compare(
        connection,
        compare_sample_ids=compare_sample_ids,
        position=position,
        alt=alt,
        af_rules=af_rules,
        metadata_filters=metadata_filters,
        statuses=statuses,
        sample_statuses=sample_statuses,
        limit=limit,
        derived_samples=derived_samples,
    )
    source_labels = fetch_sample_labels(
        connection,
        compare_sample_ids,
        derived_samples=derived_samples,
    )
    metadata = comparison_metadata(
        compare_sample_ids,
        source_labels,
        statuses=statuses,
        sample_statuses=sample_statuses or {},
        position=position,
        alt=alt,
    )
    source_description = metadata["expression"]

    calls = []
    mutations = []
    seen = set()
    for row in rows:
        allele = MutationAllele(
            position=row["pos"],
            ref=row["ref"],
            alt=row["alt"],
        )
        if allele in seen:
            continue
        seen.add(allele)

        present_calls = row.get("present", [])
        af_text = ",".join(unique_nonempty(call.get("af_text") for call in present_calls))
        filters = unique_nonempty(call.get("filter") for call in present_calls)
        vcf_refs = unique_nonempty(call.get("vcf_ref") for call in present_calls)
        numeric_afs = [
            call.get("af")
            for call in present_calls
            if isinstance(call.get("af"), (int, float))
        ]
        filter_text = ",".join(filters)
        allele_metadata = {}
        for call in present_calls:
            for key, value in (call.get("metadata") or {}).items():
                if key not in allele_metadata and value not in (None, ""):
                    allele_metadata[key] = value
        calls.append(SampleAlleleCall(
            allele=allele,
            sample_id=derived_id,
            label=f"Derived {label}",
            af=numeric_afs[0] if len(numeric_afs) == 1 else None,
            af_text=af_text,
            filter=filter_text,
            vcf_ref=vcf_refs[0] if vcf_refs else allele.ref,
            metadata=allele_metadata,
        ))
        mutations.append(vcf_parser.VCFMutation(
            position=allele.position,
            alt=allele.alt,
            metadata={
                **allele_metadata,
                "AF": af_text,
                "DERIVED_EXPRESSION": metadata["expression"],
                "DERIVED_FROM": source_description,
                "DERIVED_LABEL": label,
                "DERIVED_SET_STATUSES": ",".join(metadata["statuses"]),
                "DERIVED_SAMPLE_FILTERS": ";".join(metadata["sample_filters"]),
                "DERIVED_PRESENT_SAMPLES": ";".join(metadata["present_samples"]),
                "DERIVED_ABSENT_SAMPLES": ";".join(metadata["absent_samples"]),
                "DERIVED_UNIQUE_SAMPLES": ";".join(metadata["unique_samples"]),
                "DERIVED_SOURCE_SAMPLES": ";".join(metadata["source_samples"]),
            },
            ref=allele.ref,
            filter=filter_text,
        ))

    sample = DerivedSample(
        id=derived_id,
        label=label,
        calls=calls,
        mutations=mutations,
        source_description=source_description,
    )
    derived_samples[derived_id] = sample
    return sample


def page_html(db_path):
    escaped_path = html.escape(str(db_path))
    template = VIEWER_TEMPLATE_PATH.read_text(encoding="utf-8")
    return template.replace("{{ db_path }}", escaped_path)


class ViewerHandler(BaseHTTPRequestHandler):
    db_path = DEFAULT_DB_PATH
    derived_samples = {}
    next_derived_sample_id = 1

    def log_message(self, fmt, *args):
        print(f"{self.address_string()} - {fmt % args}")

    def open_database(self):
        if not self.db_path.exists():
            raise FileNotFoundError(
                f"Database not found: {self.db_path}. "
                "Create it with main.make_sql_database() first."
            )
        return connect_database(self.db_path)

    def read_json_body(self):
        content_length = int(self.headers.get("Content-Length", "0"))
        if content_length <= 0:
            return {}
        body = self.rfile.read(content_length).decode("utf-8")
        return json.loads(body or "{}")

    @classmethod
    def allocate_derived_sample_id(cls):
        derived_id = f"{DERIVED_SAMPLE_PREFIX}{cls.next_derived_sample_id}"
        cls.next_derived_sample_id += 1
        return derived_id

    def do_GET(self):
        parsed = urlparse(self.path)
        query = parse_qs(parsed.query)

        try:
            if parsed.path == "/":
                html_response(self, page_html(self.db_path))
            elif parsed.path == "/api/counts":
                with self.open_database() as connection:
                    json_response(self, database_counts(connection))
            elif parsed.path == "/api/subjects":
                with self.open_database() as connection:
                    json_response(self, fetch_subjects(connection))
            elif parsed.path == "/api/samples":
                with self.open_database() as connection:
                    json_response(
                        self,
                        fetch_samples(
                            connection,
                            subject_id=query.get("subject", [""])[0],
                            tags=parse_tags(query),
                            derived_samples=self.derived_samples,
                        ),
                    )
            elif parsed.path == "/api/mutations":
                with self.open_database() as connection:
                    json_response(
                        self,
                        fetch_mutations(
                            connection,
                            sample_id=query.get("sample_id", [""])[0],
                            position=query.get("position", [""])[0],
                            alt=query.get("alt", [""])[0],
                            af_rules=parse_af_rules(query.get("af_rule", [])),
                            metadata_filters=parse_metadata_filters(
                                query.get("metadata_filter", [])
                            ),
                            limit=int(query.get("limit", ["500"])[0]),
                            derived_samples=self.derived_samples,
                        ),
                    )
            elif parsed.path == "/api/compare":
                with self.open_database() as connection:
                    json_response(
                        self,
                        fetch_compare(
                            connection,
                            compare_sample_ids=query.get("compare_sample_id", []),
                            position=query.get("position", [""])[0],
                            alt=query.get("alt", [""])[0],
                            af_rules=parse_af_rules(query.get("af_rule", [])),
                            metadata_filters=parse_metadata_filters(
                                query.get("metadata_filter", [])
                            ),
                            statuses=query.get("status", []),
                            sample_statuses=parse_sample_statuses(
                                query.get("sample_status", [])
                            ),
                            limit=int(query.get("limit", ["2000"])[0]),
                            derived_samples=self.derived_samples,
                        ),
                    )
            else:
                json_response(self, {"error": "Not found"}, status=404)
        except Exception as exc:
            json_response(self, {"error": str(exc)}, status=500)

    def do_POST(self):
        parsed = urlparse(self.path)

        try:
            if parsed.path == "/api/derived-samples":
                payload = self.read_json_body()
                compare_sample_ids = payload.get("compare_sample_id", [])
                if len(compare_sample_ids) < 2:
                    json_response(
                        self,
                        {"error": "Select at least two samples before creating a derived sample."},
                        status=400,
                    )
                    return

                label = payload.get("label") or f"Comparison {self.next_derived_sample_id}"
                with self.open_database() as connection:
                    sample = create_derived_sample(
                        connection,
                        self.derived_samples,
                        self.allocate_derived_sample_id(),
                        label,
                        compare_sample_ids=compare_sample_ids,
                        position=payload.get("position", ""),
                        alt=payload.get("alt", ""),
                        af_rules=parse_af_rules(payload.get("af_rule", [])),
                        metadata_filters=parse_metadata_filters(
                            payload.get("metadata_filter", [])
                        ),
                        statuses=payload.get("status", []),
                        sample_statuses=parse_sample_statuses(
                            payload.get("sample_status", [])
                        ),
                        limit=int(payload.get("limit", 2000)),
                    )
                json_response(self, sample.sample_row(), status=201)
            else:
                json_response(self, {"error": "Not found"}, status=404)
        except Exception as exc:
            json_response(self, {"error": str(exc)}, status=500)

    def do_DELETE(self):
        parsed = urlparse(self.path)

        try:
            prefix = "/api/derived-samples/"
            if parsed.path.startswith(prefix):
                sample_id = unquote(parsed.path[len(prefix):])
                if not is_derived_sample_id(sample_id):
                    json_response(
                        self,
                        {"error": "Only derived samples can be deleted."},
                        status=400,
                    )
                    return
                removed = self.derived_samples.pop(sample_id, None)
                if removed is None:
                    json_response(self, {"error": "Derived sample not found."}, status=404)
                    return
                json_response(self, removed.sample_row())
            else:
                json_response(self, {"error": "Not found"}, status=404)
        except Exception as exc:
            json_response(self, {"error": str(exc)}, status=500)


def run_server(db_path=DEFAULT_DB_PATH, host="127.0.0.1", port=8000):
    ViewerHandler.db_path = Path(db_path)
    server = ThreadingHTTPServer((host, port), ViewerHandler)
    print(f"Serving {ViewerHandler.db_path} at http://{host}:{port}")
    server.serve_forever()


def main():
    parser = argparse.ArgumentParser(description="View the mito SQLite database.")
    parser.add_argument("--db", type=Path, default=DEFAULT_DB_PATH)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=18000)
    args = parser.parse_args()
    run_server(args.db, args.host, args.port)


if __name__ == "__main__":
    main()
