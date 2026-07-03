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
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Mito SQL Viewer</title>
  <style>
    :root {{
      color-scheme: light;
      font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: #f6f7f9;
      color: #1d2733;
    }}
    * {{ box-sizing: border-box; }}
    body {{ margin: 0; }}
    header {{
      padding: 18px 22px;
      background: #ffffff;
      border-bottom: 1px solid #d9dee7;
      position: sticky;
      top: 0;
      z-index: 2;
    }}
    h1 {{ margin: 0 0 4px; font-size: 22px; font-weight: 700; }}
    .path {{ color: #607086; font-size: 13px; }}
    main {{ padding: 18px 22px 32px; }}
    .workspace {{
      display: grid;
      grid-template-columns: minmax(0, 1fr) 320px;
      gap: 18px;
      align-items: start;
    }}
    .content {{ display: grid; gap: 18px; min-width: 0; }}
    .compare-panel {{
      position: sticky;
      top: 78px;
      max-height: calc(100vh - 96px);
      display: grid;
      grid-template-rows: auto auto minmax(0, 1fr);
      background: #fff;
      border: 1px solid #d9dee7;
      border-radius: 8px;
      overflow: hidden;
    }}
    .toolbar {{
      display: grid;
      grid-template-columns:
        minmax(120px, 0.8fr)
        minmax(260px, 1.7fr)
        minmax(110px, 0.7fr)
        minmax(110px, 0.7fr)
        minmax(110px, 0.7fr);
      grid-template-areas: "subject tags position alt apply";
      gap: 10px 10px;
      align-items: end;
      padding: 12px;
      min-width: 0;
    }}
    .filter-subject {{ grid-area: subject; }}
    .filter-tags {{ grid-area: tags; }}
    .filter-position {{ grid-area: position; }}
    .filter-alt {{ grid-area: alt; }}
    .filter-apply {{ grid-area: apply; width: 100%; }}
    label {{ display: grid; gap: 5px; font-size: 12px; color: #4a596d; font-weight: 650; min-width: 0; }}
    input, select, button {{
      width: 100%;
      min-height: 36px;
      border: 1px solid #c7cfdb;
      border-radius: 6px;
      background: #fff;
      color: #1d2733;
      font: inherit;
      padding: 7px 9px;
    }}
    button {{
      cursor: pointer;
      background: #2457c5;
      border-color: #2457c5;
      color: #fff;
      font-weight: 700;
    }}
    .delete-derived-sample {{
      width: auto;
      min-height: 26px;
      margin-left: 8px;
      padding: 3px 8px;
      background: #ffffff;
      border-color: #c94f4f;
      color: #9f2f2f;
      font-size: 12px;
    }}
    .order-cell {{ width: 92px; white-space: nowrap; }}
    .drag-handle {{
      width: auto;
      min-height: 26px;
      padding: 3px 8px;
      background: #f7f9fc;
      border-color: #c7cfdb;
      color: #35445a;
      font-size: 12px;
      cursor: grab;
    }}
    .sample-order-button {{
      width: auto;
      min-height: 26px;
      margin-left: 4px;
      padding: 3px 7px;
      background: #ffffff;
      border-color: #c7cfdb;
      color: #35445a;
      font-size: 12px;
    }}
    tr.reorder-target td {{ box-shadow: inset 0 2px 0 #2457c5; }}
    section {{
      background: #fff;
      border: 1px solid #d9dee7;
      border-radius: 8px;
      overflow: hidden;
    }}
    .section-title {{
      padding: 10px 12px;
      border-bottom: 1px solid #e4e8ef;
      font-weight: 700;
      display: flex;
      justify-content: space-between;
      gap: 12px;
    }}
    .counts {{ color: #607086; font-weight: 600; font-size: 13px; }}
    .table-wrap {{ overflow: auto; max-height: 58vh; }}
    table {{ width: 100%; border-collapse: collapse; font-size: 13px; }}
    th, td {{ border-bottom: 1px solid #edf0f4; padding: 7px 9px; text-align: left; vertical-align: top; }}
    th {{ position: sticky; top: 0; background: #f9fafc; z-index: 1; color: #3f4d60; }}
    tr:hover td {{ background: #f7faff; }}
    tr.selected-sample td {{
      background: #eef5ff;
      box-shadow: inset 4px 0 0 #2457c5;
    }}
    .mono {{ font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; }}
    .muted {{ color: #607086; }}
    .metadata-cell {{ max-width: 520px; cursor: help; }}
    .derived-metadata {{
      display: inline-block;
      max-width: 520px;
      white-space: nowrap;
      overflow: hidden;
      text-overflow: ellipsis;
      border-bottom: 1px dotted #607086;
    }}
    .metadata-popover {{
      position: fixed;
      width: min(560px, calc(100vw - 28px));
      max-height: min(560px, calc(100vh - 28px));
      overflow: auto;
      display: none;
      z-index: 20;
      padding: 14px;
      background: #ffffff;
      border: 1px solid #b8c2d2;
      border-radius: 8px;
      box-shadow: 0 16px 40px rgba(28, 39, 55, 0.18);
      color: #1d2733;
    }}
    .metadata-popover.open {{ display: grid; gap: 12px; }}
    .metadata-popover h3 {{ margin: 0; font-size: 15px; }}
    .metadata-popover .expression {{
      padding: 8px 10px;
      border-radius: 6px;
      background: #f7f9fc;
      border: 1px solid #e0e6ef;
      font-size: 12px;
      white-space: normal;
      overflow-wrap: anywhere;
    }}
    .venn-wrap {{ display: grid; grid-template-columns: 220px minmax(0, 1fr); gap: 14px; align-items: center; }}
    .venn-diagram {{ position: relative; width: 220px; height: 170px; }}
    .venn-circle {{
      position: absolute;
      width: 112px;
      height: 112px;
      display: grid;
      place-items: center;
      padding: 8px;
      border: 2px solid #6f88a7;
      border-radius: 50%;
      background: rgba(121, 167, 236, 0.22);
      text-align: center;
      font-size: 10px;
      font-weight: 700;
      color: #26364b;
      overflow: hidden;
    }}
    .venn-circle:nth-child(1) {{ left: 24px; top: 18px; }}
    .venn-circle:nth-child(2) {{ left: 84px; top: 18px; background: rgba(107, 184, 131, 0.22); }}
    .venn-circle:nth-child(3) {{ left: 54px; top: 58px; background: rgba(222, 183, 100, 0.22); }}
    .venn-circle:nth-child(4) {{ left: 54px; top: 0; background: rgba(196, 136, 210, 0.18); }}
    .venn-circle.excluded {{
      background: repeating-linear-gradient(135deg, rgba(201, 79, 79, 0.08), rgba(201, 79, 79, 0.08) 6px, rgba(255, 255, 255, 0.6) 6px, rgba(255, 255, 255, 0.6) 12px);
      border-color: #c94f4f;
    }}
    .venn-circle.excluded::before,
    .venn-circle.excluded::after,
    .intersection-core.crossed::before,
    .intersection-core.crossed::after {{
      content: "";
      position: absolute;
      left: 18%;
      right: 18%;
      top: 50%;
      height: 3px;
      background: #c94f4f;
      border-radius: 999px;
    }}
    .venn-circle.excluded::before,
    .intersection-core.crossed::before {{ transform: rotate(38deg); }}
    .venn-circle.excluded::after,
    .intersection-core.crossed::after {{ transform: rotate(-38deg); }}
    .intersection-core {{
      position: absolute;
      left: 88px;
      top: 70px;
      width: 46px;
      height: 34px;
      border-radius: 999px;
      border: 2px solid #2457c5;
      background: rgba(36, 87, 197, 0.2);
      box-shadow: 0 0 0 4px rgba(255, 255, 255, 0.64);
    }}
    .intersection-core.partial {{
      width: 92px;
      left: 64px;
      border-color: #79a7ec;
      background: rgba(121, 167, 236, 0.18);
    }}
    .venn-notes {{ display: grid; gap: 7px; font-size: 12px; color: #4a596d; }}
    .venn-notes strong {{ color: #1d2733; }}
    .header-filter-button {{
      width: auto;
      min-height: 28px;
      padding: 3px 8px;
      background: #ffffff;
      border-color: #b9c5d6;
      color: #35445a;
      font-size: 12px;
      font-weight: 700;
    }}
    .modal-backdrop {{
      position: fixed;
      inset: 0;
      display: none;
      place-items: center;
      padding: 18px;
      background: rgba(29, 39, 51, 0.28);
      z-index: 30;
    }}
    .modal-backdrop.open {{ display: grid; }}
    .modal {{
      width: min(440px, 100%);
      display: grid;
      gap: 12px;
      padding: 14px;
      background: #ffffff;
      border: 1px solid #b8c2d2;
      border-radius: 8px;
      box-shadow: 0 18px 45px rgba(28, 39, 55, 0.22);
    }}
    .modal h2 {{ margin: 0; font-size: 16px; }}
    .af-rule-editor {{
      display: grid;
      grid-template-columns: minmax(100px, 0.8fr) minmax(120px, 1fr) auto;
      gap: 8px;
      align-items: end;
    }}
    .af-rule-list {{ display: grid; gap: 6px; }}
    .af-rule-item {{
      display: grid;
      grid-template-columns: 1fr auto;
      gap: 8px;
      align-items: center;
      padding: 7px 8px;
      border: 1px solid #e0e6ef;
      border-radius: 6px;
      background: #f9fbfe;
    }}
    .af-rule-remove {{
      width: auto;
      min-height: 28px;
      padding: 3px 8px;
      background: #ffffff;
      border-color: #c94f4f;
      color: #9f2f2f;
      font-size: 12px;
    }}
    .metadata-filter-grid {{
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 10px;
    }}
    .metadata-filter-grid .full {{ grid-column: 1 / -1; }}
    .modal-actions {{ display: grid; grid-template-columns: 1fr 1fr; gap: 8px; }}
    .compare-options {{ padding: 12px; display: grid; gap: 10px; }}
    .compare-actions {{ display: grid; grid-template-columns: 1fr 1fr; gap: 8px; }}
    .secondary-button {{
      background: #ffffff;
      border-color: #c7cfdb;
      color: #35445a;
    }}
    .compare-summary {{
      border: 1px solid #e0e6ef;
      border-radius: 6px;
      background: #f9fbfe;
      overflow: hidden;
      min-width: 0;
    }}
    .compare-summary summary {{
      cursor: pointer;
      padding: 8px 10px;
      color: #35445a;
      font-size: 12px;
      font-weight: 700;
    }}
    .compare-summary-body {{
      display: grid;
      gap: 6px;
      padding: 0 10px 10px;
      min-width: 0;
    }}
    .compare-summary-expression {{
      overflow-wrap: anywhere;
      line-height: 1.35;
      font-size: 12px;
      color: #1d2733;
    }}
    .compare-summary-sentence {{
      overflow-wrap: anywhere;
      line-height: 1.4;
      font-size: 12px;
      color: #35445a;
    }}
    .compare-summary-hint {{
      font-size: 11px;
      color: #607086;
    }}
    .checkline {{
      display: flex;
      align-items: center;
      gap: 8px;
      font-size: 13px;
      color: #35445a;
      font-weight: 600;
    }}
    .checkline input {{ min-height: 0; }}
    .tag-picker {{
      min-height: 36px;
      display: flex;
      flex-wrap: wrap;
      gap: 6px 12px;
      align-items: center;
      align-content: center;
      border: 1px solid #c7cfdb;
      border-radius: 6px;
      background: #fff;
      padding: 7px 10px;
      width: 100%;
    }}
    .tag-picker label {{
      display: flex;
      flex-direction: row;
      align-items: center;
      gap: 5px;
      font-size: 12px;
      font-weight: 600;
      color: #35445a;
      white-space: nowrap;
    }}
    .tag-picker input {{ min-height: 0; padding: 0; }}
    .compare-list {{
      overflow-y: auto;
      overflow-x: hidden;
      border-top: 1px solid #e4e8ef;
      padding: 8px 12px 12px;
      display: grid;
      gap: 7px;
      min-width: 0;
    }}
    .compare-item {{
      display: grid;
      grid-template-columns: 1fr;
      gap: 8px;
      align-items: start;
      border: 1px solid #e4e8ef;
      border-radius: 6px;
      padding: 9px;
      background: #fbfcfe;
      font-size: 13px;
      min-width: 0;
    }}
    .compare-item-header {{
      display: grid;
      grid-template-columns: auto minmax(0, 1fr);
      gap: 8px;
      align-items: start;
      min-width: 0;
    }}
    .compare-item-title {{
      display: grid;
      gap: 3px;
      min-width: 0;
    }}
    .compare-item-title strong {{
      overflow-wrap: anywhere;
      line-height: 1.25;
    }}
    .compare-meta {{
      display: grid;
      gap: 2px;
      color: #607086;
      font-size: 11px;
      line-height: 1.25;
      min-width: 0;
    }}
    .compare-meta span {{
      overflow-wrap: anywhere;
    }}
    .sample-statuses {{
      display: grid;
      grid-template-columns: 1fr;
      gap: 5px;
      padding-top: 7px;
      border-top: 1px solid #e4e8ef;
      min-width: 0;
    }}
    .sample-statuses label {{
      display: flex;
      flex-direction: row;
      align-items: center;
      gap: 5px;
      min-height: 28px;
      padding: 4px 7px;
      border: 1px solid #dce2eb;
      border-radius: 6px;
      background: #ffffff;
      font-size: 12px;
      font-weight: 600;
      color: #4a596d;
      min-width: 0;
    }}
    .sample-statuses label:has(input:checked) {{
      border-color: #8aa9e6;
      background: #eef5ff;
      color: #24436e;
    }}
    .sample-statuses input {{ width: auto; min-height: 0; padding: 0; flex: 0 0 auto; }}
    .legend {{ display: grid; gap: 6px; padding: 0 12px 12px; font-size: 12px; color: #4a596d; }}
    .swatch {{ display: inline-block; width: 12px; height: 12px; border-radius: 3px; margin-right: 6px; vertical-align: -1px; }}
    .swatch.common {{ background: #dff4e6; border: 1px solid #6bb883; }}
    .swatch.partial {{ background: #e0edff; border: 1px solid #79a7ec; }}
    .swatch.unique {{ background: #fff1cf; border: 1px solid #deb764; }}
    tr.common td {{ background: #f0fbf3; }}
    tr.partial td {{ background: #eef5ff; }}
    tr.unique td {{ background: #fff7e3; }}
    tr.common.group-even td {{ background: #edf9f1; }}
    tr.common.group-odd td {{ background: #e3f5ea; }}
    tr.partial.group-even td {{ background: #eef5ff; }}
    tr.partial.group-odd td {{ background: #e3efff; }}
    tr.unique.group-even td {{ background: #fff7e3; }}
    tr.unique.group-odd td {{ background: #ffefc9; }}
    tr.group-start td {{ border-top: 2px solid #8fb1d9; }}
    tr.group-end td {{ border-bottom: 2px solid #d1dae7; }}
    tr.grouped td:first-child {{
      border-left: 5px solid #6bb883;
      font-weight: 700;
    }}
    tr.partial.grouped td:first-child {{ border-left-color: #79a7ec; }}
    tr.unique.grouped td:first-child {{ border-left-color: #deb764; }}
    tr.grouped.group-middle td:first-child, tr.grouped.group-end td:first-child {{
      color: transparent;
    }}
    tr.common:hover td {{ background: #e1f6e8; }}
    tr.partial:hover td {{ background: #deebff; }}
    tr.unique:hover td {{ background: #ffefc2; }}
    @media (max-width: 1250px) {{
      .toolbar {{
        grid-template-columns: minmax(130px, 1fr) minmax(260px, 2fr) minmax(160px, 1.2fr) minmax(110px, 0.8fr);
        grid-template-areas:
          "subject tags tags apply"
          "position alt alt apply";
      }}
      .filter-apply {{ align-self: stretch; }}
    }}
    @media (max-width: 900px) {{
      .toolbar {{
        grid-template-columns: 1fr 1fr;
        grid-template-areas:
          "subject position"
          "tags tags"
          "alt alt"
          "apply apply";
      }}
      .workspace {{ grid-template-columns: 1fr; }}
      .compare-panel {{ position: static; max-height: none; }}
    }}
    @media (max-width: 560px) {{
      header {{ padding: 14px 14px; }}
      main {{ padding: 14px; }}
      .toolbar {{
        grid-template-columns: 1fr;
        grid-template-areas:
          "subject"
          "tags"
          "position"
          "alt"
          "apply";
      }}
      .section-title {{ flex-direction: column; align-items: flex-start; }}
    }}
  </style>
</head>
<body>
  <header>
    <h1>Mito SQL Viewer</h1>
    <div class="path">Database: <span class="mono">{escaped_path}</span></div>
  </header>
  <main>
    <div class="workspace">
    <div class="content">
    <section>
      <div class="section-title">Filters <span id="counts" class="counts"></span></div>
      <div class="toolbar">
        <label class="filter-subject">Subject
          <select id="subject"><option value="">All</option></select>
        </label>
        <label class="filter-tags">Tags
          <div class="tag-picker" id="tag-picker">
            <label><input class="tag-filter" type="checkbox" value="Base"> Base</label>
            <label><input class="tag-filter" type="checkbox" value="EV"> EV</label>
            <label><input class="tag-filter" type="checkbox" value="Raw"> Raw</label>
            <label><input class="tag-filter" type="checkbox" value="PCR"> PCR</label>
          </div>
        </label>
        <label class="filter-position">Position
          <input id="position" inputmode="numeric" placeholder="e.g. 73">
        </label>
        <label class="filter-alt">ALT
          <input id="alt" placeholder="e.g. G">
        </label>
        <button class="filter-apply" id="apply">Refresh</button>
      </div>
    </section>

    <section>
      <div class="section-title">Samples <span class="counts" id="sample-count"></span></div>
      <div class="table-wrap" style="max-height: 220px;">
        <table>
          <thead><tr><th>Order</th><th>Subject</th><th>Population Tags</th><th>Mutations</th><th>Source File</th></tr></thead>
          <tbody id="samples"></tbody>
        </table>
      </div>
    </section>

    <section>
      <div class="section-title">Mutations <span class="counts" id="mutation-count"></span></div>
      <div class="table-wrap">
        <table>
          <thead id="mutation-head">
            <tr>
              <th>Subject</th><th>Population</th><th>Pos</th><th>Ref</th><th>VCF Ref</th>
              <th>ALT</th><th>AF</th><th>Filter</th><th>Metadata</th>
            </tr>
          </thead>
          <tbody id="mutations"></tbody>
        </table>
      </div>
    </section>
    </div>
    <aside class="compare-panel">
      <section style="display: contents;">
        <div class="section-title">Compare <span class="counts" id="compare-count"></span></div>
        <div class="compare-options">
          <details class="compare-summary" open>
            <summary>Comparison</summary>
            <div class="compare-summary-body">
              <div class="compare-summary-expression mono" id="compare-summary-expression">No samples selected.</div>
              <div class="compare-summary-sentence" id="compare-summary-sentence">No comparison is active.</div>
              <div class="compare-summary-hint" id="compare-summary-hint">Choose samples below to build a comparison.</div>
            </div>
          </details>
          <div class="compare-actions">
            <button id="create-derived-sample" type="button">Create sample</button>
            <button class="secondary-button" id="clear-compare" type="button">Clear</button>
          </div>
        </div>
        <div class="compare-list" id="compare-list"></div>
      </section>
    </aside>
    </div>
  </main>
  <div id="metadata-popover" class="metadata-popover" role="dialog" aria-live="polite"></div>
  <div id="af-modal" class="modal-backdrop" role="dialog" aria-modal="true">
    <div class="modal">
      <h2>AF Filter</h2>
      <div class="af-rule-editor">
        <label>Rule
          <select id="af-rule-operator">
            <option value="gt">&gt;</option>
            <option value="gte">&gt;=</option>
            <option value="lt">&lt;</option>
            <option value="lte">&lt;=</option>
            <option value="eq">=</option>
            <option value="neq">!=</option>
          </select>
        </label>
        <label>Value
          <input id="af-rule-value" type="number" min="0" max="1" step="0.0001" placeholder="0.8">
        </label>
        <button id="add-af-rule" type="button">Add</button>
      </div>
      <div class="af-rule-list" id="af-rule-list"></div>
      <div class="modal-actions">
        <button class="secondary-button" id="clear-af-rules" type="button">Clear all</button>
        <button id="close-af-modal" type="button">Done</button>
      </div>
    </div>
  </div>
  <div id="metadata-filter-modal" class="modal-backdrop" role="dialog" aria-modal="true">
    <div class="modal">
      <h2>REF / ALT Filters</h2>
      <div class="metadata-filter-grid">
        <label class="full">Polymorphism
          <select id="meta-polymorphism">
            <option value="">Any</option>
            <option value="1">Polymorphism only</option>
            <option value="0">Non-polymorphism only</option>
          </select>
        </label>
        <label>Repeat base
          <input id="meta-repeat-base" placeholder="A">
        </label>
        <label>Max repeat count
          <input id="meta-repeat-count" type="number" min="0" placeholder="e.g. 3">
        </label>
        <label>2-base repeat
          <input id="meta-repeat-2-bases" placeholder="AT">
        </label>
        <label>Max 2-base count
          <input id="meta-repeat-2-count" type="number" min="0">
        </label>
        <label>3-base repeat
          <input id="meta-repeat-3-bases" placeholder="ATC">
        </label>
        <label>Max 3-base count
          <input id="meta-repeat-3-count" type="number" min="0">
        </label>
        <label class="full">Reference context contains
          <input id="meta-reference-context" placeholder="e.g. ATGC">
        </label>
        <label class="full">Reference context single-base repeat
          <select id="meta-reference-repeat">
            <option value="">Any</option>
            <option value="either">Either side</option>
            <option value="one">One side only</option>
            <option value="before">Before only</option>
            <option value="after">After only</option>
            <option value="both">Both sides</option>
            <option value="none">No side</option>
          </select>
        </label>
      </div>
      <div class="modal-actions">
        <button class="secondary-button" id="clear-metadata-filters" type="button">Clear</button>
        <button id="close-metadata-filter-modal" type="button">Done</button>
      </div>
    </div>
  </div>
  <script>
    const state = {{
      subject: "",
      tags: new Set(),
      position: "",
      alt: "",
      compareSamples: new Set(),
      compareStatuses: new Set(["common", "partial", "unique"]),
      compareSampleStatuses: new Map(),
      sampleOrder: [],
      afRules: [],
      metadataFilters: {{}}
    }};
    const defaultSampleStatuses = ["present"];
    const afOperatorLabels = {{
      gt: ">",
      gte: ">=",
      lt: "<",
      lte: "<=",
      eq: "=",
      neq: "!="
    }};

    const subjectSelect = document.getElementById("subject");
    const tagFilterInputs = [...document.querySelectorAll(".tag-filter")];
    const positionInput = document.getElementById("position");
    const altInput = document.getElementById("alt");
    const compareList = document.getElementById("compare-list");
    const statusFilterInputs = [...document.querySelectorAll(".status-filter")];
    const afModal = document.getElementById("af-modal");
    const metadataFilterModal = document.getElementById("metadata-filter-modal");

    function esc(value) {{
      return String(value ?? "").replace(/[&<>"']/g, char => ({{
        "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"
      }}[char]));
    }}

    async function getJson(url) {{
      const response = await fetch(url);
      if (!response.ok) throw new Error(await response.text());
      return await response.json();
    }}

    async function postJson(url, payload) {{
      const response = await fetch(url, {{
        method: "POST",
        headers: {{"Content-Type": "application/json"}},
        body: JSON.stringify(payload)
      }});
      if (!response.ok) throw new Error(await response.text());
      return await response.json();
    }}

    async function deleteJson(url) {{
      const response = await fetch(url, {{ method: "DELETE" }});
      if (!response.ok) throw new Error(await response.text());
      return await response.json();
    }}

    function queryString(params) {{
      const search = new URLSearchParams();
      for (const [key, value] of Object.entries(params)) {{
        if (Array.isArray(value)) {{
          for (const item of value) {{
            if (item) search.append(key, item);
          }}
        }} else if (value) {{
          search.set(key, value);
        }}
      }}
      return search.toString();
    }}

    function sampleLabel(sample) {{
      return `${{sample.subject_id}} ${{sample.population_key.replaceAll("|", "_")}}`;
    }}

    function samplePopulationText(sample) {{
      return String(sample.population_key || "").replaceAll("|", ", ");
    }}

    function syncSampleOrder(samples) {{
      const ids = samples.map(sample => String(sample.id));
      const known = new Set(ids);
      state.sampleOrder = state.sampleOrder.filter(id => known.has(id));
      for (const id of ids) {{
        if (!state.sampleOrder.includes(id)) {{
          state.sampleOrder.push(id);
        }}
      }}
    }}

    function orderedSamples(samples) {{
      syncSampleOrder(samples);
      const position = new Map(state.sampleOrder.map((id, index) => [id, index]));
      return [...samples].sort((left, right) =>
        position.get(String(left.id)) - position.get(String(right.id))
      );
    }}

    function orderedCompareSampleIds() {{
      return state.sampleOrder.filter(sampleId => state.compareSamples.has(sampleId));
    }}

    function sampleById(samples) {{
      return new Map(samples.map(sample => [String(sample.id), sample]));
    }}

    function listSentence(items) {{
      if (!items.length) return "";
      if (items.length === 1) return items[0];
      if (items.length === 2) return `${{items[0]}} and ${{items[1]}}`;
      return `${{items.slice(0, -1).join(", ")}}, and ${{items[items.length - 1]}}`;
    }}

    function compareSummarySentence(present, unique, absent) {{
      const sentences = [];
      if (present.length) {{
        sentences.push(
          `The samples shown are mutations that are present in ${{listSentence(present)}}.`
        );
      }}
      if (unique.length) {{
        sentences.push(
          `The samples shown also include mutations unique to ${{listSentence(unique)}}.`
        );
      }}
      if (absent.length) {{
        sentences.push(
          `The samples shown exclude mutations present in ${{listSentence(absent)}}.`
        );
      }}
      return sentences.join(" ") || "No comparison is active.";
    }}

    function compareSummary(samples) {{
      const samplesById = sampleById(samples);
      const present = [];
      const unique = [];
      const absent = [];

      for (const sampleId of orderedCompareSampleIds()) {{
        const sample = samplesById.get(sampleId);
        const label = sample ? sampleLabel(sample) : sampleId;
        const statuses = state.compareSampleStatuses.get(sampleId) || new Set(defaultSampleStatuses);
        if (!statuses.size) continue;
        if (statuses.has("not_in")) {{
          absent.push(label);
          continue;
        }}
        if (statuses.has("present")) present.push(label);
        if (statuses.has("unique")) unique.push(label);
      }}

      const parts = [];
      if (present.length === 1) {{
        parts.push(present[0]);
      }} else if (present.length > 1) {{
        parts.push(`AND(${{present.join(", ")}})`);
      }}
      parts.push(...unique.map(label => `UNIQUE(${{label}})`));

      let expression = parts.length ? parts.join(" OR ") : "No active sample constraints.";
      if (absent.length) {{
        expression = `${{parts.length > 1 ? `(${{expression}})` : expression}} AND NOT(${{absent.join(", ")}})`;
      }}

      return {{
        expression,
        sentence: compareSummarySentence(present, unique, absent),
        hint: `${{state.compareSamples.size}} selected · ${{present.length}} present · ${{unique.length}} unique · ${{absent.length}} not in`
      }};
    }}

    function updateCompareSummary(samples) {{
      const summary = compareSummary(samples);
      document.getElementById("compare-summary-expression").textContent = summary.expression;
      document.getElementById("compare-summary-sentence").textContent = summary.sentence;
      document.getElementById("compare-summary-hint").textContent = summary.hint;
    }}

    function moveSample(sampleId, offset) {{
      const index = state.sampleOrder.indexOf(sampleId);
      if (index < 0) return;
      const nextIndex = Math.max(0, Math.min(state.sampleOrder.length - 1, index + offset));
      if (index === nextIndex) return;
      state.sampleOrder.splice(index, 1);
      state.sampleOrder.splice(nextIndex, 0, sampleId);
    }}

    function moveSampleBefore(sampleId, targetId) {{
      if (!sampleId || !targetId || sampleId === targetId) return;
      state.sampleOrder = state.sampleOrder.filter(id => id !== sampleId);
      const targetIndex = state.sampleOrder.indexOf(targetId);
      if (targetIndex < 0) {{
        state.sampleOrder.push(sampleId);
      }} else {{
        state.sampleOrder.splice(targetIndex, 0, sampleId);
      }}
    }}

    function parseMetadata(raw) {{
      try {{
        return JSON.parse(raw || "{{}}");
      }} catch {{
        return null;
      }}
    }}

    function metadataMode(expression) {{
      if (expression.includes(" AND NOT(")) return "and_not";
      if (expression.startsWith("AND(")) return "and";
      if (expression.startsWith("XOR(") || expression.startsWith("EXACTLY_ONE(")) return "xor";
      if (expression.startsWith("SOME_NOT_ALL(")) return "partial";
      if (expression.startsWith("OR(")) return "or";
      return "derived";
    }}

    function metadataModeLabel(mode) {{
      return {{
        and_not: "Present in the included samples; absent from the crossed-out samples",
        and: "AND: shared by every source sample",
        or: "OR: present in any included source sample",
        xor: "XOR: present in exactly one source sample; overlap is crossed out",
        partial: "SOME NOT ALL: shared by a subset, but not every source sample",
        derived: "Derived comparison set"
      }}[mode] || "Derived comparison set";
    }}

    function excludedMetadataSamples(metadata) {{
      const absentSamples = String(metadata.DERIVED_ABSENT_SAMPLES || "")
        .split(";")
        .map(item => item.trim())
        .filter(Boolean);
      if (absentSamples.length) {{
        return absentSamples;
      }}
      return String(metadata.DERIVED_SAMPLE_FILTERS || "")
        .split(";")
        .map(item => item.trim())
        .filter(Boolean)
        .filter(item => item.includes("NOT_IN") || item.includes("NONE"))
        .map(item => item.split(" IS ")[0]);
    }}

    function vennDiagramHtml(metadata) {{
      const expression = String(metadata.DERIVED_EXPRESSION || "");
      const mode = metadataMode(expression);
      const sourceSamples = String(metadata.DERIVED_SOURCE_SAMPLES || "")
        .split(";")
        .map(item => item.trim())
        .filter(Boolean);
      const excluded = new Set(excludedMetadataSamples(metadata));
      const presentSamples = String(metadata.DERIVED_PRESENT_SAMPLES || "")
        .split(";")
        .map(item => item.trim())
        .filter(Boolean);
      const absentSamples = String(metadata.DERIVED_ABSENT_SAMPLES || "")
        .split(";")
        .map(item => item.trim())
        .filter(Boolean);
      const circles = sourceSamples.slice(0, 4).map(sample => `
        <div class="venn-circle ${{excluded.has(sample) ? "excluded" : ""}}">
          ${{esc(sample)}}
        </div>
      `).join("");
      const coreClass = [
        "intersection-core",
        mode === "partial" ? "partial" : "",
        mode === "xor" ? "crossed" : ""
      ].filter(Boolean).join(" ");
      const coreHtml = mode === "or" || mode === "derived" ? "" : `<div class="${{coreClass}}"></div>`;
      const filterText = String(metadata.DERIVED_SAMPLE_FILTERS || "").trim();
      const simpleText = [
        presentSamples.length ? `Present in: ${{presentSamples.join(", ")}}` : "",
        absentSamples.length ? `Absent from: ${{absentSamples.join(", ")}}` : ""
      ].filter(Boolean).join(" / ");

      return `
        <div class="venn-wrap">
          <div class="venn-diagram" aria-hidden="true">
            ${{circles}}
            ${{coreHtml}}
          </div>
          <div class="venn-notes">
            <div><strong>${{esc(metadataModeLabel(mode))}}</strong></div>
            ${{simpleText ? `<div><strong>Simple:</strong> ${{esc(simpleText)}}</div>` : ""}}
            ${{filterText ? `<div><strong>Crossed out:</strong> ${{esc(filterText)}}</div>` : ""}}
            <div><strong>Sources:</strong> ${{sourceSamples.map(esc).join(", ") || "Unknown"}}</div>
            <div><strong>Statuses:</strong> ${{esc(metadata.DERIVED_SET_STATUSES || "")}}</div>
          </div>
        </div>
      `;
    }}

    function showMetadataPopover(cell, event) {{
      const metadata = parseMetadata(cell.dataset.metadata);
      if (!metadata || !metadata.DERIVED_EXPRESSION) {{
        return;
      }}

      const popover = document.getElementById("metadata-popover");
      popover.innerHTML = `
        <h3>${{esc(metadata.DERIVED_LABEL || "Derived comparison")}}</h3>
        ${{vennDiagramHtml(metadata)}}
        <div class="expression mono">${{esc(metadata.DERIVED_EXPRESSION)}}</div>
      `;
      popover.classList.add("open");

      const width = popover.offsetWidth;
      const height = popover.offsetHeight;
      const left = Math.min(event.clientX + 18, window.innerWidth - width - 14);
      const top = Math.min(event.clientY + 18, window.innerHeight - height - 14);
      popover.style.left = `${{Math.max(14, left)}}px`;
      popover.style.top = `${{Math.max(14, top)}}px`;
    }}

    function hideMetadataPopover() {{
      document.getElementById("metadata-popover").classList.remove("open");
    }}

    function setupMetadataPopovers() {{
      document.querySelectorAll(".metadata-cell[data-metadata]").forEach(cell => {{
        cell.addEventListener("mouseenter", event => showMetadataPopover(cell, event));
        cell.addEventListener("mousemove", event => showMetadataPopover(cell, event));
        cell.addEventListener("mouseleave", hideMetadataPopover);
      }});
    }}

    function metadataCellHtml(rawMetadata) {{
      const metadata = parseMetadata(rawMetadata);
      if (metadata && metadata.DERIVED_EXPRESSION) {{
        return `<span class="derived-metadata">${{esc(metadata.DERIVED_EXPRESSION)}}</span>`;
      }}
      return esc(rawMetadata);
    }}

    function metadataList(metadata, key) {{
      return String(metadata?.[key] || "")
        .split(";")
        .map(item => item.trim())
        .filter(Boolean);
    }}

    function derivedRowStatus(metadata) {{
      const expression = String(metadata?.DERIVED_EXPRESSION || "");
      if (metadataList(metadata, "DERIVED_UNIQUE_SAMPLES").length || expression.includes("UNIQUE(")) {{
        return "unique";
      }}
      if (metadataList(metadata, "DERIVED_ABSENT_SAMPLES").length || expression.includes("NOT(")) {{
        return "partial";
      }}
      if (expression.startsWith("AND(")) {{
        return "common";
      }}
      return "partial";
    }}

    function renderListCell(items) {{
      return items.length ? items.map(item => esc(item)).join("<br>") : '<span class="muted">-</span>';
    }}

    function renderDerivedMutationRows(mutations) {{
      return mutations.map(mutation => {{
        const metadata = parseMetadata(mutation.metadata_json) || {{}};
        const status = derivedRowStatus(metadata);
        const expression = metadata.DERIVED_EXPRESSION || "";
        return `
          <tr class="${{esc(status)}}">
            <td>${{esc(status.replaceAll("_", " "))}}</td>
            <td class="mono">${{mutation.pos}}</td>
            <td class="mono">${{esc(mutation.ref)}}</td>
            <td class="mono">${{esc(mutation.alt)}}</td>
            <td class="mono">${{esc(mutation.af)}}</td>
            <td>${{renderListCell(metadataList(metadata, "DERIVED_PRESENT_SAMPLES"))}}</td>
            <td>${{renderListCell(metadataList(metadata, "DERIVED_UNIQUE_SAMPLES"))}}</td>
            <td>${{renderListCell(metadataList(metadata, "DERIVED_ABSENT_SAMPLES"))}}</td>
            <td class="mono metadata-cell" data-metadata="${{esc(mutation.metadata_json)}}">
              ${{expression ? `<span class="derived-metadata">${{esc(expression)}}</span>` : metadataCellHtml(mutation.metadata_json)}}
            </td>
          </tr>
        `;
      }}).join("");
    }}

    function afRuleText(rule) {{
      return `${{afOperatorLabels[rule.operator] || rule.operator}} ${{rule.value}}`;
    }}

    function afHeaderLabel(baseLabel) {{
      if (!state.afRules.length) return baseLabel;
      return `${{baseLabel}} (${{state.afRules.map(afRuleText).join(", ")}})`;
    }}

    function afRuleParams() {{
      return state.afRules.map(rule => `${{rule.operator}}:${{rule.value}}`);
    }}

    function metadataFilterParams() {{
      const filters = state.metadataFilters;
      const params = [];
      if (filters.polymorphism) params.push(`polymorphism:${{filters.polymorphism}}`);
      if (filters.repeat_base) params.push(`repeat_base:${{filters.repeat_base}}`);
      if (filters.repeat_count) params.push(`repeat_count:lte|${{filters.repeat_count}}`);
      if (filters.repeat_2_bases) params.push(`repeat_2_bases:${{filters.repeat_2_bases}}`);
      if (filters.repeat_2_count) params.push(`repeat_2_count:lte|${{filters.repeat_2_count}}`);
      if (filters.repeat_3_bases) params.push(`repeat_3_bases:${{filters.repeat_3_bases}}`);
      if (filters.repeat_3_count) params.push(`repeat_3_count:lte|${{filters.repeat_3_count}}`);
      if (filters.reference_context) params.push(`reference_context:${{filters.reference_context}}`);
      if (filters.reference_repeat) params.push(`reference_repeat:${{filters.reference_repeat}}`);
      return params;
    }}

    function metadataFilterSummary() {{
      const labels = [];
      const filters = state.metadataFilters;
      if (filters.polymorphism === "1") labels.push("poly");
      if (filters.polymorphism === "0") labels.push("non-poly");
      if (filters.repeat_base) labels.push(`repeat ${{filters.repeat_base}}`);
      if (filters.repeat_count) labels.push(`repeat <= ${{filters.repeat_count}}`);
      if (filters.repeat_2_bases) labels.push(`2-repeat ${{filters.repeat_2_bases}}`);
      if (filters.repeat_2_count) labels.push(`2-count <= ${{filters.repeat_2_count}}`);
      if (filters.repeat_3_bases) labels.push(`3-repeat ${{filters.repeat_3_bases}}`);
      if (filters.repeat_3_count) labels.push(`3-count <= ${{filters.repeat_3_count}}`);
      if (filters.reference_context) labels.push(`ctx ${{
        filters.reference_context
      }}`);
      const referenceRepeatLabels = {{
        either: "ctx repeat either",
        one: "ctx repeat one",
        before: "ctx repeat before",
        after: "ctx repeat after",
        both: "ctx repeat both",
        none: "ctx repeat none"
      }};
      if (filters.reference_repeat) labels.push(referenceRepeatLabels[filters.reference_repeat] || "ctx repeat");
      return labels.join(", ");
    }}

    function refAltHeaderLabel(baseLabel) {{
      const summary = metadataFilterSummary();
      return summary ? `${{baseLabel}} (${{summary}})` : baseLabel;
    }}

    function renderAfRules() {{
      const list = document.getElementById("af-rule-list");
      list.innerHTML = state.afRules.length ? state.afRules.map((rule, index) => `
        <div class="af-rule-item">
          <span class="mono">AF ${{esc(afRuleText(rule))}}</span>
          <button class="af-rule-remove" type="button" data-index="${{index}}">Remove</button>
        </div>
      `).join("") : '<div class="muted">No AF rules active.</div>';

      list.querySelectorAll(".af-rule-remove").forEach(button => {{
        button.addEventListener("click", async () => {{
          state.afRules.splice(Number(button.dataset.index), 1);
          renderAfRules();
          await loadMutations();
        }});
      }});
    }}

    function openAfModal() {{
      renderAfRules();
      afModal.classList.add("open");
      document.getElementById("af-rule-value").focus();
    }}

    function closeAfModal() {{
      afModal.classList.remove("open");
    }}

    function metadataFilterValue(id) {{
      return document.getElementById(id).value.trim();
    }}

    function syncMetadataFilterState() {{
      state.metadataFilters = {{
        polymorphism: metadataFilterValue("meta-polymorphism"),
        repeat_base: metadataFilterValue("meta-repeat-base").toUpperCase(),
        repeat_count: metadataFilterValue("meta-repeat-count"),
        repeat_2_bases: metadataFilterValue("meta-repeat-2-bases").toUpperCase(),
        repeat_2_count: metadataFilterValue("meta-repeat-2-count"),
        repeat_3_bases: metadataFilterValue("meta-repeat-3-bases").toUpperCase(),
        repeat_3_count: metadataFilterValue("meta-repeat-3-count"),
        reference_context: metadataFilterValue("meta-reference-context").toUpperCase(),
        reference_repeat: metadataFilterValue("meta-reference-repeat")
      }};
    }}

    function populateMetadataFilterModal() {{
      const filters = state.metadataFilters;
      document.getElementById("meta-polymorphism").value = filters.polymorphism || "";
      document.getElementById("meta-repeat-base").value = filters.repeat_base || "";
      document.getElementById("meta-repeat-count").value = filters.repeat_count || "";
      document.getElementById("meta-repeat-2-bases").value = filters.repeat_2_bases || "";
      document.getElementById("meta-repeat-2-count").value = filters.repeat_2_count || "";
      document.getElementById("meta-repeat-3-bases").value = filters.repeat_3_bases || "";
      document.getElementById("meta-repeat-3-count").value = filters.repeat_3_count || "";
      document.getElementById("meta-reference-context").value = filters.reference_context || "";
      document.getElementById("meta-reference-repeat").value = filters.reference_repeat || "";
    }}

    function openMetadataFilterModal() {{
      populateMetadataFilterModal();
      metadataFilterModal.classList.add("open");
    }}

    function closeMetadataFilterModal() {{
      syncMetadataFilterState();
      metadataFilterModal.classList.remove("open");
      loadMutations();
    }}

    function renderMutationHead(mode) {{
      document.getElementById("mutation-head").innerHTML = mode === "compare" ? `
        <tr>
          <th>Status</th><th>Pos</th><th><button class="header-filter-button metadata-header-filter" type="button">${{esc(refAltHeaderLabel("Ref"))}}</button></th><th><button class="header-filter-button metadata-header-filter" type="button">${{esc(refAltHeaderLabel("ALT"))}}</button></th>
          <th>Present In</th><th><button class="header-filter-button af-header-filter" type="button">${{esc(afHeaderLabel("AFs"))}}</button></th><th>Missing From</th>
        </tr>
      ` : mode === "derived" ? `
        <tr>
          <th>Status</th><th>Pos</th><th><button class="header-filter-button metadata-header-filter" type="button">${{esc(refAltHeaderLabel("Ref"))}}</button></th><th><button class="header-filter-button metadata-header-filter" type="button">${{esc(refAltHeaderLabel("ALT"))}}</button></th><th><button class="header-filter-button af-header-filter" type="button">${{esc(afHeaderLabel("AF"))}}</button></th>
          <th>Present</th><th>Unique</th><th>Absent</th><th>Expression</th>
        </tr>
      ` : `
        <tr>
          <th>Subject</th><th>Population</th><th>Pos</th><th><button class="header-filter-button metadata-header-filter" type="button">${{esc(refAltHeaderLabel("Ref"))}}</button></th><th>VCF Ref</th>
          <th><button class="header-filter-button metadata-header-filter" type="button">${{esc(refAltHeaderLabel("ALT"))}}</button></th><th><button class="header-filter-button af-header-filter" type="button">${{esc(afHeaderLabel("AF"))}}</button></th><th>Filter</th><th>Metadata</th>
        </tr>
      `;
      document.querySelectorAll(".af-header-filter").forEach(button => {{
        button.addEventListener("click", openAfModal);
      }});
      document.querySelectorAll(".metadata-header-filter").forEach(button => {{
        button.addEventListener("click", openMetadataFilterModal);
      }});
    }}

    function comparisonParams(limit = 2000) {{
      const compareSamples = orderedCompareSampleIds();
      const selectedStatuses = [...state.compareStatuses];
      const sampleStatusParams = [];
      for (const sampleId of compareSamples) {{
        const statuses = state.compareSampleStatuses.get(sampleId)
          || new Set(defaultSampleStatuses);
        if (statuses.size === 0) {{
          sampleStatusParams.push(`${{sampleId}}:__none__`);
        }} else {{
          for (const status of statuses) {{
            sampleStatusParams.push(`${{sampleId}}:${{status}}`);
          }}
        }}
      }}

      return {{
        compare_sample_id: compareSamples,
        status: selectedStatuses.length ? selectedStatuses : ["__none__"],
        sample_status: sampleStatusParams,
        position: state.position,
        alt: state.alt,
        af_rule: afRuleParams(),
        metadata_filter: metadataFilterParams(),
        limit
      }};
    }}

    async function loadSubjects() {{
      const subjects = await getJson("/api/subjects");
      subjectSelect.innerHTML = '<option value="">All</option>' + subjects.map(subject =>
        `<option value="${{esc(subject.subject_id)}}">${{esc(subject.subject_id)}} (${{subject.sample_count}})</option>`
      ).join("");
    }}

    async function loadCounts() {{
      const counts = await getJson("/api/counts");
      document.getElementById("counts").textContent =
        `${{counts.subjects}} subjects, ${{counts.samples}} samples, ${{counts.mutations}} mutations`;
    }}

    async function loadSamples() {{
      const params = {{}};
      if (state.subject) params.subject = state.subject;
      params.tag = [...state.tags];
      const samples = orderedSamples(await getJson(`/api/samples?${{queryString(params)}}`));
      document.getElementById("sample-count").textContent = `${{samples.length}} shown`;
      document.getElementById("samples").innerHTML = samples.map(sample => `
        <tr data-sample="${{sample.id}}" class="${{state.compareSamples.has(String(sample.id)) ? "selected-sample" : ""}}">
          <td class="order-cell">
            <button class="drag-handle" type="button" draggable="true" data-sample="${{esc(sample.id)}}">Drag</button>
            <button class="sample-order-button" type="button" data-sample="${{esc(sample.id)}}" data-offset="-1">Up</button>
            <button class="sample-order-button" type="button" data-sample="${{esc(sample.id)}}" data-offset="1">Down</button>
          </td>
          <td class="mono">${{esc(sample.subject_id)}}</td>
          <td>${{esc(sample.population_key.replaceAll("|", ", "))}}</td>
          <td>${{sample.mutation_count}}</td>
          <td class="mono muted">
            ${{esc(sample.source_file)}}
            ${{sample.is_derived ? `<button class="delete-derived-sample" type="button" data-sample="${{esc(sample.id)}}">Delete</button>` : ""}}
          </td>
        </tr>
      `).join("");
      document.querySelectorAll("tr[data-sample]").forEach(row => {{
        row.addEventListener("dragover", event => {{
          if (!Array.from(event.dataTransfer.types).includes("text/plain")) return;
          event.preventDefault();
          row.classList.add("reorder-target");
        }});
        row.addEventListener("dragleave", () => {{
          row.classList.remove("reorder-target");
        }});
        row.addEventListener("drop", async event => {{
          event.preventDefault();
          row.classList.remove("reorder-target");
          moveSampleBefore(event.dataTransfer.getData("text/plain"), row.dataset.sample);
          await loadSamples();
          await loadMutations();
        }});
        row.addEventListener("click", () => {{
          const sampleId = row.dataset.sample;
          if (state.compareSamples.has(sampleId)) {{
            state.compareSamples.delete(sampleId);
            state.compareSampleStatuses.delete(sampleId);
          }} else {{
            state.compareSamples.add(sampleId);
            state.compareSampleStatuses.set(
              sampleId,
              new Set(defaultSampleStatuses)
            );
          }}
          renderCompareList(samples);
          loadMutations();
        }});
      }});
      document.querySelectorAll(".drag-handle").forEach(button => {{
        button.addEventListener("click", event => event.stopPropagation());
        button.addEventListener("dragstart", event => {{
          event.stopPropagation();
          event.dataTransfer.setData("text/plain", button.dataset.sample);
          event.dataTransfer.effectAllowed = "move";
        }});
      }});
      document.querySelectorAll(".sample-order-button").forEach(button => {{
        button.addEventListener("click", async event => {{
          event.stopPropagation();
          moveSample(button.dataset.sample, Number(button.dataset.offset));
          await loadSamples();
          await loadMutations();
        }});
      }});
      document.querySelectorAll(".delete-derived-sample").forEach(button => {{
        button.addEventListener("click", async event => {{
          event.stopPropagation();
          const sampleId = button.dataset.sample;
          await deleteJson(`/api/derived-samples/${{encodeURIComponent(sampleId)}}`);
          state.compareSamples.delete(sampleId);
          state.compareSampleStatuses.delete(sampleId);
          state.sampleOrder = state.sampleOrder.filter(id => id !== sampleId);
          await loadSamples();
          await loadMutations();
        }});
      }});
      renderCompareList(samples);
    }}

    function renderCompareList(samples) {{
      samples = orderedSamples(samples);
      state.compareSamples = new Set(
        [...state.compareSamples].filter(sampleId =>
          samples.some(sample => String(sample.id) === sampleId)
        )
      );
      for (const sampleId of [...state.compareSampleStatuses.keys()]) {{
        if (!state.compareSamples.has(sampleId)) {{
          state.compareSampleStatuses.delete(sampleId);
        }}
      }}
      const available = samples;
      document.getElementById("compare-count").textContent = `${{state.compareSamples.size}} selected`;
      updateCompareSummary(samples);

      compareList.innerHTML = available.length ? available.map(sample => {{
        const id = String(sample.id);
        const checked = state.compareSamples.has(id) ? "checked" : "";
        const sampleStatuses = state.compareSampleStatuses.get(id)
          || new Set(defaultSampleStatuses);
        const showStatuses = state.compareSamples.has(id);
        const statusControls = showStatuses ? `
          <div class="sample-statuses" data-sample-statuses="${{id}}">
            <label><input class="sample-status-filter" data-sample="${{id}}" type="checkbox" value="present" ${{sampleStatuses.has("present") ? "checked" : ""}}> Present</label>
            <label><input class="sample-status-filter" data-sample="${{id}}" type="checkbox" value="unique" ${{sampleStatuses.has("unique") ? "checked" : ""}}> Unique</label>
            <label><input class="sample-status-filter exclusive" data-sample="${{id}}" type="checkbox" value="not_in" ${{sampleStatuses.has("not_in") ? "checked" : ""}}> Not In</label>
            <label><input class="sample-status-filter exclusive" data-sample="${{id}}" type="checkbox" value="none" ${{sampleStatuses.size === 0 ? "checked" : ""}}> None</label>
          </div>
        ` : "";
        return `
          <div class="compare-item">
            <div class="compare-item-header">
              <input class="compare-sample-toggle" type="checkbox" value="${{id}}" ${{checked}} aria-label="Compare ${{esc(sampleLabel(sample))}}">
              <div class="compare-item-title">
                <strong>${{esc(sample.subject_id)}} · ${{esc(samplePopulationText(sample))}}</strong>
                <div class="compare-meta">
                  <span>${{sample.mutation_count}} mutations${{sample.is_derived ? " · derived sample" : ""}}</span>
                  <span class="mono">${{esc(sample.source_file)}}</span>
                </div>
              </div>
            </div>
            ${{statusControls}}
          </div>
        `;
      }}).join("") : '<div class="muted">No samples match the current filters.</div>';

      compareList.querySelectorAll(".compare-sample-toggle").forEach(input => {{
        input.addEventListener("change", () => {{
          if (input.checked) {{
            state.compareSamples.add(input.value);
            if (!state.compareSampleStatuses.has(input.value)) {{
              state.compareSampleStatuses.set(
                input.value,
                new Set(defaultSampleStatuses)
              );
            }}
          }} else {{
            state.compareSamples.delete(input.value);
            state.compareSampleStatuses.delete(input.value);
          }}
          document.getElementById("compare-count").textContent = `${{state.compareSamples.size}} selected`;
          renderCompareList(samples);
          loadMutations();
        }});
      }});
      compareList.querySelectorAll(".sample-status-filter").forEach(input => {{
        input.addEventListener("change", () => {{
          const sampleId = input.dataset.sample;
          const statuses = new Set(state.compareSampleStatuses.get(sampleId) || defaultSampleStatuses);
          if (input.value === "none" && input.checked) {{
            state.compareSampleStatuses.set(sampleId, new Set());
          }} else if (input.value === "not_in" && input.checked) {{
            state.compareSampleStatuses.set(sampleId, new Set(["not_in"]));
          }} else if (input.value === "not_in") {{
            state.compareSampleStatuses.set(sampleId, new Set());
          }} else {{
            statuses.delete("not_in");
            if (input.checked) {{
              statuses.add(input.value);
            }} else {{
              statuses.delete(input.value);
            }}
            state.compareSampleStatuses.set(sampleId, statuses);
          }}
          renderCompareList(samples);
          loadMutations();
        }});
      }});
    }}

    async function loadMutations() {{
      const compareSamples = [...state.compareSamples];
      const orderedCompareSamples = orderedCompareSampleIds();
      if (orderedCompareSamples.length >= 2) {{
        hideMetadataPopover();
        renderMutationHead("compare");
        const params = comparisonParams(2000);
        const rows = await getJson(`/api/compare?${{queryString(params)}}`);
        document.getElementById("mutation-count").textContent = `${{rows.length}} comparison rows`;
        const groupNumbers = new Map();
        let nextGroupNumber = 0;
        document.getElementById("mutations").innerHTML = rows.map(row => {{
          const present = row.present || [];
          const missing = row.missing || [];
          if (!groupNumbers.has(row.group_key)) {{
            groupNumbers.set(row.group_key, nextGroupNumber++);
          }}
          const groupNumber = groupNumbers.get(row.group_key);
          const groupClasses = [
            row.status,
            row.group_size > 1 ? "grouped" : "single",
            row.group_start ? "group-start" : "group-middle",
            row.group_end ? "group-end" : "",
            groupNumber % 2 === 0 ? "group-even" : "group-odd"
          ].filter(Boolean).join(" ");
          return `
            <tr class="${{esc(groupClasses)}}">
              <td>${{esc(row.status.replaceAll("_", " "))}}</td>
              <td class="mono">${{row.pos}}</td>
              <td class="mono">${{esc(row.ref)}}</td>
              <td class="mono">${{esc(row.alt)}}</td>
              <td>${{present.map(item => esc(item.label)).join("<br>")}}</td>
              <td class="mono">${{present.map(item => esc(item.af_text || "")).join("<br>")}}</td>
              <td>${{missing.map(item => esc(item.label)).join("<br>")}}</td>
            </tr>
          `;
        }}).join("");
        return;
      }}

      hideMetadataPopover();
      const params = {{
        sample_id: orderedCompareSamples.length === 1 ? orderedCompareSamples[0] : "",
        position: state.position,
        alt: state.alt,
        af_rule: afRuleParams(),
        metadata_filter: metadataFilterParams(),
        limit: 1000
      }};
      const mutations = await getJson(`/api/mutations?${{queryString(params)}}`);
      document.getElementById("mutation-count").textContent = `${{mutations.length}} shown`;
      const showingDerivedSample = orderedCompareSamples.length === 1
        && orderedCompareSamples[0].startsWith("derived:");
      if (showingDerivedSample) {{
        renderMutationHead("derived");
        document.getElementById("mutations").innerHTML = renderDerivedMutationRows(mutations);
        setupMetadataPopovers();
        return;
      }}

      renderMutationHead("normal");
      document.getElementById("mutations").innerHTML = mutations.map(mutation => `
        <tr>
          <td class="mono">${{esc(mutation.subject_id)}}</td>
          <td>${{esc(mutation.population_key.replaceAll("|", ", "))}}</td>
          <td class="mono">${{mutation.pos}}</td>
          <td class="mono">${{esc(mutation.ref)}}</td>
          <td class="mono">${{esc(mutation.vcf_ref)}}</td>
          <td class="mono">${{esc(mutation.alt)}}</td>
          <td class="mono">${{esc(mutation.af)}}</td>
          <td>${{esc(mutation.filter)}}</td>
          <td class="mono metadata-cell" data-metadata="${{esc(mutation.metadata_json)}}">
            ${{metadataCellHtml(mutation.metadata_json)}}
          </td>
        </tr>
      `).join("");
      setupMetadataPopovers();
    }}

    async function createDerivedSampleFromComparison() {{
      syncState();
      const compareSamples = [...state.compareSamples];
      if (compareSamples.length < 2) {{
        return;
      }}

      const label = window.prompt("Name this comparison sample:", `Comparison ${{Date.now()}}`);
      if (label === null) {{
        return;
      }}

      const sample = await postJson("/api/derived-samples", {{
        ...comparisonParams(2000),
        label: label.trim()
      }});
      const sampleId = String(sample.id);
      if (!state.sampleOrder.includes(sampleId)) {{
        state.sampleOrder.push(sampleId);
      }}
      state.compareSamples.add(sampleId);
      state.compareSampleStatuses.set(sampleId, new Set(defaultSampleStatuses));
      await loadSamples();
      await loadMutations();
    }}

    async function clearComparison() {{
      state.compareSamples.clear();
      state.compareSampleStatuses.clear();
      await loadSamples();
      await loadMutations();
    }}

    function syncState() {{
      state.subject = subjectSelect.value;
      state.tags = new Set(
        tagFilterInputs
          .filter(input => input.checked)
          .map(input => input.value)
      );
      state.position = positionInput.value.trim();
      state.alt = altInput.value.trim();
      if (statusFilterInputs.length) {{
        state.compareStatuses = new Set(
          statusFilterInputs
            .filter(input => input.checked)
            .map(input => input.value)
        );
      }}
    }}

    async function refreshFromFilters() {{
      syncState();
      await loadSamples();
      await loadMutations();
    }}

    function debounce(fn, delay = 250) {{
      let timeoutId;
      return (...args) => {{
        window.clearTimeout(timeoutId);
        timeoutId = window.setTimeout(() => fn(...args), delay);
      }};
    }}

    const refreshFromTextFilters = debounce(refreshFromFilters, 250);

    document.getElementById("apply").addEventListener("click", refreshFromFilters);

    subjectSelect.addEventListener("change", refreshFromFilters);
    tagFilterInputs.forEach(input => {{
      input.addEventListener("change", refreshFromFilters);
    }});
    [positionInput, altInput].forEach(input => {{
      input.addEventListener("input", refreshFromTextFilters);
      input.addEventListener("change", refreshFromFilters);
    }});
    statusFilterInputs.forEach(input => {{
      input.addEventListener("change", async () => {{
        syncState();
        await loadMutations();
      }});
    }});
    document.getElementById("create-derived-sample").addEventListener("click", async () => {{
      await createDerivedSampleFromComparison();
    }});
    document.getElementById("clear-compare").addEventListener("click", async () => {{
      await clearComparison();
    }});
    document.getElementById("add-af-rule").addEventListener("click", async () => {{
      const operator = document.getElementById("af-rule-operator").value;
      const valueInput = document.getElementById("af-rule-value");
      const value = Number(valueInput.value);
      if (!Number.isFinite(value) || value < 0 || value > 1) {{
        return;
      }}
      state.afRules.push({{ operator, value }});
      valueInput.value = "";
      renderAfRules();
      await loadMutations();
    }});
    document.getElementById("clear-af-rules").addEventListener("click", async () => {{
      state.afRules = [];
      renderAfRules();
      await loadMutations();
    }});
    document.getElementById("close-af-modal").addEventListener("click", closeAfModal);
    afModal.addEventListener("click", event => {{
      if (event.target === afModal) {{
        closeAfModal();
      }}
    }});
    document.getElementById("close-metadata-filter-modal").addEventListener("click", closeMetadataFilterModal);
    document.getElementById("clear-metadata-filters").addEventListener("click", async () => {{
      state.metadataFilters = {{}};
      populateMetadataFilterModal();
      metadataFilterModal.classList.remove("open");
      await loadMutations();
    }});
    metadataFilterModal.addEventListener("click", event => {{
      if (event.target === metadataFilterModal) {{
        closeMetadataFilterModal();
      }}
    }});

    loadSubjects().then(loadCounts).then(loadSamples).then(loadMutations);
  </script>
</body>
</html>"""


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
    parser.add_argument("--host", default="0.0.0.0")
    parser.add_argument("--port", type=int, default=18000)
    args = parser.parse_args()
    run_server(args.db, args.host, args.port)


if __name__ == "__main__":
    main()
