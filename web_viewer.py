#!/usr/bin/env python3
"""Small web viewer for the SQLite mutation database."""

import argparse
import email.policy
import html
import json
from email.parser import BytesParser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, unquote, urlparse

import vcf_parser
from mito_viewer.domain import (
    AFRule,
    AlleleKey,
    DerivedSample,
    MetadataFilter,
    MutationFilters,
    SampleAlleleCall,
    SampleAlleleSet,
)
from mito_viewer.catalog import CatalogRepository, CatalogViewerService
from mito_viewer.domain.filters import AF_OPERATORS
from mito_viewer.repositories import (
    DATABASE_EXTENSIONS,
    NO_TAGS_FILTER,  # noqa: F401 - compatibility export
    AnnotationRepository,
    StudyRepository,
    discover_study_databases,
    inspect_study_database,
)
from mito_viewer.repositories.schema import read_only_connection
from mito_viewer.reference_features import mitochondrial_gene_locus


# Compatibility name for callers that imported web_viewer.MutationAllele.
MutationAllele = AlleleKey


DEFAULT_DATABASE_DIR = Path(".")
DEFAULT_ANNOTATION_DATABASE_PATH = (
    Path(__file__).resolve().parent / "mutation_annotations.sqlite"
)
DEFAULT_CATALOG_DATABASE_PATH = (
    Path(__file__).resolve().parent / "analysis_catalog.sqlite"
)
VIEWER_TEMPLATE_PATH = Path(__file__).resolve().parent / "templates" / "viewer.html"
ROADMAP_TEMPLATE_PATH = Path(__file__).resolve().parent / "templates" / "roadmap.html"
ROADMAP_DATA_PATH = Path(__file__).resolve().parent / "roadmap.json"
DEFAULT_COMPARE_STATUSES = {"common", "partial", "unique"}
DEFAULT_SAMPLE_COMPARE_STATUSES = {"present"}
DERIVED_SAMPLE_PREFIX = "derived:"
ROADMAP_STATUSES = {"now", "next", "later", "done"}
ROADMAP_PRIORITIES = {"None", "Low", "Medium", "High"}
REFERENCE_REPEAT_BASES = ("A", "C", "G", "T", "N")


def connect_database(db_path):
    """Compatibility helper returning a read-only study connection."""
    inspect_study_database(db_path).require_valid()
    return read_only_connection(db_path)


def fetch_cached_annotations(annotation_db_path, position, ref, alt):
    """Compatibility wrapper around the annotation repository."""
    return AnnotationRepository(annotation_db_path).fetch_cached(
        position,
        ref,
        alt,
    )


def fetch_variant_annotations(annotation_db_path, position, ref, alt):
    """Combine cached provider data with the authoritative local gene locus."""
    position = int(position)
    ref = str(ref).strip().upper()
    alt = str(alt).strip().upper()
    annotation = fetch_cached_annotations(
        annotation_db_path,
        position,
        ref,
        alt,
    )
    if annotation is None:
        annotation = {
            "variant": {
                "pos": position,
                "ref": ref,
                "alt": alt,
            },
            "cache": {
                "database": Path(annotation_db_path).name,
                "hit": False,
            },
        }
    else:
        annotation.setdefault("cache", {})["hit"] = True
    annotation["locus"] = mitochondrial_gene_locus(position, ref)
    return annotation


def fetch_annotation_vocabulary(annotation_db_path):
    """Compatibility wrapper around the annotation repository."""
    return AnnotationRepository(annotation_db_path).vocabulary()


def fetch_mutation_samples(connection, position, ref, alt):
    """Return every study sample containing an exact VCF allele."""
    return StudyRepository(connection).mutation_samples(position, ref, alt)


def discover_databases(database_dir):
    """Compatibility wrapper returning schema-compatible study databases."""
    return discover_study_databases(database_dir)


def database_id_for_path(db_path):
    return Path(db_path).name


def database_options(database_dir, default_db_id, selected_db_id=None):
    databases = discover_databases(database_dir)
    selected_db_id = selected_db_id or default_db_id
    return [
        {
            "id": database_id,
            "label": database_id,
            "path": str(path),
            "selected": database_id == selected_db_id,
            "schema_version": inspect_study_database(path).version_label,
        }
        for database_id, path in databases.items()
    ]


def is_probable_sqlite_database(content):
    return content.startswith(b"SQLite format 3\x00")


def parse_multipart_form(headers, body):
    content_type = headers.get("Content-Type", "")
    message_bytes = (
        f"Content-Type: {content_type}\r\nMIME-Version: 1.0\r\n\r\n".encode("utf-8")
        + body
    )
    message = BytesParser(policy=email.policy.default).parsebytes(message_bytes)
    fields = {}
    files = {}
    if not message.is_multipart():
        return fields, files

    for part in message.iter_parts():
        name = part.get_param("name", header="content-disposition")
        if not name:
            continue
        filename = part.get_filename()
        payload = part.get_payload(decode=True) or b""
        if filename:
            files[name] = {
                "filename": filename,
                "content": payload,
            }
        else:
            fields[name] = payload.decode(part.get_content_charset() or "utf-8")
    return fields, files


def save_uploaded_database(database_dir, filename, content, replace=False):
    safe_name = Path(filename or "").name
    if not safe_name:
        raise ValueError("Choose a database file to upload.")
    if Path(safe_name).suffix.lower() not in DATABASE_EXTENSIONS:
        extensions = ", ".join(sorted(DATABASE_EXTENSIONS))
        raise ValueError(f"Database uploads must use one of these extensions: {extensions}.")
    if not is_probable_sqlite_database(content):
        raise ValueError("Uploaded file does not look like a SQLite database.")

    database_dir = Path(database_dir).resolve()
    target_path = (database_dir / safe_name).resolve()
    if target_path.parent != database_dir:
        raise ValueError("Invalid database filename.")
    if target_path.exists() and not replace:
        raise FileExistsError(f"Database already exists: {safe_name}")

    temp_path = target_path.with_name(f".{target_path.name}.uploading")
    temp_path.write_bytes(content)
    try:
        inspect_study_database(temp_path).require_valid()
    except Exception:
        temp_path.unlink(missing_ok=True)
        raise
    temp_path.replace(target_path)
    return target_path


def is_derived_sample_id(sample_id):
    return str(sample_id).startswith(DERIVED_SAMPLE_PREFIX)


def parse_af_rules(values):
    """Parse valid AF query values such as ``gt:0.8``."""
    rules = []
    for value in values or []:
        try:
            rules.append(AFRule.parse(value))
        except ValueError:
            continue
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
    """Parse valid metadata-filter query values."""
    filters = []
    for value in values or []:
        try:
            filters.append(MetadataFilter.parse(value))
        except ValueError:
            continue
    return filters


def single_base_repeat_seen(sequence):
    """Return True when a reference context contains an adjacent single-base run."""
    sequence = str(sequence or "").upper()
    return any(base * 2 in sequence for base in REFERENCE_REPEAT_BASES)


def metadata_filters_match(metadata, metadata_filters, alt=""):
    if not metadata_filters:
        return True
    for field, raw_value in metadata_filters:
        if field == "polymorphism":
            if str(metadata.get("POLYMORPHISM", "")) != raw_value:
                return False
        elif field == "reference_contains_alt":
            alt_text = str(alt).upper()
            before = str(metadata.get("REFERENCE_6_BEFORE", "")).upper()
            after = str(metadata.get("REFERENCE_6_AFTER", "")).upper()
            contains_alt = bool(alt_text) and (alt_text in before or alt_text in after)
            if raw_value == "contains" and not contains_alt:
                return False
            if raw_value == "not_contains" and contains_alt:
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
    return StudyRepository(connection).subjects()


def fetch_population_tags(connection):
    return StudyRepository(connection).population_tags()


def fetch_samples(connection, subject_id=None, tags=None, derived_samples=None):
    samples = StudyRepository(connection).samples(
        subject_id=subject_id,
        tags=tags,
    )

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
        if not metadata_filters_match(metadata, metadata_filters or [], allele.alt):
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

    return StudyRepository(connection).mutation_rows(
        sample_id=sample_id,
        position=position,
        alt=alt,
        af_rules=af_rules or (),
        metadata_filters=metadata_filters or (),
        limit=limit,
    )


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
            call.allele.alt,
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

    calls.extend(
        StudyRepository(connection).allele_calls(
            real_sample_ids,
            position=position,
            alt=alt,
            af_rules=af_rules or (),
            metadata_filters=metadata_filters or (),
        )
    )
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
    labels.update(StudyRepository(connection).sample_labels(real_sample_ids))
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
    return StudyRepository(connection).counts()


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


def normalize_roadmap_card(card):
    if not isinstance(card, dict):
        raise ValueError("Each roadmap card must be an object.")

    card_id = str(card.get("id") or "").strip()
    title = str(card.get("title") or "").strip()
    notes = str(card.get("notes") or "").strip()
    system_notes = str(card.get("system_notes") or "").strip()
    status = str(card.get("status") or "next")
    priority = str(card.get("priority") or "None")

    if not card_id:
        raise ValueError("Roadmap cards must include an id.")
    if not title:
        raise ValueError("Roadmap cards must include a title.")
    if status not in ROADMAP_STATUSES:
        raise ValueError(f"Invalid roadmap status: {status}")
    if priority not in ROADMAP_PRIORITIES:
        raise ValueError(f"Invalid roadmap priority: {priority}")

    return {
        "id": card_id,
        "title": title,
        "notes": notes,
        "system_notes": system_notes,
        "status": status,
        "priority": priority,
    }


def normalize_roadmap_payload(payload):
    if isinstance(payload, list):
        cards = payload
    elif isinstance(payload, dict):
        cards = payload.get("cards", [])
    else:
        raise ValueError("Roadmap payload must be an object or list.")

    if not isinstance(cards, list):
        raise ValueError("Roadmap cards must be a list.")

    seen_ids = set()
    normalized_cards = []
    for card in cards:
        normalized_card = normalize_roadmap_card(card)
        if normalized_card["id"] in seen_ids:
            raise ValueError(f"Duplicate roadmap card id: {normalized_card['id']}")
        seen_ids.add(normalized_card["id"])
        normalized_cards.append(normalized_card)

    return {"cards": normalized_cards}


def read_roadmap():
    if not ROADMAP_DATA_PATH.exists():
        return {"cards": []}
    payload = json.loads(ROADMAP_DATA_PATH.read_text(encoding="utf-8"))
    return normalize_roadmap_payload(payload)


def write_roadmap(payload):
    normalized_payload = normalize_roadmap_payload(payload)
    temp_path = ROADMAP_DATA_PATH.with_name(f".{ROADMAP_DATA_PATH.name}.saving")
    temp_path.write_text(
        json.dumps(normalized_payload, indent=2) + "\n",
        encoding="utf-8",
    )
    temp_path.replace(ROADMAP_DATA_PATH)
    return normalized_payload


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


def roadmap_html():
    return ROADMAP_TEMPLATE_PATH.read_text(encoding="utf-8")


class ViewerHandler(BaseHTTPRequestHandler):
    default_db_id = ""
    database_dir = None
    derived_samples_by_db = {}
    next_derived_sample_ids = {}
    annotation_db_path = DEFAULT_ANNOTATION_DATABASE_PATH
    catalog_db_path = DEFAULT_CATALOG_DATABASE_PATH

    def log_message(self, fmt, *args):
        print(f"{self.address_string()} - {fmt % args}")

    @classmethod
    def default_database_id(cls):
        return cls.default_db_id

    @classmethod
    def available_databases(cls):
        return discover_databases(cls.database_dir)

    @classmethod
    def resolve_database(cls, db_id=None):
        db_id = db_id or cls.default_database_id()
        databases = cls.available_databases()
        if db_id not in databases:
            available = ", ".join(databases) or "none"
            raise ValueError(f"Unknown database: {db_id}. Available databases: {available}.")
        return db_id, databases[db_id]

    @staticmethod
    def database_id_from_query(query):
        return query.get("db", [""])[0] or None

    @classmethod
    def derived_samples(cls, db_id):
        return cls.derived_samples_by_db.setdefault(db_id, {})

    @staticmethod
    def catalog_context(query=None, payload=None):
        query = query or {}
        payload = payload or {}
        perspective_id = (
            payload.get("perspective_id")
            or query.get("perspective_id", [""])[0]
        )
        dataset_id = (
            payload.get("dataset_id")
            or query.get("dataset_id", [""])[0]
        )
        return str(perspective_id or ""), str(dataset_id or "")

    def request_derived_samples(self, db_id, query):
        perspective_id, dataset_id = self.catalog_context(query)
        if not perspective_id or not dataset_id:
            return dict(self.derived_samples(db_id))
        with CatalogRepository(self.catalog_db_path) as catalog:
            return CatalogViewerService(catalog).durable_derived_samples(
                perspective_id,
                dataset_id,
            )

    def open_database(self, db_path):
        if not db_path.exists():
            raise FileNotFoundError(
                f"Database not found: {db_path}. "
                "Create it with main.make_sql_database() first."
            )
        return connect_database(db_path)

    def read_json_body(self):
        content_length = int(self.headers.get("Content-Length", "0"))
        if content_length <= 0:
            return {}
        body = self.rfile.read(content_length).decode("utf-8")
        return json.loads(body or "{}")

    def read_body(self):
        content_length = int(self.headers.get("Content-Length", "0"))
        if content_length <= 0:
            return b""
        return self.rfile.read(content_length)

    @classmethod
    def allocate_derived_sample_id(cls, db_id):
        next_id = cls.next_derived_sample_ids.get(db_id, 1)
        derived_id = f"{DERIVED_SAMPLE_PREFIX}{next_id}"
        cls.next_derived_sample_ids[db_id] = next_id + 1
        return derived_id

    def do_GET(self):
        parsed = urlparse(self.path)
        query = parse_qs(parsed.query)

        try:
            if parsed.path == "/roadmap":
                html_response(self, roadmap_html())
                return
            if parsed.path == "/api/roadmap":
                json_response(self, read_roadmap())
                return
            if parsed.path == "/api/annotations":
                annotation = fetch_variant_annotations(
                    self.annotation_db_path,
                    query.get("position", [""])[0],
                    query.get("ref", [""])[0],
                    query.get("alt", [""])[0],
                )
                json_response(self, annotation)
                return
            if parsed.path == "/api/annotation-vocabulary":
                json_response(
                    self,
                    fetch_annotation_vocabulary(self.annotation_db_path),
                )
                return
            if parsed.path == "/api/catalog/perspectives":
                with CatalogRepository(self.catalog_db_path) as catalog:
                    json_response(
                        self,
                        CatalogViewerService(catalog).list_perspectives(),
                    )
                return
            if parsed.path == "/api/catalog/datasets":
                perspective_id, _dataset_id = self.catalog_context(query)
                if not perspective_id:
                    raise ValueError("Select a Study Perspective.")
                with CatalogRepository(self.catalog_db_path) as catalog:
                    json_response(
                        self,
                        CatalogViewerService(catalog).list_datasets(
                            perspective_id
                        ),
                    )
                return
            if parsed.path == "/api/catalog/workspace":
                perspective_id, dataset_id = self.catalog_context(query)
                if not perspective_id or not dataset_id:
                    raise ValueError(
                        "Select a Study Perspective and Dataset."
                    )
                db_id, db_path = self.resolve_database(
                    self.database_id_from_query(query)
                )
                with CatalogRepository(self.catalog_db_path) as catalog:
                    json_response(
                        self,
                        CatalogViewerService(catalog).workspace(
                            perspective_id,
                            dataset_id,
                            database_id=db_id,
                            database_path=db_path,
                        ),
                    )
                return
            if parsed.path == "/api/mutation-samples":
                study_db_id = query.get("study_db", [""])[0]
                _study_db_id, study_db_path = self.resolve_database(study_db_id)
                with self.open_database(study_db_path) as connection:
                    samples = fetch_mutation_samples(
                        connection,
                        query.get("position", [""])[0],
                        query.get("ref", [""])[0],
                        query.get("alt", [""])[0],
                    )
                json_response(
                    self,
                    {
                        "database": study_db_id,
                        "samples": samples,
                    },
                )
                return

            db_id, db_path = self.resolve_database(self.database_id_from_query(query))
            derived_samples = self.request_derived_samples(
                db_id,
                query,
            )

            if parsed.path == "/":
                html_response(self, page_html(db_path))
            elif parsed.path == "/api/databases":
                json_response(
                    self,
                    database_options(
                        self.database_dir,
                        self.default_database_id(),
                        selected_db_id=db_id,
                    ),
                )
            elif parsed.path == "/api/counts":
                with self.open_database(db_path) as connection:
                    json_response(self, database_counts(connection))
            elif parsed.path == "/api/subjects":
                with self.open_database(db_path) as connection:
                    json_response(self, fetch_subjects(connection))
            elif parsed.path == "/api/tags":
                with self.open_database(db_path) as connection:
                    json_response(self, fetch_population_tags(connection))
            elif parsed.path == "/api/samples":
                with self.open_database(db_path) as connection:
                    json_response(
                        self,
                        fetch_samples(
                            connection,
                            subject_id=query.get("subject", [""])[0],
                            tags=parse_tags(query),
                            derived_samples=derived_samples,
                        ),
                    )
            elif parsed.path == "/api/mutations":
                with self.open_database(db_path) as connection:
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
                            derived_samples=derived_samples,
                        ),
                    )
            elif parsed.path == "/api/compare":
                with self.open_database(db_path) as connection:
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
                            derived_samples=derived_samples,
                        ),
                    )
            else:
                json_response(self, {"error": "Not found"}, status=404)
        except ValueError as exc:
            json_response(self, {"error": str(exc)}, status=400)
        except Exception as exc:
            json_response(self, {"error": str(exc)}, status=500)

    def do_POST(self):
        parsed = urlparse(self.path)
        query = parse_qs(parsed.query)

        try:
            if parsed.path == "/api/roadmap":
                json_response(self, write_roadmap(self.read_json_body()))
            elif parsed.path == "/api/catalog/perspectives":
                payload = self.read_json_body()
                with CatalogRepository(self.catalog_db_path) as catalog:
                    result = CatalogViewerService(
                        catalog
                    ).create_perspective(payload.get("name", ""))
                json_response(self, result, status=201)
            elif parsed.path == "/api/catalog/datasets":
                payload = self.read_json_body()
                perspective_id, _dataset_id = self.catalog_context(
                    query,
                    payload,
                )
                if not perspective_id:
                    raise ValueError("Select a Study Perspective.")
                db_id, db_path = self.resolve_database(
                    payload.get("db") or self.database_id_from_query(query)
                )
                with CatalogRepository(self.catalog_db_path) as catalog:
                    result = CatalogViewerService(catalog).create_dataset(
                        perspective_id,
                        payload.get("name", ""),
                        database_id=db_id,
                        database_path=db_path,
                    )
                json_response(self, result, status=201)
            elif parsed.path == "/api/catalog/groups":
                payload = self.read_json_body()
                perspective_id, dataset_id = self.catalog_context(
                    query,
                    payload,
                )
                if not perspective_id or not dataset_id:
                    raise ValueError(
                        "Select a Study Perspective and Dataset."
                    )
                db_id, db_path = self.resolve_database(
                    payload.get("db") or self.database_id_from_query(query)
                )
                with CatalogRepository(self.catalog_db_path) as catalog:
                    result = CatalogViewerService(catalog).create_group(
                        perspective_id,
                        dataset_id,
                        payload.get("name", ""),
                        payload.get("sample_ids", []),
                        database_id=db_id,
                        database_path=db_path,
                    )
                json_response(self, result, status=201)
            elif parsed.path == "/api/catalog/dataset-cohorts":
                payload = self.read_json_body()
                perspective_id, dataset_id = self.catalog_context(
                    query,
                    payload,
                )
                if not perspective_id or not dataset_id:
                    raise ValueError(
                        "Select a Study Perspective and Dataset."
                    )
                db_id, db_path = self.resolve_database(
                    payload.get("db") or self.database_id_from_query(query)
                )
                with CatalogRepository(self.catalog_db_path) as catalog:
                    result = CatalogViewerService(catalog).attach_database(
                        perspective_id,
                        dataset_id,
                        database_id=db_id,
                        database_path=db_path,
                    )
                json_response(self, result, status=201)
            elif parsed.path == "/api/upload-database":
                fields, files = parse_multipart_form(self.headers, self.read_body())
                upload = files.get("database")
                if not upload:
                    json_response(self, {"error": "Choose a database file to upload."}, status=400)
                    return
                replace = fields.get("replace") == "1"
                db_path = save_uploaded_database(
                    self.database_dir,
                    upload["filename"],
                    upload["content"],
                    replace=replace,
                )
                json_response(
                    self,
                    {
                        "database": database_id_for_path(db_path),
                        "path": str(db_path),
                    },
                    status=201,
                )
            elif parsed.path == "/api/derived-samples":
                payload = self.read_json_body()
                db_id, db_path = self.resolve_database(
                    payload.get("db") or self.database_id_from_query(query)
                )
                compare_sample_ids = payload.get("compare_sample_id", [])
                if len(compare_sample_ids) < 2:
                    json_response(
                        self,
                        {"error": "Select at least two samples before creating a derived sample."},
                        status=400,
                    )
                    return

                perspective_id, dataset_id = self.catalog_context(
                    query,
                    payload,
                )
                next_id = self.next_derived_sample_ids.get(db_id, 1)
                label = payload.get("label") or f"Comparison {next_id}"
                if perspective_id and dataset_id:
                    with CatalogRepository(self.catalog_db_path) as catalog:
                        result = CatalogViewerService(
                            catalog
                        ).save_comparison(
                            perspective_id,
                            dataset_id,
                            label,
                            compare_sample_ids,
                            database_id=db_id,
                            database_path=db_path,
                            sample_statuses=parse_sample_statuses(
                                payload.get("sample_status", [])
                            ),
                            global_statuses=payload.get("status", []),
                            filters=MutationFilters(
                                position=payload.get("position", ""),
                                alt=payload.get("alt", ""),
                                af_rules=tuple(
                                    parse_af_rules(
                                        payload.get("af_rule", [])
                                    )
                                ),
                                metadata_filters=tuple(
                                    parse_metadata_filters(
                                        payload.get(
                                            "metadata_filter",
                                            [],
                                        )
                                    )
                                ),
                            ),
                        )
                    json_response(self, result, status=201)
                else:
                    derived_samples = self.derived_samples(db_id)
                    with self.open_database(db_path) as connection:
                        sample = create_derived_sample(
                            connection,
                            derived_samples,
                            self.allocate_derived_sample_id(db_id),
                            label,
                            compare_sample_ids=compare_sample_ids,
                            position=payload.get("position", ""),
                            alt=payload.get("alt", ""),
                            af_rules=parse_af_rules(
                                payload.get("af_rule", [])
                            ),
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
        except ValueError as exc:
            json_response(self, {"error": str(exc)}, status=400)
        except FileExistsError as exc:
            json_response(self, {"error": str(exc)}, status=409)
        except Exception as exc:
            json_response(self, {"error": str(exc)}, status=500)

    def do_DELETE(self):
        parsed = urlparse(self.path)
        query = parse_qs(parsed.query)

        try:
            db_id, _db_path = self.resolve_database(self.database_id_from_query(query))
            derived_samples = self.derived_samples(db_id)

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
                perspective_id, dataset_id = self.catalog_context(query)
                if perspective_id and dataset_id:
                    json_response(
                        self,
                        {
                            "error": (
                                "Durable derived samples and their run "
                                "history are immutable."
                            )
                        },
                        status=405,
                    )
                    return
                removed = derived_samples.pop(sample_id, None)
                if removed is None:
                    json_response(self, {"error": "Derived sample not found."}, status=404)
                    return
                json_response(self, removed.sample_row())
            else:
                json_response(self, {"error": "Not found"}, status=404)
        except ValueError as exc:
            json_response(self, {"error": str(exc)}, status=400)
        except Exception as exc:
            json_response(self, {"error": str(exc)}, status=500)


def first_database_id(databases):
    return next(iter(databases), "")


def configure_databases(db_path=None, database_dir=None):
    if db_path is not None:
        db_path = Path(db_path).resolve()
        database_dir = db_path.parent
        default_db_id = database_id_for_path(db_path)
    else:
        database_dir = Path(database_dir).resolve()
        default_db_id = ""

    databases = discover_databases(database_dir)
    if not databases:
        raise FileNotFoundError(f"No SQLite databases found in: {database_dir}")

    if not default_db_id:
        default_db_id = first_database_id(databases)
    if default_db_id not in databases:
        raise FileNotFoundError(f"Database not found: {database_dir / default_db_id}")

    ViewerHandler.database_dir = database_dir
    ViewerHandler.default_db_id = default_db_id
    return databases[default_db_id], databases


def run_server(
    db_path=None,
    database_dir=None,
    annotation_db_path=DEFAULT_ANNOTATION_DATABASE_PATH,
    catalog_db_path=DEFAULT_CATALOG_DATABASE_PATH,
    host="127.0.0.1",
    port=8000,
):
    selected_db_path, databases = configure_databases(
        db_path=db_path,
        database_dir=database_dir,
    )
    ViewerHandler.annotation_db_path = Path(annotation_db_path).resolve()
    ViewerHandler.catalog_db_path = Path(catalog_db_path).resolve()
    with CatalogRepository(ViewerHandler.catalog_db_path):
        pass
    server = ThreadingHTTPServer((host, port), ViewerHandler)
    database_names = ", ".join(databases)
    print(f"Serving {selected_db_path} at http://{host}:{port}")
    print(f"Available databases: {database_names}")
    server.serve_forever()


def main():
    parser = argparse.ArgumentParser(description="View the mito SQLite database.")
    source = parser.add_mutually_exclusive_group()
    source.add_argument(
        "--db",
        type=Path,
        help="SQLite database to open initially.",
    )
    source.add_argument(
        "--db-dir",
        type=Path,
        default=DEFAULT_DATABASE_DIR,
        help="Directory containing SQLite databases for the dropdown.",
    )
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=18000)
    parser.add_argument(
        "--annotations-db",
        type=Path,
        default=DEFAULT_ANNOTATION_DATABASE_PATH,
        help="Persistent annotation cache used by the Annotations tab.",
    )
    parser.add_argument(
        "--catalog-db",
        type=Path,
        default=DEFAULT_CATALOG_DATABASE_PATH,
        help=(
            "Writable catalog used for perspectives, datasets, groups, "
            "and durable derived samples."
        ),
    )
    args = parser.parse_args()
    run_server(
        db_path=args.db,
        database_dir=None if args.db else args.db_dir,
        annotation_db_path=args.annotations_db,
        catalog_db_path=args.catalog_db,
        host=args.host,
        port=args.port,
    )


if __name__ == "__main__":
    main()
