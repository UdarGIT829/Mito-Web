#!/usr/bin/env python3
"""Print one sample's SQLite database entry to stdout."""

import argparse
import json
import sqlite3
import sys
from pathlib import Path


DEFAULT_DB_PATH = Path(__file__).resolve().parent / "master.sqlite"


def connect_database(db_path):
    """Open the mutation database with row access by column name."""
    connection = sqlite3.connect(db_path)
    connection.row_factory = sqlite3.Row
    return connection


def normalize_population_key(values):
    """Return a database population_key from CLI tag values."""
    if not values:
        return None
    if len(values) == 1:
        return values[0].replace("_", "|")
    return "|".join(values)


def sample_lookup_parts(sample):
    """Split a compact sample label into subject and optional tags."""
    if "|" in sample:
        subject_id, _, tags = sample.partition("|")
        return subject_id, normalize_population_key([tags])

    parts = sample.split("_")
    if len(parts) > 1:
        return parts[0], normalize_population_key(parts[1:])

    return sample, None


def fetch_sample_options(connection, subject_id):
    """Return all database sample rows for a subject."""
    rows = connection.execute(
        """
        SELECT
            samples.id,
            subjects.subject_id,
            samples.population_key,
            samples.source_file,
            COUNT(mutations.id) AS mutation_count
        FROM samples
        JOIN subjects ON subjects.id = samples.subject_id
        LEFT JOIN mutations ON mutations.sample_id = samples.id
        WHERE subjects.subject_id = ?
        GROUP BY samples.id
        ORDER BY samples.population_key
        """,
        (subject_id,),
    ).fetchall()
    return [dict(row) for row in rows]


def fetch_sample(connection, sample, population_key=None):
    """Return one sample row resolved by database id or subject/population."""
    if str(sample).isdigit() and population_key is None:
        row = connection.execute(
            """
            SELECT
                samples.id,
                subjects.subject_id,
                samples.population_key,
                samples.source_file,
                COUNT(mutations.id) AS mutation_count
            FROM samples
            JOIN subjects ON subjects.id = samples.subject_id
            LEFT JOIN mutations ON mutations.sample_id = samples.id
            WHERE samples.id = ?
            GROUP BY samples.id
            """,
            (int(sample),),
        ).fetchone()
        return dict(row) if row else None

    subject_id, compact_population_key = sample_lookup_parts(str(sample))
    population_key = population_key or compact_population_key

    if population_key is None:
        options = fetch_sample_options(connection, subject_id)
        if len(options) == 1:
            return options[0]
        if not options:
            return None
        choices = ", ".join(
            f"{row['id']}={row['subject_id']}_{row['population_key'].replace('|', '_')}"
            for row in options
        )
        raise ValueError(
            f"Subject {subject_id!r} has multiple samples. "
            f"Choose one with tags or --population. Options: {choices}"
        )

    row = connection.execute(
        """
        SELECT
            samples.id,
            subjects.subject_id,
            samples.population_key,
            samples.source_file,
            COUNT(mutations.id) AS mutation_count
        FROM samples
        JOIN subjects ON subjects.id = samples.subject_id
        LEFT JOIN mutations ON mutations.sample_id = samples.id
        WHERE subjects.subject_id = ? AND samples.population_key = ?
        GROUP BY samples.id
        """,
        (subject_id, population_key),
    ).fetchone()
    return dict(row) if row else None


def fetch_population_tags(connection, sample_id):
    """Return ordered population tags for a sample."""
    rows = connection.execute(
        """
        SELECT tag
        FROM sample_population_tags
        WHERE sample_id = ?
        ORDER BY tag_order
        """,
        (sample_id,),
    ).fetchall()
    return [row["tag"] for row in rows]


def fetch_mutations(connection, sample_id):
    """Return mutation rows for a sample."""
    rows = connection.execute(
        """
        SELECT
            id,
            pos,
            ref,
            vcf_ref,
            alt,
            af,
            filter,
            metadata_json
        FROM mutations
        WHERE sample_id = ?
        ORDER BY pos, id
        """,
        (sample_id,),
    ).fetchall()
    return [dict(row) for row in rows]


def fetch_mutation_alts(connection, mutation_ids):
    """Return per-ALT rows grouped by mutation id."""
    if not mutation_ids:
        return {}

    placeholders = ",".join("?" for _ in mutation_ids)
    rows = connection.execute(
        f"""
        SELECT mutation_id, alt_index, alt, af, af_text
        FROM mutation_alts
        WHERE mutation_id IN ({placeholders})
        ORDER BY mutation_id, alt_index
        """,
        mutation_ids,
    ).fetchall()

    grouped = {}
    for row in rows:
        grouped.setdefault(row["mutation_id"], []).append(dict(row))
    return grouped


def metadata_text(metadata_json):
    """Return metadata_json as stable key=value text."""
    try:
        metadata = json.loads(metadata_json)
    except json.JSONDecodeError:
        return metadata_json
    return "; ".join(f"{key}={metadata[key]}" for key in sorted(metadata))


def print_sample_database_entry(sample, db_path=DEFAULT_DB_PATH, output=None):
    """Print a sample's full database entry to a text stream.

    Args:
        sample: Database sample id, subject id, or compact label such as
            AP4_EV_PCR. Subject-only values are accepted when they resolve to
            exactly one database sample.
        db_path: Path to master.sqlite.
        output: Text stream to write to. Defaults to stdout.
    """
    output = output or sys.stdout

    with connect_database(db_path) as connection:
        sample_row = fetch_sample(connection, sample)
        if sample_row is None:
            raise ValueError(f"No sample found for {sample!r} in {db_path}")

        tags = fetch_population_tags(connection, sample_row["id"])
        mutations = fetch_mutations(connection, sample_row["id"])
        alts_by_mutation = fetch_mutation_alts(
            connection,
            [mutation["id"] for mutation in mutations],
        )

    print("Sample", file=output)
    print(f"  database id: {sample_row['id']}", file=output)
    print(f"  subject id: {sample_row['subject_id']}", file=output)
    print(f"  population key: {sample_row['population_key']}", file=output)
    print(f"  population tags: {', '.join(tags)}", file=output)
    print(f"  source file: {sample_row['source_file']}", file=output)
    print(f"  mutation count: {sample_row['mutation_count']}", file=output)
    print("", file=output)
    print("Mutations", file=output)

    for mutation in mutations:
        print(f"  mutation id: {mutation['id']}", file=output)
        print(f"    pos: {mutation['pos']}", file=output)
        print(f"    ref: {mutation['ref']}", file=output)
        print(f"    vcf_ref: {mutation['vcf_ref']}", file=output)
        print(f"    alt: {mutation['alt']}", file=output)
        print(f"    af: {mutation['af']}", file=output)
        print(f"    filter: {mutation['filter']}", file=output)
        print(f"    metadata: {metadata_text(mutation['metadata_json'])}", file=output)
        print("    alt rows:", file=output)
        for alt in alts_by_mutation.get(mutation["id"], []):
            print(
                "      "
                f"{alt['alt_index']}: "
                f"alt={alt['alt']} "
                f"af={'' if alt['af'] is None else alt['af']} "
                f"af_text={alt['af_text'] or ''}",
                file=output,
            )


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Print one sample's full SQLite database entry to stdout.",
    )
    parser.add_argument(
        "sample",
        help="Sample database id, subject id, or compact label like AP4_EV_PCR.",
    )
    parser.add_argument(
        "tags",
        nargs="*",
        help="Optional population tags, for example EV PCR.",
    )
    parser.add_argument(
        "--population",
        help="Population key, for example EV|PCR or EV_PCR.",
    )
    parser.add_argument(
        "--db",
        default=DEFAULT_DB_PATH,
        type=Path,
        help=f"SQLite database path. Default: {DEFAULT_DB_PATH}",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    population_key = args.population or normalize_population_key(args.tags)

    with connect_database(args.db) as connection:
        sample_row = fetch_sample(connection, args.sample, population_key)
        if sample_row is None:
            print(f"No sample found for {args.sample!r}", file=sys.stderr)
            return 1

    print_sample_database_entry(sample_row["id"], db_path=args.db)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
