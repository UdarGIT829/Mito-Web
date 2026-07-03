#!/usr/bin/env python3
"""Print the SQLite database schema and table links to stdout."""

import argparse
import sqlite3
import sys
from pathlib import Path


DEFAULT_DB_PATH = Path(__file__).resolve().parent / "master.sqlite"


def connect_database(db_path):
    """Open a SQLite database with row access by column name."""
    connection = sqlite3.connect(db_path)
    connection.row_factory = sqlite3.Row
    return connection


def fetch_tables(connection):
    """Return user-defined table names in dependency-friendly order."""
    rows = connection.execute(
        """
        SELECT name
        FROM sqlite_master
        WHERE type = 'table'
          AND name NOT LIKE 'sqlite_%'
        ORDER BY name
        """
    ).fetchall()
    return [row["name"] for row in rows]


def fetch_columns(connection, table):
    """Return column metadata from PRAGMA table_info."""
    return [
        dict(row)
        for row in connection.execute(f"PRAGMA table_info({quote_identifier(table)})")
    ]


def fetch_foreign_keys(connection, table):
    """Return foreign-key metadata from PRAGMA foreign_key_list."""
    return [
        dict(row)
        for row in connection.execute(
            f"PRAGMA foreign_key_list({quote_identifier(table)})"
        )
    ]


def fetch_indexes(connection, table):
    """Return index metadata from PRAGMA index_list and index_info."""
    indexes = []
    rows = connection.execute(f"PRAGMA index_list({quote_identifier(table)})").fetchall()
    for row in rows:
        index = dict(row)
        index["columns"] = [
            column["name"]
            for column in connection.execute(
                f"PRAGMA index_info({quote_identifier(index['name'])})"
            )
        ]
        indexes.append(index)
    return indexes


def fetch_create_sql(connection, table):
    """Return the CREATE TABLE statement for a table."""
    row = connection.execute(
        """
        SELECT sql
        FROM sqlite_master
        WHERE type = 'table' AND name = ?
        """,
        (table,),
    ).fetchone()
    return row["sql"] if row else ""


def fetch_row_count(connection, table):
    """Return a table's row count."""
    row = connection.execute(
        f"SELECT COUNT(*) AS count FROM {quote_identifier(table)}"
    ).fetchone()
    return row["count"]


def quote_identifier(identifier):
    """Quote a SQLite identifier."""
    return '"' + str(identifier).replace('"', '""') + '"'


def column_flags(column):
    """Return a compact description of column flags."""
    flags = []
    if column["pk"]:
        flags.append("PRIMARY KEY")
    if column["notnull"]:
        flags.append("NOT NULL")
    if column["dflt_value"] is not None:
        flags.append(f"DEFAULT {column['dflt_value']}")
    return ", ".join(flags)


def print_database_schema(db_path=DEFAULT_DB_PATH, output=None, include_sql=False):
    """Print schema details and table links to a text stream."""
    output = output or sys.stdout

    with connect_database(db_path) as connection:
        tables = fetch_tables(connection)
        table_details = {
            table: {
                "columns": fetch_columns(connection, table),
                "foreign_keys": fetch_foreign_keys(connection, table),
                "indexes": fetch_indexes(connection, table),
                "create_sql": fetch_create_sql(connection, table),
                "row_count": fetch_row_count(connection, table),
            }
            for table in tables
        }

    print(f"Database: {Path(db_path)}", file=output)
    print("", file=output)
    print("Tables", file=output)
    for table in tables:
        details = table_details[table]
        print(f"  {table} ({details['row_count']} rows)", file=output)
    print("", file=output)

    print("Table Details", file=output)
    for table in tables:
        details = table_details[table]
        print(f"  {table}", file=output)
        print("    columns:", file=output)
        for column in details["columns"]:
            flags = column_flags(column)
            suffix = f" [{flags}]" if flags else ""
            print(
                f"      {column['name']}: {column['type'] or 'ANY'}{suffix}",
                file=output,
            )

        if details["indexes"]:
            print("    indexes:", file=output)
            for index in details["indexes"]:
                unique = " UNIQUE" if index["unique"] else ""
                origin = f" origin={index['origin']}" if index["origin"] else ""
                columns = ", ".join(index["columns"])
                print(
                    f"      {index['name']}:{unique} ({columns}){origin}",
                    file=output,
                )

        if details["foreign_keys"]:
            print("    links:", file=output)
            for foreign_key in details["foreign_keys"]:
                print(
                    "      "
                    f"{table}.{foreign_key['from']} -> "
                    f"{foreign_key['table']}.{foreign_key['to']}",
                    file=output,
                )

        if include_sql and details["create_sql"]:
            print("    create sql:", file=output)
            for line in details["create_sql"].splitlines():
                print(f"      {line}", file=output)
        print("", file=output)

    print("Relationship Map", file=output)
    links = []
    for table in tables:
        for foreign_key in table_details[table]["foreign_keys"]:
            links.append(
                (
                    table,
                    foreign_key["from"],
                    foreign_key["table"],
                    foreign_key["to"],
                )
            )

    if not links:
        print("  No foreign-key links found.", file=output)
    else:
        for source_table, source_column, target_table, target_column in links:
            print(
                f"  {source_table}.{source_column} -> "
                f"{target_table}.{target_column}",
                file=output,
            )


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Print SQLite table schema, indexes, and foreign-key links.",
    )
    parser.add_argument(
        "--db",
        default=DEFAULT_DB_PATH,
        type=Path,
        help=f"SQLite database path. Default: {DEFAULT_DB_PATH}",
    )
    parser.add_argument(
        "--sql",
        action="store_true",
        help="Include raw CREATE TABLE SQL for each table.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    print_database_schema(args.db, include_sql=args.sql)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
