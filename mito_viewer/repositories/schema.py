"""Database schema inspection shared by repository implementations."""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping


@dataclass(frozen=True)
class DatabaseSchemaReport:
    """Result of checking a SQLite file against one repository contract."""

    path: Path
    database_kind: str
    user_version: int = 0
    missing_tables: tuple[str, ...] = ()
    missing_columns: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    error: str = ""

    @property
    def valid(self) -> bool:
        return not self.error and not self.missing_tables and not self.missing_columns

    @property
    def version_label(self) -> str:
        if self.user_version:
            return str(self.user_version)
        return "legacy-unversioned"

    def require_valid(self) -> "DatabaseSchemaReport":
        """Return this report or raise a concrete compatibility error."""
        if self.valid:
            return self

        details = []
        if self.error:
            details.append(self.error)
        if self.missing_tables:
            details.append("missing tables: " + ", ".join(self.missing_tables))
        for table, columns in sorted(self.missing_columns.items()):
            details.append(f"{table} missing columns: {', '.join(columns)}")
        joined = "; ".join(details) or "unknown schema error"
        raise ValueError(
            f"{self.path.name} is not a compatible {self.database_kind} "
            f"database: {joined}."
        )


def read_only_connection(path: str | Path) -> sqlite3.Connection:
    """Open an existing SQLite database in enforced read-only mode."""
    database_path = Path(path).resolve()
    if not database_path.is_file():
        raise FileNotFoundError(f"Database not found: {database_path}")
    connection = sqlite3.connect(
        f"file:{database_path}?mode=ro",
        uri=True,
    )
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA query_only = ON")
    return connection


def inspect_database_schema(
    path: str | Path,
    *,
    database_kind: str,
    required_schema: Mapping[str, frozenset[str]],
) -> DatabaseSchemaReport:
    """Inspect a SQLite file without modifying it."""
    database_path = Path(path).resolve()
    if not database_path.is_file():
        return DatabaseSchemaReport(
            path=database_path,
            database_kind=database_kind,
            error="file does not exist",
        )

    try:
        connection = read_only_connection(database_path)
        try:
            user_version = int(
                connection.execute("PRAGMA user_version").fetchone()[0]
            )
            table_names = {
                row["name"]
                for row in connection.execute(
                    "SELECT name FROM sqlite_master WHERE type = 'table'"
                )
            }
            missing_tables = tuple(
                sorted(set(required_schema) - table_names)
            )
            missing_columns = {}
            for table in sorted(set(required_schema) & table_names):
                columns = {
                    row["name"]
                    for row in connection.execute(
                        f"PRAGMA table_info({table})"
                    )
                }
                missing = tuple(sorted(required_schema[table] - columns))
                if missing:
                    missing_columns[table] = missing
        finally:
            connection.close()
    except sqlite3.Error as exc:
        return DatabaseSchemaReport(
            path=database_path,
            database_kind=database_kind,
            error=str(exc),
        )

    return DatabaseSchemaReport(
        path=database_path,
        database_kind=database_kind,
        user_version=user_version,
        missing_tables=missing_tables,
        missing_columns=missing_columns,
    )
