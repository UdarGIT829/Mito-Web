"""Versioned SQLite schema and migrations for the writable catalog."""

from __future__ import annotations

import sqlite3
from pathlib import Path

from mito_viewer.repositories.schema import (
    DatabaseSchemaReport,
    inspect_database_schema,
)


CATALOG_SCHEMA_VERSION = 3


class CatalogMigrationError(RuntimeError):
    """Raised when a catalog cannot be migrated safely."""


CATALOG_REQUIRED_SCHEMA = {
    "study_perspectives": frozenset(
        {"id", "name", "description", "visibility", "created_at", "updated_at"}
    ),
    "perspective_credentials": frozenset(
        {
            "perspective_id",
            "password_hash",
            "credential_version",
            "updated_at",
        }
    ),
    "cohorts": frozenset(
        {
            "id",
            "owner_perspective_id",
            "name",
            "description",
            "cohort_type",
            "source_database_identifier",
            "source_database_fingerprint",
            "reference_id",
            "normalization_version",
            "provenance_json",
            "visibility",
            "created_at",
            "updated_at",
        }
    ),
    "perspective_cohort_access": frozenset(
        {"perspective_id", "cohort_id", "access_level", "granted_at"}
    ),
    "datasets": frozenset(
        {
            "id",
            "perspective_id",
            "name",
            "description",
            "derived_results_cohort_id",
            "visibility",
            "created_at",
            "updated_at",
        }
    ),
    "dataset_cohorts": frozenset(
        {"dataset_id", "cohort_id", "display_order", "added_at"}
    ),
    "catalog_samples": frozenset(
        {
            "id",
            "cohort_id",
            "sample_type",
            "source_sample_id",
            "display_label",
            "source_fingerprint",
            "metadata_json",
            "created_at",
        }
    ),
    "sample_groups": frozenset(
        {
            "id",
            "perspective_id",
            "dataset_id",
            "name",
            "description",
            "membership_fingerprint",
            "visibility",
            "created_at",
            "updated_at",
        }
    ),
    "sample_group_members": frozenset(
        {"group_id", "catalog_sample_id", "display_order", "added_at"}
    ),
    "derived_samples": frozenset(
        {
            "id",
            "catalog_sample_id",
            "perspective_id",
            "dataset_id",
            "name",
            "description",
            "visibility",
            "current_run_id",
            "created_at",
            "updated_at",
        }
    ),
    "derived_runs": frozenset(
        {
            "id",
            "derived_sample_id",
            "definition_json",
            "input_snapshot_fingerprint",
            "reference_id",
            "normalization_version",
            "status",
            "output_count",
            "output_fingerprint",
            "error_text",
            "created_at",
            "completed_at",
        }
    ),
    "derived_run_input_samples": frozenset(
        {
            "run_id",
            "input_clause_index",
            "catalog_sample_id",
            "source_group_id",
            "group_membership_fingerprint",
            "sample_fingerprint",
            "source_database_fingerprint",
            "input_role",
            "display_order",
        }
    ),
    "derived_run_parent_runs": frozenset(
        {
            "run_id",
            "input_clause_index",
            "parent_run_id",
            "source_group_id",
            "group_membership_fingerprint",
            "input_role",
            "display_order",
        }
    ),
    "derived_run_alleles": frozenset(
        {"id", "run_id", "position", "ref", "alt", "result_metadata_json"}
    ),
    "derived_allele_sources": frozenset(
        {
            "id",
            "derived_allele_id",
            "source_catalog_sample_id",
            "parent_derived_allele_id",
            "original_catalog_sample_id",
            "source_mutation_id",
            "source_alt_index",
            "input_role",
            "evaluation_status",
            "af",
            "af_text",
            "filter_text",
            "metadata_json",
            "lineage_json",
        }
    ),
}


MIGRATION_1 = """
CREATE TABLE study_perspectives (
    id TEXT PRIMARY KEY,
    name TEXT NOT NULL COLLATE NOCASE UNIQUE,
    description TEXT NOT NULL DEFAULT '',
    visibility TEXT NOT NULL DEFAULT 'private'
        CHECK (visibility IN ('private', 'perspective', 'shared', 'public')),
    created_at TEXT NOT NULL,
    updated_at TEXT NOT NULL
);

CREATE TABLE perspective_credentials (
    perspective_id TEXT PRIMARY KEY,
    password_hash TEXT NOT NULL,
    credential_version INTEGER NOT NULL DEFAULT 1
        CHECK (credential_version > 0),
    updated_at TEXT NOT NULL,
    FOREIGN KEY (perspective_id)
        REFERENCES study_perspectives(id) ON DELETE CASCADE
);

CREATE TABLE cohorts (
    id TEXT PRIMARY KEY,
    owner_perspective_id TEXT,
    name TEXT NOT NULL,
    description TEXT NOT NULL DEFAULT '',
    cohort_type TEXT NOT NULL
        CHECK (cohort_type IN ('source', 'derived')),
    source_database_identifier TEXT NOT NULL DEFAULT '',
    source_database_fingerprint TEXT NOT NULL DEFAULT '',
    reference_id TEXT NOT NULL DEFAULT '',
    normalization_version TEXT NOT NULL DEFAULT '',
    provenance_json TEXT NOT NULL DEFAULT '{}',
    visibility TEXT NOT NULL DEFAULT 'private'
        CHECK (visibility IN ('private', 'perspective', 'shared', 'public')),
    created_at TEXT NOT NULL,
    updated_at TEXT NOT NULL,
    CHECK (
        cohort_type != 'source'
        OR (
            source_database_identifier != ''
            AND source_database_fingerprint != ''
        )
    ),
    CHECK (cohort_type != 'derived' OR owner_perspective_id IS NOT NULL),
    FOREIGN KEY (owner_perspective_id)
        REFERENCES study_perspectives(id) ON DELETE RESTRICT
);

CREATE TABLE perspective_cohort_access (
    perspective_id TEXT NOT NULL,
    cohort_id TEXT NOT NULL,
    access_level TEXT NOT NULL
        CHECK (access_level IN ('view', 'use', 'manage')),
    granted_at TEXT NOT NULL,
    PRIMARY KEY (perspective_id, cohort_id),
    FOREIGN KEY (perspective_id)
        REFERENCES study_perspectives(id) ON DELETE CASCADE,
    FOREIGN KEY (cohort_id)
        REFERENCES cohorts(id) ON DELETE CASCADE
);

CREATE TABLE datasets (
    id TEXT PRIMARY KEY,
    perspective_id TEXT NOT NULL,
    name TEXT NOT NULL COLLATE NOCASE,
    description TEXT NOT NULL DEFAULT '',
    derived_results_cohort_id TEXT UNIQUE,
    visibility TEXT NOT NULL DEFAULT 'private'
        CHECK (visibility IN ('private', 'perspective', 'shared', 'public')),
    created_at TEXT NOT NULL,
    updated_at TEXT NOT NULL,
    UNIQUE (perspective_id, name),
    UNIQUE (id, perspective_id),
    FOREIGN KEY (perspective_id)
        REFERENCES study_perspectives(id) ON DELETE CASCADE,
    FOREIGN KEY (derived_results_cohort_id)
        REFERENCES cohorts(id) ON DELETE RESTRICT
);

CREATE TABLE dataset_cohorts (
    dataset_id TEXT NOT NULL,
    cohort_id TEXT NOT NULL,
    display_order INTEGER NOT NULL DEFAULT 0,
    added_at TEXT NOT NULL,
    PRIMARY KEY (dataset_id, cohort_id),
    FOREIGN KEY (dataset_id) REFERENCES datasets(id) ON DELETE CASCADE,
    FOREIGN KEY (cohort_id) REFERENCES cohorts(id) ON DELETE RESTRICT
);

CREATE TABLE catalog_samples (
    id TEXT PRIMARY KEY,
    cohort_id TEXT NOT NULL,
    sample_type TEXT NOT NULL
        CHECK (sample_type IN ('observed', 'derived')),
    source_sample_id TEXT,
    display_label TEXT NOT NULL,
    source_fingerprint TEXT NOT NULL DEFAULT '',
    metadata_json TEXT NOT NULL DEFAULT '{}',
    created_at TEXT NOT NULL,
    CHECK (
        (sample_type = 'observed' AND source_sample_id IS NOT NULL)
        OR (sample_type = 'derived' AND source_sample_id IS NULL)
    ),
    UNIQUE (cohort_id, source_sample_id),
    FOREIGN KEY (cohort_id) REFERENCES cohorts(id) ON DELETE RESTRICT
);

CREATE TABLE sample_groups (
    id TEXT PRIMARY KEY,
    perspective_id TEXT NOT NULL,
    dataset_id TEXT NOT NULL,
    name TEXT NOT NULL COLLATE NOCASE,
    description TEXT NOT NULL DEFAULT '',
    visibility TEXT NOT NULL DEFAULT 'private'
        CHECK (visibility IN ('private', 'perspective', 'shared', 'public')),
    membership_fingerprint TEXT NOT NULL DEFAULT '',
    created_at TEXT NOT NULL,
    updated_at TEXT NOT NULL,
    UNIQUE (dataset_id, name),
    FOREIGN KEY (dataset_id, perspective_id)
        REFERENCES datasets(id, perspective_id) ON DELETE CASCADE
);

CREATE TABLE sample_group_members (
    group_id TEXT NOT NULL,
    catalog_sample_id TEXT NOT NULL,
    display_order INTEGER NOT NULL DEFAULT 0,
    added_at TEXT NOT NULL,
    PRIMARY KEY (group_id, catalog_sample_id),
    FOREIGN KEY (group_id) REFERENCES sample_groups(id) ON DELETE CASCADE,
    FOREIGN KEY (catalog_sample_id)
        REFERENCES catalog_samples(id) ON DELETE RESTRICT
);

CREATE TABLE derived_samples (
    id TEXT PRIMARY KEY,
    catalog_sample_id TEXT NOT NULL UNIQUE,
    perspective_id TEXT NOT NULL,
    dataset_id TEXT NOT NULL,
    name TEXT NOT NULL COLLATE NOCASE,
    description TEXT NOT NULL DEFAULT '',
    visibility TEXT NOT NULL DEFAULT 'private'
        CHECK (visibility IN ('private', 'perspective', 'shared', 'public')),
    current_run_id TEXT,
    created_at TEXT NOT NULL,
    updated_at TEXT NOT NULL,
    UNIQUE (dataset_id, name),
    FOREIGN KEY (catalog_sample_id)
        REFERENCES catalog_samples(id) ON DELETE RESTRICT,
    FOREIGN KEY (dataset_id, perspective_id)
        REFERENCES datasets(id, perspective_id) ON DELETE CASCADE,
    FOREIGN KEY (current_run_id)
        REFERENCES derived_runs(id) ON DELETE SET NULL
);

CREATE TABLE derived_runs (
    id TEXT PRIMARY KEY,
    derived_sample_id TEXT NOT NULL,
    definition_json TEXT NOT NULL,
    input_snapshot_fingerprint TEXT NOT NULL,
    reference_id TEXT NOT NULL DEFAULT '',
    normalization_version TEXT NOT NULL DEFAULT '',
    status TEXT NOT NULL
        CHECK (status IN ('pending', 'running', 'completed', 'failed')),
    output_count INTEGER NOT NULL DEFAULT 0 CHECK (output_count >= 0),
    error_text TEXT NOT NULL DEFAULT '',
    created_at TEXT NOT NULL,
    completed_at TEXT,
    FOREIGN KEY (derived_sample_id)
        REFERENCES derived_samples(id) ON DELETE CASCADE
);

CREATE TABLE derived_run_input_samples (
    run_id TEXT NOT NULL,
    catalog_sample_id TEXT NOT NULL,
    source_group_id TEXT,
    input_role TEXT NOT NULL,
    display_order INTEGER NOT NULL DEFAULT 0,
    PRIMARY KEY (run_id, catalog_sample_id, input_role),
    FOREIGN KEY (run_id) REFERENCES derived_runs(id) ON DELETE CASCADE,
    FOREIGN KEY (catalog_sample_id)
        REFERENCES catalog_samples(id) ON DELETE RESTRICT,
    FOREIGN KEY (source_group_id)
        REFERENCES sample_groups(id) ON DELETE SET NULL
);

CREATE TABLE derived_run_parent_runs (
    run_id TEXT NOT NULL,
    parent_run_id TEXT NOT NULL,
    input_role TEXT NOT NULL,
    display_order INTEGER NOT NULL DEFAULT 0,
    CHECK (run_id != parent_run_id),
    PRIMARY KEY (run_id, parent_run_id, input_role),
    FOREIGN KEY (run_id) REFERENCES derived_runs(id) ON DELETE CASCADE,
    FOREIGN KEY (parent_run_id)
        REFERENCES derived_runs(id) ON DELETE RESTRICT
);

CREATE TABLE derived_run_alleles (
    id TEXT PRIMARY KEY,
    run_id TEXT NOT NULL,
    position INTEGER NOT NULL CHECK (position > 0),
    ref TEXT NOT NULL,
    alt TEXT NOT NULL,
    result_metadata_json TEXT NOT NULL DEFAULT '{}',
    UNIQUE (run_id, position, ref, alt),
    FOREIGN KEY (run_id) REFERENCES derived_runs(id) ON DELETE CASCADE
);

CREATE TABLE derived_allele_sources (
    id TEXT PRIMARY KEY,
    derived_allele_id TEXT NOT NULL,
    source_catalog_sample_id TEXT,
    parent_derived_allele_id TEXT,
    original_catalog_sample_id TEXT,
    source_mutation_id TEXT,
    source_alt_index INTEGER,
    input_role TEXT NOT NULL,
    evaluation_status TEXT NOT NULL
        CHECK (
            evaluation_status IN (
                'qualifying_present',
                'filtered_out',
                'not_observed'
            )
        ),
    af REAL,
    af_text TEXT NOT NULL DEFAULT '',
    filter_text TEXT NOT NULL DEFAULT '',
    metadata_json TEXT NOT NULL DEFAULT '{}',
    lineage_json TEXT NOT NULL DEFAULT '{}',
    CHECK (
        source_catalog_sample_id IS NOT NULL
        OR parent_derived_allele_id IS NOT NULL
    ),
    FOREIGN KEY (derived_allele_id)
        REFERENCES derived_run_alleles(id) ON DELETE CASCADE,
    FOREIGN KEY (source_catalog_sample_id)
        REFERENCES catalog_samples(id) ON DELETE RESTRICT,
    FOREIGN KEY (parent_derived_allele_id)
        REFERENCES derived_run_alleles(id) ON DELETE RESTRICT,
    FOREIGN KEY (original_catalog_sample_id)
        REFERENCES catalog_samples(id) ON DELETE RESTRICT
);

CREATE INDEX idx_cohorts_owner ON cohorts(owner_perspective_id);
CREATE INDEX idx_cohorts_source
    ON cohorts(source_database_identifier, source_database_fingerprint);
CREATE INDEX idx_cohort_access_cohort
    ON perspective_cohort_access(cohort_id);
CREATE INDEX idx_datasets_perspective ON datasets(perspective_id);
CREATE INDEX idx_dataset_cohorts_cohort ON dataset_cohorts(cohort_id);
CREATE INDEX idx_catalog_samples_cohort ON catalog_samples(cohort_id);
CREATE INDEX idx_sample_groups_dataset ON sample_groups(dataset_id);
CREATE INDEX idx_group_members_sample
    ON sample_group_members(catalog_sample_id);
CREATE INDEX idx_derived_samples_dataset ON derived_samples(dataset_id);
CREATE INDEX idx_derived_runs_sample ON derived_runs(derived_sample_id);
CREATE INDEX idx_derived_input_samples_sample
    ON derived_run_input_samples(catalog_sample_id);
CREATE INDEX idx_derived_parent_runs_parent
    ON derived_run_parent_runs(parent_run_id);
CREATE INDEX idx_derived_alleles_identity
    ON derived_run_alleles(position, ref, alt);
CREATE INDEX idx_derived_sources_allele
    ON derived_allele_sources(derived_allele_id);
CREATE INDEX idx_derived_sources_observed
    ON derived_allele_sources(original_catalog_sample_id);
"""

MIGRATION_2 = """
CREATE UNIQUE INDEX idx_cohorts_source_identifier_unique
    ON cohorts(source_database_identifier)
    WHERE cohort_type = 'source';
"""

MIGRATION_3 = """
ALTER TABLE derived_runs
    ADD COLUMN output_fingerprint TEXT NOT NULL DEFAULT '';

ALTER TABLE derived_run_input_samples
    RENAME TO derived_run_input_samples_v2;

CREATE TABLE derived_run_input_samples (
    run_id TEXT NOT NULL,
    input_clause_index INTEGER NOT NULL CHECK (input_clause_index >= 0),
    catalog_sample_id TEXT NOT NULL,
    source_group_id TEXT,
    group_membership_fingerprint TEXT NOT NULL DEFAULT '',
    sample_fingerprint TEXT NOT NULL DEFAULT '',
    source_database_fingerprint TEXT NOT NULL DEFAULT '',
    input_role TEXT NOT NULL,
    display_order INTEGER NOT NULL DEFAULT 0,
    PRIMARY KEY (run_id, input_clause_index, catalog_sample_id),
    FOREIGN KEY (run_id) REFERENCES derived_runs(id) ON DELETE CASCADE,
    FOREIGN KEY (catalog_sample_id)
        REFERENCES catalog_samples(id) ON DELETE RESTRICT,
    FOREIGN KEY (source_group_id)
        REFERENCES sample_groups(id) ON DELETE SET NULL
);

INSERT INTO derived_run_input_samples(
    run_id,
    input_clause_index,
    catalog_sample_id,
    source_group_id,
    input_role,
    display_order
)
SELECT
    run_id,
    rowid,
    catalog_sample_id,
    source_group_id,
    input_role,
    display_order
FROM derived_run_input_samples_v2;

DROP TABLE derived_run_input_samples_v2;

ALTER TABLE derived_run_parent_runs
    RENAME TO derived_run_parent_runs_v2;

CREATE TABLE derived_run_parent_runs (
    run_id TEXT NOT NULL,
    input_clause_index INTEGER NOT NULL CHECK (input_clause_index >= 0),
    parent_run_id TEXT NOT NULL,
    source_group_id TEXT,
    group_membership_fingerprint TEXT NOT NULL DEFAULT '',
    input_role TEXT NOT NULL,
    display_order INTEGER NOT NULL DEFAULT 0,
    CHECK (run_id != parent_run_id),
    PRIMARY KEY (run_id, input_clause_index, parent_run_id),
    FOREIGN KEY (run_id) REFERENCES derived_runs(id) ON DELETE CASCADE,
    FOREIGN KEY (parent_run_id)
        REFERENCES derived_runs(id) ON DELETE RESTRICT,
    FOREIGN KEY (source_group_id)
        REFERENCES sample_groups(id) ON DELETE SET NULL
);

INSERT INTO derived_run_parent_runs(
    run_id,
    input_clause_index,
    parent_run_id,
    input_role,
    display_order
)
SELECT
    run_id,
    rowid,
    parent_run_id,
    input_role,
    display_order
FROM derived_run_parent_runs_v2;

DROP TABLE derived_run_parent_runs_v2;

CREATE INDEX idx_derived_input_samples_sample_v3
    ON derived_run_input_samples(catalog_sample_id);
CREATE INDEX idx_derived_input_samples_group_v3
    ON derived_run_input_samples(source_group_id);
CREATE INDEX idx_derived_parent_runs_parent_v3
    ON derived_run_parent_runs(parent_run_id);
CREATE INDEX idx_derived_parent_runs_group_v3
    ON derived_run_parent_runs(source_group_id);
"""


MIGRATIONS = {
    1: MIGRATION_1,
    2: MIGRATION_2,
    3: MIGRATION_3,
}


def catalog_user_version(connection: sqlite3.Connection) -> int:
    return int(connection.execute("PRAGMA user_version").fetchone()[0])


def migrate_catalog(connection: sqlite3.Connection) -> int:
    """Migrate a writable catalog connection to the current schema."""
    connection.execute("PRAGMA foreign_keys = ON")
    current_version = catalog_user_version(connection)
    if current_version > CATALOG_SCHEMA_VERSION:
        raise CatalogMigrationError(
            "Catalog schema version "
            f"{current_version} is newer than supported version "
            f"{CATALOG_SCHEMA_VERSION}."
        )

    for target_version in range(
        current_version + 1,
        CATALOG_SCHEMA_VERSION + 1,
    ):
        migration = MIGRATIONS.get(target_version)
        if migration is None:
            raise CatalogMigrationError(
                f"No migration is available for catalog version {target_version}."
            )
        try:
            connection.executescript(
                "BEGIN IMMEDIATE;\n"
                + migration
                + f"\nPRAGMA user_version = {target_version};\n"
                + "COMMIT;\n"
            )
        except sqlite3.Error as exc:
            if connection.in_transaction:
                connection.execute("ROLLBACK")
            raise CatalogMigrationError(
                f"Catalog migration {target_version} failed: {exc}"
            ) from exc

    return catalog_user_version(connection)


def inspect_catalog_database(path: str | Path) -> DatabaseSchemaReport:
    """Return the structural schema report for an existing catalog."""
    return inspect_database_schema(
        path,
        database_kind="analysis catalog",
        required_schema=CATALOG_REQUIRED_SCHEMA,
    )
