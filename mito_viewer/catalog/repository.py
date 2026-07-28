"""Transactional repository for the writable analysis catalog."""

from __future__ import annotations

import json
import sqlite3
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator

from .models import (
    AccessLevel,
    CatalogSample,
    Cohort,
    CohortType,
    Dataset,
    SampleGroup,
    SampleGroupSnapshot,
    SampleType,
    StudyPerspective,
    Visibility,
    sample_membership_fingerprint,
)
from .schema import (
    CATALOG_SCHEMA_VERSION,
    catalog_user_version,
    inspect_catalog_database,
    migrate_catalog,
)


class CatalogNotFoundError(LookupError):
    """Raised when a requested catalog entity does not exist."""


class CatalogAccessError(PermissionError):
    """Raised when one perspective cannot perform a catalog operation."""


class CatalogRegistrationConflictError(RuntimeError):
    """Raised when a source registration conflicts with catalog identity."""


class CatalogConcurrencyError(RuntimeError):
    """Raised when editable catalog state changed since it was loaded."""


def _new_id() -> str:
    return str(uuid.uuid4())


def _utc_now() -> str:
    return (
        datetime.now(timezone.utc)
        .isoformat(timespec="microseconds")
        .replace("+00:00", "Z")
    )


def _json_dump(value: dict) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    )


def _json_load(value: str) -> dict:
    loaded = json.loads(value or "{}")
    if not isinstance(loaded, dict):
        raise ValueError("Catalog JSON fields must contain objects.")
    return loaded


class CatalogRepository:
    """Own a schema-validated, writable SQLite catalog connection."""

    def __init__(
        self,
        path: str | Path,
        *,
        create: bool = True,
        validate: bool = True,
    ) -> None:
        self.path = Path(path).resolve()
        if not create and not self.path.is_file():
            raise FileNotFoundError(f"Catalog database not found: {self.path}")
        if create:
            self.path.parent.mkdir(parents=True, exist_ok=True)

        self.connection = sqlite3.connect(
            self.path,
            isolation_level=None,
            timeout=5,
        )
        self.connection.row_factory = sqlite3.Row
        self._savepoint_counter = 0
        self._closed = False

        try:
            self.connection.execute("PRAGMA foreign_keys = ON")
            self.connection.execute("PRAGMA busy_timeout = 5000")
            migrate_catalog(self.connection)
            if validate:
                report = inspect_catalog_database(self.path)
                report.require_valid()
                if report.user_version != CATALOG_SCHEMA_VERSION:
                    raise ValueError(
                        "Catalog schema version does not match the "
                        "repository version."
                    )
            self.connection.execute("PRAGMA journal_mode = WAL")
            self.connection.execute("PRAGMA synchronous = NORMAL")
        except Exception:
            self.connection.close()
            self._closed = True
            raise

    @classmethod
    def open(
        cls,
        path: str | Path,
        *,
        create: bool = True,
        validate: bool = True,
    ) -> "CatalogRepository":
        return cls(path, create=create, validate=validate)

    def __enter__(self) -> "CatalogRepository":
        return self

    def __exit__(self, exc_type, exc, traceback) -> None:
        self.close()

    def close(self) -> None:
        if not self._closed:
            self.connection.close()
            self._closed = True

    @property
    def schema_version(self) -> int:
        return catalog_user_version(self.connection)

    @contextmanager
    def transaction(self) -> Iterator["CatalogRepository"]:
        """Commit a unit of work or roll it back, including when nested."""
        if self._closed:
            raise RuntimeError("Catalog repository is closed.")

        nested = self.connection.in_transaction
        savepoint = ""
        if nested:
            self._savepoint_counter += 1
            savepoint = f"catalog_savepoint_{self._savepoint_counter}"
            self.connection.execute(f"SAVEPOINT {savepoint}")
        else:
            self.connection.execute("BEGIN IMMEDIATE")

        try:
            yield self
        except Exception:
            if nested:
                self.connection.execute(f"ROLLBACK TO {savepoint}")
                self.connection.execute(f"RELEASE {savepoint}")
            else:
                self.connection.execute("ROLLBACK")
            raise
        else:
            if nested:
                self.connection.execute(f"RELEASE {savepoint}")
            else:
                self.connection.execute("COMMIT")

    def create_perspective(
        self,
        name: str,
        *,
        perspective_id: str | None = None,
        description: str = "",
        visibility: Visibility = Visibility.PRIVATE,
    ) -> StudyPerspective:
        now = _utc_now()
        perspective = StudyPerspective(
            id=perspective_id or _new_id(),
            name=name,
            description=description,
            visibility=visibility,
            created_at=now,
            updated_at=now,
        )
        with self.transaction():
            self.connection.execute(
                """
                INSERT INTO study_perspectives(
                    id,
                    name,
                    description,
                    visibility,
                    created_at,
                    updated_at
                )
                VALUES (?, ?, ?, ?, ?, ?)
                """,
                (
                    perspective.id,
                    perspective.name,
                    perspective.description,
                    perspective.visibility.value,
                    perspective.created_at,
                    perspective.updated_at,
                ),
            )
        return perspective

    def get_perspective(
        self,
        perspective_id: str,
    ) -> StudyPerspective | None:
        row = self.connection.execute(
            """
            SELECT id, name, description, visibility, created_at, updated_at
            FROM study_perspectives
            WHERE id = ?
            """,
            (perspective_id,),
        ).fetchone()
        return self._perspective_from_row(row) if row else None

    def list_perspectives(self) -> list[StudyPerspective]:
        rows = self.connection.execute(
            """
            SELECT id, name, description, visibility, created_at, updated_at
            FROM study_perspectives
            ORDER BY name COLLATE NOCASE, id
            """
        ).fetchall()
        return [self._perspective_from_row(row) for row in rows]

    def set_perspective_credential(
        self,
        perspective_id: str,
        password_hash: str,
        *,
        credential_version: int = 1,
    ) -> None:
        """Store an already-hashed credential; plaintext is never accepted."""
        password_hash = str(password_hash or "").strip()
        if not password_hash:
            raise ValueError("Password hash cannot be empty.")
        if credential_version < 1:
            raise ValueError("Credential version must be positive.")
        self._require_perspective(perspective_id)
        with self.transaction():
            self.connection.execute(
                """
                INSERT INTO perspective_credentials(
                    perspective_id,
                    password_hash,
                    credential_version,
                    updated_at
                )
                VALUES (?, ?, ?, ?)
                ON CONFLICT(perspective_id) DO UPDATE SET
                    password_hash = excluded.password_hash,
                    credential_version = excluded.credential_version,
                    updated_at = excluded.updated_at
                """,
                (
                    perspective_id,
                    password_hash,
                    credential_version,
                    _utc_now(),
                ),
            )

    def create_cohort(
        self,
        name: str,
        cohort_type: CohortType,
        *,
        cohort_id: str | None = None,
        owner_perspective_id: str | None = None,
        description: str = "",
        source_database_identifier: str = "",
        source_database_fingerprint: str = "",
        reference_id: str = "",
        normalization_version: str = "",
        provenance: dict | None = None,
        visibility: Visibility = Visibility.PRIVATE,
    ) -> Cohort:
        now = _utc_now()
        cohort = Cohort(
            id=cohort_id or _new_id(),
            name=name,
            cohort_type=cohort_type,
            owner_perspective_id=owner_perspective_id,
            description=description,
            source_database_identifier=source_database_identifier,
            source_database_fingerprint=source_database_fingerprint,
            reference_id=reference_id,
            normalization_version=normalization_version,
            provenance=provenance or {},
            visibility=visibility,
            created_at=now,
            updated_at=now,
        )
        if cohort.owner_perspective_id is not None:
            self._require_perspective(cohort.owner_perspective_id)

        with self.transaction():
            self.connection.execute(
                """
                INSERT INTO cohorts(
                    id,
                    owner_perspective_id,
                    name,
                    description,
                    cohort_type,
                    source_database_identifier,
                    source_database_fingerprint,
                    reference_id,
                    normalization_version,
                    provenance_json,
                    visibility,
                    created_at,
                    updated_at
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    cohort.id,
                    cohort.owner_perspective_id,
                    cohort.name,
                    cohort.description,
                    cohort.cohort_type.value,
                    cohort.source_database_identifier,
                    cohort.source_database_fingerprint,
                    cohort.reference_id,
                    cohort.normalization_version,
                    _json_dump(cohort.provenance),
                    cohort.visibility.value,
                    cohort.created_at,
                    cohort.updated_at,
                ),
            )
            if cohort.owner_perspective_id is not None:
                self._upsert_cohort_access(
                    cohort.owner_perspective_id,
                    cohort.id,
                    AccessLevel.MANAGE,
                )
        return cohort

    def get_cohort(self, cohort_id: str) -> Cohort | None:
        row = self.connection.execute(
            """
            SELECT
                id,
                owner_perspective_id,
                name,
                description,
                cohort_type,
                source_database_identifier,
                source_database_fingerprint,
                reference_id,
                normalization_version,
                provenance_json,
                visibility,
                created_at,
                updated_at
            FROM cohorts
            WHERE id = ?
            """,
            (cohort_id,),
        ).fetchone()
        return self._cohort_from_row(row) if row else None

    def get_source_cohort(
        self,
        source_database_identifier: str,
    ) -> Cohort | None:
        row = self.connection.execute(
            """
            SELECT
                id,
                owner_perspective_id,
                name,
                description,
                cohort_type,
                source_database_identifier,
                source_database_fingerprint,
                reference_id,
                normalization_version,
                provenance_json,
                visibility,
                created_at,
                updated_at
            FROM cohorts
            WHERE
                cohort_type = 'source'
                AND source_database_identifier = ?
            """,
            (str(source_database_identifier).strip(),),
        ).fetchone()
        return self._cohort_from_row(row) if row else None

    def register_source_snapshot(
        self,
        cohort: Cohort,
        samples: list[CatalogSample] | tuple[CatalogSample, ...],
    ) -> bool:
        """Atomically register one immutable view of a source database.

        Returns ``True`` for a new registration and ``False`` when the exact
        source snapshot was already present. A changed source is rejected
        rather than silently rewriting identities used by later analyses.
        """
        if cohort.cohort_type is not CohortType.SOURCE:
            raise ValueError("Source registration requires a source cohort.")

        samples = tuple(samples)
        sample_ids: set[str] = set()
        source_sample_ids: set[str] = set()
        for sample in samples:
            if sample.sample_type is not SampleType.OBSERVED:
                raise ValueError(
                    "Source registration accepts only observed samples."
                )
            if sample.cohort_id != cohort.id:
                raise ValueError(
                    "Registered sample cohort does not match source cohort."
                )
            if sample.id in sample_ids:
                raise ValueError(f"Duplicate catalog sample ID: {sample.id}")
            if sample.source_sample_id in source_sample_ids:
                raise ValueError(
                    "Duplicate source sample ID: "
                    f"{sample.source_sample_id}"
                )
            sample_ids.add(sample.id)
            source_sample_ids.add(sample.source_sample_id)

        with self.transaction():
            existing_by_id = self.get_cohort(cohort.id)
            existing_by_source = self.get_source_cohort(
                cohort.source_database_identifier
            )
            if (
                existing_by_id is not None
                and existing_by_id.cohort_type is not CohortType.SOURCE
            ):
                raise CatalogRegistrationConflictError(
                    f"Catalog ID {cohort.id!r} already belongs to a "
                    "non-source cohort."
                )
            if (
                existing_by_source is not None
                and existing_by_source.id != cohort.id
            ):
                raise CatalogRegistrationConflictError(
                    "Source database identifier is already registered under "
                    f"cohort {existing_by_source.id!r}."
                )
            if existing_by_id is not None:
                self._verify_source_snapshot(
                    existing_by_id,
                    cohort,
                    samples,
                )
                return False

            self.connection.execute(
                """
                INSERT INTO cohorts(
                    id,
                    owner_perspective_id,
                    name,
                    description,
                    cohort_type,
                    source_database_identifier,
                    source_database_fingerprint,
                    reference_id,
                    normalization_version,
                    provenance_json,
                    visibility,
                    created_at,
                    updated_at
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    cohort.id,
                    cohort.owner_perspective_id,
                    cohort.name,
                    cohort.description,
                    cohort.cohort_type.value,
                    cohort.source_database_identifier,
                    cohort.source_database_fingerprint,
                    cohort.reference_id,
                    cohort.normalization_version,
                    _json_dump(cohort.provenance),
                    cohort.visibility.value,
                    cohort.created_at,
                    cohort.updated_at,
                ),
            )
            self.connection.executemany(
                """
                INSERT INTO catalog_samples(
                    id,
                    cohort_id,
                    sample_type,
                    source_sample_id,
                    display_label,
                    source_fingerprint,
                    metadata_json,
                    created_at
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                [
                    (
                        sample.id,
                        sample.cohort_id,
                        sample.sample_type.value,
                        sample.source_sample_id,
                        sample.display_label,
                        sample.source_fingerprint,
                        _json_dump(sample.metadata),
                        sample.created_at,
                    )
                    for sample in samples
                ],
            )
        return True

    def list_cohorts(
        self,
        *,
        perspective_id: str | None = None,
        usable_only: bool = False,
    ) -> list[Cohort]:
        rows = self.connection.execute(
            """
            SELECT
                id,
                owner_perspective_id,
                name,
                description,
                cohort_type,
                source_database_identifier,
                source_database_fingerprint,
                reference_id,
                normalization_version,
                provenance_json,
                visibility,
                created_at,
                updated_at
            FROM cohorts
            ORDER BY name COLLATE NOCASE, id
            """
        ).fetchall()
        cohorts = [self._cohort_from_row(row) for row in rows]
        if perspective_id is None:
            return cohorts
        self._require_perspective(perspective_id)
        required = AccessLevel.USE if usable_only else AccessLevel.VIEW
        return [
            cohort
            for cohort in cohorts
            if self.has_cohort_access(
                perspective_id,
                cohort.id,
                required=required,
            )
        ]

    def grant_cohort_access(
        self,
        cohort_id: str,
        perspective_id: str,
        access_level: AccessLevel,
        *,
        acting_perspective_id: str,
    ) -> None:
        self._require_perspective(perspective_id)
        self._require_cohort(cohort_id)
        self.require_cohort_access(
            acting_perspective_id,
            cohort_id,
            required=AccessLevel.MANAGE,
        )
        with self.transaction():
            self._upsert_cohort_access(
                perspective_id,
                cohort_id,
                AccessLevel(access_level),
            )

    def cohort_access_level(
        self,
        perspective_id: str,
        cohort_id: str,
    ) -> AccessLevel | None:
        cohort = self._require_cohort(cohort_id)
        if cohort.owner_perspective_id == perspective_id:
            return AccessLevel.MANAGE
        row = self.connection.execute(
            """
            SELECT access_level
            FROM perspective_cohort_access
            WHERE perspective_id = ? AND cohort_id = ?
            """,
            (perspective_id, cohort_id),
        ).fetchone()
        if row:
            return AccessLevel(row["access_level"])
        if cohort.visibility is Visibility.PUBLIC:
            return AccessLevel.USE
        return None

    def has_cohort_access(
        self,
        perspective_id: str,
        cohort_id: str,
        *,
        required: AccessLevel = AccessLevel.VIEW,
    ) -> bool:
        self._require_perspective(perspective_id)
        level = self.cohort_access_level(perspective_id, cohort_id)
        if level is None:
            return False
        ranks = {
            AccessLevel.VIEW: 1,
            AccessLevel.USE: 2,
            AccessLevel.MANAGE: 3,
        }
        return ranks[level] >= ranks[AccessLevel(required)]

    def require_cohort_access(
        self,
        perspective_id: str,
        cohort_id: str,
        *,
        required: AccessLevel = AccessLevel.VIEW,
    ) -> None:
        if not self.has_cohort_access(
            perspective_id,
            cohort_id,
            required=required,
        ):
            raise CatalogAccessError(
                f"Perspective {perspective_id!r} does not have "
                f"{AccessLevel(required).value} access to cohort "
                f"{cohort_id!r}."
            )

    def create_dataset(
        self,
        perspective_id: str,
        name: str,
        *,
        dataset_id: str | None = None,
        description: str = "",
        derived_results_cohort_id: str | None = None,
        visibility: Visibility = Visibility.PRIVATE,
    ) -> Dataset:
        self._require_perspective(perspective_id)
        if derived_results_cohort_id is not None:
            self._validate_derived_results_cohort(
                perspective_id,
                derived_results_cohort_id,
            )
        now = _utc_now()
        dataset = Dataset(
            id=dataset_id or _new_id(),
            perspective_id=perspective_id,
            name=name,
            description=description,
            derived_results_cohort_id=derived_results_cohort_id,
            visibility=visibility,
            created_at=now,
            updated_at=now,
        )
        with self.transaction():
            self.connection.execute(
                """
                INSERT INTO datasets(
                    id,
                    perspective_id,
                    name,
                    description,
                    derived_results_cohort_id,
                    visibility,
                    created_at,
                    updated_at
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    dataset.id,
                    dataset.perspective_id,
                    dataset.name,
                    dataset.description,
                    dataset.derived_results_cohort_id,
                    dataset.visibility.value,
                    dataset.created_at,
                    dataset.updated_at,
                ),
            )
        return dataset

    def get_dataset(self, dataset_id: str) -> Dataset | None:
        row = self.connection.execute(
            """
            SELECT
                id,
                perspective_id,
                name,
                description,
                derived_results_cohort_id,
                visibility,
                created_at,
                updated_at
            FROM datasets
            WHERE id = ?
            """,
            (dataset_id,),
        ).fetchone()
        return self._dataset_from_row(row) if row else None

    def list_datasets(self, perspective_id: str) -> list[Dataset]:
        self._require_perspective(perspective_id)
        rows = self.connection.execute(
            """
            SELECT
                id,
                perspective_id,
                name,
                description,
                derived_results_cohort_id,
                visibility,
                created_at,
                updated_at
            FROM datasets
            WHERE perspective_id = ?
            ORDER BY name COLLATE NOCASE, id
            """,
            (perspective_id,),
        ).fetchall()
        return [self._dataset_from_row(row) for row in rows]

    def set_derived_results_cohort(
        self,
        dataset_id: str,
        cohort_id: str,
        *,
        acting_perspective_id: str,
    ) -> Dataset:
        dataset = self._require_dataset(dataset_id)
        self._require_dataset_owner(dataset, acting_perspective_id)
        self._validate_derived_results_cohort(
            dataset.perspective_id,
            cohort_id,
        )
        if (
            dataset.derived_results_cohort_id is not None
            and dataset.derived_results_cohort_id != cohort_id
        ):
            raise ValueError(
                f"Dataset {dataset_id!r} already has a Derived Results "
                "cohort."
            )
        with self.transaction():
            self.connection.execute(
                """
                UPDATE datasets
                SET derived_results_cohort_id = ?, updated_at = ?
                WHERE id = ?
                """,
                (cohort_id, _utc_now(), dataset_id),
            )
        return self._require_dataset(dataset_id)

    def add_cohort_to_dataset(
        self,
        dataset_id: str,
        cohort_id: str,
        *,
        acting_perspective_id: str,
        display_order: int = 0,
    ) -> None:
        dataset = self._require_dataset(dataset_id)
        self._require_dataset_owner(dataset, acting_perspective_id)
        self.require_cohort_access(
            dataset.perspective_id,
            cohort_id,
            required=AccessLevel.USE,
        )
        with self.transaction():
            self.connection.execute(
                """
                INSERT INTO dataset_cohorts(
                    dataset_id,
                    cohort_id,
                    display_order,
                    added_at
                )
                VALUES (?, ?, ?, ?)
                ON CONFLICT(dataset_id, cohort_id) DO UPDATE SET
                    display_order = excluded.display_order
                """,
                (dataset_id, cohort_id, int(display_order), _utc_now()),
            )

    def list_dataset_cohorts(self, dataset_id: str) -> list[Cohort]:
        self._require_dataset(dataset_id)
        rows = self.connection.execute(
            """
            SELECT
                cohorts.id,
                cohorts.owner_perspective_id,
                cohorts.name,
                cohorts.description,
                cohorts.cohort_type,
                cohorts.source_database_identifier,
                cohorts.source_database_fingerprint,
                cohorts.reference_id,
                cohorts.normalization_version,
                cohorts.provenance_json,
                cohorts.visibility,
                cohorts.created_at,
                cohorts.updated_at
            FROM dataset_cohorts
            JOIN cohorts ON cohorts.id = dataset_cohorts.cohort_id
            WHERE dataset_cohorts.dataset_id = ?
            ORDER BY
                dataset_cohorts.display_order,
                cohorts.name COLLATE NOCASE,
                cohorts.id
            """,
            (dataset_id,),
        ).fetchall()
        return [self._cohort_from_row(row) for row in rows]

    def create_catalog_sample(
        self,
        cohort_id: str,
        sample_type: SampleType,
        display_label: str,
        *,
        acting_perspective_id: str,
        catalog_sample_id: str | None = None,
        source_sample_id: str | None = None,
        source_fingerprint: str = "",
        metadata: dict | None = None,
    ) -> CatalogSample:
        cohort = self._require_cohort(cohort_id)
        self.require_cohort_access(
            acting_perspective_id,
            cohort_id,
            required=AccessLevel.MANAGE,
        )
        sample_type = SampleType(sample_type)
        expected_cohort_type = (
            CohortType.SOURCE
            if sample_type is SampleType.OBSERVED
            else CohortType.DERIVED
        )
        if cohort.cohort_type is not expected_cohort_type:
            raise ValueError(
                f"{sample_type.value.title()} samples require a "
                f"{expected_cohort_type.value} cohort."
            )

        sample = CatalogSample(
            id=catalog_sample_id or _new_id(),
            cohort_id=cohort_id,
            sample_type=sample_type,
            display_label=display_label,
            source_sample_id=source_sample_id,
            source_fingerprint=source_fingerprint,
            metadata=metadata or {},
            created_at=_utc_now(),
        )
        with self.transaction():
            self.connection.execute(
                """
                INSERT INTO catalog_samples(
                    id,
                    cohort_id,
                    sample_type,
                    source_sample_id,
                    display_label,
                    source_fingerprint,
                    metadata_json,
                    created_at
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    sample.id,
                    sample.cohort_id,
                    sample.sample_type.value,
                    sample.source_sample_id,
                    sample.display_label,
                    sample.source_fingerprint,
                    _json_dump(sample.metadata),
                    sample.created_at,
                ),
            )
        return sample

    def get_catalog_sample(
        self,
        catalog_sample_id: str,
    ) -> CatalogSample | None:
        row = self.connection.execute(
            """
            SELECT
                id,
                cohort_id,
                sample_type,
                source_sample_id,
                display_label,
                source_fingerprint,
                metadata_json,
                created_at
            FROM catalog_samples
            WHERE id = ?
            """,
            (catalog_sample_id,),
        ).fetchone()
        return self._sample_from_row(row) if row else None

    def list_catalog_samples(self, cohort_id: str) -> list[CatalogSample]:
        self._require_cohort(cohort_id)
        rows = self.connection.execute(
            """
            SELECT
                id,
                cohort_id,
                sample_type,
                source_sample_id,
                display_label,
                source_fingerprint,
                metadata_json,
                created_at
            FROM catalog_samples
            WHERE cohort_id = ?
            ORDER BY display_label COLLATE NOCASE, id
            """,
            (cohort_id,),
        ).fetchall()
        return [self._sample_from_row(row) for row in rows]

    def list_dataset_samples(
        self,
        dataset_id: str,
        *,
        acting_perspective_id: str,
    ) -> list[CatalogSample]:
        dataset = self._require_dataset(dataset_id)
        self._require_dataset_owner(dataset, acting_perspective_id)
        rows = self.connection.execute(
            """
            SELECT
                catalog_samples.id,
                catalog_samples.cohort_id,
                catalog_samples.sample_type,
                catalog_samples.source_sample_id,
                catalog_samples.display_label,
                catalog_samples.source_fingerprint,
                catalog_samples.metadata_json,
                catalog_samples.created_at
            FROM dataset_cohorts
            JOIN catalog_samples
                ON catalog_samples.cohort_id = dataset_cohorts.cohort_id
            WHERE dataset_cohorts.dataset_id = ?
            ORDER BY
                dataset_cohorts.display_order,
                catalog_samples.display_label COLLATE NOCASE,
                catalog_samples.id
            """,
            (dataset_id,),
        ).fetchall()
        return [self._sample_from_row(row) for row in rows]

    def create_sample_group(
        self,
        perspective_id: str,
        dataset_id: str,
        name: str,
        sample_ids: list[str] | tuple[str, ...],
        *,
        group_id: str | None = None,
        description: str = "",
        visibility: Visibility = Visibility.PRIVATE,
    ) -> SampleGroup:
        sample_ids = self._normalize_sample_ids(sample_ids)
        now = _utc_now()
        group = SampleGroup(
            id=group_id or _new_id(),
            perspective_id=perspective_id,
            dataset_id=dataset_id,
            name=name,
            description=description,
            visibility=visibility,
            membership_fingerprint=sample_membership_fingerprint(sample_ids),
            created_at=now,
            updated_at=now,
        )

        with self.transaction():
            dataset = self._require_dataset(dataset_id)
            self._require_dataset_owner(dataset, perspective_id)
            self._validate_samples_available_in_dataset(
                dataset_id,
                sample_ids,
            )
            self.connection.execute(
                """
                INSERT INTO sample_groups(
                    id,
                    perspective_id,
                    dataset_id,
                    name,
                    description,
                    visibility,
                    membership_fingerprint,
                    created_at,
                    updated_at
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    group.id,
                    group.perspective_id,
                    group.dataset_id,
                    group.name,
                    group.description,
                    group.visibility.value,
                    group.membership_fingerprint,
                    group.created_at,
                    group.updated_at,
                ),
            )
            self._insert_sample_group_members(group.id, sample_ids, now)
        return group

    def get_sample_group(
        self,
        group_id: str,
        *,
        acting_perspective_id: str,
    ) -> SampleGroup | None:
        group = self._find_sample_group(group_id)
        if group is None:
            return None
        self._require_sample_group_owner(group, acting_perspective_id)
        return group

    def list_sample_groups(
        self,
        dataset_id: str,
        *,
        acting_perspective_id: str,
    ) -> list[SampleGroup]:
        dataset = self._require_dataset(dataset_id)
        self._require_dataset_owner(dataset, acting_perspective_id)
        rows = self.connection.execute(
            """
            SELECT
                id,
                perspective_id,
                dataset_id,
                name,
                description,
                visibility,
                membership_fingerprint,
                created_at,
                updated_at
            FROM sample_groups
            WHERE dataset_id = ?
            ORDER BY name COLLATE NOCASE, id
            """,
            (dataset_id,),
        ).fetchall()
        return [self._sample_group_from_row(row) for row in rows]

    def sample_group_members(
        self,
        group_id: str,
        *,
        acting_perspective_id: str,
    ) -> list[CatalogSample]:
        group = self._require_sample_group(group_id)
        self._require_sample_group_owner(group, acting_perspective_id)
        rows = self.connection.execute(
            """
            SELECT
                catalog_samples.id,
                catalog_samples.cohort_id,
                catalog_samples.sample_type,
                catalog_samples.source_sample_id,
                catalog_samples.display_label,
                catalog_samples.source_fingerprint,
                catalog_samples.metadata_json,
                catalog_samples.created_at
            FROM sample_group_members
            JOIN catalog_samples
                ON catalog_samples.id =
                    sample_group_members.catalog_sample_id
            WHERE sample_group_members.group_id = ?
            ORDER BY
                sample_group_members.display_order,
                catalog_samples.id
            """,
            (group_id,),
        ).fetchall()
        return [self._sample_from_row(row) for row in rows]

    def update_sample_group(
        self,
        group_id: str,
        *,
        acting_perspective_id: str,
        name: str | None = None,
        description: str | None = None,
        visibility: Visibility | None = None,
        sample_ids: list[str] | tuple[str, ...] | None = None,
        expected_membership_fingerprint: str | None = None,
    ) -> SampleGroup:
        normalized_sample_ids = (
            None
            if sample_ids is None
            else self._normalize_sample_ids(sample_ids)
        )
        with self.transaction():
            current = self._require_sample_group(group_id)
            self._require_sample_group_owner(
                current,
                acting_perspective_id,
            )
            if (
                expected_membership_fingerprint is not None
                and current.membership_fingerprint
                != expected_membership_fingerprint
            ):
                raise CatalogConcurrencyError(
                    f"Sample group {group_id!r} changed after it was loaded."
                )

            if normalized_sample_ids is None:
                membership_fingerprint = current.membership_fingerprint
            else:
                self._validate_samples_available_in_dataset(
                    current.dataset_id,
                    normalized_sample_ids,
                )
                membership_fingerprint = sample_membership_fingerprint(
                    normalized_sample_ids
                )
            updated_at = _utc_now()
            updated = SampleGroup(
                id=current.id,
                perspective_id=current.perspective_id,
                dataset_id=current.dataset_id,
                name=current.name if name is None else name,
                description=(
                    current.description
                    if description is None
                    else description
                ),
                visibility=(
                    current.visibility
                    if visibility is None
                    else visibility
                ),
                membership_fingerprint=membership_fingerprint,
                created_at=current.created_at,
                updated_at=updated_at,
            )
            self.connection.execute(
                """
                UPDATE sample_groups
                SET
                    name = ?,
                    description = ?,
                    visibility = ?,
                    membership_fingerprint = ?,
                    updated_at = ?
                WHERE id = ?
                """,
                (
                    updated.name,
                    updated.description,
                    updated.visibility.value,
                    updated.membership_fingerprint,
                    updated.updated_at,
                    updated.id,
                ),
            )
            if normalized_sample_ids is not None:
                self.connection.execute(
                    "DELETE FROM sample_group_members WHERE group_id = ?",
                    (group_id,),
                )
                self._insert_sample_group_members(
                    group_id,
                    normalized_sample_ids,
                    updated_at,
                )
        return updated

    def snapshot_sample_group(
        self,
        group_id: str,
        *,
        acting_perspective_id: str,
    ) -> SampleGroupSnapshot:
        with self.transaction():
            group = self._require_sample_group(group_id)
            self._require_sample_group_owner(
                group,
                acting_perspective_id,
            )
            sample_ids = tuple(
                sample.id
                for sample in self.sample_group_members(
                    group_id,
                    acting_perspective_id=acting_perspective_id,
                )
            )
            actual_fingerprint = sample_membership_fingerprint(sample_ids)
            if actual_fingerprint != group.membership_fingerprint:
                raise RuntimeError(
                    f"Sample group {group_id!r} membership is inconsistent."
                )
            return SampleGroupSnapshot(
                group_id=group.id,
                perspective_id=group.perspective_id,
                dataset_id=group.dataset_id,
                group_name=group.name,
                sample_ids=sample_ids,
                membership_fingerprint=group.membership_fingerprint,
                captured_at=_utc_now(),
            )

    def _require_perspective(
        self,
        perspective_id: str,
    ) -> StudyPerspective:
        perspective = self.get_perspective(perspective_id)
        if perspective is None:
            raise CatalogNotFoundError(
                f"Perspective not found: {perspective_id}"
            )
        return perspective

    def _require_cohort(self, cohort_id: str) -> Cohort:
        cohort = self.get_cohort(cohort_id)
        if cohort is None:
            raise CatalogNotFoundError(f"Cohort not found: {cohort_id}")
        return cohort

    def _require_dataset(self, dataset_id: str) -> Dataset:
        dataset = self.get_dataset(dataset_id)
        if dataset is None:
            raise CatalogNotFoundError(f"Dataset not found: {dataset_id}")
        return dataset

    def _find_sample_group(self, group_id: str) -> SampleGroup | None:
        row = self.connection.execute(
            """
            SELECT
                id,
                perspective_id,
                dataset_id,
                name,
                description,
                visibility,
                membership_fingerprint,
                created_at,
                updated_at
            FROM sample_groups
            WHERE id = ?
            """,
            (group_id,),
        ).fetchone()
        return self._sample_group_from_row(row) if row else None

    def _require_sample_group(self, group_id: str) -> SampleGroup:
        group = self._find_sample_group(group_id)
        if group is None:
            raise CatalogNotFoundError(
                f"Sample group not found: {group_id}"
            )
        return group

    @staticmethod
    def _require_dataset_owner(
        dataset: Dataset,
        acting_perspective_id: str,
    ) -> None:
        if dataset.perspective_id != acting_perspective_id:
            raise CatalogAccessError(
                f"Perspective {acting_perspective_id!r} does not own "
                f"dataset {dataset.id!r}."
            )

    @staticmethod
    def _require_sample_group_owner(
        group: SampleGroup,
        acting_perspective_id: str,
    ) -> None:
        if group.perspective_id != acting_perspective_id:
            raise CatalogAccessError(
                f"Perspective {acting_perspective_id!r} does not own "
                f"sample group {group.id!r}."
            )

    def _validate_samples_available_in_dataset(
        self,
        dataset_id: str,
        sample_ids: tuple[str, ...],
    ) -> None:
        if not sample_ids:
            return
        placeholders = ",".join("?" for _ in sample_ids)
        rows = self.connection.execute(
            f"""
            SELECT catalog_samples.id
            FROM catalog_samples
            JOIN dataset_cohorts
                ON dataset_cohorts.cohort_id =
                    catalog_samples.cohort_id
            WHERE
                dataset_cohorts.dataset_id = ?
                AND catalog_samples.id IN ({placeholders})
            """,
            (dataset_id, *sample_ids),
        ).fetchall()
        available = {row["id"] for row in rows}
        missing = [
            sample_id
            for sample_id in sample_ids
            if sample_id not in available
        ]
        if missing:
            raise ValueError(
                "Samples are not available in dataset "
                f"{dataset_id!r}: {', '.join(missing)}"
            )

    @staticmethod
    def _normalize_sample_ids(
        sample_ids: list[str] | tuple[str, ...],
    ) -> tuple[str, ...]:
        normalized = tuple(
            str(sample_id or "").strip()
            for sample_id in sample_ids
        )
        if any(not sample_id for sample_id in normalized):
            raise ValueError("Catalog sample IDs cannot be empty.")
        if len(set(normalized)) != len(normalized):
            raise ValueError("Sample groups cannot contain duplicate samples.")
        return normalized

    def _insert_sample_group_members(
        self,
        group_id: str,
        sample_ids: tuple[str, ...],
        added_at: str,
    ) -> None:
        self.connection.executemany(
            """
            INSERT INTO sample_group_members(
                group_id,
                catalog_sample_id,
                display_order,
                added_at
            )
            VALUES (?, ?, ?, ?)
            """,
            [
                (group_id, sample_id, index, added_at)
                for index, sample_id in enumerate(sample_ids)
            ],
        )

    def _validate_derived_results_cohort(
        self,
        perspective_id: str,
        cohort_id: str,
    ) -> Cohort:
        cohort = self._require_cohort(cohort_id)
        if cohort.cohort_type is not CohortType.DERIVED:
            raise ValueError("A dataset results cohort must be derived.")
        if cohort.owner_perspective_id != perspective_id:
            raise CatalogAccessError(
                "A dataset's Derived Results cohort must be owned by "
                "the same perspective."
            )
        return cohort

    def _verify_source_snapshot(
        self,
        existing: Cohort,
        requested: Cohort,
        samples: tuple[CatalogSample, ...],
    ) -> None:
        identity_fields = (
            "source_database_identifier",
            "source_database_fingerprint",
            "reference_id",
            "normalization_version",
        )
        for field_name in identity_fields:
            if getattr(existing, field_name) != getattr(requested, field_name):
                raise CatalogRegistrationConflictError(
                    f"Registered source {existing.id!r} has changed "
                    f"{field_name.replace('_', ' ')}."
                )

        existing_samples = self.list_catalog_samples(existing.id)
        existing_by_source = {
            sample.source_sample_id: sample
            for sample in existing_samples
        }
        requested_by_source = {
            sample.source_sample_id: sample
            for sample in samples
        }
        if existing_by_source.keys() != requested_by_source.keys():
            raise CatalogRegistrationConflictError(
                f"Registered source {existing.id!r} has a changed sample set."
            )

        for source_sample_id, requested_sample in requested_by_source.items():
            existing_sample = existing_by_source[source_sample_id]
            if existing_sample.id != requested_sample.id:
                raise CatalogRegistrationConflictError(
                    f"Source sample {source_sample_id!r} resolved to a "
                    "different catalog identity."
                )
            if (
                existing_sample.source_fingerprint
                != requested_sample.source_fingerprint
            ):
                raise CatalogRegistrationConflictError(
                    f"Source sample {source_sample_id!r} has changed content."
                )

    def _upsert_cohort_access(
        self,
        perspective_id: str,
        cohort_id: str,
        access_level: AccessLevel,
    ) -> None:
        self.connection.execute(
            """
            INSERT INTO perspective_cohort_access(
                perspective_id,
                cohort_id,
                access_level,
                granted_at
            )
            VALUES (?, ?, ?, ?)
            ON CONFLICT(perspective_id, cohort_id) DO UPDATE SET
                access_level = excluded.access_level,
                granted_at = excluded.granted_at
            """,
            (
                perspective_id,
                cohort_id,
                access_level.value,
                _utc_now(),
            ),
        )

    @staticmethod
    def _perspective_from_row(row: sqlite3.Row) -> StudyPerspective:
        return StudyPerspective(
            id=row["id"],
            name=row["name"],
            description=row["description"],
            visibility=row["visibility"],
            created_at=row["created_at"],
            updated_at=row["updated_at"],
        )

    @staticmethod
    def _cohort_from_row(row: sqlite3.Row) -> Cohort:
        return Cohort(
            id=row["id"],
            owner_perspective_id=row["owner_perspective_id"],
            name=row["name"],
            description=row["description"],
            cohort_type=row["cohort_type"],
            source_database_identifier=row["source_database_identifier"],
            source_database_fingerprint=row["source_database_fingerprint"],
            reference_id=row["reference_id"],
            normalization_version=row["normalization_version"],
            provenance=_json_load(row["provenance_json"]),
            visibility=row["visibility"],
            created_at=row["created_at"],
            updated_at=row["updated_at"],
        )

    @staticmethod
    def _dataset_from_row(row: sqlite3.Row) -> Dataset:
        return Dataset(
            id=row["id"],
            perspective_id=row["perspective_id"],
            name=row["name"],
            description=row["description"],
            derived_results_cohort_id=row["derived_results_cohort_id"],
            visibility=row["visibility"],
            created_at=row["created_at"],
            updated_at=row["updated_at"],
        )

    @staticmethod
    def _sample_from_row(row: sqlite3.Row) -> CatalogSample:
        return CatalogSample(
            id=row["id"],
            cohort_id=row["cohort_id"],
            sample_type=row["sample_type"],
            source_sample_id=row["source_sample_id"],
            display_label=row["display_label"],
            source_fingerprint=row["source_fingerprint"],
            metadata=_json_load(row["metadata_json"]),
            created_at=row["created_at"],
        )

    @staticmethod
    def _sample_group_from_row(row: sqlite3.Row) -> SampleGroup:
        return SampleGroup(
            id=row["id"],
            perspective_id=row["perspective_id"],
            dataset_id=row["dataset_id"],
            name=row["name"],
            description=row["description"],
            visibility=row["visibility"],
            membership_fingerprint=row["membership_fingerprint"],
            created_at=row["created_at"],
            updated_at=row["updated_at"],
        )
