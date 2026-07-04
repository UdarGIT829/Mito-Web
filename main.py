import argparse
import json
import re
import sqlite3
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

import vcf_parser


MITOCHONDRIAL_LENGTH = 16569
DEFAULT_MITO_REFERENCE_PATH = (
    Path(__file__).resolve().parent / "reference" / "hg38_chrM.fa"
)


class mitochondrial_base_positions(Iterator):
    def __init__(self):
        self._positions = iter(range(1, MITOCHONDRIAL_LENGTH + 1))

    def __iter__(self):
        return self

    def __next__(self):
        return next(self._positions)


@dataclass(frozen=True)
class MutationAllele:
    """Comparable mutation allele used for sample set operations."""

    position: int
    alt: str
    ref: str = field(compare=False)
    af: float | str | None = field(default=None, compare=False)
    filter: str = field(default="", compare=False)
    mutation: vcf_parser.VCFMutation | None = field(
        default=None,
        compare=False,
        repr=False,
    )


@dataclass
class Sample:
    sample_id: str
    population: list[str]
    source_path: Path
    mutations: list[vcf_parser.VCFMutation] = field(default_factory=list)

    def __contains__(self, tags):
        """Return True when this sample has the given population tag(s)."""
        if isinstance(tags, str):
            return tags in self.population

        return all(tag in self.population for tag in tags)

    def __iter__(self):
        """Iterate over population tags."""
        return iter(self.population)

    @property
    def population_key(self):
        """Return a stable serialized population tag key."""
        return "|".join(self.population)

    @property
    def label(self):
        """Return a compact subject/population label."""
        if not self.population:
            return self.sample_id
        return f"{self.sample_id}_{'_'.join(self.population)}"

    def has_any(self, tags):
        """Return True when this sample has at least one of the given tags."""
        return any(tag in self for tag in tags)

    def has_all(self, tags):
        """Return True when this sample has every given tag."""
        return tags in self

    def is_subject(self, sample_id):
        """Return True when this sample belongs to the given subject."""
        return self.sample_id == sample_id

    @property
    def mutation_alleles(self):
        """Return comparable mutation alleles for this sample."""
        alleles = set()
        for mutation in self.mutations:
            for index, alt in enumerate(mutation.alts):
                af = mutation.afs[index] if index < len(mutation.afs) else None
                alleles.add(MutationAllele(
                    position=mutation.position,
                    ref=mutation.ref,
                    alt=alt,
                    af=af,
                    filter=mutation.filter,
                    mutation=mutation,
                ))
        return alleles

    def __and__(self, other):
        """Return mutation alleles common to both samples."""
        return self.mutation_alleles & other.mutation_alleles

    def __or__(self, other):
        """Return mutation alleles present in either sample."""
        return self.mutation_alleles | other.mutation_alleles

    def __sub__(self, other):
        """Return mutation alleles present in this sample but not the other."""
        return self.mutation_alleles - other.mutation_alleles

    def __xor__(self, other):
        """Return mutation alleles that differ between two samples."""
        return self.mutation_alleles ^ other.mutation_alleles

    def common_mutations(self, *others):
        """Return mutation alleles common to this sample and all others."""
        common = self.mutation_alleles
        for other in others:
            common &= other.mutation_alleles
        return common

    def different_mutations(self, *others):
        """Return mutation alleles not shared by every provided sample."""
        samples = (self, *others)
        if not samples:
            return set()

        union = set()
        common = samples[0].mutation_alleles
        for sample in samples:
            union |= sample.mutation_alleles
            common &= sample.mutation_alleles

        return union - common


def is_vcf_file(path):
    """Return True for plain or gzipped VCF files."""
    name = path.name
    return name.endswith(".vcf") or name.endswith(".vcf.gz")


def strip_vcf_suffix(path):
    """Return a VCF filename without .vcf or .vcf.gz."""
    name = Path(path).name
    if name.endswith(".vcf.gz"):
        return name[:-len(".vcf.gz")]
    if name.endswith(".vcf"):
        return name[:-len(".vcf")]
    return Path(path).stem


def default_sample_id(path):
    """Return a generic sample ID inferred from a VCF filename."""
    stem = strip_vcf_suffix(path)
    return re.sub(r"_S\d+$", "", stem)


def parse_sample_filename(path, input_dir=None):
    """Return generic sample metadata inferred from a VCF path."""
    sample_id = default_sample_id(path)
    population = []
    if input_dir is not None:
        try:
            relative_parent = Path(path).parent.relative_to(Path(input_dir))
        except ValueError:
            relative_parent = Path()
        population = [
            part
            for part in relative_parent.parts
            if part and part != "."
        ]

    return sample_id, population


def is_dna_motif(motif):
    """Return True when all motif characters are DNA bases."""
    return bool(motif) and all(base in {"A", "C", "G", "T", "N"} for base in motif)


def is_minimal_motif(motif):
    """Return True when a motif is not made from a smaller repeated motif."""
    for size in range(1, len(motif)):
        if len(motif) % size:
            continue
        smaller = motif[:size]
        if smaller * (len(motif) // size) == motif:
            return False
    return True


def longest_repeat(sequence, motif_length):
    """Return the longest repeated motif run for a fixed motif length."""
    sequence = sequence.upper()
    best_motif = ""
    best_count = 0
    limit = len(sequence) - motif_length + 1

    for start in range(max(limit, 0)):
        motif = sequence[start:start + motif_length]
        if not is_dna_motif(motif) or not is_minimal_motif(motif):
            continue

        count = 1
        next_start = start + motif_length
        while sequence[next_start:next_start + motif_length] == motif:
            count += 1
            next_start += motif_length

        if count > best_count and count > 1:
            best_motif = motif
            best_count = count

    return best_motif, best_count


def repeat_profile(sequence):
    """Return 1-, 2-, and 3-base repeat motifs and counts for a sequence."""
    return {
        motif_length: longest_repeat(sequence, motif_length)
        for motif_length in (1, 2, 3)
    }


def strongest_alt_repeat_profile(alts):
    """Return the strongest 1-, 2-, and 3-base repeats across ALT alleles."""
    strongest = {
        1: ("", 0),
        2: ("", 0),
        3: ("", 0),
    }

    for alt in alts:
        for motif_length, (motif, count) in repeat_profile(alt).items():
            if count > strongest[motif_length][1]:
                strongest[motif_length] = (motif, count)

    return strongest


def reference_context(reference, position, flank_size=6):
    """Return reference bases before and after a 1-indexed circular position."""
    reference_length = len(reference)

    def base_at(wrapped_position):
        index = (wrapped_position - 1) % reference_length + 1
        return reference[index]

    before = "".join(
        base_at(position - offset)
        for offset in range(flank_size, 0, -1)
    )
    after = "".join(
        base_at(position + offset)
        for offset in range(1, flank_size + 1)
    )

    return before, after


def iter_variant_files(input_dir):
    """Yield one VCF per sample from an input directory."""
    base_dir = Path(input_dir)
    if not base_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {base_dir}")

    selected_paths = {}
    for path in sorted(base_dir.rglob("*")):
        if not path.is_file() or not is_vcf_file(path):
            continue

        sample_id, population = parse_sample_filename(path, base_dir)
        sample_key = (sample_id, tuple(population))
        existing_path = selected_paths.get(sample_key)
        if existing_path is None or existing_path.suffix == ".gz":
            selected_paths[sample_key] = path

    for path in selected_paths.values():
        yield path


def load_sample(path, input_dir=None):
    """Build a Sample from a filtered VCF path and its mutations."""
    sample_id, population = parse_sample_filename(path, input_dir)
    sample = Sample(
        sample_id=sample_id,
        population=population,
        source_path=Path(path),
    )

    with vcf_parser.VCFIterator(str(path)) as rows:
        for mutation in rows:
            sample.mutations.append(mutation)

    return sample


def load_samples(input_dir):
    """Load all logical samples from a variant call set."""
    return [
        load_sample(path, input_dir)
        for path in iter_variant_files(input_dir)
    ]


def load_mito_reference(path=DEFAULT_MITO_REFERENCE_PATH):
    """Load a mitochondrial FASTA sequence as a 1-indexed position map."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"Missing mitochondrial reference FASTA: {path}. "
            "Run scripts/download_mito_ref.py first."
        )

    sequence_parts = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            sequence_parts.append(line.upper())

    sequence = "".join(sequence_parts)
    if len(sequence) != MITOCHONDRIAL_LENGTH:
        raise ValueError(
            f"Expected mitochondrial reference length {MITOCHONDRIAL_LENGTH}, "
            f"got {len(sequence)} from {path}"
        )

    return {
        position: base
        for position, base in enumerate(sequence, start=1)
    }


def summarize_sample(sample):
    """Return basic mutation counts for a loaded sample."""
    total_variants = len(sample.mutations)
    pass_variants = sum(
        1 for mutation in sample.mutations
        if mutation.filter == "PASS"
    )
    multi_alt_variants = sum(
        1 for mutation in sample.mutations
        if mutation.has_multiple_alt_alleles
    )
    example_metadata = sample.mutations[-1].metadata if sample.mutations else {}

    return total_variants, pass_variants, multi_alt_variants, example_metadata


def format_mutation(mutation):
    """Return one readable line for a parsed VCF mutation."""
    return (
        f"  POS {mutation.position}: "
        f"{mutation.ref}>{mutation.alt} "
        f"FILTER={mutation.filter} "
        f"ALTS={mutation.alts} "
        f"AFS={mutation.afs} "
        f"ALT_AFS={mutation.alt_afs} "
        f"METADATA={mutation.metadata}"
    )


def sample_label(sample):
    """Return a compact label for sample-specific output columns."""
    return sample.label


def append_unique(values, value):
    """Append value if present and not already in values."""
    if value and value not in values:
        values.append(value)


def mutation_alt_text(mutations):
    """Return comma-packed ALT values for one sample at one position."""
    return ",".join(
        mutation.alt
        for mutation in mutations
        if mutation.alt
    )


def mutation_af_text(mutations):
    """Return comma-packed AF values for one sample at one position."""
    return ",".join(
        mutation.allele_fraction_text
        for mutation in mutations
        if mutation.allele_fraction_text
    )


def mutation_metadata_text(mutation):
    """Return mutation metadata as semicolon-delimited key=value pairs."""
    return ";".join(
        f"{key}={value}"
        for key, value in mutation.metadata.items()
    )


def mutation_alt_af_text(mutation):
    """Return ALT/AF pairs as semicolon-delimited key=value pairs."""
    return ";".join(
        f"{alt}={af}"
        for alt, af in mutation.alt_afs.items()
    )


def safe_sheet_title(title, existing_titles):
    """Return an Excel-safe sheet title, unique within a workbook."""
    invalid_chars = "[]:*?/\\"
    cleaned = "".join("_" if char in invalid_chars else char for char in title)
    cleaned = cleaned[:31] or "Sheet"

    if cleaned not in existing_titles:
        return cleaned

    counter = 2
    while True:
        suffix = f"_{counter}"
        candidate = f"{cleaned[:31 - len(suffix)]}{suffix}"
        if candidate not in existing_titles:
            return candidate
        counter += 1


def append_sample_sheet(workbook, sample, reference):
    """Append one worksheet containing all mutations from one VCF sample."""
    title = safe_sheet_title(sample_label(sample), set(workbook.sheetnames))
    worksheet = workbook.create_sheet(title)
    worksheet.append([
        "pos",
        "ref",
        "vcf ref",
        "alt",
        "af",
        "filter",
        "alts",
        "afs",
        "alt_afs",
        "metadata",
        "source file",
    ])

    for mutation in sample.mutations:
        worksheet.append([
            mutation.position,
            reference.get(mutation.position, ""),
            mutation.ref,
            mutation.alt,
            mutation.allele_fraction_text,
            mutation.filter,
            ",".join(mutation.alts),
            ",".join(str(af) for af in mutation.afs),
            mutation_alt_af_text(mutation),
            mutation_metadata_text(mutation),
            sample.source_path.name,
        ])


def make_master_xlsx(
    output_path,
    input_dir,
    reference_path=DEFAULT_MITO_REFERENCE_PATH,
):
    """Write a master mutation table as an XLSX workbook."""
    try:
        from openpyxl import Workbook
    except ImportError as exc:
        raise ImportError(
            "openpyxl is required to write master.xlsx. "
            "Install it with: python3 -m pip install openpyxl"
        ) from exc

    samples = load_samples(input_dir)
    reference = load_mito_reference(reference_path)
    mutation_lookup = {}

    for sample in samples:
        label = sample_label(sample)
        for mutation in sample.mutations:
            mutation_lookup.setdefault((label, mutation.position), []).append(mutation)

    headers = ["pos", "ref"]
    for sample in samples:
        label = sample_label(sample)
        headers.extend([f"alt of {label}", f"af of {label}"])

    workbook = Workbook(write_only=True)
    worksheet = workbook.create_sheet("Master")
    worksheet.append(headers)

    for position in mitochondrial_base_positions():
        row = [position, reference[position]]
        for sample in samples:
            label = sample_label(sample)
            mutations = mutation_lookup.get((label, position), [])
            row.extend([
                mutation_alt_text(mutations),
                mutation_af_text(mutations),
            ])
        worksheet.append(row)

    for sample in samples:
        append_sample_sheet(workbook, sample, reference)

    output_path = Path(output_path)
    workbook.save(output_path)
    return output_path


def create_database_schema(connection):
    """Create the subject/sample/mutation database schema."""
    connection.executescript("""
        PRAGMA foreign_keys = ON;

        DROP TABLE IF EXISTS mutation_alts;
        DROP TABLE IF EXISTS mutations;
        DROP TABLE IF EXISTS sample_population_tags;
        DROP TABLE IF EXISTS samples;
        DROP TABLE IF EXISTS subjects;

        CREATE TABLE subjects (
            id INTEGER PRIMARY KEY,
            subject_id TEXT NOT NULL UNIQUE
        );

        CREATE TABLE samples (
            id INTEGER PRIMARY KEY,
            subject_id INTEGER NOT NULL,
            population_key TEXT NOT NULL,
            source_file TEXT NOT NULL,
            UNIQUE(subject_id, population_key),
            FOREIGN KEY(subject_id) REFERENCES subjects(id)
        );

        CREATE TABLE sample_population_tags (
            sample_id INTEGER NOT NULL,
            tag TEXT NOT NULL,
            tag_order INTEGER NOT NULL,
            PRIMARY KEY(sample_id, tag),
            FOREIGN KEY(sample_id) REFERENCES samples(id)
        );

        CREATE TABLE mutations (
            id INTEGER PRIMARY KEY,
            sample_id INTEGER NOT NULL,
            pos INTEGER NOT NULL,
            ref TEXT NOT NULL,
            vcf_ref TEXT NOT NULL,
            alt TEXT NOT NULL,
            af TEXT NOT NULL,
            polymorphism INTEGER NOT NULL,
            repeat_base TEXT NOT NULL,
            repeat_count INTEGER NOT NULL,
            repeat_2_bases TEXT NOT NULL,
            repeat_2_count INTEGER NOT NULL,
            repeat_3_bases TEXT NOT NULL,
            repeat_3_count INTEGER NOT NULL,
            filter TEXT NOT NULL,
            metadata_json TEXT NOT NULL,
            FOREIGN KEY(sample_id) REFERENCES samples(id)
        );

        CREATE TABLE mutation_alts (
            mutation_id INTEGER NOT NULL,
            alt_index INTEGER NOT NULL,
            alt TEXT NOT NULL,
            af REAL,
            af_text TEXT,
            polymorphism INTEGER NOT NULL,
            repeat_base TEXT NOT NULL,
            repeat_count INTEGER NOT NULL,
            repeat_2_bases TEXT NOT NULL,
            repeat_2_count INTEGER NOT NULL,
            repeat_3_bases TEXT NOT NULL,
            repeat_3_count INTEGER NOT NULL,
            PRIMARY KEY(mutation_id, alt_index),
            FOREIGN KEY(mutation_id) REFERENCES mutations(id)
        );

        CREATE INDEX idx_samples_subject_id ON samples(subject_id);
        CREATE INDEX idx_sample_population_tags_tag ON sample_population_tags(tag);
        CREATE INDEX idx_mutations_sample_pos ON mutations(sample_id, pos);
        CREATE INDEX idx_mutations_pos ON mutations(pos);
        CREATE INDEX idx_mutation_alts_alt ON mutation_alts(alt);
    """)


def insert_subject(connection, subject_id):
    """Insert or fetch a subject row."""
    connection.execute(
        "INSERT OR IGNORE INTO subjects(subject_id) VALUES (?)",
        (subject_id,),
    )
    cursor = connection.execute(
        "SELECT id FROM subjects WHERE subject_id = ?",
        (subject_id,),
    )
    return cursor.fetchone()[0]


def insert_sample(connection, sample, subject_db_id):
    """Insert or fetch a sample row for a subject/population combination."""
    connection.execute(
        """
        INSERT OR IGNORE INTO samples(subject_id, population_key, source_file)
        VALUES (?, ?, ?)
        """,
        (subject_db_id, sample.population_key, sample.source_path.name),
    )
    cursor = connection.execute(
        """
        SELECT id FROM samples
        WHERE subject_id = ? AND population_key = ?
        """,
        (subject_db_id, sample.population_key),
    )
    sample_db_id = cursor.fetchone()[0]

    for tag_order, tag in enumerate(sample.population):
        connection.execute(
            """
            INSERT OR IGNORE INTO sample_population_tags(sample_id, tag, tag_order)
            VALUES (?, ?, ?)
            """,
            (sample_db_id, tag, tag_order),
        )

    return sample_db_id


def insert_mutation(connection, sample_db_id, mutation, reference):
    """Insert one mutation and its individual ALT/AF rows."""
    polymorphism = int(mutation.has_multiple_alt_alleles)
    repeats = strongest_alt_repeat_profile(mutation.alts)
    repeat_base, repeat_count = repeats[1]
    repeat_2_bases, repeat_2_count = repeats[2]
    repeat_3_bases, repeat_3_count = repeats[3]
    reference_6_before, reference_6_after = reference_context(
        reference,
        mutation.position,
    )
    metadata = {
        **mutation.metadata,
        "POLYMORPHISM": str(polymorphism),
        "REPEAT_1_BASE": repeat_base,
        "REPEAT_1_BASE_COUNT": str(repeat_count),
        "REPEAT_2_BASES": repeat_2_bases,
        "REPEAT_2_BASES_COUNT": str(repeat_2_count),
        "REPEAT_3_BASES": repeat_3_bases,
        "REPEAT_3_BASES_COUNT": str(repeat_3_count),
        "REFERENCE_6_BEFORE": reference_6_before,
        "REFERENCE_6_AFTER": reference_6_after,
    }
    cursor = connection.execute(
        """
        INSERT INTO mutations(
            sample_id, pos, ref, vcf_ref, alt, af, polymorphism, repeat_base,
            repeat_count, repeat_2_bases, repeat_2_count, repeat_3_bases,
            repeat_3_count, filter, metadata_json
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            sample_db_id,
            mutation.position,
            reference[mutation.position],
            mutation.ref,
            mutation.alt,
            mutation.allele_fraction_text,
            polymorphism,
            repeat_base,
            repeat_count,
            repeat_2_bases,
            repeat_2_count,
            repeat_3_bases,
            repeat_3_count,
            mutation.filter,
            json.dumps(metadata, sort_keys=True),
        ),
    )
    mutation_db_id = cursor.lastrowid

    for alt_index, alt in enumerate(mutation.alts):
        af_value = mutation.afs[alt_index] if alt_index < len(mutation.afs) else None
        numeric_af = af_value if isinstance(af_value, float) else None
        af_text = "" if af_value is None else str(af_value)
        alt_repeats = repeat_profile(alt)
        alt_repeat_base, alt_repeat_count = alt_repeats[1]
        alt_repeat_2_bases, alt_repeat_2_count = alt_repeats[2]
        alt_repeat_3_bases, alt_repeat_3_count = alt_repeats[3]
        connection.execute(
            """
            INSERT INTO mutation_alts(
                mutation_id, alt_index, alt, af, af_text, polymorphism,
                repeat_base, repeat_count, repeat_2_bases, repeat_2_count,
                repeat_3_bases, repeat_3_count
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                mutation_db_id,
                alt_index,
                alt,
                numeric_af,
                af_text,
                polymorphism,
                alt_repeat_base,
                alt_repeat_count,
                alt_repeat_2_bases,
                alt_repeat_2_count,
                alt_repeat_3_bases,
                alt_repeat_3_count,
            ),
        )


def make_sql_database(
    output_path,
    input_dir,
    reference_path=DEFAULT_MITO_REFERENCE_PATH,
):
    """Create a SQLite database for subjects, samples, populations, and mutations."""
    samples = load_samples(input_dir)
    reference = load_mito_reference(reference_path)
    output_path = Path(output_path)

    with sqlite3.connect(output_path) as connection:
        create_database_schema(connection)
        skipped_out_of_reference = 0

        for sample in samples:
            subject_db_id = insert_subject(connection, sample.sample_id)
            sample_db_id = insert_sample(connection, sample, subject_db_id)
            for mutation in sample.mutations:
                if mutation.position not in reference:
                    skipped_out_of_reference += 1
                    continue
                insert_mutation(connection, sample_db_id, mutation, reference)

    if skipped_out_of_reference:
        print(
            f"Skipped {skipped_out_of_reference} mutation(s) outside the "
            f"{len(reference)} bp reference."
        )

    return output_path


def print_sample(sample, path):
    """Print a loaded sample and all of its parsed mutations."""
    total_variants, pass_variants, multi_alt_variants, _ = summarize_sample(sample)
    print(
        f"{path.name}: "
        f"{sample.label}: "
        f"{total_variants} variants "
        f"({pass_variants} PASS, {multi_alt_variants} multi-ALT)"
    )

    for mutation in sample.mutations:
        print(format_mutation(mutation))


def list_variant_files(input_dir):
    """Parse filtered VCFs and print each file's parsed mutations."""
    for path in iter_variant_files(input_dir):
        sample = load_sample(path, input_dir)
        print_sample(sample, path)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Build mutation SQLite/XLSX outputs from VCF files.",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing VCF files for the study.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="SQLite output path.",
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=DEFAULT_MITO_REFERENCE_PATH,
        help=f"Mitochondrial reference FASTA. Default: {DEFAULT_MITO_REFERENCE_PATH}",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="Print parsed VCF summaries instead of creating a SQLite database.",
    )
    return parser.parse_args(argv)

def file_in_folder_matching(path, match):
    """Return True if any parent folder name contains the given string."""
    dir_path = Path(path).resolve().parent
    while dir_path != dir_path.parent:
        if match in dir_path.name:
            return True
        dir_path = dir_path.parent
    return match in dir_path.name


if __name__ == "__main__":
    args = parse_args()
    if args.list:
        list_variant_files(input_dir=args.input_dir)
    else:
        output_path = make_sql_database(
            output_path=args.output,
            input_dir=args.input_dir,
            reference_path=args.reference,
        )
        print(f"Wrote {output_path}")
