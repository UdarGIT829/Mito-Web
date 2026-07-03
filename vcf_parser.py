import gzip
from dataclasses import dataclass
from itertools import zip_longest
from typing import Dict, Iterator, List


@dataclass
class VCFMutation:
    """A richer representation of one VCF mutation row."""

    position: int
    alt: str
    metadata: Dict[str, str]
    ref: str = ""
    filter: str = ""
    raw_row: Dict[str, str] | None = None

    @property
    def alts(self):
        """Return the ALT field split into individual alternate alleles."""
        if not self.alt:
            return []
        return self.alt.split(",")

    @property
    def alt_alleles(self):
        """Return the ALT field split into individual alternate alleles."""
        return self.alts

    @property
    def afs(self):
        """Return allele fractions from metadata AF as floats when possible."""
        af = self.metadata.get("AF", "")
        if not af:
            return []

        values = []
        for value in af.split(","):
            try:
                values.append(float(value))
            except ValueError:
                values.append(value)
        return values

    @property
    def alt_afs(self):
        """Return alternate alleles paired with their AF values."""
        return dict(zip_longest(self.alts, self.afs, fillvalue=""))

    @property
    def has_multiple_alt_alleles(self):
        """Return True when this row contains more than one ALT allele."""
        return len(self.alts) > 1

    @classmethod
    def from_row(cls, row: Dict[str, str]):
        format_fields = row.get("FORMAT", "").split(":")
        metadata_values = row.get("metadata", "").split(":")
        metadata = {
            key: value
            for key, value in zip_longest(format_fields, metadata_values, fillvalue="")
            if key
        }

        return cls(
            position=int(row["POS"]),
            alt=row.get("ALT", ""),
            metadata=metadata,
            ref=row.get("REF", ""),
            filter=row.get("FILTER", ""),
            raw_row=row,
        )


class VCFIterator:
    """Simple VCF iterator that yields VCFMutation rows.

    Usage:
        it = VCFIterator(path)
        for mutation in it:
            # mutation.position, mutation.alt, mutation.metadata
    """

    def __init__(self, path: str):
        self.path = path
        self._fh = None
        self._mutations = None
        self.columns: List[str] = []
        self.sample_columns: List[str] = []
        self._open()

    @classmethod
    def from_mutations(cls, mutations, path: str = "<memory>"):
        """Build a VCFIterator-compatible object from in-memory mutations."""
        iterator = cls.__new__(cls)
        iterator.path = path
        iterator._fh = None
        iterator._mutations = iter(mutations)
        iterator.columns = []
        iterator.sample_columns = []
        return iterator

    def _open(self):
        opener = gzip.open if self.path.endswith('.gz') else open
        self._fh = opener(self.path, 'rt', encoding='utf-8', errors='replace')

        # advance until we find the column header (starts with #CHROM)
        for raw in self._fh:
            if raw.startswith('#CHROM'):
                self.columns = raw.lstrip('#').rstrip('\n').split('\t')
                self._normalize_sample_columns()
                break

    def _normalize_sample_columns(self):
        """Rename sample data columns after FORMAT to stable metadata keys."""
        if 'FORMAT' not in self.columns:
            return

        format_index = self.columns.index('FORMAT')
        self.sample_columns = self.columns[format_index + 1:]
        if not self.sample_columns:
            return

        self.columns[format_index + 1] = 'metadata'
        for index in range(format_index + 2, len(self.columns)):
            metadata_index = index - format_index
            self.columns[index] = f'metadata_{metadata_index}'

    def __iter__(self) -> Iterator[VCFMutation]:
        return self

    def __next__(self) -> VCFMutation:
        if self._mutations is not None:
            return next(self._mutations)

        if self._fh is None:
            raise StopIteration

        for raw in self._fh:
            if not raw:
                continue
            if raw.startswith('#'):
                continue
            line = raw.rstrip('\n')
            if not line:
                continue
            fields = line.split('\t')
            # Map columns to values; unmatched columns get empty string
            row = {col: val for col, val in zip_longest(self.columns, fields, fillvalue='')}
            return VCFMutation.from_row(row)

        # EOF
        self._fh.close()
        self._fh = None
        raise StopIteration

    def close(self):
        if self._fh is not None:
            try:
                self._fh.close()
            finally:
                self._fh = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()


if __name__ == '__main__':
    import sys

    if len(sys.argv) < 2:
        print('Usage: python vcf_parser.py path/to/file.vcf[.gz]')
        sys.exit(1)

    path = sys.argv[1]
    with VCFIterator(path) as it:
        # print first 5 rows as a quick smoke check
        for i, row in enumerate(it):
            print(row)
            if i >= 4:
                break
