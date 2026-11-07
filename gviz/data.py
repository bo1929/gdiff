"""Read and validate named gdiff TSV columns without loading a dataframe."""

import csv
import gzip
import math
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

ALIASES = {
    "query_id": "seq", "qid": "seq", "ref_id": "reference",
    "interval_start": "start", "interval_end": "end", "dist": "d",
}
MISSING = {"", ".", "na", "n/a", "nan", "-nan", "inf", "-inf"}
MAX_SERIES = 12
MATRIX_LIMIT = 60
RAW_INTERVAL_LIMIT = 256


def detect_kind(columns):
    cols = set(columns)
    if {"genome_a", "genome_b", "d"} <= cols:
        if "direction" in cols:
            if not {"seq", "start", "end"} <= cols:
                raise ValueError("Sampled-window output needs seq, start, end, and direction columns")
            return "samples"
        return "dist"
    if {"seq", "start", "end", "reference"} <= cols:
        if "d" in cols and cols & {"side", "alpha", "threshold", "mask", "d_interval"}:
            return "map"
        if "d" in cols or {"d_fw", "d_rc"} <= cols:
            return "roll"
    raise ValueError("Unrecognized header; expected map, roll, dist, or dist --output-samples columns")


def number(value):
    if value.strip().lower() in MISSING:
        return math.nan
    result = float(value)
    if not math.isfinite(result):
        return math.nan
    if result < 0:
        raise ValueError(f"negative distance: {value}")
    return result


class Table:
    def __init__(self, path):
        self.path = Path(path)
        opener = gzip.open if self.path.suffix == ".gz" else open
        self.opener = opener
        with opener(self.path, "rt", encoding="utf-8-sig") as handle:
            self.provenance = []
            for line in handle:
                if line.startswith("#"):
                    self.provenance.append(line[1:].strip())
                elif line.strip():
                    self.columns = [ALIASES.get(c.strip().lower(), c.strip().lower())
                                    for c in next(csv.reader([line], delimiter="\t"))]
                    break
            else:
                raise ValueError("Input has no TSV header")
        if len(set(self.columns)) != len(self.columns):
            raise ValueError("Duplicate column names after normalization")
        self.kind = detect_kind(self.columns)
        if self.kind == "roll" and (("d_fw" in self.columns) != ("d_rc" in self.columns)):
            raise ValueError("Strand-aware roll output needs both d_fw and d_rc")

    def rows(self):
        """Reopen for each pass. Memory is bounded by a single input row."""
        with self.opener(self.path, "rt", encoding="utf-8-sig") as handle:
            header = False
            for line_no, line in enumerate(handle, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                if not header:
                    header = True
                    continue
                try:
                    fields = next(csv.reader([line], delimiter="\t"))
                    if len(fields) != len(self.columns):
                        raise ValueError(f"expected {len(self.columns)} fields, found {len(fields)}")
                    row = dict(zip(self.columns, fields))
                    for name in ("d", "d_fw", "d_rc"):
                        if name in row:
                            row[name] = number(row[name])
                    if self.kind != "dist":
                        for name in ("start", "end"):
                            row[name] = int(row[name])
                        if row["start"] < 1 or row["end"] < row["start"]:
                            raise ValueError("coordinates must be 1-based with end >= start")
                        if "seq_len" in row:
                            row["seq_len"] = int(row["seq_len"])
                            if row["seq_len"] < row["end"]:
                                raise ValueError("end exceeds seq_len")
                        if not row["seq"]:
                            raise ValueError("empty sequence name")
                    names = ("genome_a", "genome_b") if self.kind in {"dist", "samples"} else ("reference",)
                    if any(not row[name] for name in names):
                        raise ValueError("empty genome/reference name")
                    if self.kind == "samples" and row["direction"] not in {"ab", "ba"}:
                        raise ValueError("sample direction must be ab or ba")
                    yield row
                except (ValueError, csv.Error) as error:
                    raise ValueError(f"{self.path}:{line_no}: {error}") from error


@dataclass
class Selection:
    sequence: str | None = None
    region: tuple | None = None
    references: tuple = ()

    def matches(self, row, kind):
        if self.sequence is not None and row.get("seq") != self.sequence:
            return False
        if self.references:
            names = (row["genome_a"], row["genome_b"]) if kind in {"dist", "samples"} else (row["reference"],)
            if not any(name in self.references for name in names):
                return False
        if self.region:
            start, end = self.region
            if kind == "map":
                return row["end"] >= start and row["start"] <= end
            center = (row["start"] - 1 + row["end"]) / 2
            return start - 1 <= center < end
        return True


def series(row, kind):
    if kind == "roll" and "d_fw" in row:
        return [((row["reference"], "forward"), row["d_fw"]),
                ((row["reference"], "reverse complement"), row["d_rc"])]
    if kind == "samples":
        a, b = row["genome_a"], row["genome_b"]
        if row["direction"] == "ba":
            a, b = b, a
        return [((a, b, row.get("config", "gdiff")), row["d"])]
    key = (row["reference"], row.get("side", "intervals") if kind == "map" else "")
    if kind == "map" and "is_rc" in row:
        key += (f"is_rc={row['is_rc']}",)
    return [(key, row["d"])]


@dataclass
class Overview:
    total: int = 0
    selected: int = 0
    missing: int = 0
    low: float = math.inf
    high: float = -math.inf
    lengths: dict = field(default_factory=dict)
    groups: dict = field(default_factory=lambda: defaultdict(Counter))
    genomes: set = field(default_factory=set)
    genomes_a: set = field(default_factory=set)
    genomes_b: set = field(default_factory=set)
    pairs: set = field(default_factory=set)


def inspect(table, selection):
    overview = Overview()
    for row in table.rows():
        overview.total += 1
        if not selection.matches(row, table.kind):
            continue
        overview.selected += 1
        if table.kind == "dist":
            overview.genomes.update((row["genome_a"], row["genome_b"]))
            overview.genomes_a.add(row["genome_a"])
            overview.genomes_b.add(row["genome_b"])
            if len(overview.genomes) <= MATRIX_LIMIT:
                if row["genome_a"] != row["genome_b"]:
                    overview.pairs.add(tuple(sorted((row["genome_a"], row["genome_b"]))))
            else:
                overview.pairs.clear()
            values = [(None, row["d"])]
        else:
            seq = row["seq"]
            overview.lengths[seq] = max(overview.lengths.get(seq, 0), row.get("seq_len", row["end"]))
            values = series(row, table.kind)
            scope = "all" if table.kind == "samples" else seq
            for key, _ in values:
                overview.groups[scope][key] += 1
        for _, value in values:
            if math.isfinite(value):
                overview.low = min(overview.low, value)
                overview.high = max(overview.high, value)
            else:
                overview.missing += 1
    if not overview.selected:
        raise ValueError("No rows match the selection (use exact sequence/reference names)")
    return overview
