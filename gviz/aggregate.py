"""Bounded arrays for positional summaries, distance matrices, and distributions."""

import math
import heapq
from dataclasses import dataclass

import numpy as np

from .data import MATRIX_LIMIT, MAX_SERIES, RAW_INTERVAL_LIMIT, inspect, series


class WindowBins:
    def __init__(self, edges):
        self.edges = edges
        n = len(edges) - 1
        self.count = np.zeros(n, dtype=np.int64)
        self.positions = np.zeros(n)
        self.valid = np.zeros(n, dtype=np.int64)
        self.total = np.zeros(n)
        self.minimum = np.full(n, np.inf)
        self.maximum = np.full(n, -np.inf)

    def add(self, start, end, value):
        # Bin centers in 0-based base coordinates; input is 1-based inclusive.
        position = (start - 1 + end) / 2
        i = min(len(self.count) - 1, int((position - self.edges[0]) /
                                        (self.edges[-1] - self.edges[0]) * len(self.count)))
        if position < self.edges[0] or position >= self.edges[-1]:
            return
        self.count[i] += 1
        self.positions[i] += position
        if math.isfinite(value):
            self.valid[i] += 1
            self.total[i] += value
            self.minimum[i] = min(self.minimum[i], value)
            self.maximum[i] = max(self.maximum[i], value)

    def finish(self):
        mean = np.divide(self.total, self.valid, out=np.full_like(self.total, np.nan), where=self.valid > 0)
        missing = np.divide(self.count - self.valid, self.count, out=np.full_like(self.total, np.nan), where=self.count > 0)
        return {"mean": mean, "minimum": np.where(self.valid, self.minimum, np.nan),
                "maximum": np.where(self.valid, self.maximum, np.nan), "missing": missing,
                "positions": np.divide(self.positions, self.count,
                                       out=np.full_like(self.positions, np.nan), where=self.count > 0),
                "observed": self.count > 0, "rows": int(self.count.sum()),
                "missing_rows": int((self.count - self.valid).sum())}


class IntervalBins:
    def __init__(self, edges):
        self.edges = edges
        n = len(edges)
        self.touches = np.zeros(n, dtype=np.int64)
        self.valid_touches = np.zeros(n, dtype=np.int64)
        self.weight = np.zeros(n)
        self.total = np.zeros(n)
        self.partial_weight = np.zeros(n - 1)
        self.partial_total = np.zeros(n - 1)
        self.coverage = np.zeros(n)
        self.partial_coverage = np.zeros(n - 1)
        self.intervals = []
        self.rows = 0
        self.missing_rows = 0

    def add(self, start, end, value):
        """Range additions are O(1), even for intervals spanning whole chromosomes."""
        left, right = max(start - 1, self.edges[0]), min(end, self.edges[-1])
        if right <= left:
            return
        self.rows += 1
        self.missing_rows += not math.isfinite(value)
        if self.intervals is not None:
            if len(self.intervals) < RAW_INTERVAL_LIMIT:
                self.intervals.append((start, end, value))
            else:
                self.intervals = None
        n = len(self.edges) - 1
        width = (self.edges[-1] - self.edges[0]) / n
        first = min(n - 1, max(0, int((left - self.edges[0]) / width)))
        last = min(n - 1, max(0, math.ceil((right - self.edges[0]) / width) - 1))
        self.touches[first] += 1
        self.touches[last + 1] -= 1
        if first == last:
            self.partial_coverage[first] += right - left
        else:
            self.partial_coverage[first] += self.edges[first + 1] - left
            self.partial_coverage[last] += right - self.edges[last]
            self.coverage[first + 1] += width
            self.coverage[last] -= width
        if not math.isfinite(value):
            return
        self.valid_touches[first] += 1
        self.valid_touches[last + 1] -= 1
        if first == last:
            self.partial_weight[first] += right - left
            self.partial_total[first] += value * (right - left)
            return
        for index, weight in ((first, self.edges[first + 1] - left), (last, right - self.edges[last])):
            self.partial_weight[index] += weight
            self.partial_total[index] += value * weight
        self.weight[first + 1] += width
        self.weight[last] -= width
        self.total[first + 1] += width * value
        self.total[last] -= width * value

    def finish(self):
        weight = np.cumsum(self.weight[:-1]) + self.partial_weight
        total = np.cumsum(self.total[:-1]) + self.partial_total
        # Roundoff in the range sums must not turn uncovered bins into distances.
        count = np.cumsum(self.touches[:-1])
        valid = np.cumsum(self.valid_touches[:-1])
        mean = np.divide(total, weight, out=np.full_like(total, np.nan), where=(weight > 0) & (valid > 0))
        coverage = np.cumsum(self.coverage[:-1]) + self.partial_coverage
        overlap = np.where(count > 0, np.maximum(coverage / np.diff(self.edges), 0), 0)
        return {"mean": mean, "count": count, "overlap": overlap,
                "intervals": self.intervals,
                "rows": self.rows, "missing_rows": self.missing_rows}


@dataclass
class Report:
    kind: str
    title: str
    notes: list
    data: dict
    rows: int
    total: int


def histogram_edges(overview, bins):
    low = overview.low if math.isfinite(overview.low) else 0.0
    high = overview.high if math.isfinite(overview.high) else 1.0
    if high == low:
        low, high = max(0, low - 0.005), high + 0.005
    return np.linspace(low, high, min(bins, 200) + 1)


def summarize(table, selection, bins=1000):
    overview = inspect(table, selection)
    notes = []
    if table.kind in {"roll", "map"}:
        seq = selection.sequence or max(overview.lengths, key=lambda s: (overview.lengths[s], s))
        start, end = selection.region or (1, overview.lengths[seq])
        groups = overview.groups[seq]
        resolution = min(bins, max(groups.values())) if table.kind == "roll" else bins
        edges = np.linspace(start - 1, end, min(resolution, end - start + 1) + 1)
        keys = sorted(groups, key=lambda key: (-groups[key], key))[:MAX_SERIES]
        keys.sort(key=lambda key: (key[0], {"low": 0, "high": 1}.get(key[1], 0), key[1:]))
        accumulator = IntervalBins if table.kind == "map" else WindowBins
        accumulators = {key: accumulator(edges) for key in keys}
        rows = 0
        for row in table.rows():
            if row["seq"] != seq or not selection.matches(row, table.kind):
                continue
            rows += 1
            for key, value in series(row, table.kind):
                if key in accumulators:
                    accumulators[key].add(row["start"], row["end"], value)
        notes.append(f"Sequence: {seq}; region: {start:,}–{end:,} (1-based, inclusive).")
        notes.append(f"{len(edges) - 1:,} equal-width bins; use --seq and --region START:END for finer detail.")
        notes.append("Sequence extent comes from seq_len when available, otherwise the largest reported end coordinate.")
        if len(groups) > MAX_SERIES:
            notes.append(f"Showing the {MAX_SERIES} series with the most rows out of {len(groups):,}; select others with --reference.")
        if table.kind == "map":
            notes.append("Each interval track shows reported start and end coordinates; gray segments have no distance estimate. Crowded tracks show mean overlap per base in each bin instead. No interval boundaries are inferred by the plot.")
        else:
            notes.append("Each reference has its own panel on a shared distance scale. Lines join observed window centers; the band shows the within-bin range. Gray ticks mark windows without a distance estimate.")
        summaries = {k: v.finish() for k, v in accumulators.items()}
        return Report(table.kind, seq, notes,
                      {"edges": edges, "series": summaries}, rows, overview.total)

    histogram_bins = min(bins, max(10, math.ceil(math.sqrt(overview.selected)))) if table.kind == "dist" else bins
    edges = histogram_edges(overview, histogram_bins)
    if table.kind == "samples":
        groups = overview.groups["all"]
        keys = sorted(groups, key=lambda key: (-groups[key], key))[:MAX_SERIES]
        histograms = {key: np.zeros(len(edges) - 1, dtype=np.int64) for key in keys}
        missing = dict.fromkeys(keys, 0)
        for row in table.rows():
            if not selection.matches(row, table.kind):
                continue
            for key, value in series(row, table.kind):
                if key in histograms:
                    if math.isfinite(value):
                        i = min(len(edges) - 2, int(np.searchsorted(edges, value, side="right") - 1))
                        histograms[key][i] += 1
                    else:
                        missing[key] += 1
        notes.append("Directional window distributions over all sequences unless --seq is supplied. ECDFs are approximated by histogram bins; missing distances are excluded.")
        if len(groups) > MAX_SERIES:
            notes.append(f"Showing {MAX_SERIES} of {len(groups):,} directional series by row count; use --reference to filter.")
        return Report("samples", "Sampled window distances", notes,
                      {"edges": edges, "histograms": histograms, "missing": missing}, overview.selected, overview.total)

    names = sorted(overview.genomes)
    symmetric = len(names) > 1 and len(overview.pairs) == len(names) * (len(names) - 1) // 2
    y_names = names if symmetric else sorted(overview.genomes_a)
    x_names = names if symmetric else sorted(overview.genomes_b)
    matrix = np.zeros((len(y_names), len(x_names))) if len(names) <= MATRIX_LIMIT else None
    counts = np.zeros_like(matrix, dtype=np.int64) if matrix is not None else None
    y_indices = {name: i for i, name in enumerate(y_names)}
    x_indices = {name: i for i, name in enumerate(x_names)}
    histogram = np.zeros(len(edges) - 1, dtype=np.int64)
    closest = []
    for row in table.rows():
        if not selection.matches(row, table.kind) or not math.isfinite(row["d"]):
            continue
        value = row["d"]
        histogram[min(len(edges) - 2, int(np.searchsorted(edges, value, side="right") - 1))] += 1
        item = (-value, row["genome_a"], row["genome_b"])
        if len(closest) < 20:
            heapq.heappush(closest, item)
        else:
            heapq.heappushpop(closest, item)
        if matrix is not None:
            a, b = y_indices[row["genome_a"]], x_indices[row["genome_b"]]
            matrix[a, b] += value
            counts[a, b] += 1
            if symmetric and a != b:
                matrix[b, a] += value
                counts[b, a] += 1
    if matrix is not None:
        matrix = np.divide(matrix, counts, out=np.full_like(matrix, np.nan), where=counts > 0)
        notes.append("Complete all-pairs inputs use a symmetric matrix; cross-comparisons use the actual A × B layout. Gray cells have no finite estimate. Repeated pairs are averaged.")
    else:
        notes.append(f"{len(names):,} genomes: showing all pair distances as a distribution and the 20 closest reported comparisons instead of a dense matrix (limit {MATRIX_LIMIT}). Use --reference to filter.")
    notes.append(f"{overview.missing:,} rows have no finite distance estimate.")
    return Report("dist", "Genome-wide distances", notes,
                  {"names": names, "x_names": x_names, "y_names": y_names, "symmetric": symmetric,
                   "matrix": matrix, "edges": edges, "histogram": histogram,
                   "closest": sorted([(-d, a, b) for d, a, b in closest])}, overview.selected, overview.total)
