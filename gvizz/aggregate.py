"""Summaries for gvizz: exact-coordinate interval tracks, empirical CDFs, pair matrices.

Nothing is binned or counted per position: every reported interval or window is
kept at its exact coordinates, and distribution views use the raw values.
"""

import math
from dataclasses import dataclass, field

import numpy as np

from .data import MATRIX_LIMIT, MAX_SERIES, inspect, series, short_name

MAX_SEGMENTS = 20000   # per track; longer tracks keep their longest intervals
ECDF_CAP = 20000       # per direction; longer ECDFs are stride-subsampled
CLOSEST_LIMIT = 20


@dataclass
class Report:
    kind: str
    title: str
    notes: list
    data: dict
    rows: int
    total: int
    heading: str | None = None
    stats: list = field(default_factory=list)


def _trim(segments):
    """Keep the longest intervals when a track exceeds MAX_SEGMENTS."""
    if len(segments) <= MAX_SEGMENTS:
        return segments, 0
    segments.sort(key=lambda seg: seg[1] - seg[0], reverse=True)
    kept = segments[:MAX_SEGMENTS]
    kept.sort(key=lambda seg: (seg[0], seg[1]))
    return kept, len(segments) - MAX_SEGMENTS


def _rank(keys, spans):
    """Greatest covered span first: an interval-based ranking, not a row count."""
    ordered = sorted(keys, key=lambda key: (-spans[key], key))[:MAX_SERIES]
    ordered.sort(key=lambda key: (key[0], {"low": 0, "high": 1}.get(key[1], 0), key[1:]))
    return ordered


def summarize(table, selection, overlay=None, annotations=None, heading=None):
    overview = inspect(table, selection)
    notes = []
    if table.kind in {"map", "roll"}:
        # No implicit sequence choice: every sequence becomes a traversable
        # section of the report; --seq restricts the report to one.
        seqs = [selection.sequence] if selection.sequence else list(overview.lengths)
        tracks_by_seq = {s: {key: [] for key in overview.groups[s]} for s in seqs}
        rows = 0
        for row in table.rows():
            buckets = tracks_by_seq.get(row["seq"])
            if buckets is None or not selection.matches(row, table.kind):
                continue
            rows += 1
            for key, series_value in series(row, table.kind):
                if key in buckets:
                    buckets[key].append((row["start"], row["end"], series_value))
        overlay_by_seq, oview, orows = {}, None, 0
        if overlay is not None:
            if overlay.kind == table.kind:
                raise ValueError("overlay input must be a different kind (map over roll, or roll over map)")
            oview = inspect(overlay, selection)
            overlay_by_seq = {s: {key: [] for key in oview.groups[s]} for s in oview.lengths}
            orows = 0
            for row in overlay.rows():
                buckets = overlay_by_seq.get(row["seq"])
                if buckets is None or not selection.matches(row, overlay.kind):
                    continue
                orows += 1
                for key, series_value in series(row, overlay.kind):
                    if key in buckets:
                        buckets[key].append((row["start"], row["end"], series_value))
        sequences, trimmed = [], 0
        truncated = False
        if annotations:
            _, truncated = annotations
        ospans = oview.spans if oview is not None else {}
        for s in seqs:
            start, end = 1, overview.lengths[s]
            keys = _rank(tracks_by_seq[s], overview.spans[s])
            dropped_series = len(tracks_by_seq[s]) - len(keys)
            kept = {}
            for key in keys:
                kept[key], dropped = _trim(tracks_by_seq[s][key])
                trimmed += dropped
            overlays = {}
            if s in overlay_by_seq:
                primary_by_short = {}
                for ref in {key[0] for key in keys}:
                    primary_by_short.setdefault(short_name(ref), []).append(ref)
                for key in _rank(overlay_by_seq[s], ospans[s]):
                    candidates = primary_by_short.get(short_name(key[0]), [])
                    if len(candidates) != 1:
                        continue
                    segments, dropped = _trim(overlay_by_seq[s][key])
                    trimmed += dropped
                    overlays[(candidates[0], *key[1:])] = segments
            genes = []
            if annotations:
                features, _ = annotations
                genes = [f for f in features if f.seq == s and f.end >= start and f.start <= end]
            sequences.append({"seq": s, "extent": (start, end), "tracks": kept,
                              "overlay": overlays, "genes": genes,
                              "dropped_series": dropped_series})
        notes.append("Every reported interval and window is drawn at its exact coordinates on one "
                     "level per track; where intervals overlap, the lower distance is drawn on top. "
                     "Gray segments have no distance estimate. No boundaries are inferred or binned.")
        if overlay is not None:
            notes.append(f"Overlay: {orows:,} {overlay.kind} rows from {overlay.path.name}; tracks join the "
                         "matching reference by name (common file suffixes ignored when unambiguous).")
        if trimmed:
            notes.append(f"{trimmed:,} shortest intervals omitted to keep tracks under "
                         f"{MAX_SEGMENTS:,} segments.")
        capped = [(len(s["tracks"]), len(s["tracks"]) + s["dropped_series"])
                  for s in sequences if s["dropped_series"]]
        if capped:
            kept, total = max(capped, key=lambda item: item[1])
            notes.append(f"Showing the {kept} tracks with the greatest covered span out of "
                         f"{total}; the rest are omitted.")
        if annotations:
            notes.append("Annotation track: features on each plotted sequence"
                         + (" (truncated to the first 2,000 features in the file)" if truncated else "") + ".")
        if not selection.sequence and len(seqs) > 1:
            notes.append(f"{len(seqs)} sequences: the HTML report traverses them with the arrow links; "
                         "--seq NAME restricts the report to one sequence.")
        stats = [("Selected rows", f"{rows:,}"),
                 ("Sequences", f"{len(seqs):,}")]
        if math.isfinite(overview.low):
            stats.append(("d range", f"{overview.low:.3g} – {overview.high:.3g}"))
        title = selection.sequence or f"{len(seqs)} sequences"
        return Report(table.kind, title, notes,
                      {"sequences": sequences, "overlay_kind": overlay.kind if overlay else None,
                       "overlay_name": overlay.path.name if overlay else None},
                      rows, overview.total, heading, stats)

    if table.kind == "samples":
        values, missing = {}, {}
        for row in table.rows():
            if not selection.matches(row, table.kind):
                continue
            for key, series_value in series(row, table.kind):
                if math.isfinite(series_value):
                    values.setdefault(key, []).append(series_value)
                else:
                    missing[key] = missing.get(key, 0) + 1
        keys = _rank(values, overview.spans["all"])
        dropped_series = len(values) - len(keys)
        ecdfs, dropped = {}, 0
        for key in keys:
            array = np.sort(np.asarray(values[key], dtype=float))
            if array.size > ECDF_CAP:
                stride = math.ceil(array.size / ECDF_CAP)
                dropped += array.size - len(array[::stride])
                array = array[::stride]
            ecdfs[key] = array
        notes.append("Each ECDF is the empirical cumulative distribution of the raw window "
                     "distances (no binning); directions are kept separate.")
        if dropped_series > 0:
            notes.append(f"Showing the {MAX_SERIES} directions with the greatest covered span out of "
                         f"{len(values)}; the rest are omitted.")
        if dropped:
            notes.append(f"{dropped:,} values omitted by subsampling to keep ECDFs bounded.")
        pairs = len({tuple(sorted(key[:2])) for key in keys})
        stats = [("Windows", f"{sum(a.size for a in ecdfs.values()):,}"), ("Pairs", f"{pairs:,}")]
        if overview.missing:
            stats.append(("Missing", f"{overview.missing:,}"))
        if math.isfinite(overview.low):
            stats.append(("d range", f"{overview.low:.3g} – {overview.high:.3g}"))
        return Report("samples", "Sampled window distances", notes,
                      {"ecdfs": ecdfs, "missing": missing},
                      overview.selected, overview.total, heading, stats)

    # dist
    cells, rows, missing_rows = {}, 0, 0
    for row in table.rows():
        if not selection.matches(row, table.kind):
            continue
        rows += 1
        cell = row["d"]
        if not math.isfinite(cell):
            missing_rows += 1
            continue
        total, count = cells.get((row["genome_a"], row["genome_b"]), (0.0, 0))
        cells[(row["genome_a"], row["genome_b"])] = (total + cell, count + 1)
    means = {pair: total / count for pair, (total, count) in cells.items()}
    names = sorted(overview.genomes)
    symmetric = len(names) > 1 and len(overview.pairs) == len(names) * (len(names) - 1) // 2
    closest = sorted((d, a, b) for (a, b), d in means.items() if a < b)[:CLOSEST_LIMIT]
    if len(names) > MATRIX_LIMIT:
        notes.append(f"{len(names)} genomes exceed the {MATRIX_LIMIT}-genome matrix limit; showing the "
                     f"{len(closest)} closest reported comparisons instead.")
        matrix, x_names, y_names = None, names, names
    else:
        matrix = np.full((len(names), len(names)), np.nan)
        index = {name: i for i, name in enumerate(names)}
        for (a, b), d in means.items():
            matrix[index[a], index[b]] = d
            if symmetric:
                matrix[index[b], index[a]] = d
        if symmetric:
            with np.errstate(invalid="ignore"):
                row_means = np.nanmean(matrix, axis=1)
            order = np.argsort(np.where(np.isfinite(row_means), row_means, np.inf))
            matrix = matrix[np.ix_(order, order)]
            names = [names[i] for i in order]
            notes.append("Genomes are ordered by their mean distance, so similar genomes sit adjacent.")
        x_names = y_names = names
        notes.append("Matrix cells are the reported distances per genome pair"
                     + ("; repeated pairs are averaged" if rows > len(cells) else "") + ".")
    if missing_rows:
        notes.append(f"{missing_rows:,} rows have no finite distance estimate.")
    stats = [("Pairs", f"{len(means):,}"), ("Genomes", f"{len(names):,}")]
    if closest:
        d0, a0, b0 = closest[0]
        stats.append(("Closest pair", f"{short_name(a0)} ↔ {short_name(b0)} ({d0:.4g})"))
    return Report("dist", "Genome-wide distances", notes,
                  {"names": names, "x_names": x_names, "y_names": y_names, "symmetric": symmetric,
                   "matrix": matrix, "closest": closest},
                  rows, overview.total, heading, stats)
