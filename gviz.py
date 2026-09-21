#!/usr/bin/env python3
"""gviz.py - one interactive HTML report for a gdiff output file.

Reads map, roll, dist or dist --output-samples results and writes a small static
HTML page: two or three position-centred panels and a filter bar. The panels are
not interactive - a filter change redraws them, which also keeps the page fast.

    python gviz.py map.tsv
    python gviz.py roll.tsv --seq chr7 --bins 8000 --open
    python gviz.py dist.tsv --annotations genes.gff3 --inline
"""

import argparse
import base64
import json
import math
import re
import shlex
import sys
import webbrowser
from pathlib import Path

import numpy as np
import pandas as pd

D_EPS = 1e-5
# muted lichen tones, used whenever colour encodes a series rather than a value
PALETTE = ["#6E8B7B", "#B08D57", "#7E9BB5", "#A9714B", "#8A8F5C", "#9C7C93",
           "#4F7A72", "#C2A26B", "#7B8FA1", "#8C8375"]

# legacy column names -> canonical ones
ALIASES = {
    "qid": "seq", "query_id": "seq", "refseq_accn": "seq", "l": "seq_len",
    "length": "seq_len", "interval_start": "start", "interval_end": "end",
    "ref_id": "reference", "dist": "d", "dist_th": "d", "d_interval": "d_bin",
    "dist_contig": "d_q", "dist_genome": "d_acc", "strand_diff": "d_diff",
    "i": "info", "distance": "d", "median": "d_median", "alternative_mean": "d_mean",
    "max_unfiltered_distance": "d_upper", "max_distance": "d_highest",
    "num_na": "n_na", "num_filtered": "n_filtered", "n_lr_zero": "n_ub",
}
NUM = {"seq_len", "start", "end", "d", "mask", "d_q", "d_acc", "d_diff", "percentile",
       "fold", "qvalue", "info", "lr_ub", "d_fw", "d_rc", "is_rc", "d_median",
       "d_mean", "d_upper", "d_highest", "d_ab", "d_ba", "n_ab", "n_ba", "n_ub",
       "n_na", "n_filtered"}
CAT = {"seq", "reference", "config", "direction", "strand", "genome_a", "genome_b"}


def kind_of(cols):
    if "config" in cols or ("direction" in cols and "lr_ub" in cols):
        return "samples"
    if "d_ab" in cols and ("n_na" in cols or "n_filtered" in cols):
        return "dist"
    if "mask" in cols or "d_bin" in cols:
        return "map"
    if "d_fw" in cols and "d_rc" in cols:
        return "roll"
    if "reference" in cols and "d" in cols:
        return "roll"
    raise SystemExit(f"gviz: unrecognised column set: {' '.join(sorted(cols))}")


def read_table(path, kind=None):
    """Read a gdiff TSV into (frame, kind, run info)."""
    with open(path, errors="replace") as fh:
        head = []
        for line in fh:
            head.append(line)
            if len(head) >= 200:
                break
    run = {"cmd": "", "version": ""}
    for line in head:
        if line.startswith("# invocation:"):
            run["cmd"] = line.split(":", 1)[1].strip()
        elif line.startswith("# version:"):
            run["version"] = line.split(":", 1)[1].strip()
    body = [ln for ln in head if not ln.startswith("#") and ln.strip()]
    sep = "\t" if "\t" in body[0] else r"\s+"

    def isnum(t):
        try:
            float(t)
            return True
        except ValueError:
            return False

    first = re.split(sep, body[0].strip())
    if any(isnum(f) for f in first[1:4]):        # old `detect` output: 14 unnamed columns
        names = ["seq", "seq_len", "start", "end", "reference", "d", "legacy_dir", "_xl",
                 "dist_th", "d_bin", "d_q", "info", "_n_kmers", "_n_hits"]
        df = pd.read_csv(path, sep=sep, comment="#", header=None, names=names,
                         engine="python" if sep != "\t" else "c")
        kind = "map"
    else:
        cols = [c.strip().lower() for c in first]
        kind = kind or kind_of(set(cols))
        df = pd.read_csv(path, sep=sep, comment="#", engine="python" if sep != "\t" else "c")
        df.columns = [c.strip().lower() for c in df.columns]
        df = df.rename(columns={c: ALIASES.get(c, c) for c in df.columns})
    for col in df.columns:
        if col in NUM:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("float64")
        elif col in CAT:
            df[col] = df[col].astype(str)
    run["opts"] = parse_flags(run["cmd"])
    return df, kind, run


def parse_flags(cmd):
    """The few gdiff options the report wants to quote."""
    out, toks = {}, shlex.split(cmd) if cmd else []
    for i, tok in enumerate(toks):
        if tok in {"-l", "-s", "-b", "-d", "--levels"} and i + 1 < len(toks):
            j, vals = i + 1, []
            while j < len(toks) and not toks[j].startswith("-"):
                vals.append(toks[j])
                j += 1
                if tok not in {"-d", "--levels"}:
                    break
            out.setdefault(tok, []).extend(vals)
        elif tok in ("--enum-only", "--per-sequence"):
            out[tok.lstrip("-")] = True
    return out


# ------------------------------------------------------------------- metadata

def read_report(path):
    """assembly_report.txt -> [(accession, label, role, length)] in file order."""
    rows = []
    for line in open(path, errors="replace"):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) >= 10 and p[0].strip().lower() != "sequence-name":
            accn = p[4] if p[4] not in ("", "na", "<>") else p[0]
            label = f"chr{p[2]}" if p[3] == "Chromosome" and p[2] not in ("", "na") else accn
            rows.append((accn, label, p[1], int(p[8]) if p[8].isdigit() else 0))
    return rows


def find_report(path, cmd):
    for token in (cmd.split()[1:4] if cmd else []):
        stem = Path(token).name.split(".")[0]
        if not stem or stem in {"map", "roll", "dist"}:
            continue
        for d in (path.parent, Path.cwd()):
            hits = sorted(d.glob(f"{stem}*assembly_report.txt"))
            hits += sorted(d.glob(f"*/{stem}*assembly_report.txt"))
            if len(hits) == 1:
                return hits[0]
    return None


def read_genes(path):
    keep = []
    for line in open(path, errors="replace"):
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) >= 9 and p[3].isdigit():
            m = re.search(r"(?:^|;)(?:Name|gene_name|gene|ID)=([^;]+)", p[8])
            keep.append((p[0], int(p[3]), int(p[4]), m.group(1)[:40] if m else p[2]))
    return pd.DataFrame(keep, columns=["seq", "start", "end", "name"])


def add_ref(df):
    """Short reference labels; mapping the unique values keeps this O(1) per row."""
    if "reference" in df:
        df["ref"] = df["reference"].map({v: short(v) for v in df["reference"].unique()})
    return df


def short(name):
    """DOM_WSB.GCA_001624835.1.fna.gz -> DOM_WSB"""
    base = str(name)
    for pattern in (r"\.(gz|bgz|zst)$", r"\.(fna|fa|fasta|fas|fq|fastq)$",
                    r"\.(GC[AF]_\d+(\.\d+)?)$", r"\.(skc|gs|gdsk)$"):
        base = re.sub(pattern, "", base, flags=re.I)
    return base or str(name)


def label_seqs(df, report, order):
    """Replace accessions with chromosome names, keeping the assembly lane order."""
    if not report or "seq" not in df:
        return
    names = {r[0]: r[1] for r in report}
    df["seq"] = df["seq"].map(names).fillna(df["seq"])
    seen = set(df["seq"])
    order.extend([r[1] for r in report if r[1] in seen and r[1] not in order])


# -------------------------------------------------------------------- prepare

def prepare_map(df, run, args):
    df["length"] = (df["end"] - df["start"] + 1).clip(lower=0)
    df["pos"] = (df["start"] + df["end"]) / 2e6                  # Mb inside the sequence
    if "dist_th" in df:
        df["d"] = df["d"].fillna(df["dist_th"])
    lo = df["d_bin"].astype(str).str.extract(r"[\(\[]\s*([-\d.eE+]+)\s*,")[0] if "d_bin" in df else None
    hi = df["d_bin"].astype(str).str.extract(r",\s*([-\d.eE+]+)\s*[\)\]]")[0] if "d_bin" in df else None
    df["d_lo"] = pd.to_numeric(lo, errors="coerce") if lo is not None else np.nan
    df["d_hi"] = pd.to_numeric(hi, errors="coerce") if hi is not None else np.nan
    df["d_plot"] = df["d"].fillna(np.sqrt(df["d_lo"] * df["d_hi"]))
    if "d_q" in df:
        df["d_plot"] = df["d_plot"].fillna(df["d_q"])

    if "strand" in df and (df["strand"] != ".").any():
        df["dir"] = df["strand"].map({"+": "closer", "-": "farther"}).fillna("unknown")
    elif "d_q" in df:
        df["dir"] = np.where(df["d_plot"] < df["d_q"], "closer", "farther")
    else:
        df["dir"] = "unknown"
    if "dir" in df and "ref" in df:              # one row per reference and direction
        df["lane"] = df["ref"] + " \u00b7 " + df["dir"]
    df["tested"] = df["percentile"].notna() if "percentile" in df else False
    if "qvalue" in df:
        df["nq"] = -np.log10(df["qvalue"].clip(lower=1e-300))

    # MASK bits, and the threshold ladder they encode (the values are not in the output)
    mask = df["mask"].fillna(0).to_numpy("int64") if "mask" in df else np.zeros(len(df), "int64")
    nbits = int(mask.max()).bit_length() if len(mask) else 0
    bit = np.full(len(df), -1, dtype="int64")
    for b in range(nbits):
        bit[((mask >> b) & 1).astype(bool) & (bit < 0)] = b
    df["bit"] = bit
    thr = {}
    for ref, g in df.groupby("ref"):
        vals = []
        for b in range(nbits):
            rows = g[g["bit"] == b]
            if rows.empty:
                continue
            edge = rows["d_hi"] if (rows["dir"] == "closer").any() else rows["d_lo"]
            edge = edge[np.isfinite(edge)]
            if len(edge):
                vals.append(float(edge.median()))
        thr[ref] = sorted(v for v in set(vals) if D_EPS < v < 1.0)
    level = {}
    for b in range(nbits):
        rows = df[df["bit"] == b]
        ref = rows["ref"].iloc[0] if len(rows) else ""
        val = next((v for i, v in enumerate(thr.get(ref, [])) if i == b), np.nan)
        side = "farther" if len(rows) and (rows["dir"] == "farther").mean() > 0.5 else "closer"
        sign = "\u2264" if side == "closer" else ">"
        level[b] = (f"{b + 1} \u00b7 {side} " + (f"d {sign} {val:.3g}" if np.isfinite(val)
                                                  else "than background"))
    df["level"] = df["bit"].map(level).fillna("no threshold")
    args.note("the thresholds themselves are not in the output, so every interval is labelled "
              "with the value recovered from the bracket it fell into: the lower half is closer "
              "to the reference than the background, the upper half farther")

    if run["opts"].get("enum-only"):
        args.note("this run enumerated intervals per threshold, so intervals may overlap")
    if nbits and not df["tested"].any():
        args.note("no interval carries a percentile: significance was not tested in this run")
    return df, thr, nbits


def prepare_roll(df, bins, args):
    if "d_fw" in df and "d_rc" in df:
        df["d"] = df[["d_fw", "d_rc"]].min(axis=1)
    df["pos"] = (df["start"] + df["end"]) / 2e6
    df["length"] = df["end"] - df["start"] + 1
    if bins and len(df) > bins:
        # one sequence is shown at a time, so bin each sequence on its own scale
        seq = pd.factorize(df["seq"])[0]
        lo = df["pos"].groupby(seq).transform("min")
        span = df["pos"].groupby(seq).transform("max") - lo
        width = np.maximum(max(float(df["length"].median()) / 1e6, 1e-6), span / bins)
        slot = ((df["pos"] - lo) / width).astype("int64").clip(upper=99_999)
        nref = df["reference"].nunique()
        cell = (seq * nref + pd.factorize(df["reference"])[0]) * 100_000 + slot
        df = df.assign(nan=df["d"].isna())
        g = df.groupby(cell, sort=True)
        df = pd.DataFrame({
            "seq": g["seq"].first(), "reference": g["reference"].first(),
            "pos": g["pos"].mean(), "d": g["d"].mean(), "nan_frac": g["nan"].mean(),
        }).reset_index(drop=True)
        df = add_ref(df)
        args.note(f"each sequence is averaged into at most {bins:,} bins "
                  f"({float(width.median()) * 1e3:,.0f} kb median width); "
                  "--bins 0 keeps every single window")
    else:
        bins = 0
    return df, bins


def prepare_dist(df, args):
    df["a"] = df["genome_a"].map(short)
    df["b"] = df["genome_b"].map(short)
    df["pair"] = df["a"] + " vs " + df["b"]
    df["asym"] = df["d_ab"] - df["d_ba"] if {"d_ab", "d_ba"} <= set(df) else np.nan
    counts = [c for c in ("n_ab", "n_ba", "n_ub", "n_na", "n_filtered") if c in df]
    total = df[counts].sum(axis=1) if counts else pd.Series(1.0, index=df.index)
    for c in counts:
        df["f_" + c[2:]] = df[c] / total
    keys = {tuple(sorted(p)) for p in zip(df["genome_a"], df["genome_b"])}
    n = len(set(df["genome_a"]) | set(df["genome_b"]))
    form = "all-by-all" if len(keys) == n * (n - 1) // 2 else "cross product"
    # rows and columns share one order, nearest genomes adjacent, so the matrix reads
    names = sorted(set(df["a"]) | set(df["b"]),
                   key=lambda g: float(np.nanmedian(df.loc[(df["a"] == g) | (df["b"] == g), "d"])))
    args.note(f"this file is a {form}"
              + (", each unordered pair written once" if form == "all-by-all" else "")
              + "; genomes are ordered by their median distance")
    return df, form, names


def prepare_samples(df, args):
    # the two directions of a pair are the same comparison, so label them together
    df["a"] = df["genome_a"].map(short)
    df["b"] = df["genome_b"].map(short)
    lo = df[["a", "b"]].min(axis=1)
    df["pair"] = lo + " vs " + df[["a", "b"]].max(axis=1)
    df["pos"] = (df["start"] + df["end"]) / 2e6
    args.note("GENOME_A is the side whose windows are scored, GENOME_B the other side")
    return df


# -------------------------------------------------------------------- filters

FILTERS = {
    "map": [("seq", "seq", "Sequence", "one"), ("ref", "ref", "Reference", "cat"),
            ("dir", "dir", "Direction", "cat"),
            ("level", "level", "Threshold", "cat"),
            ("d", "d_plot", "Distance d", "log"), ("len", "length", "Length (bp)", "log"),
            ("q", "qvalue", "q-value", "lin")],
    "roll": [("seq", "seq", "Sequence", "one"), ("ref", "ref", "Reference", "cat"),
             ("d", "d", "Distance d", "log")],
    "dist": [("pair", "pair", "Pair", "cat"), ("d", "d", "Distance d", "log"),
             ("na", "f_na", "Unmapped fraction", "lin")],
    "samples": [("seq", "seq", "Sequence", "one"), ("pair", "pair", "Pair", "cat"),
                ("dir", "direction", "Direction", "cat"), ("d", "d", "Distance d", "log")],
}


def make_filters(df, kind, seq_order, sizes):
    out = []
    for fid, col, label, typ in FILTERS[kind]:
        if col not in df:
            continue
        if typ == "one":
            vals = [str(v) for v in pd.unique(df[col].dropna())]
            vals.sort(key=lambda v: -float(sizes.get(v, 0)))
            if len(vals) < 2:
                continue
            out.append({"id": fid, "col": col, "label": label, "type": "one",
                        "values": vals, "default": vals[0],
                        "hint": "one sequence at a time"})
        elif typ == "cat":
            vals = sorted(str(v) for v in pd.unique(df[col].dropna()))
            if len(vals) < 2:
                continue
            out.append({"id": fid, "col": col, "label": label, "type": "cat",
                        "values": vals, "default": vals})
        else:
            v = df[col].to_numpy("float64")
            v = v[np.isfinite(v) & ((v > 0) if typ == "log" else True)]
            if v.size == 0:
                continue
            out.append({"id": fid, "col": col, "label": label, "type": "range",
                        "lo": float(v.min()), "hi": float(v.max()), "log": typ == "log"})
    return out


# ------------------------------------------------------------------- encoding

def b64(a):
    return base64.b64encode(np.ascontiguousarray(a).tobytes()).decode()


def enc(series, col, cats=None, n=0):
    a = series.to_numpy()
    if cats is not None:
        code = pd.Categorical(a, categories=cats).codes.astype("<i2")
        return {"t": "cat", "dt": "i2", "d": b64(code), "cats": cats}
    if col in ("start", "end") or n < 25_000:     # float32 costs nothing visible below that
        return {"t": "num", "dt": "f8", "d": b64(a.astype("<f8"))}
    return {"t": "num", "dt": "f4", "d": b64(a.astype("<f4"))}


def columns_of(panels, filters):
    cols = set()
    for pnl in panels:
        for lay in pnl["layers"]:
            for key in ("x", "x2", "y", "lane", "split", "color", "size", "z"):
                if lay.get(key) and lay[key] not in ("none", "split"):
                    cols.add(lay[key])
            cols.update(lay.get("cols") or [])
        cols.update(s.split(":")[-1] for s in pnl["stats"] if ":" in s)
    cols.update(f["col"] for f in filters)
    return cols


def table_json(df, cols, order):
    out = {}
    for col in sorted(cols):
        if col not in df:
            continue
        cats = order.get(col)
        if cats is None and not pd.api.types.is_numeric_dtype(df[col]):
            cats = sorted(map(str, pd.unique(df[col].dropna())))
        out[col] = enc(df[col], col, cats, len(df))
    return {"n": int(len(df)), "cols": out}


# --------------------------------------------------------------------- panels

def P(pid, title, desc, layers, x, y, h=380, stats=("n",), controls=(), notes="", table="rows",
      width=None):
    x = dict(x, col=layers[0].get("x") if layers else None)
    return {"id": pid, "title": title, "desc": desc, "layers": list(layers), "x": x, "y": y,
            "h": h, "stats": list(stats), "controls": list(controls), "notes": notes,
            "table": table, "width": width}


def C(title, value, hint=""):
    return {"title": title, "value": value, "hint": hint}


def map_panels(d, thr, nbits, refs):
    """What the map output says about where: position on x, one sequence at a time."""
    lanes = [x for x in (f"{r} \u00b7 {s}" for r in refs for s in ("closer", "farther"))
             if (d["lane"] == x).any()] or sorted(set(d["lane"]))
    return [P("track", "Interval track",
              "One row per reference and direction: the closer rows hold intervals that are "
              "closer to that reference than its genome background, the farther rows the "
              "opposite. Marker area grows with interval length.",
              [{"k": "points", "x": "pos", "lane": "lane", "split": "ref", "color": "ref",
                "size": "length", "opacity": 0.8}],
              {"t": "Position (Mb)"}, {"t": "reference and direction", "cats": lanes},
              46 + 21 * max(len(lanes), 1), stats=("n", "median:d_plot", "max:length")),
            P("hits", "Significance along the sequence",
              "The same intervals with -log10(q) on y; the dashed line is q = 0.05.",
              [{"k": "points", "x": "pos", "y": "nq", "split": "ref", "color": "ref",
                "size": "length", "opacity": 0.8, "hline": 1.301}],
              {"t": "Position (Mb)"}, {"t": "-log10 q-value"}, 280,
              stats=("n", "median:nq", "max:nq"))]


def gene_panel(genes):
    return P("genes", "Annotation track",
             "Features from the annotation file at their positions, marker area by feature "
             "length. This panel ignores the row filters because it is not gdiff output.",
             [{"k": "points", "x": "pos", "lane": "seq", "color": "seq", "size": "length",
               "opacity": 0.9, "table": "genes"}],
             {"t": "Position (Mb)"},
             {"t": "sequence", "cats": sorted(map(str, genes["seq"].unique()))},
             46 + 22 * genes["seq"].nunique(), stats=(), table="genes")


def roll_panels(d, refs, binned):
    """Position on x, plain distance on y, one line per reference."""
    out = [P("profile", "Distance profile",
             "Window distance along the sequence, one line per reference"
             + (f", averaged into {binned:,} bins and smoothed with a running median"
                if binned else ", smoothed with a running median")
             + ". Windows with no k-mer hit are left out instead of being drawn as gaps.",
             [{"k": "lines", "x": "pos", "y": "d", "split": "ref"}],
             {"t": "Position (Mb)"}, {"t": "distance d", "log": True}, 420,
             stats=("n", "median:d", "max:d"),
             controls=[{"label": "smoothing", "target": "smooth", "layer": 0, "num": 5},
                       {"label": "log y", "target": "ylog", "layer": 0, "default": True}]),
           P("landscape", "Reference x position",
             "The same profile as a heat map: rows are references, colour is the window distance, "
             "so a single-row stripe is a segment specific to one reference.",
             [{"k": "heatmap", "x": "pos", "lane": "ref", "z": "d", "agg": "mean",
               "bins": 600, "scale": "Viridis"}],
             {"t": "Position (Mb)"}, {"t": "reference", "cats": refs},
             70 + 40 * len(refs), stats=("n", "median:d"))]
    if "nan_frac" in d and d["nan_frac"].max() > 0:
        out.append(P("unmapped", "Windows with no k-mer hit",
                     "Share of windows per bin with no hit at all: gaps, N runs or sequence too "
                     "divergent to match. Those stretches are invisible above.",
                     [{"k": "heatmap", "x": "pos", "lane": "ref", "z": "nan_frac", "agg": "mean",
                       "bins": 600, "scale": "Reds", "zmin": 0, "zmax": 1}],
                     {"t": "Position (Mb)"}, {"t": "reference", "cats": refs},
                     70 + 40 * len(refs), stats=("n", "mean:nan_frac")))
    return out


def dist_panels(d, form, order):
    a = b = order
    pairs = sorted(d["pair"].unique(), key=lambda s: d.loc[d["pair"] == s, "d"].median())
    stats = [c for c in ("d", "d_median", "d_mean", "d_ab", "d_ba", "d_upper", "d_highest")
             if c in d]
    out = [P("matrix", "Pairwise distance matrix",
             "Every genome against every other one; each tile holds the median of the selected "
             "statistic, and empty tiles are pairs this file does not contain"
             + (" (each unordered pair is written once, so the matrix is mirrored)."
                if form == "all-by-all" else "."),
             [{"k": "heatmap", "x": "b", "y": "a", "z": "d", "agg": "median", "scale": "Viridis",
               "text": len(a) * len(b) <= 400, "sym": form == "all-by-all", "tiles": True}],
             {"t": "genome" if form == "all-by-all" else "genome B"},
             {"t": "genome" if form == "all-by-all" else "genome A", "cats": a},
             54 + 76 * max(len(a), 1), width=120 + 78 * max(len(a), 1),
             stats=("n", "min:d", "median:d", "max:d"),
             controls=[{"label": "value", "target": "z", "layer": 0, "default": "d",
                        "opts": stats}])]
    if len(stats) > 1:
        out.append(P("spread", "Estimator spread per pair",
                     "Every distance statistic on its own row per pair, sorted by d. Points that "
                     "sit far apart mean the pair's windows disagree, so the point estimate carries "
                     "little information.",
                     [{"k": "melt", "lane": "pair", "cols": stats}],
                     {"t": "distance", "log": True}, {"t": "pair", "cats": pairs},
                     120 + 26 * len(pairs), stats=("n", "median:d")))
    counts = [c for c in ("n_ab", "n_ba", "n_ub", "n_na", "n_filtered") if c in d]
    if counts:
        out.append(P("accounting", "Window accounting per pair",
                     "What happened to every sampled window of a pair: valid on either side, "
                     "unmapped, at the upper bound, or rejected by the likelihood-ratio filter.",
                     [{"k": "bars", "x": "pair", "cols": counts, "agg": "sum",
                       "labels": {"n_ab": "valid a>b", "n_ba": "valid b>a", "n_ub": "upper bound",
                                  "n_na": "unmapped", "n_filtered": "filtered by LR"}}],
                     {"t": "pair"}, {"t": "windows"}, 380, stats=("n",),
                     controls=[{"label": "fractions", "target": "frac", "layer": 0,
                                "default": False}]))
    return out


def samples_panels(d):
    return [P("profile", "Window distances along their sequence",
              "Every stored window at its position in the sequence it came from; a row belongs to "
              "the genome in GENOME_A, so filter to one pair to read it as a profile.",
              [{"k": "points", "x": "pos", "y": "d", "split": "pair", "color": "pair",
                "opacity": 0.45}],
              {"t": "Position (Mb)"}, {"t": "distance d", "log": True},
              400, stats=("n", "median:d", "max:d"),
              controls=[{"label": "split", "target": "split", "layer": 0, "default": "pair",
                         "opts": ["pair", "direction", "none"]},
                        {"label": "log y", "target": "ylog", "layer": 0, "default": True}]),
            P("ecdf", "Per-pair distribution",
              "ECDF of the window distances of each pair: the median is the pair's estimate, the "
              "tail is what the reconciliation filter decides about.",
              [{"k": "ecdf", "x": "d", "split": "pair"}],
              {"t": "distance d", "log": True}, {"t": "cumulative share"}, 380,
              stats=("n", "median:d", "max:d"))]


# --------------------------------------------------------------------- report

class Args:
    """Collects the notes that end up in the report header."""

    def __init__(self, **kw):
        self.__dict__.update(kw)
        self._notes = []

    def note(self, text):
        self._notes.append(text)

    @property
    def notes(self):
        return self._notes


def cards_for(df, kind, form="", n_raw=0):
    if kind == "map":
        out = [C("Intervals", f"{len(df):,}", f"{df['ref'].nunique()} references"),
               C("Tested", f"{100 * df['tested'].mean():.1f}%" if "tested" in df else "n/a",
                 "the rest are intact or fallback rows")]
        if "qvalue" in df and df["qvalue"].notna().any():
            out.append(C("q \u2264 0.05", f"{int((df['qvalue'] <= 0.05).sum()):,}",
                         "Benjamini-Hochberg, per strand"))
        out.append(C("Median interval", f"{df['length'].median() / 1e3:,.1f} kb",
                     f"longest {df['length'].max() / 1e3:,.0f} kb"))
    elif kind == "roll":
        out = [C("Windows", f"{n_raw or len(df):,}",
                 f"{df['ref'].nunique()} references"
                 + (f" · shown as {len(df):,} bins" if n_raw and n_raw != len(df) else "")),
               C("Median d", f"{df['d'].median():.4g}", "over all windows")]
        nan = float(df["d"].isna().mean())
        if nan:
            out.append(C("Unmapped", f"{100 * nan:.1f}%", "no k-mer hit in the window"))
    elif kind == "dist":
        out = [C("Pairs", f"{len(df):,}", form),
               C("Median d", f"{df['d'].median():.4g}",
                 f"range {df['d'].min():.4g} - {df['d'].max():.4g}")]
        if "d_mean" in df:
            out.append(C("Filter effect", f"{(df['d_mean'] - df['d']).median():+.4g}",
                         "median d_mean - d, what the LR filter removes"))
        if "f_na" in df:
            out.append(C("Worst unmapped", f"{100 * df['f_na'].max():.1f}%",
                         "share of windows without a hit for one pair"))
    else:
        out = [C("Windows", f"{len(df):,}", "one row per stored window"),
               C("Pairs", f"{df['pair'].nunique():,}", "as written in the file"),
               C("Median d", f"{df['d'].median():.4g}", "over all windows")]
    return out


def chips(run, kind, n):
    o = run["opts"]
    out = [(kind, f"{n:,} rows")] if not kind == "roll" else []
    if o.get("--levels"):
        out.append(("levels", " ".join(o["--levels"])))
    if o.get("-d"):
        out.append(("-d", " ".join(o["-d"])))
    if o.get("enum-only"):
        out.append(("enum-only", "per-threshold enumeration"))
    if o.get("per-sequence"):
        out.append(("per-sequence", "background per query"))
    for flag, label in (("-l", "window (k-mers)"), ("-s", "step"), ("-b", "bin shift")):
        if o.get(flag):
            out.append((label, o[flag][0]))
    if run["version"]:
        out.append(("version", " ".join(run["version"].split()[:2])))
    return out


def add_ranges(panels, df, extents):
    """Pin every axis to the full data extent, so filtering never rescales a plot."""
    def span(vals, log):
        v = np.asarray(vals, dtype="float64")
        v = v[np.isfinite(v)]
        if log:
            v = v[v > 0]
        if not v.size:
            return None
        lo, hi = (float(np.quantile(v, 0.0005)), float(np.quantile(v, 0.9995))) if not log \
            else (float(v.min()), float(v.max()))
        if log:
            return [10 ** math.floor(math.log10(lo)), 10 ** math.ceil(math.log10(hi))]
        pad = (hi - lo) * 0.04 or max(abs(hi) * 0.04, 1e-9)
        return [min(lo, 0.0) - pad if lo < 0 else lo - pad, hi + pad]

    for pnl in panels:
        for key, axis in (("x", pnl["x"]), ("y", pnl["y"])):
            if axis.get("cats"):
                continue
            for lay in pnl["layers"]:
                c = lay.get(key)
                if not c or c not in df or not pd.api.types.is_numeric_dtype(df[c]):
                    continue
                if key == "x" and c == "pos":
                    axis["seq_extent"] = True          # follows the selected sequence only
                    axis["max"] = float(np.nanmax(df[c])) if len(df) else 1.0
                else:
                    r = span(df[c], bool(axis.get("log")))
                    if r:
                        axis["range"] = r
                break
        for lay in pnl["layers"]:
            if lay.get("k") != "heatmap":
                continue
            if lay.get("x") == "pos":
                hi = float(np.nanmax(df["pos"])) if len(df) else 1.0
                lay["x0"], lay["x1"] = 0.0, round(hi * 1.001, 4)
            zc = lay.get("z")
            if zc and zc in df and pd.api.types.is_numeric_dtype(df[zc]):
                v = df[zc].to_numpy("float64")
                v = v[np.isfinite(v)]
                if v.size:
                    lay["zmin"] = float(v.min())
                    lay["zmax"] = float(np.quantile(v, 0.995))
    return panels


def build(path, kind, opts):
    df, kind, run = read_table(path, kind)
    if len(df) == 0:
        raise SystemExit(f"gviz: {path} has no data rows")
    rep_path = opts.assembly_report or (None if opts.no_report else find_report(path, run["cmd"]))
    report = read_report(rep_path) if rep_path else None
    seq_order = []
    label_seqs(df, report, seq_order)
    if opts.seq and "seq" in df:
        keep = [s for s in opts.seq for s in ([s] if s in set(df["seq"]) else [])]
        keep += [s for s in df["seq"].unique() if s.startswith(tuple(opts.seq)) and s not in keep]
        if not keep:
            raise SystemExit(f"gviz: no sequence matches {', '.join(opts.seq)}")
        df = df[df["seq"].isin(keep)].reset_index(drop=True)

    df = add_ref(df)                            # the short label every prepare step uses
    sizes = {}
    if "seq" in df:
        end = (df.groupby("seq")["seq_len"].max() if "seq_len" in df
               else df.groupby("seq")["end"].max() if "end" in df
               else df.groupby("seq")["pos"].max() * 1e6)
        sizes = {str(k): float(v) for k, v in end.items()}
    if "seq" in df and not opts.keep_all:
        size = (df.groupby("seq")["seq_len"].max() if "seq_len" in df
                else df.groupby("seq")["end"].max()).sort_values(ascending=False)
        keep, acc = [], 0.0
        for name, bp in size.items():
            keep.append(name)
            acc += float(bp)
            if acc >= 0.998 * float(size.sum()):
                break
        if len(keep) < len(size):
            opts.note(f"{len(size) - len(keep)} short sequences (<0.2% of the assembly) are left "
                      "out; --keep-all to include them")
            df = df[df["seq"].isin(keep)].reset_index(drop=True)
    n_raw, thr, nbits, form = len(df), {}, 0, ""
    if kind == "map":
        df, thr, nbits = prepare_map(df, run, opts)
        panels = map_panels(df, thr, nbits, sorted(df["ref"].unique()))
    elif kind == "roll":
        df, bins = prepare_roll(df, opts.bins, opts)
        panels = roll_panels(df, sorted(df["ref"].unique()), bins)
    elif kind == "dist":
        df, form, order = prepare_dist(df, opts)
        panels = dist_panels(df, form, order)
    else:
        df = prepare_samples(df, opts)
        panels = samples_panels(df)
    if opts.annotations and kind in ("map", "roll", "samples"):
        genes = read_genes(opts.annotations)
        if report:
            genes["seq"] = genes["seq"].map({r[0]: r[1] for r in report}).fillna(genes["seq"])
        genes["pos"] = (genes["start"] + genes["end"]) / 2e6
        genes["length"] = genes["end"] - genes["start"] + 1
        genes = genes[genes["seq"].isin(set(df["seq"]) if "seq" in df else set())]
        if len(genes) > 3000:                      # a whole chromosome of genes is a smear
            genes = genes.nlargest(3000, "length")
            opts.note("the annotation track shows the 3,000 longest features; whole-chromosome "
                      "annotation is dense by nature, so use --seq to narrow it")
        if len(genes):
            panels.insert(2, gene_panel(genes))
            opts.tables = {"genes": table_json(genes, {"pos", "seq", "length", "name"},
                                               {"seq": sorted(map(str, genes["seq"].unique()))})}

    if "seq" in df:
        span = (df.groupby("seq")["end"].max() if "end" in df
                else df.groupby("seq")["pos"].max() * 1e6)
        extents = {str(k): round(float(v) / 1e6, 4) for k, v in span.items()}
    else:
        extents = {}
    panels = add_ranges(panels, df, extents)

    have = set(df.columns)

    def missing(pnl):
        cols = set()
        for lay in pnl["layers"]:
            if lay.get("table", "rows") != "rows":
                continue
            for key in ("x", "x2", "y", "lane", "split", "color", "size", "z"):
                if lay.get(key) and lay[key] not in ("none", "split"):
                    cols.add(lay[key])
            cols.update(lay.get("cols") or [])
        return cols - have

    skipped = [pnl for pnl in panels if missing(pnl)]
    panels = [pnl for pnl in panels if not missing(pnl)]
    filters = make_filters(df, kind, seq_order, sizes)
    notes = list(opts.notes)
    if rep_path:
        notes.append(f"sequence names and order come from {rep_path.name}")
    notes += [f"view skipped: {s['title']} (no {', '.join(sorted(missing(s)))})" for s in skipped]

    hints = {"seq": ([str(v) for v in df["seq"].unique()] if "seq" in df else [])}
    if kind == "dist":
        hints["a"] = hints["b"] = order
    for col in ("ref", "pair", "level", "dir"):
        if col in df:
            hints[col] = sorted(map(str, pd.unique(df[col])))
    out = {"title": opts.title or path.name, "kind": kind, "cmd": run["cmd"],
           "extents": extents,
           "chips": chips(run, kind, n_raw), "cards": cards_for(df, kind, form, n_raw),
           "filters": filters, "panels": panels, "notes": notes,
           "meta": f"{len(df):,} rows · {df['seq'].nunique() if 'seq' in df else 0:,} sequences",
           "palette": PALETTE}
    out["tables"] = {"rows": table_json(df, columns_of(panels, filters), hints)}
    out["tables"].update(getattr(opts, "tables", {}))
    return out


# ----------------------------------------------------------------------- html

CSS = """
/* plain figure document: one serif heading, a sans body, monospace only for the invocation */
html{font:14px/1.55 "Helvetica Neue",Helvetica,Arial,sans-serif;color:#222;background:#fff;
-webkit-font-smoothing:antialiased}
body{max-width:1040px;margin:0 auto;padding:40px 20px 64px}
h1{font:600 23px/1.25 Georgia,"Iowan Old Style","Times New Roman",serif;margin:0 0 6px;
letter-spacing:-.01em}
.sub{margin:0;color:#555;font-size:13.5px}
.cmd{margin:10px 0 0;font:11.5px/1.5 ui-monospace,SFMono-Regular,Menlo,monospace;color:#8c8c8c;
word-break:break-all}
#filters{display:flex;flex-wrap:wrap;gap:8px 26px;align-items:baseline;
padding:16px 0;border-top:1px solid #EAEAEA;border-bottom:1px solid #EAEAEA;margin:24px 0 20px}
.f{display:flex;align-items:baseline;gap:6px;font-size:13px;position:relative}
.f .nm{color:#666}
.f summary{cursor:pointer;list-style:none;color:#222}
.f summary::-webkit-details-marker{display:none}
select{font:inherit;font-size:13px;padding:2px 4px;border:1px solid #DCDCDC;border-radius:2px;
background:#fff;color:#222}
.sl{display:flex;align-items:center;gap:7px}
.sl input[type=range]{width:100px;margin:0;accent-color:#6E8B7B}
.sl b{font-weight:400;color:#666;font-size:12.5px;font-variant-numeric:tabular-nums;min-width:54px;
text-align:right}
.pop{position:absolute;z-index:12;top:100%;left:0;margin-top:6px;background:#fff;border:1px solid #E4E4E4;
border-radius:3px;padding:8px 10px;box-shadow:0 6px 20px rgba(0,0,0,.09);max-height:260px;
overflow:auto;font-size:13px;min-width:120px}
.pop label{display:block;white-space:nowrap;padding:1px 0}
.pop .row{display:flex;gap:10px;padding-bottom:6px}
button{font:inherit;font-size:13px;background:none;border:none;color:#3A6B8F;cursor:pointer;
padding:0;margin:0}
button:hover{text-decoration:underline}
.pop button{font-size:12.5px}
#summary{margin:0 0 32px;color:#555;font-size:13.5px}
#summary b{font-weight:400;color:#222}
figure{margin:0 0 36px}
figcaption{margin:0 0 10px;font-size:13.5px;color:#555;max-width:92ch;line-height:1.5}
figcaption b{display:block;font-size:15px;font-weight:600;color:#111;margin-bottom:3px}
.pctl{display:flex;gap:18px;align-items:baseline;margin:0 0 8px;font-size:13px;color:#555}
.pctl label{margin-right:2px}
.pctl input[type=number]{font:inherit;font-size:13px;width:54px;padding:1px 3px;border:1px solid #DCDCDC;
border-radius:2px}
.pstat{margin:8px 0 0;font-size:12.5px;color:#8c8c8c}
footer{border-top:1px solid #EAEAEA;margin-top:12px;padding-top:18px;color:#666;font-size:13px}
footer ul{margin:8px 0 0 18px;padding:0}footer li{margin:4px 0}
"""

JS = r"""
"use strict";
const R = __REPORT__, PAL = R.palette;
const GG = {panel:"#EBEBEB", grid:"#FFFFFF", ink:"#333333", tick:"#4D4D4D", font:"Helvetica Neue,Helvetica,Arial,sans-serif"};
// static panels: no zoom, no pan, no hover - the filters are the only interaction
const CFG = {responsive:true, displaylogo:false, staticPlot:false, scrollZoom:false,
  doubleClick:false, displayModeBar:true, modeBarButtons:[["toImage"]],
  toImageButtonOptions:{format:"png", scale:2}};
const TB = {}, ST = {f:{}, c:{}, mask:null, rows:null, timer:null};
const CAP = 12000;                      // points per trace; more gets strided

function tb(name){
  if (TB[name]) return TB[name];
  const spec = R.tables[name], out = {n:spec.n, c:{}};
  for (const k in spec.cols){
    const cs = spec.cols[k], bin = atob(cs.d), u = new Uint8Array(bin.length);
    for (let i=0;i<bin.length;i++) u[i] = bin.charCodeAt(i);
    out.c[k] = {a: cs.dt==="f8" ? new Float64Array(u.buffer) : cs.dt==="f4" ? new Float32Array(u.buffer)
      : cs.dt==="i2" ? new Int16Array(u.buffer) : new Int32Array(u.buffer), cats: cs.cats || null};
  }
  return TB[name] = out;
}
const col = (n,t) => (n && tb(t||"rows").c[n]) || null;
const isCat = n => { const c = col(n); return !!(c && c.cats); };
function fmt(v){ if(!Number.isFinite(v)) return "n/a"; const a = Math.abs(v);
  return a && (a<1e-3 || a>=1e5) ? v.toExponential(2) : String(Math.round(v*1e5)/1e5); }
function lab(n,i,t){ const c = col(n,t); if(!c) return "";
  if(!c.cats) return fmt(c.a[i]); const k = c.a[i];
  if (!(k>=0 && k<c.cats.length)) return "(none)";
  const v = c.cats[k]; return v.length > 26 ? v.slice(0, 25) + "\u2026" : v; }
function pick(n,key,t){ const c = col(n,t); if(!c||!c.cats) return PAL[0];
  const i = c.cats.indexOf(key); return PAL[(i<0?0:i) % PAL.length]; }
function colour(L,key,t){
  if (L.color && isCat(L.color)) return pick(L.color, key, t);
  // no colour column given: split series still need distinct hues
  if (L.split && L.split !== "none" && isCat(L.split)) return pick(L.split, key, t);
  return PAL[0];
}
function groups(split, r, t){
  const g = new Map();
  if (!split || split === "none"){ g.set("", r); return g; }
  const c = col(split, t);
  if (!c || !c.cats){ g.set("", r); return g; }
  for (const v of c.cats) g.set(v, []);
  for (const i of r){ const k = c.a[i]; if (k>=0 && k<c.cats.length) g.get(c.cats[k]).push(i); }
  for (const [k,v] of [...g]) if (!v.length) g.delete(k);
  return g;
}
function thin(r){
  if (r.length <= CAP) return r;
  const n = Math.ceil(r.length / Math.ceil(r.length / CAP));
  const out = new Int32Array(n), step = Math.ceil(r.length / n);
  for (let i=0,j=0;i<r.length && j<n;i+=step) out[j++] = r[i];
  return out;
}

/* ------------------------------------------------------------------ filters */
function filtUI(){
  const host = document.getElementById("filters");
  host.innerHTML = "";
  for (const f of R.filters){
    const box = document.createElement("div"); box.className = "f";
    if (f.type === "one"){
      ST.f[f.id] = ST.f[f.id] || f.default;
      const sel = document.createElement("select");
      for (const v of f.values){ const o = document.createElement("option");
        o.value = v; o.textContent = v; o.selected = v === ST.f[f.id]; sel.appendChild(o); }
      sel.onchange = () => { ST.f[f.id] = sel.value; schedule(); };
      const step = d => {
        const i = Math.min(Math.max(sel.selectedIndex + d, 0), sel.options.length - 1);
        if (i !== sel.selectedIndex){ sel.selectedIndex = i; sel.onchange(); }
      };
      const prev = document.createElement("button"), next = document.createElement("button");
      prev.textContent = "\u2039"; next.textContent = "\u203a";
      prev.title = "previous sequence"; next.title = "next sequence";
      prev.onclick = () => step(-1); next.onclick = () => step(1);
      const nm = document.createElement("span"); nm.className = "nm"; nm.textContent = f.label;
      box.append(nm, prev, sel, next);
      if (f.hint) box.title = f.hint;
    } else if (f.type === "cat"){
      const st = ST.f[f.id] = new Set(ST.f[f.id] || f.default);
      const det = document.createElement("details");
      const sum = document.createElement("summary");
      const label = () => f.label + " " + st.size + "/" + f.values.length;
      sum.textContent = label(); det.appendChild(sum);
      const pop = document.createElement("div"); pop.className = "pop";
      const row = document.createElement("div"); row.className = "row";
      const all = document.createElement("button"); all.textContent = "all";
      const none = document.createElement("button"); none.textContent = "none";
      row.append(all, none); pop.appendChild(row);
      const boxes = [];
      for (const v of f.values){
        const l = document.createElement("label"), cb = document.createElement("input");
        cb.type = "checkbox"; cb.checked = st.has(v);
        cb.onchange = () => { cb.checked ? st.add(v) : st.delete(v);
          sum.textContent = label(); schedule(); };
        l.append(cb, document.createTextNode(" " + v)); pop.appendChild(l); boxes.push(cb);
      }
      all.onclick = e => { e.preventDefault(); f.values.forEach(v=>st.add(v));
        boxes.forEach(b=>b.checked=true); sum.textContent = label(); schedule(); };
      none.onclick = e => { e.preventDefault(); st.clear();
        boxes.forEach(b=>b.checked=false); sum.textContent = label(); schedule(); };
      det.appendChild(pop); box.appendChild(det);
      if (f.hint) box.title = f.hint;
    } else {
      const isLog = f.log, t = v => isLog ? Math.log10(Math.max(v,1e-12)) : v,
            inv = v => isLog ? Math.pow(10,v) : v, a = t(f.lo), b = t(f.hi);
      const cur = ST.f[f.id];
      const s1 = document.createElement("input"), s2 = document.createElement("input");
      for (const s of [s1,s2]){ s.type="range"; s.min=a; s.max=b; s.step=(b-a)/200 || 1e-9; }
      s1.value = cur ? t(cur[0]) : a;
      s2.value = cur ? t(cur[1]) : b;
      const read = document.createElement("b");
      read.textContent = cur ? fmt(cur[0]) + "\u2013" + fmt(cur[1]) : "all";
      const sync = () => {
        let lo = +s1.value, hi = +s2.value;
        if (lo > hi){ const m = lo; lo = hi; hi = m; }
        s1.value = lo; s2.value = hi;
        const full = lo <= a + 1e-12 && hi >= b - 1e-12;
        ST.f[f.id] = full ? null : [inv(lo), inv(hi)];
        read.textContent = full ? "all" : fmt(ST.f[f.id][0]) + "\u2013" + fmt(ST.f[f.id][1]);
        schedule();
      };
      s1.oninput = sync; s2.oninput = sync;
      const wrap = document.createElement("div"); wrap.className = "sl";
      const nm = document.createElement("span"); nm.className = "nm"; nm.textContent = f.label;
      wrap.append(nm, s1, s2, read); box.appendChild(wrap);
    }
    host.appendChild(box);
  }
}
function mask(){
  const t = tb("rows"); let m = null;
  for (const f of R.filters){
    const st = ST.f[f.id];
    if (!st || (st instanceof Set && st.size === f.values.length)) continue;
    const c = col(f.col); if (!c) continue;
    if (!m) m = new Uint8Array(t.n).fill(1);
    if (typeof st === "string"){
      const want = c.cats ? c.cats.indexOf(st) : -2;
      for (let i=0;i<t.n;i++) if (m[i] && c.a[i] !== want) m[i] = 0;
    } else if (st instanceof Set){
      const keep = new Uint8Array(c.cats.length);
      c.cats.forEach((v,i)=>{ if (st.has(v)) keep[i] = 1; });
      for (let i=0;i<t.n;i++) if (m[i] && !(c.a[i]>=0 && keep[c.a[i]])) m[i] = 0;
    } else {
      for (let i=0;i<t.n;i++) if (m[i] && !(c.a[i]>=st[0] && c.a[i]<=st[1])) m[i] = 0;
    }
  }
  return m;
}
function rowsOf(m, name){
  const n = tb(name||"rows").n;
  if (!m) { const a = new Int32Array(n); for (let i=0;i<n;i++) a[i]=i; return a; }
  const o = []; for (let i=0;i<n;i++) if (m[i]) o.push(i);
  return Int32Array.from(o);
}

/* ------------------------------------------------------------- trace builders */
function shapes(L, ex){
  if (L.k === "vlines") (L.pos||[]).forEach(x => ex.shapes.push({type:"line", xref:"x", yref:"paper",
    x0:x, x1:x, y0:0, y1:1, line:{color:"#9a9a9a", width:1, dash:"dot"}}));
  if (L.k === "hlines") (L.values||[]).forEach(v => {
    ex.shapes.push({type:"line", xref:"paper", yref:"y", x0:0, x1:1, y0:v.y, y1:v.y,
      line:{color:"#7a7a7a", width:1, dash:"dash"}});
    if (v.label) ex.ann.push({x:1, y:v.y, xref:"paper", yref:"y", text:v.label, showarrow:false,
      font:{size:9, color:"#7a7a7a"}, xanchor:"right", yanchor:"bottom"});
  });
  if (L.diag) ex.shapes.push({type:"line", xref:"paper", yref:"paper", x0:0, y0:0, x1:1, y1:1,
    line:{color:"#9a9a9a", width:1, dash:"dash"}});
  if (L.hline !== undefined) ex.shapes.push({type:"line", xref:"paper", yref:"y", x0:0, x1:1,
    y0:L.hline, y1:L.hline, line:{color:"#7a7a7a", width:1, dash:"dash"}});
  if (L.vline !== undefined) ex.shapes.push({type:"line", xref:"x", yref:"paper", x0:L.vline,
    x1:L.vline, y0:0, y1:1, line:{color:"#7a7a7a", width:1, dash:"dash"}});
}
function pts(L, r, t, ex){
  shapes(L, ex);
  const X = col(L.x,t), Y = col(L.y,t), SZ = col(L.size,t), C = col(L.color,t);
  const out = [];
  let lo = Infinity, hi = -Infinity;
  if (SZ) for (const i of r){ const v = SZ.a[i]; if (v>0 && Number.isFinite(v)){ lo=Math.min(lo,v); hi=Math.max(hi,v); } }
  for (const [key, ix0] of groups(L.split, r, t)){
    const ix = thin(ix0);
    const x = new Array(ix.length), y = new Array(ix.length);
    for (let j=0;j<ix.length;j++){ const i = ix[j];
      x[j] = X ? X.a[i] : null; y[j] = Y ? Y.a[i] : (L.lane ? lab(L.lane,i,t) : null); }
    const mk = {size:5.5, opacity:L.opacity===undefined?0.8:L.opacity, line:{width:0}};
    if (C && !C.cats){
      mk.color = Array.from(ix, i => C.a[i]); mk.colorscale = L.scale || "Viridis";
      mk.showscale = true;
      mk.colorbar = {thickness:9, len:0.75, outlinewidth:0, ticks:"outside",
        tickfont:{size:9.5, color:"#4D4D4D"}, title:{text:L.color, side:"right", font:{size:9.5}}};
    } else if (C && C.cats){
      // a categorical colour column varies point by point (e.g. closer/farther on a lane)
      mk.color = Array.from(ix, i => { const k = C.a[i];
        return (k >= 0 && k < C.cats.length) ? colour(L, C.cats[k], t) : "#B3B3B3"; });
    } else mk.color = colour(L, key, t);
    if (SZ && hi > lo){
      const root = Math.sqrt(hi);
      mk.size = Array.from(ix, i => { const v = SZ.a[i];
        return (v > 0 && Number.isFinite(v)) ? 2.5 + 11 * Math.sqrt(Math.min(v, hi)) / root : 2.5; });
    }
    out.push({type:"scatter", mode:"markers", x:x, y:y, hoverinfo:"skip",
      name: key || L.name || "points", marker: mk, legendgroup: String(key),
      showlegend: !!(L.split && L.split !== "none")});
  }
  return out;
}
function smooth(x, y, w){
  // running median: robust to the spikes a mean would smear across the profile
  if (!w || w < 2) return y;
  const n = y.length, out = new Float64Array(n), half = Math.floor(w / 2), buf = [];
  for (let i=0;i<n;i++){
    if (!Number.isFinite(y[i])){ out[i] = NaN; continue; }
    buf.length = 0;
    for (let k = Math.max(0, i-half); k <= Math.min(n-1, i+half); k++)
      if (Number.isFinite(y[k])) buf.push(y[k]);
    buf.sort((a,b) => a-b);
    out[i] = buf.length ? buf[buf.length >> 1] : NaN;
  }
  return out;
}
function lines(L, r, t, ex){
  shapes(L, ex);
  const X = col(L.x,t), Y = col(L.y,t), Y2 = col(L.y2,t);
  if (!X || !Y) return [];
  const out = [];
  for (const [key, ix0] of groups(L.split, r, t)){
    const s = Array.from(ix0).filter(i => Number.isFinite(X.a[i]) && Number.isFinite(Y.a[i]))
      .sort((a,b) => X.a[a]-X.a[b]);
    if (!s.length) continue;
    const x = s.map(i => X.a[i]), y = s.map(i => Y.a[i]), c = colour(L, key, t);
    out.push({type:"scatter", mode:"lines", x:x, y:Array.from(smooth(x, y, L.smooth)),
      line:{color:c, width:1.4}, hoverinfo:"skip", name: key || L.name || "value",
      legendgroup:String(key), showlegend: !!(L.split && L.split !== "none")});
  }
  return out;
}
function ecdf(L, r, t, ex){
  const X = col(L.x,t); if (!X) return [];
  const out = [];
  for (const [key, ix] of groups(L.split, r, t)){
    const v = []; for (const i of ix) if (Number.isFinite(X.a[i])) v.push(X.a[i]);
    if (!v.length) continue;
    v.sort((a,b) => a-b);
    const step = Math.max(1, Math.floor(v.length/1200)), x = [], y = [];
    for (let k=0;k<v.length;k+=step){ x.push(v[k]); y.push((k+1)/v.length); }
    x.push(v[v.length-1]); y.push(1);
    const c = key ? colour(L, key, t) : PAL[0];
    out.push({type:"scatter", mode:"lines", x:x, y:y, name: key || "ECDF", hoverinfo:"skip",
      line:{color:c, width:1.4}, legendgroup:String(key),
      showlegend: !!(L.split && L.split !== "none")});
  }
  return out;
}
function binsOf(v, n, log){
  let lo = Infinity, hi = -Infinity;
  for (const x of v){ if (!Number.isFinite(x) || (log && x <= 0)) continue;
    if (x<lo) lo=x; if (x>hi) hi=x; }
  if (!Number.isFinite(lo)) return null;
  if (hi <= lo) hi = lo*1.0001 + 1e-12;
  if (log){ lo = Math.log10(lo); hi = Math.log10(hi); }
  const w = (hi-lo)/n;
  return {lo:lo, w:w, n:n, log:log,
    edge: i => log ? Math.pow(10, lo+i*w) : lo+i*w,
    idx: x => { if (!Number.isFinite(x) || (log && x<=0)) return -1;
      const k = Math.floor(((log ? Math.log10(x) : x) - lo)/w); return k<0||k>=n ? -1 : k; }};
}
function heat(L, r, t, ex){
  const X = col(L.x,t), Y = col(L.y || L.lane,t), Z = col(L.z,t);
  if (!X || !Y) return [];
  const xc = X.cats, yc = Y.cats;
  let xb = null, yb = null, nx = xc ? xc.length : 0, ny = yc ? yc.length : 0;
  if (!xc){
    // fixed bin edges keep the heat map cells in place when filters change
    xb = (L.x0 !== undefined && L.x1 !== undefined)
      ? binsOf([L.x0, L.x1], L.bins||300, !!L.logx) : binsOf(Array.from(r, i => X.a[i]), L.bins||300, !!L.logx);
    if (!xb) return [];
    nx = xb.n;
  }
  if (!yc){ yb = binsOf(Array.from(r, i => Y.a[i]), L.bins||60, !!L.logy); if (!yb) return []; ny = yb.n; }
  const sum = new Float64Array(nx*ny), cnt = new Float64Array(nx*ny),
        mx = new Float64Array(nx*ny).fill(-Infinity);
  const acc = L.agg === "median" ? new Map() : null;
  for (const i of r){
    const ix = xc ? X.a[i] : xb.idx(X.a[i]), iy = yc ? Y.a[i] : yb.idx(Y.a[i]);
    if (ix<0 || iy<0 || ix>=nx || iy>=ny) continue;
    const k = iy*nx + ix;
    if (L.agg === "count"){ cnt[k]++; continue; }
    const v = Z ? Z.a[i] : 1;
    if (!Number.isFinite(v)) continue;
    cnt[k]++; sum[k] += v;
    if (v > mx[k]) mx[k] = v;
    if (acc){ if(!acc.has(k)) acc.set(k, []); acc.get(k).push(v); }
  }
  const z = [];
  for (let j=0;j<ny;j++){ const row = new Array(nx);
    for (let i=0;i<nx;i++){ const k = j*nx+i; let v = NaN;
      if (cnt[k]){
        if (L.agg === "count") v = cnt[k];
        else if (L.agg === "mean") v = sum[k]/cnt[k];
        else if (L.agg === "sum") v = sum[k];
        else if (L.agg === "max") v = mx[k];
        else { const a = acc.get(k) || []; a.sort((p,q)=>p-q); v = a.length ? a[a.length>>1] : NaN; }
      }
      row[i] = v; }
    z.push(row); }
  if (L.sym) for (let j=0;j<ny;j++) for (let i=0;i<j;i++){
    const up = z[i][j], lo = z[j][i] || null;
    if (Number.isFinite(up) && !Number.isFinite(z[j][i])) z[j][i] = up;
    else if (Number.isFinite(lo) && !Number.isFinite(z[i][j])) z[i][j] = lo;
  }
  const xv = xc ? xc.slice() : Array.from({length:nx}, (_,i) => (xb.edge(i)+xb.edge(i+1))/2);
  const yv = yc ? yc.slice() : Array.from({length:ny}, (_,i) => (yb.edge(i)+yb.edge(i+1))/2);
  const tr = {type:"heatmap", z:z, x:xv, y:yv, colorscale:L.scale||"Viridis", hoverinfo:"skip",
    xgap: (xc && yc && L.tiles) ? 2 : 0, ygap: (xc && yc && L.tiles) ? 2 : 0,
    colorbar:{thickness:9, len:0.75, outlinewidth:0, ticks:"outside",
      tickfont:{size:9.5, color:"#4D4D4D"}, title:{text:L.z||"", side:"right", font:{size:9.5}}}};
  if (L.zmin !== undefined) tr.zmin = L.zmin;
  if (L.zmax !== undefined) tr.zmax = L.zmax;
  if (L.text){ tr.text = z.map(row => row.map(v => Number.isFinite(v) ? fmt(v) : ""));
    tr.texttemplate = "%{text}"; tr.textfont = {size:10, color:"#333"}; }
  return [tr];
}
function bars(L, r, t, ex){
  const X = col(L.x,t), Y = col(L.y,t);
  const xs = X.cats || [];
  const keys = L.cols ? L.cols.slice() : [...groups(L.split, r, t).keys()];
  const g = L.cols ? null : groups(L.split, r, t);
  const val = new Map(keys.map(k => [k, new Float64Array(xs.length)]));
  const num = new Map(keys.map(k => [k, new Float64Array(xs.length)]));
  for (const k of keys){
    const ix = g ? g.get(k) : r;
    for (const i of (ix || r)){
      const xi = X.a[i]; if (!(xi >= 0 && xi < xs.length)) continue;
      if (L.cols){ const v = col(k,t).a[i]; if (Number.isFinite(v)) val.get(k)[xi] += v; continue; }
      if (L.agg === "count"){ val.get(k)[xi]++; continue; }
      const v = Y ? Y.a[i] : NaN; if (!Number.isFinite(v)) continue;
      val.get(k)[xi] += v; num.get(k)[xi]++;
    }
  }
  if (L.frac) for (let i=0;i<xs.length;i++){ let tot = 0;
    for (const k of keys) tot += val.get(k)[i];
    if (tot) for (const k of keys) val.get(k)[i] /= tot; }
  const label = k => (L.labels && L.labels[k]) || k || L.name || "value";
  return keys.map((k, ki) => {
    const v = val.get(k);
    if (L.agg === "mean") for (let i=0;i<v.length;i++) if (num.get(k)[i]) v[i] /= num.get(k)[i];
    return {type:"bar", x:xs, y:Array.from(v), name:label(k), hoverinfo:"skip",
      marker:{color: (g && L.split) ? colour(L, k, t) : (keys.length > 1 ? PAL[ki % PAL.length] : "#595959"),
        line:{width:0}},
      legendgroup:String(k), showlegend:keys.length>1};
  });
}
function melt(L, r, t, ex){
  const lane = col(L.lane,t), out = [];
  L.cols.forEach((name, ki) => {
    const C = col(name,t); if (!C) return;
    const x = [], y = [];
    for (const i of r){ if (!Number.isFinite(C.a[i])) continue;
      x.push(C.a[i]); y.push(lane ? lab(L.lane, i, t) : ""); }
    out.push({type:"scatter", mode:"markers", x:x, y:y, name:name, hoverinfo:"skip",
      marker:{size:6, color:PAL[ki % PAL.length], symbol:"diamond", line:{width:0}},
      legendgroup:name, showlegend:true});
  });
  return out;
}
const BUILD = {points:pts, lines:lines, heatmap:heat, bars:bars, ecdf:ecdf, melt:melt};

/* ---------------------------------------------------------------- panel draw */
function ctl(panel, host){
  const st = ST.c[panel.id] = ST.c[panel.id] || {};
  for (const c of panel.controls || []){
    const w = document.createElement("span");
    const l = document.createElement("label"); l.textContent = c.label; w.appendChild(l);
    const v = (c.label in st) ? st[c.label] : c.default;
    if (c.opts){
      const s = document.createElement("select");
      for (const o of c.opts){ const op = document.createElement("option");
        op.value = o; op.textContent = o; op.selected = o === v; s.appendChild(op); }
      s.onchange = () => { st[c.label] = s.value; draw(panel); };
      w.appendChild(s);
    } else if (typeof v === "number"){
      const i = document.createElement("input"); i.type = "number"; i.value = v; i.min = 1;
      i.style.width = "52px"; i.style.font = "inherit"; i.style.fontSize = "12px";
      i.onchange = () => { st[c.label] = +i.value; draw(panel); };
      w.appendChild(i);
    } else {
      const i = document.createElement("input"); i.type = "checkbox"; i.checked = !!v;
      i.onchange = () => { st[c.label] = i.checked; draw(panel); };
      w.appendChild(i);
    }
    host.appendChild(w);
  }
}
function patch(panel, L, li){
  const st = ST.c[panel.id] || {};
  for (const c of panel.controls || []){
    if (c.layer !== li) continue;
    const v = (c.label in st) ? st[c.label] : c.default;
    if (c.target === "smooth") L.smooth = v;
    else if (c.target === "frac") L.frac = v;
    else if (c.target === "show") L.show = v;
    else L[c.target] = (v === "none" ? null : v);
  }
  return L;
}
// white panel with light grid, no axis lines: theme_minimal rather than theme_grey
function axis(o, log, cats, range){
  const a = {title:{text:o.t, font:{size:11.5, color:"#333"}}, showline:false, zeroline:false,
    ticks:"outside", ticklen:3, tickwidth:1, tickcolor:"#B3B3B3",
    tickfont:{size:10.5, color:"#4D4D4D"}, gridcolor:"#E8E8E8", gridwidth:1,
    fixedrange:true, automargin:true, autorange:false,
    exponentformat:"power", showexponent:"all", separatethousands:false,
    type: log ? "log" : (cats ? "category" : "linear")};
  if (cats){ a.categoryorder = "array"; a.categoryarray = cats; }
  else if (range){ a.range = range; }
  return a;
}
function xrange(x){
  if (!x.seq_extent) return x.range;
  const f = R.filters.find(f => f.type === "one");
  const seq = f && ST.f[f.id];
  const hi = (R.extents && R.extents[seq]) || x.max || 1;
  return [0, hi * 1.01];
}
function layout(panel, ex){
  const x = panel.x, y = panel.y;
  return {
    height: panel.h, margin:{l:64, r:14, t:10, b:44}, hovermode:false, dragmode:false,
    barmode:"stack", bargap:0.18, paper_bgcolor:"#FFFFFF", plot_bgcolor:"#FFFFFF",
    font:{family:GG.font, size:11.5, color:"#222"}, colorway:PAL,
    legend:{orientation:"v", x:1.02, y:1, xanchor:"left", yanchor:"top",
      bgcolor:"rgba(0,0,0,0)", borderwidth:0, font:{size:10.5}, itemsizing:"constant"},
    showlegend:false, shapes:ex.shapes, annotations:ex.ann,
    xaxis: axis(x, x.log, null, xrange(x)),
    yaxis: axis(y, y.log, y.cats, y.range)
  };
}
function draw(panel){
  const host = document.getElementById("p-" + panel.id);
  if (!host) return;
  try {
    const t = panel.table || "rows";
    const r = rowsOf(t === "rows" ? ST.mask : null, t);
    const ex = {shapes:[], ann:[]};
    let tr = [];
    panel.layers.forEach((L0, li) => {
      const L = patch(panel, JSON.parse(JSON.stringify(L0)), li);
      const fn = BUILD[L.k];
      if (fn) for (const x of fn(L, r, t, ex)) tr.push(x);
      else shapes(L, ex);
    });
    const lay = layout(panel, ex);
    lay.showlegend = tr.some(x => x.showlegend);
    Plotly.react(host, tr, lay, CFG);
    stat(panel, r, t);
  } catch (err){
    host.innerHTML = "<p style='color:#b91c1c;font-size:12.5px'>" + err.message + "</p>";
    console.error(err);
  }
}
function stat(panel, r, t){
  const host = document.getElementById("s-" + panel.id);
  if (!host) return;
  if (t !== "rows"){ host.textContent = tb(t).n + " features (annotation table, not filtered)"; return; }
  const parts = [r.length.toLocaleString() + " of " + tb("rows").n.toLocaleString() + " rows"];
  for (const spec of panel.stats || []){
    const i = spec.indexOf(":"), op = i < 0 ? spec : spec.slice(0,i), name = i < 0 ? "" : spec.slice(i+1);
    if (op === "n") continue;
    const c = col(name); if (!c) continue;
    const v = []; const step = Math.max(1, Math.floor(r.length/30000));
    for (let k=0;k<r.length;k+=step) if (Number.isFinite(c.a[r[k]])) v.push(c.a[r[k]]);
    if (!v.length) continue;
    let out;
    if (op === "sum") out = fmt(v.reduce((a,b)=>a+b,0));
    else if (op === "mean") out = fmt(v.reduce((a,b)=>a+b,0)/v.length);
    else { v.sort((a,b)=>a-b);
      out = op === "min" ? fmt(v[0]) : op === "max" ? fmt(v[v.length-1]) : fmt(v[v.length>>1]); }
    parts.push(op + "(" + name + ") " + out);
  }
  host.textContent = parts.join("   ·   ");
}
function redraw(){
  ST.mask = mask();
  const n = rowsOf(ST.mask).length;
  const btn = document.getElementById("count");
  if (btn) btn.textContent = n.toLocaleString();
  for (const panel of R.panels) draw(panel);
}
function schedule(){
  if (ST.timer) clearTimeout(ST.timer);
  ST.timer = setTimeout(redraw, 120);
}

/* ------------------------------------------------------------------- boot */
function boot(){
  document.getElementById("title").textContent = R.title;
  document.getElementById("sub").textContent = R.meta;
  document.getElementById("cmd").textContent =
    (R.chips || []).map(c => c[0] + " " + c[1]).join("  ·  ") + (R.cmd ? "  ·  " + R.cmd : "");
  const sum = document.getElementById("summary");
  for (const c of R.cards){
    const sp = document.createElement("span");
    sp.innerHTML = "<b>" + c.title + "</b> " + c.value;
    sum.appendChild(sp);
    sum.appendChild(document.createTextNode(c.hint ? " (" + c.hint + ")  ·  " : "  ·  "));
  }
  const main = document.getElementById("panels");
  for (const panel of R.panels){
    const fig = document.createElement("figure");
    const cap = document.createElement("figcaption");
    const b = document.createElement("b"); b.textContent = panel.title;
    cap.appendChild(b); cap.appendChild(document.createTextNode(panel.desc));
    fig.appendChild(cap);
    if ((panel.controls || []).length){
      const c = document.createElement("div"); c.className = "pctl";
      ctl(panel, c); fig.appendChild(c);
    }
    const plot = document.createElement("div"); plot.className = "plot"; plot.id = "p-" + panel.id;
    plot.style.width = "100%";
    if (panel.width) plot.style.maxWidth = panel.width + "px";
    fig.appendChild(plot);
    const st = document.createElement("div"); st.className = "pstat"; st.id = "s-" + panel.id;
    fig.appendChild(st);
    main.appendChild(fig);
  }
  const notes = document.getElementById("notes");
  if (R.notes.length) notes.innerHTML = "<ul>" + R.notes.map(n => "<li>" + n + "</li>").join("") + "</ul>";
  document.getElementById("reset").onclick = () => { ST.f = {}; ST.c = {}; filtUI(); redraw(); };
  filtUI();
  redraw();
}
window.gviz = {draw, redraw, mask, rowsOf, tb, col, ST, R, groups, BUILD, layout};
document.addEventListener("DOMContentLoaded", boot);
"""

HTML = """<!DOCTYPE html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>__TITLE__</title>
__PLOTLY__<style>__CSS__</style></head><body>
<h1 id="title"></h1><p class="sub" id="sub"></p><p class="cmd" id="cmd"></p>
<div id="filters"></div>
<p id="summary"></p>
<div id="panels"></div>
<footer><span>rows shown <b id="count">-</b> · static panels, change a filter to redraw ·
<button id="reset">reset filters</button></span><div id="notes"></div></footer>
<script>__JS__</script></body></html>"""


def write_html(report, path, inline):
    payload = json.dumps(report, separators=(",", ":"), default=str).replace("</", "<\\/")
    if inline:
        from plotly.offline import get_plotlyjs
        lib = "<script>" + get_plotlyjs() + "</script>"
    else:
        lib = '<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>'
    html = (HTML.replace("__TITLE__", report["title"]).replace("__PLOTLY__", lib)
            .replace("__CSS__", CSS).replace("__JS__", JS.replace("__REPORT__", payload)))
    Path(path).write_text(html, encoding="utf-8")
    return Path(path).stat().st_size


def main(argv=None):
    ap = argparse.ArgumentParser(description="interactive HTML report for one gdiff output")
    ap.add_argument("input", type=Path)
    ap.add_argument("-o", "--output", type=Path)
    ap.add_argument("--kind", choices=["map", "roll", "dist", "samples"])
    ap.add_argument("--seq", action="append", default=[],
                    help="only these sequences (name or prefix, repeatable)")
    ap.add_argument("--assembly-report", type=Path)
    ap.add_argument("--no-report", action="store_true", help="do not look for an assembly report")
    ap.add_argument("--annotations", type=Path, help="GFF3 file for a gene track")
    ap.add_argument("--title")
    ap.add_argument("--bins", type=int, default=1000,
                    help="roll bins per sequence (~one per screen column) [1000]")
    ap.add_argument("--keep-all", action="store_true",
                    help="keep every sequence, including unplaced scaffolds")
    ap.add_argument("--inline", action="store_true", help="embed plotly.js (offline, +4.8 MB)")
    ap.add_argument("--open", action="store_true")
    ap.add_argument("--print-schema", action="store_true")
    args = ap.parse_args(argv)

    if args.print_schema:
        df, kind, run = read_table(args.input, args.kind)
        print(f"kind:   {kind}\nrows:   {len(df):,}\ncols:   {' '.join(df.columns)}\ncmd:    {run['cmd']}")
        return 0
    opts = Args(**vars(args))
    report = build(args.input, args.kind, opts)
    out = args.output or args.input.with_suffix(args.input.suffix + ".gviz.html")
    size = write_html(report, out, args.inline)
    print(f"{out}  {size/1e6:.1f} MB  {len(report['panels'])} panels  "
          f"{report['tables']['rows']['n']:,} points"
          + ("" if args.inline else "  (plotly.js from the CDN; --inline to embed it)"))
    if args.open:
        webbrowser.open(out.resolve().as_uri())
    return 0


if __name__ == "__main__":
    sys.exit(main())
