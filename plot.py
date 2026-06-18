import argparse
import base64
import importlib.metadata
import traceback
from collections import defaultdict
from functools import lru_cache
from typing import Optional, Union

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from dash import Dash, ctx, dcc, html, Input, Output, State, no_update
from ete3 import Tree
from plotly.subplots import make_subplots

from plotrc import *  # noqa: F403

# =============================================================================
# GDIFF TSV SCHEMA
# =============================================================================

GDIFF_COLUMNS = [
    "QUERY_ID",
    "SEQ_LEN",
    "INTERVAL_START",
    "INTERVAL_END",
    "STRAND",
    "IS_RC",
    "REF_ID",
    "DIST",
    "MASK",
    "D_INTERVAL",
    "DIST_CONTIG",
    "STRAND_DIFF",
    "DIST_GENOME",
    "PERCENTILE",
    "FOLD",
    "QVALUE",
]

PLOT_DIST_COL = "_PLOT_DIST"

# STRAND: '+' = closer (lower-distance) strand, '-' = farther, '.' = unknown.
STRAND_LABELS = {"+": "closer", "-": "farther", ".": "unknown"}


def clamp_distance(val) -> float:
    """Clamp finite distances to the valid MLE open interval [D_EPS, D_UB)."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return float("nan")
    try:
        d = float(val)
    except (TypeError, ValueError):
        return float("nan")
    if not np.isfinite(d) or d < D_EPS or d >= D_UB:
        return float("nan")
    return d


def format_strand_diff(val) -> Optional[str]:
    """Format STRAND_DIFF for hover text (finite diff or sentinel labels)."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return None
    try:
        d = float(val)
    except (TypeError, ValueError):
        return None
    if np.isfinite(d):
        return f"{d:.4g}"
    if d == float("-inf"):
        return "fw only"
    if d == float("inf"):
        return "rc only"
    return None


def has_informative_strand(df: pd.DataFrame) -> bool:
    """True when STRAND has at least one '+' or '-' row (not only '.')."""
    if "STRAND" not in df.columns:
        return False
    chars = df["STRAND"].astype(str).str.strip().str[:1]
    return chars.isin(["+", "-"]).any()


def parse_is_rc(val) -> float:
    """Normalize IS_RC to 0 (forward) or 1 (reverse-complement)."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return float("nan")
    s = str(val).strip().lower()
    if s in ("1", "true", "rc", "yes"):
        return 1.0
    if s in ("0", "false", "fw", "no"):
        return 0.0
    try:
        return 1.0 if int(float(s)) else 0.0
    except (ValueError, TypeError):
        return float("nan")


def is_rc_label(val) -> str:
    """Return 'fw' or 'rc' for hover text."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return "?"
    return "rc" if int(val) else "fw"


def has_is_rc_column(df: pd.DataFrame) -> bool:
    """True when IS_RC is present with at least one parsed value."""
    if "IS_RC" not in df.columns:
        return False
    parsed = df["IS_RC"].map(parse_is_rc)
    return parsed.notna().any()


def is_enum_lite_dataframe(df: pd.DataFrame) -> bool:
    """Enum lite: 15-column output with no segment MLE or significance testing."""
    if "DIST" not in df.columns or "MASK" not in df.columns:
        return False
    dist_ok = df["DIST"].notna() & np.isfinite(df["DIST"])
    if dist_ok.any():
        return False
    if "PERCENTILE" in df.columns:
        p_ok = df["PERCENTILE"].notna() & np.isfinite(df["PERCENTILE"])
        if p_ok.any():
            return False
    return True


def parse_d_interval_upper(val) -> float:
    """Return the upper bound from a D_INTERVAL string like '(0, 0.1)'."""
    if not isinstance(val, str) or not val.startswith("("):
        return float("nan")
    try:
        parts = val.strip("()").split(",")
        return float(parts[1].strip())
    except (ValueError, IndexError):
        return float("nan")


def normalize_gdiff_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Add optional columns and a unified distance column for plotting."""
    df = df.copy()
    for col in ("STRAND_DIFF", "QVALUE", "DIST_CONTIG", "DIST_GENOME", "FOLD"):
        if col not in df.columns:
            df[col] = np.nan
    if "STRAND" in df.columns:
        df["STRAND"] = df["STRAND"].astype(str).str.strip().str[:1]
        df.loc[~df["STRAND"].isin(["+", "-", "."]), "STRAND"] = "."
    if "IS_RC" in df.columns:
        df["IS_RC"] = df["IS_RC"].map(parse_is_rc)
    if "DIST" in df.columns:
        df["DIST"] = pd.to_numeric(df["DIST"], errors="coerce").map(clamp_distance)
    if "DIST_CONTIG" in df.columns:
        df["DIST_CONTIG"] = pd.to_numeric(df["DIST_CONTIG"], errors="coerce").map(
            clamp_distance
        )
    if "DIST_GENOME" in df.columns:
        df["DIST_GENOME"] = pd.to_numeric(df["DIST_GENOME"], errors="coerce").map(
            clamp_distance
        )
    if "DIST_TH" in df.columns:
        df[PLOT_DIST_COL] = pd.to_numeric(df["DIST_TH"], errors="coerce")
    else:
        df[PLOT_DIST_COL] = pd.to_numeric(
            df["DIST"] if "DIST" in df.columns else np.nan, errors="coerce"
        )
        if "D_INTERVAL" in df.columns:
            missing = df[PLOT_DIST_COL].isna()
            if missing.any():
                df.loc[missing, PLOT_DIST_COL] = df.loc[missing, "D_INTERVAL"].map(
                    parse_d_interval_upper
                )
    if PLOT_DIST_COL in df.columns:
        df[PLOT_DIST_COL] = pd.to_numeric(df[PLOT_DIST_COL], errors="coerce")
    return df


def detect_plot_mode(
    df: pd.DataFrame,
    enum_only: Optional[bool] = None,
) -> str:
    """Return 'legacy_enum', 'enum', or 'continuous'.

    legacy_enum: 7-column DIST_TH ground truth.
    enum: 15-column enum lite (no segment MLE / no significance).
    continuous: segment MLE output (continuous mode or enum + significance test).
    """
    if enum_only is True:
        return "legacy_enum" if "DIST_TH" in df.columns else "enum"
    if "DIST_TH" in df.columns and "PERCENTILE" not in df.columns:
        return "legacy_enum"
    if "DIST" in df.columns and "MASK" in df.columns:
        if is_enum_lite_dataframe(df):
            return "enum"
        return "continuous"
    raise ValueError(
        "Unrecognized TSV format. Expected legacy enum (DIST_TH), "
        "15-column enum lite (DIST/MASK, NaN MLE), or continuous/enum+test "
        "(finite DIST and/or PERCENTILE)."
    )


def is_enum_plot_mode(mode: str) -> bool:
    return mode in ("legacy_enum", "enum")


def load_gdiff_tsv(path: str) -> pd.DataFrame:
    return normalize_gdiff_dataframe(pd.read_csv(path, sep="\t"))


# =============================================================================
# DATA LOADING & FILTERING
# =============================================================================


def load_tree(path: str) -> Tree:
    """Load Newick tree and ladderize.

    Args:
        path: Path to Newick format tree file.

    Returns:
        Ladderized ete3 Tree object.
    """
    tree = Tree(path, format=1)
    tree.ladderize()
    return tree


def prune_tree(tree: Tree, retained: set[str]) -> Tree:
    """Return pruned tree keeping only retained leaves.

    Args:
        tree: Source tree to prune.
        retained: Set of leaf names to keep.

    Returns:
        Pruned copy of the tree.
    """
    if not retained:
        return tree
    t = tree.copy()
    for leaf in t.iter_leaves():
        if leaf.name not in retained:
            leaf.delete()
    return t


def get_tip_order(tree: Tree) -> list[str]:
    """Return leaf names in tree order (top to bottom).

    Args:
        tree: Input tree.

    Returns:
        List of leaf names in traversal order.
    """
    return [leaf.name for leaf in tree.iter_leaves()]


def get_retained_leaves(df: pd.DataFrame) -> set[str]:
    """Return REF_IDs that should be retained (have non-trivial intervals).

    Args:
        df: Input dataframe with REF_ID, INTERVAL_START, INTERVAL_END columns.

    Returns:
        Set of REF_IDs with at least one non-empty interval.
    """
    return set(df[df["INTERVAL_START"] != df["INTERVAL_END"]]["REF_ID"].unique())


def get_query_retained_leaves(df: pd.DataFrame, query: str) -> set[str]:
    """Return REF_IDs with non-trivial intervals for a specific query.

    Args:
        df: Input dataframe.
        query: Query identifier to filter on.

    Returns:
        Set of REF_IDs matching the query with non-empty intervals.
    """
    df_q = df[df["QUERY_ID"] == query]
    return set(df_q[df_q["INTERVAL_START"] != df_q["INTERVAL_END"]]["REF_ID"].unique())


def compute_distance_range(df: pd.DataFrame) -> tuple[float, float]:
    """Adaptive [min, max] distance range clamped to valid MLE bounds."""
    if df.empty or "DIST" not in df.columns:
        return D_EPS, min(0.5, D_UB - D_EPS)
    finite = df["DIST"].replace([np.inf, -np.inf], np.nan).dropna()
    if finite.empty:
        return D_EPS, min(0.5, D_UB - D_EPS)
    d_min = max(float(finite.min()), D_EPS)
    d_max = min(float(finite.max()), D_UB - D_EPS)
    if d_max <= d_min:
        d_max = min(d_min + 0.01, D_UB - D_EPS)
    return d_min, d_max


def get_distance_thresholds(df: pd.DataFrame) -> list[float]:
    """Return sorted unique plotting distance values."""
    if PLOT_DIST_COL not in df.columns:
        return []
    vals = df[PLOT_DIST_COL].dropna().unique()
    return sorted(float(v) for v in vals)


def get_sequence_identifiers(df: pd.DataFrame) -> list[str]:
    """Return sorted unique QUERY_IDs.

    Args:
        df: Input dataframe with QUERY_ID column.

    Returns:
        Sorted list of unique query identifiers.
    """
    return sorted(df["QUERY_ID"].unique())


def get_seq_len(df: pd.DataFrame) -> Optional[int]:
    """Extract sequence length from dataframe (SEQ_LEN column).

    Args:
        df: Input dataframe with SEQ_LEN column.

    Returns:
        Sequence length if uniform, None otherwise.
    """
    if df.empty or "SEQ_LEN" not in df.columns:
        return None
    vals = df["SEQ_LEN"].dropna().unique()
    return int(vals[0]) if len(vals) == 1 else None


def add_tip_order(df: pd.DataFrame, tip_order: list[str]) -> pd.DataFrame:
    ref_to_y = {name: i for i, name in enumerate(tip_order)}
    df = df.copy()
    df["y"] = df["REF_ID"].map(ref_to_y)
    return df


def filter_intervals(
    df: pd.DataFrame,
    seq_id: str,
    tip_order: list[str],
    dist_hi: float,
    dist_lo: Optional[float] = None,
    direction: Optional[str] = None,
) -> pd.DataFrame:
    """Filter intervals by query, distance threshold, and forward/rc direction."""
    df_q = df[df["QUERY_ID"] == seq_id].copy()

    if direction and direction != "both" and "IS_RC" in df_q.columns:
        want_rc = direction == "rc"
        df_q = df_q[df_q["IS_RC"] == (1.0 if want_rc else 0.0)]

    grp_cols = ["REF_ID", "y", "QUERY_ID", "INTERVAL_START", "INTERVAL_END", "SEQ_LEN"]
    if "IS_RC" in df_q.columns:
        grp_cols.append("IS_RC")
    elif "STRAND" in df_q.columns:
        grp_cols.append("STRAND")
    df_q = df_q.groupby(grp_cols, as_index=False)[PLOT_DIST_COL].min()

    df_q = df_q[df_q[PLOT_DIST_COL] <= dist_hi]
    if dist_lo is not None:
        df_q = df_q[df_q[PLOT_DIST_COL] >= dist_lo]

    return df_q.dropna(subset=["y"]).sort_values("y")


def filter_intervals_continuous(
    df: pd.DataFrame,
    seq_id: str,
    tip_order: list[str],
    strand: Optional[str] = None,
    sig_th: Optional[float] = None,
    sig_col: str = "PERCENTILE",
    fold: Optional[str] = None,
) -> pd.DataFrame:
    """Filter intervals by query, strand, significance cutoff, and fold change.

    strand: '+', '-', or 'both'. '+' = closer strand, '-' = farther strand.
    sig_th: if set, keep rows with sig_col <= sig_th (NaN rows dropped).
    sig_col: 'PERCENTILE' (p-value) or 'QVALUE' (BH-adjusted).
    fold: '<1' keeps FOLD < 1 (closer), '>1' keeps FOLD > 1 (more distant).
    """
    df_q = df[df["QUERY_ID"] == seq_id].copy()

    if sig_th is not None and sig_col in df_q.columns:
        sig = pd.to_numeric(df_q[sig_col], errors="coerce")
        df_q = df_q[sig.notna() & np.isfinite(sig) & (sig <= sig_th)]

    if fold is not None and "FOLD" in df_q.columns:
        fold_v = pd.to_numeric(df_q["FOLD"], errors="coerce")
        if fold == "<1":
            df_q = df_q[fold_v.notna() & np.isfinite(fold_v) & (fold_v < 1.0)]
        elif fold == ">1":
            df_q = df_q[fold_v.notna() & np.isfinite(fold_v) & (fold_v > 1.0)]

    if strand and strand != "both" and "STRAND" in df_q.columns:
        df_q = df_q[df_q["STRAND"] == strand]

    return df_q.dropna(subset=["y"]).sort_values("y")


def make_cached_filter(df: pd.DataFrame, enum_only=True):
    if enum_only:

        @lru_cache(maxsize=CACHE_SIZE)
        def cached(seq_id, tip_order_tuple, dist_hi, dist_lo, direction):
            return filter_intervals(
                df, seq_id, list(tip_order_tuple), dist_hi, dist_lo, direction
            )

        return cached
    else:

        @lru_cache(maxsize=CACHE_SIZE)
        def cached_cont(
            seq_id,
            tip_order_tuple,
            strand,
            sig_log=None,
            sig_col="PERCENTILE",
            fold=None,
        ):
            sig_th = 10.0**sig_log if sig_log is not None and sig_log < 0 else None
            return filter_intervals_continuous(
                df,
                seq_id,
                list(tip_order_tuple),
                strand,
                sig_th,
                sig_col,
                fold,
            )

        return cached_cont


# =============================================================================
# TREE LAYOUT
# =============================================================================


def compute_path_distances(tree: Tree, query_name: str) -> dict:
    """Compute tree distance from query leaf to all nodes."""
    query_leaf = next(
        (leaf for leaf in tree.iter_leaves() if leaf.name == query_name), None
    )
    if query_leaf is None:
        return {}
    return {node: query_leaf.get_distance(node) for node in tree.traverse()}


def compute_tree_layout(tree: Tree, mode: str = "phylogeny", distances: dict = None):
    """Compute plotly-ready tree data, tip order, and max coordinates."""
    data = defaultdict(list)
    tip_order = get_tip_order(tree)
    y_pos = {name: i for i, name in enumerate(tip_order)}

    for leaf in tree.iter_leaves():
        leaf.add_feature("y", y_pos[leaf.name])
    for node in tree.traverse("postorder"):
        if not node.is_leaf():
            node.add_feature("y", np.mean([c.y for c in node.children]))

    if mode == "cladogram":
        max_depth = 0
        for node in tree.traverse("postorder"):
            depth = 0 if node.is_leaf() else max(c.depth for c in node.children) + 1
            node.add_feature("depth", depth)
            max_depth = max(max_depth, depth)
        for node in tree.traverse("preorder"):
            node.add_feature("x", float(max_depth - node.depth))
        max_x = float(max_depth)
        count_offset = 1
    else:
        max_x = 0.0
        for node in tree.traverse("preorder"):
            x = 0.0 if node.is_root() else node.up.x + (node.dist or 0)
            node.add_feature("x", x)
            max_x = max(max_x, x)
        count_offset = 0

    for node in tree.traverse():
        if node.is_root():
            continue
        hover_text = _node_hover(node, distances, count_offset)
        data["x"].extend([node.up.x, node.x, None])
        data["y"].extend([node.y, node.y, None])
        data["text"].extend([hover_text, hover_text, None])

    for node in tree.traverse():
        if len(node.children) < 2:
            continue
        child_ys = [c.y for c in node.children]
        imin, imax = int(np.argmin(child_ys)), int(np.argmax(child_ys))
        data["x"].extend([node.x, node.x, None])
        data["y"].extend([node.children[imin].y, node.children[imax].y, None])
        data["text"].extend(
            [
                _node_hover(node.children[imin], distances, count_offset),
                _node_hover(node.children[imax], distances, count_offset),
                None,
            ]
        )
    return data, tip_order, (max_x, len(tip_order) - 1)


def _node_hover(node, distances, count_offset=0) -> str:
    """Build hover text for a tree node."""
    branch_len = node.dist or 0.0
    if node.is_leaf():
        return f"{node.name}<br>Branch length: {branch_len:.4f}"

    name = node.name or "(internal)"
    size = len(node.get_leaves()) + count_offset
    text = f"{name}<br>Subtree size: {size}<br>Branch length: {branch_len:.4f}"
    if distances and node in distances:
        text += f"<br>Distance to query: {distances[node]:.4f}"
    return text


# =============================================================================
# ZOOM & RANGE ENFORCEMENT
# =============================================================================


def enforce_min_span(
    v_range: Optional[Union[list, tuple]],
    min_span: float,
    bounds: Optional[tuple] = None,
) -> Optional[list]:
    """Clamp zoom range to enforce minimum span and respect bounds.

    Ensures that:
    1. The span is at least min_span
    2. The range stays within the specified bounds
    3. Extreme zoom values are gracefully handled
    """
    if v_range is None:
        return None

    lo, hi = sorted(v_range)
    span = hi - lo

    # Handle edge cases
    if span <= 0 or np.isnan(span) or np.isinf(span):
        if bounds:
            return [int(bounds[0]), int(bounds[1])]
        return None

    # Enforce minimum span
    if span < min_span:
        center = (lo + hi) / 2
        lo = center - min_span / 2
        hi = center + min_span / 2

    # Enforce bounds
    if bounds is not None:
        b_lo, b_hi = bounds
        b_span = b_hi - b_lo

        # If zoomed range is larger than bounds, reset to bounds
        if hi - lo > b_span:
            lo, hi = float(b_lo), float(b_hi)
        else:
            # Clamp to bounds while respecting min_span
            if lo < b_lo:
                lo = float(b_lo)
                hi = min(float(b_hi), lo + (hi - sorted(v_range)[0]))
            if hi > b_hi:
                hi = float(b_hi)
                lo = max(float(b_lo), hi - (sorted(v_range)[1] - sorted(v_range)[0]))

            # Final enforcement of minimum span within bounds
            if hi - lo < min_span:
                mid = (lo + hi) / 2
                lo = max(float(b_lo), mid - min_span / 2)
                hi = min(float(b_hi), mid + min_span / 2)

    # Ensure values are valid floats
    lo = float(np.clip(lo, -1e10, 1e10))
    hi = float(np.clip(hi, -1e10, 1e10))

    return [int(np.floor(lo)), int(np.ceil(hi))]


def visible_rows(y_range, total):
    return max(1, abs(int(y_range[1] - y_range[0]))) if y_range else total


def scaled_interval_width(y_range, total):
    visible = visible_rows(y_range, total)
    thickness = (FIG_HEIGHT / visible) * INTERVAL_ROW_FILL
    thickness = min(thickness, (FIG_HEIGHT / visible) * INTERVAL_LINE_MAX_RATIO)
    return min(max(INTERVAL_LINE_MIN, thickness), INTERVAL_LINE_MAX)


def scaled_tree_width(y_range, total):
    visible = visible_rows(y_range, total)
    width = (FIG_HEIGHT / visible) * TREE_ROW_FILL
    return min(max(TREE_LINE_MIN, width), TREE_LINE_MAX)


# =============================================================================
# COLOR MAPPING
# =============================================================================


@lru_cache(maxsize=CACHE_SIZE)
def get_binned_colors(thresholds: tuple, scheme: str):
    """Return bin edges and sampled colors for the colorscale."""
    n_bins = len(thresholds)
    frac = [i / max(n_bins - 1, 1) for i in range(n_bins)]
    colors = tuple(px.colors.sample_colorscale(scheme, frac))
    return (0.0,) + thresholds, colors


def assign_color_indices(distances, bin_edges, n_colors):
    indices = np.searchsorted(bin_edges, distances, side="left") - 1
    return np.clip(indices, 0, n_colors - 1)


def _build_interval_traces(starts, ends, ys, hovers, bin_idx, bin_colors,
                           descending=False):
    """Build trace dicts keyed by color from pre-computed arrays.

    Each interval is rendered as 4 points (start, mid, end, None) so the
    midpoint acts as a reliable hover target even for short intervals.
    Plotly renders later traces on top.  Use descending=True when lower bin
    indices are the most significant (enum mode).
    """
    mids = (starts + ends) / 2.0
    hovers_arr = np.asarray(hovers, dtype=object)
    traces = {}
    order = np.unique(bin_idx)
    if descending:
        order = order[::-1]
    for ci in order:
        sel = bin_idx == ci
        n = sel.sum()
        x = np.empty(n * 4, dtype=object)
        x[0::4], x[1::4], x[2::4], x[3::4] = (
            starts[sel], mids[sel], ends[sel], None,
        )
        y = np.empty(n * 4, dtype=object)
        y[0::4], y[1::4], y[2::4], y[3::4] = (
            ys[sel], ys[sel], ys[sel], None,
        )
        sh = hovers_arr[sel]
        text = np.empty(n * 4, dtype=object)
        text[0::4], text[1::4], text[2::4], text[3::4] = sh, sh, sh, None
        traces[bin_colors[ci]] = {
            "x": x.tolist(), "y": y.tolist(), "text": text.tolist(),
        }
    return traces


def batch_by_color(df, bin_edges, colors, leaf_distances=None, flip=False):
    """Group intervals by color bin for efficient trace creation (enum mode)."""
    if df.empty:
        return {}
    mask = df["INTERVAL_START"].values != df["INTERVAL_END"].values
    if not mask.any():
        return {}

    if flip:
        colors = list(reversed(colors))

    starts = df["INTERVAL_START"].values[mask]
    ends = df["INTERVAL_END"].values[mask]
    ys = df["y"].values[mask]
    dists = df[PLOT_DIST_COL].values[mask]
    refs = df["REF_ID"].values[mask]
    cidx = assign_color_indices(dists, bin_edges, len(colors))

    hovers = []
    masks = df["MASK"].values[mask] if "MASK" in df.columns else None
    strands = df["STRAND"].values[mask] if "STRAND" in df.columns else None
    is_rcs = df["IS_RC"].values[mask] if "IS_RC" in df.columns else None
    d_intervals = (
        df["D_INTERVAL"].values[mask] if "D_INTERVAL" in df.columns else None
    )
    for i, (rid, s, e, d) in enumerate(zip(refs, starts, ends, dists)):
        strand_str = ""
        if strands is not None:
            st = strands[i]
            if st in STRAND_LABELS:
                strand_str = f" ({st} {STRAND_LABELS[st]})"
        h = f"{rid}{strand_str}<br>Pos: {s:,}-{e:,}<br>Dist: {d:.4f}"
        if is_rcs is not None:
            h += f"  dir: {is_rc_label(is_rcs[i])}"
        if masks is not None:
            h += f"  mask: {int(masks[i])}"
        if d_intervals is not None:
            h += f"  D-int: {d_intervals[i]}"
        if leaf_distances and rid in leaf_distances:
            h += f"<br>Tree dist: {leaf_distances[rid]:.4f}"
        hovers.append(h)

    return _build_interval_traces(starts, ends, ys, hovers, cidx, colors,
                                   descending=True)


# Descending edges so lower-index bins = less significant (higher p-values).
# _pval_bin_index negates both edges and values for ascending searchsorted.
PVAL_BINS = (1.0, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001, 0.0)


def _pval_bin_colors(scheme, flip=False):
    """Build colors for p-value bins: first bin whitish, rest from colorscale.

    The gray bin (index 0, for p > 0.1 and NaN) is never affected by flip.
    Only the significant-color bins (indices 1..n-1) are reversed when flip=True.
    """
    n = len(PVAL_BINS) - 1  # number of bins
    fracs = [i / max(n - 2, 1) for i in range(n - 1)]
    sampled = px.colors.sample_colorscale(scheme, fracs)
    if flip:
        sampled = list(reversed(sampled))
    return ["rgba(230,230,230,0.6)"] + list(sampled)


def _pval_bin_index(pvals):
    """Assign p-values to bins. Bin 0 = [PVAL_BINS[1], PVAL_BINS[0]], etc.

    NaN and p > 0.1 are always mapped to bin 0 (the gray/whitish bin).
    """
    pvals = np.asarray(pvals, dtype=float)
    edges = np.array(PVAL_BINS)  # descending
    nan_mask = np.isnan(pvals)
    safe_pvals = np.where(nan_mask, 1.0, pvals)
    idx = np.searchsorted(-edges, -safe_pvals, side="left") - 1
    idx = np.clip(idx, 0, len(PVAL_BINS) - 2)
    idx[nan_mask] = 0
    return idx


def compute_color_values(df, color_by, dist_range):
    """Compute normalized color values and hover text for continuous coloring."""
    mask = df["INTERVAL_START"].values != df["INTERVAL_END"].values
    if not mask.any():
        return None

    starts = df["INTERVAL_START"].values[mask]
    ends = df["INTERVAL_END"].values[mask]
    ys = df["y"].values[mask]
    dists = df["DIST"].values[mask]
    pvals = df["PERCENTILE"].values[mask]
    refs = df["REF_ID"].values[mask]
    dist_contigs = (
        df["DIST_CONTIG"].values[mask] if "DIST_CONTIG" in df.columns else dists
    )
    dist_genomes = (
        df["DIST_GENOME"].values[mask] if "DIST_GENOME" in df.columns else dists
    )
    d_intervals = (
        df["D_INTERVAL"].values[mask] if "D_INTERVAL" in df.columns else None
    )
    folds = df["FOLD"].values[mask] if "FOLD" in df.columns else None
    qvalues = (
        df["QVALUE"].values[mask] if "QVALUE" in df.columns else None
    )
    strand_diffs = (
        df["STRAND_DIFF"].values[mask] if "STRAND_DIFF" in df.columns else None
    )
    strands = df["STRAND"].values[mask] if "STRAND" in df.columns else None
    is_rcs = df["IS_RC"].values[mask] if "IS_RC" in df.columns else None

    # Select color column
    if color_by in ("pval", "qval"):
        sig_src = qvalues if color_by == "qval" else pvals
        if sig_src is None:
            sig_src = pvals
        raw_colors = sig_src
        cmin, cmax = 0.0, 1.0
    else:  # "dist"
        raw_colors = dists
        cmin = max(dist_range[0], D_EPS)
        cmax = min(dist_range[1], D_UB - D_EPS)
        if cmax <= cmin:
            cmax = min(cmin + 0.01, D_UB - D_EPS)

    return {
        "starts": starts,
        "ends": ends,
        "ys": ys,
        "refs": refs,
        "dists": dists,
        "pvals": pvals,
        "dist_contigs": dist_contigs,
        "dist_genomes": dist_genomes,
        "d_intervals": d_intervals,
        "folds": folds,
        "qvalues": qvalues,
        "strand_diffs": strand_diffs,
        "strands": strands,
        "is_rcs": is_rcs,
        "raw_colors": raw_colors,
        "cmin": cmin,
        "cmax": cmax,
    }


def batch_by_continuous_color(
    df, color_by, dist_range, scheme, leaf_distances=None, flip=False
):
    """Create traces with continuous or binned coloring for non-enum mode.

    For 'pval' / 'qval' color_by: uses fixed significance bins with whitish first bin.
    For 'dist' color_by: 65-bin continuous discretization (bin 0 = NaN gray).

    Each interval is rendered as 4 points (start, mid, end, None) so the
    midpoint acts as a reliable hover target even for short intervals.
    Bins are iterated from least to most significant so that the most
    significant intervals are added last and rendered on top.
    """
    cv = compute_color_values(df, color_by, dist_range)
    if cv is None:
        return {}, 0.0, 1.0

    if color_by in ("pval", "qval"):
        sig_arr = cv["qvalues"] if color_by == "qval" else cv["pvals"]
        if sig_arr is None:
            sig_arr = cv["pvals"]
        bin_idx = _pval_bin_index(sig_arr)
        bin_colors = _pval_bin_colors(scheme, flip=flip)
    else:
        # 64 color bins + 1 gray NaN bin at index 0
        n_bins = 65
        cmin, cmax = cv["cmin"], cv["cmax"]
        span = cmax - cmin if cmax > cmin else 1.0
        raw = np.asarray(cv["raw_colors"], dtype=float)
        nan_mask = np.isnan(raw)
        fracs = np.clip((np.where(nan_mask, cmin, raw) - cmin) / span, 0.0, 1.0)
        bin_idx = np.clip((fracs * (n_bins - 2)).astype(int) + 1, 1, n_bins - 1)
        bin_idx[nan_mask] = 0
        sample_fracs = [i / max(n_bins - 2, 1) for i in range(n_bins - 1)]
        sampled = px.colors.sample_colorscale(scheme, sample_fracs)
        if flip:
            sampled = list(reversed(sampled))
        bin_colors = ["rgba(230,230,230,0.6)"] + list(sampled)

    def _fmt(val, fmt=".4f"):
        return f"{val:{fmt}}" if np.isfinite(val) else "nan"

    d_intervals = cv["d_intervals"]
    folds = cv["folds"]
    qvalues = cv["qvalues"]
    strand_diffs = cv["strand_diffs"]
    strands = cv["strands"]
    is_rcs = cv["is_rcs"]
    n_items = len(cv["refs"])
    hovers = []
    for i in range(n_items):
        rid = cv["refs"][i]
        s, e = int(cv["starts"][i]), int(cv["ends"][i])
        d, p = cv["dists"][i], cv["pvals"][i]
        dc, dg = cv["dist_contigs"][i], cv["dist_genomes"][i]
        strand_str = ""
        if strands is not None:
            st = strands[i]
            if st in STRAND_LABELS:
                strand_str = f" ({st} {STRAND_LABELS[st]})"
        h = f"{rid}{strand_str}<br>Pos: {s:,}-{e:,}<br>Dist: {_fmt(d)}"
        if is_rcs is not None:
            h += f"  dir: {is_rc_label(is_rcs[i])}"
        h += f"  p: {_fmt(p, '.3e')}"
        if qvalues is not None:
            h += f"  q: {_fmt(qvalues[i], '.3e')}"
        if folds is not None:
            h += f"  fold: {_fmt(folds[i], '.3f')}"
        h += f"<br>Contig: {_fmt(dc)}  Genome: {_fmt(dg)}"
        sd = format_strand_diff(strand_diffs[i]) if strand_diffs is not None else None
        if sd is not None:
            h += f"  d_diff: {sd}"
        if d_intervals is not None:
            h += f"  D-int: {d_intervals[i]}"
        if leaf_distances and rid in leaf_distances:
            h += f"<br>Tree dist: {leaf_distances[rid]:.4f}"
        hovers.append(h)

    traces = _build_interval_traces(
        cv["starts"], cv["ends"], cv["ys"], hovers, bin_idx, bin_colors,
    )
    return traces, cv["cmin"], cv["cmax"]


def make_continuous_colorbar(
    cmin, cmax, color_by, scheme, y_center=0.5, length=None, flip=False
):
    """Create a colorbar trace for non-enum mode."""
    if color_by in ("pval", "qval"):
        # Binned significance colorbar matching PVAL_BINS
        bins = PVAL_BINS
        n = len(bins) - 1
        colors = _pval_bin_colors(scheme, flip=flip)
        edges = np.linspace(0, 1, n + 1)
        scale = [[edges[i + j], colors[i]] for i in range(n) for j in (0, 1)]
        tickvals = [edges[i] for i in range(n + 1)]
        ticktext = [f"{v:g}" for v in bins]
        sig_title = "q-value" if color_by == "qval" else "p-value"
        return go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(
                colorscale=scale,
                showscale=True,
                cmin=0,
                cmax=1,
                color=[0],
                colorbar=dict(
                    title=sig_title,
                    len=length or COLORBAR_LEN,
                    x=COLORBAR_X,
                    y=y_center,
                    yanchor="middle",
                    thickness=COLORBAR_THICKNESS,
                    title_font=dict(size=COLORBAR_TITLE_SIZE),
                    tickfont=dict(size=COLORBAR_TICK_SIZE),
                    tickvals=tickvals,
                    ticktext=ticktext,
                ),
            ),
            hoverinfo="skip",
            showlegend=False,
        )
    else:
        title = "Distance"
        cmax_plot = min(cmax, D_UB - D_EPS) if np.isfinite(cmax) else D_UB - D_EPS
        cmin_plot = max(cmin, D_EPS) if np.isfinite(cmin) else D_EPS
        # Reverse colorscale if flipped (swap positions 0↔1)
        if flip:
            orig = px.colors.get_colorscale(scheme)
            cs = [[1 - pos, col] for pos, col in reversed(orig)]
        else:
            cs = scheme
        return go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(
                colorscale=cs,
                showscale=True,
                cmin=cmin_plot,
                cmax=cmax_plot,
                color=[cmin_plot],
                colorbar=dict(
                    title=title,
                    len=length or COLORBAR_LEN,
                    x=COLORBAR_X,
                    y=y_center,
                    yanchor="middle",
                    thickness=COLORBAR_THICKNESS,
                    title_font=dict(size=COLORBAR_TITLE_SIZE),
                    tickfont=dict(size=COLORBAR_TICK_SIZE),
                ),
            ),
            hoverinfo="skip",
            showlegend=False,
        )


def make_colorbar(bin_edges, colors, y_center=0.5, length=None, flip=False):
    """Create a colorbar trace for the interval panel."""
    # Reverse colors if flipped
    if flip:
        colors = list(reversed(colors))
    n = len(colors)
    edges = np.linspace(0, 1, n + 1)
    scale = [[edges[i + j], colors[i]] for i in range(n) for j in (0, 1)]
    ticks = [f"≤{bin_edges[i+1]:.3f}" for i in range(n)]
    dummy_y = np.linspace(0, 1, n)

    return go.Scatter(
        x=[None],
        y=[None],
        mode="markers",
        marker=dict(
            colorscale=scale,
            showscale=True,
            cmin=0,
            cmax=1,
            color=dummy_y,
            colorbar=dict(
                title="Distance",
                tickvals=[(edges[i] + edges[i + 1]) / 2 for i in range(n)],
                ticktext=ticks,
                len=length or COLORBAR_LEN,
                x=COLORBAR_X,
                y=y_center,
                yanchor="middle",
                thickness=COLORBAR_THICKNESS,
                title_font=dict(size=COLORBAR_TITLE_SIZE),
                tickfont=dict(size=COLORBAR_TICK_SIZE),
            ),
        ),
        hoverinfo="skip",
        showlegend=False,
    )


def compute_y_ticks(
    tip_order: list[str],
    y_range: Optional[tuple] = None,
    max_ticks: int = AXIS_Y_MAX_TICKS,
) -> tuple:
    """Compute y-axis tick positions and labels."""
    n = len(tip_order)
    if y_range is None:
        step = max(1, n // max_ticks)
        indices = list(range(0, n, step))
        if indices[-1] != n - 1:
            indices.append(n - 1)
    else:
        lo, hi = sorted(y_range)
        i_lo = max(0, min(int(round(lo)), n - 1))
        i_hi = max(0, min(int(round(hi)), n - 1))
        indices = list(range(i_lo, i_hi + 1))
        if len(indices) > max_ticks:
            indices = indices[:: len(indices) // max_ticks]
    return indices, [tip_order[i] for i in indices]


# =============================================================================
# ANNOTATIONS (GFF/GTF/TSV)
# =============================================================================


def _parse_gff_attrs(attr_str: str) -> dict:
    """Parse GFF3 (key=value) or GTF (key "value") attributes."""
    attrs = {}
    if pd.isna(attr_str) or not str(attr_str).strip():
        return attrs
    for field in str(attr_str).split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            key, _, val = field.partition("=")
            attrs[key.strip()] = val.strip()
        elif " " in field:
            key, _, val = field.partition(" ")
            attrs[key.strip()] = val.strip().strip('"')
    return attrs


def _is_gff_format(path: str) -> bool:
    """Heuristic: file is GFF3/GTF if the first non-comment, non-header line
    has ≥ 9 tab-separated columns and columns 3 and 4 (0-indexed) are integers.

    A plain TSV with a text header is rejected correctly because the header
    row's 4th and 5th fields ("start" / "stop") are not integers.
    """
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            # GFF has at least 9 columns; reject TSVs with text headers early
            if len(parts) < 9:
                return False
            try:
                int(parts[3])
                int(parts[4])
                return True
            except ValueError:
                return False
    return False


def _load_gff(path: str) -> Optional[pd.DataFrame]:
    """Load and normalize GFF3/GTF to internal schema."""
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 8:
                continue

            seqname, source, feature, start, end, score, strand, frame = cols[:8]
            try:
                start_i, end_i = int(start), int(end)
            except ValueError:
                continue

            attrs = _parse_gff_attrs(cols[8]) if len(cols) > 8 else {}

            locus_tag = (
                attrs.get("ID")
                or attrs.get("locus_tag")
                or attrs.get("Name")
                or attrs.get("ann_id")
                or f"{seqname}_{start_i}_{end_i}"
            )

            rows.append(
                {
                    "contig_id": seqname,
                    "locus_tag": locus_tag,
                    "ftype": feature if feature != "." else "misc",
                    "start": min(start_i, end_i),
                    "stop": max(start_i, end_i),
                    "strand": strand if strand in ("+", "-") else "+",
                    "ann_name": attrs.get("ann")
                    or attrs.get("ann_name")
                    or attrs.get("Name"),
                    "product": attrs.get("product") or attrs.get("description"),
                    "ec_number": attrs.get("ec_number"),
                    "source": source,
                    "score": score if score != "." else None,
                    "frame": frame if frame != "." else None,
                }
            )

    return pd.DataFrame(rows) if rows else None


def _load_gff_with_gffutils(path: str) -> Optional[pd.DataFrame]:
    """Load GFF3/GTF via gffutils (optional dependency).

    gffutils is imported lazily so its absence doesn't crash the script.
    Returns None if gffutils is not installed or parsing fails.
    """
    try:
        import gffutils  # noqa: PLC0415 — intentional lazy import
    except ImportError:
        return None

    try:
        db = gffutils.create_db(
            path,
            dbfn=":memory:",
            force=True,
            keep_order=True,
            merge_strategy="merge",
            sort_attribute_values=True,
        )

        rows = []
        for feature in db.all_features():
            if feature.start is None or feature.end is None:
                continue

            attrs = dict(feature.attributes)

            def _first(key):
                v = attrs.get(key, [None])
                return v[0] if v else None

            locus_tag = (
                _first("ID")
                or _first("locus_tag")
                or _first("Name")
                or _first("gene")
                or f"{feature.chrom}_{feature.start}_{feature.end}"
            )

            rows.append(
                {
                    "contig_id": feature.chrom,
                    "locus_tag": locus_tag,
                    "ftype": (
                        feature.featuretype if feature.featuretype != "." else "misc"
                    ),
                    "start": int(feature.start),
                    "stop": int(feature.end),
                    "strand": feature.strand if feature.strand in ("+", "-") else "+",
                    "ann_name": _first("gene")
                    or _first("ann")
                    or _first("ann_name")
                    or _first("Name"),
                    "product": _first("product")
                    or _first("description")
                    or _first("note"),
                    "ec_number": _first("ec_number"),
                    "source": feature.source,
                    "score": feature.score if feature.score != "." else None,
                    "frame": feature.frame if feature.frame != "." else None,
                }
            )

        return pd.DataFrame(rows) if rows else None

    except Exception as e:
        print(f"Warning: gffutils parsing failed: {e}")
        return None


def _load_custom_tsv(path: str) -> Optional[pd.DataFrame]:
    """Load the custom TSV annotation format (e.g. Prokka/SwissProt merged tables).

    Required columns: contig_id, locus_tag, ftype, start, stop, strand.
    Extra columns (prokka_gene, prokka_EC_number, prokka_product, swissprot_*)
    are kept as-is; _ann_hover knows how to read them directly.

    Additionally synthesises a unified ``ann_name`` column from whichever
    source columns are available so downstream code always has one place to look.
    """
    try:
        df = pd.read_csv(path, sep="\t", na_values=[".", ""], keep_default_na=True)

        # Drop unnamed numeric index column sometimes prepended by shell tools
        first_col = df.columns[0]
        if str(first_col) not in {
            "contig_id",
            "locus_tag",
            "ftype",
            "start",
            "stop",
            "strand",
        }:
            try:
                pd.to_numeric(df[first_col])
                df = df.drop(columns=[first_col])
            except (ValueError, TypeError):
                pass

        required = {"contig_id", "locus_tag", "ftype", "start", "stop", "strand"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")

        df = df.copy()
        # Normalise start/stop order
        df["start"] = df[["start", "stop"]].min(axis=1)
        df["stop"] = df[["start", "stop"]].max(axis=1)

        # Normalise strand
        df["strand"] = df["strand"].apply(lambda x: x if x in ("+", "-") else "+")

        # Synthesise a unified ann_name so rendering code has one canonical column
        if "ann_name" not in df.columns:
            gene_candidates = ["prokka_gene", "swissprot_gene", "gene"]
            existing = [c for c in gene_candidates if c in df.columns]
            if existing:
                df["ann_name"] = df[existing].apply(
                    lambda row: next(
                        (v for v in row if pd.notna(v) and str(v).strip()), None
                    ),
                    axis=1,
                )
            else:
                df["ann_name"] = None

        # Synthesise unified ec_number
        if "ec_number" not in df.columns:
            ec_candidates = ["prokka_EC_number", "swissprot_EC_number"]
            existing = [c for c in ec_candidates if c in df.columns]
            if existing:
                df["ec_number"] = df[existing].apply(
                    lambda row: next(
                        (v for v in row if pd.notna(v) and str(v).strip()), None
                    ),
                    axis=1,
                )
            else:
                df["ec_number"] = None

        # Synthesise unified product
        if "product" not in df.columns:
            prod_candidates = ["prokka_product", "swissprot_product"]
            existing = [c for c in prod_candidates if c in df.columns]
            if existing:
                df["product"] = df[existing].apply(
                    lambda row: next(
                        (v for v in row if pd.notna(v) and str(v).strip()), None
                    ),
                    axis=1,
                )
            else:
                df["product"] = None

        return df

    except Exception as e:
        print(f"Warning: custom TSV parsing failed: {e}")
        traceback.print_exc()
        return None


def load_annotations(path: Optional[str]) -> Optional[pd.DataFrame]:
    """Load annotations from GFF3, GTF, or custom TSV.

    Routing logic (in order):
      1. If the file looks like GFF/GTF → try gffutils first (better attribute
         handling), then fall back to the built-in line parser.
      2. Otherwise → treat as custom TSV directly; never run gffutils on it.

    gffutils is optional: if it is not installed the GFF path still works via
    the built-in parser; TSV files are unaffected.
    """
    if path is None:
        return None

    try:
        if _is_gff_format(path):
            # Try gffutils (no-op if not installed), then fall back to built-in
            result = _load_gff_with_gffutils(path)
            if result is not None:
                return result
            return _load_gff(path)
        else:
            return _load_custom_tsv(path)

    except Exception as e:
        print(f"Warning: could not load annotation file '{path}': {e}")
        traceback.print_exc()
        return None


def filter_annotations(
    df: Optional[pd.DataFrame], query_id: str, x_range: Optional[tuple] = None
) -> Optional[pd.DataFrame]:
    """Filter annotations by contig and optionally by genomic position."""
    if df is None:
        return None

    df_filt = df[df["contig_id"] == query_id].copy()
    if df_filt.empty:
        return None

    if x_range is not None:
        x_min, x_max = sorted(x_range)
        df_filt = df_filt[
            (df_filt["stop"] >= x_min) & (df_filt["start"] <= x_max)
        ].copy()

    return df_filt if not df_filt.empty else None


# =============================================================================
# GENE ARROWS (IGV-style)
# =============================================================================


def create_ann_arrow(start, end, strand, y_center, height=None):
    """Create polygon vertices for a direction-aware gene arrow.

    Arrow height is capped at ARROW_HEIGHT_MAX to prevent overflow into
    adjacent tracks.  The body/head ratio follows ARROW_BODY_HEIGHT_RATIO.
    Feature-specific styling (line_width, opacity) is applied by the caller
    via FEATURE_STYLES so this function stays pure-geometry.
    """
    if height is None:
        height = ANNOTATION_TRACK_HEIGHT

    length = max(end - start, 1)
    arrow_size = max(ARROW_SIZE_MIN, length * ARROW_HEAD_RATIO)

    # Body half-height and arrowhead half-height, both capped at ARROW_HEIGHT_MAX
    h = min(height * ARROW_BODY_HEIGHT_RATIO, ARROW_HEIGHT_MAX)  # body
    H = min(height * 0.80, ARROW_HEIGHT_MAX)  # head

    if strand == "+":
        body_end = end - arrow_size
        if body_end <= start:
            # Gene too short for a body — draw a simple triangle
            x = [start, end, start, start]
            y = [y_center - H, y_center, y_center + H, y_center - H]
        else:
            x = [start, body_end, body_end, end, body_end, body_end, start, start]
            y = [
                y_center - h,
                y_center - h,
                y_center - H,
                y_center,
                y_center + H,
                y_center + h,
                y_center + h,
                y_center - h,
            ]
    else:
        body_start = start + arrow_size
        if body_start >= end:
            x = [end, start, end, end]
            y = [y_center - H, y_center, y_center + H, y_center - H]
        else:
            x = [end, body_start, body_start, start, body_start, body_start, end, end]
            y = [
                y_center - h,
                y_center - h,
                y_center - H,
                y_center,
                y_center + H,
                y_center + h,
                y_center + h,
                y_center - h,
            ]

    return {"x": x, "y": y}


def _ann_hover(ann_row) -> str:
    """Build rich HTML hover text for a feature annotation."""
    length = int(ann_row["stop"]) - int(ann_row["start"])
    strand_symbol = "→" if ann_row["strand"] == "+" else "←"

    hover_text = f"<b>{ann_row['locus_tag']}</b>"
    hover_text += f"<br><i>{ann_row['ftype']}</i>"
    hover_text += f"<br>Position: {int(ann_row['start']):,} – {int(ann_row['stop']):,}"
    hover_text += f"<br>Length: {length:,} bp  {ann_row['strand']} {strand_symbol}"

    # Gene name — first non-null across candidate columns
    ann_name = None
    for col in ["ann_name", "prokka_gene", "swissprot_gene", "gene"]:
        if col in ann_row and pd.notna(ann_row[col]):
            ann_name = ann_row[col]
            break
    if ann_name:
        hover_text += f"<br><b>Gene:</b> {ann_name}"

    # EC number — prefer prokka, then swissprot, then generic
    ec_number = None
    for col in ["ec_number", "prokka_EC_number", "swissprot_EC_number"]:
        if col in ann_row and pd.notna(ann_row[col]):
            ec_number = ann_row[col]
            break
    if ec_number:
        hover_text += f"<br><b>EC:</b> {ec_number}"

    # Product description
    product = None
    for col in ["product", "prokka_product", "swissprot_product", "description"]:
        if col in ann_row and pd.notna(ann_row[col]):
            product = ann_row[col]
            break
    if product:
        product = str(product)
        if len(product) > 120:
            product = product[:117] + "…"
        hover_text += f"<br><b>Product:</b> {product}"

    # Functional annotations (COG/eggNOG, KEGG KO, Pfam) — shown if present
    extras = []
    for col in ["swissprot_eggNOG", "eggNOG", "cog"]:
        if col in ann_row and pd.notna(ann_row[col]):
            extras.append(f"COG: {ann_row[col]}")
            break
    for col in ["swissprot_KO", "KO", "kegg_ko"]:
        if col in ann_row and pd.notna(ann_row[col]):
            extras.append(f"KO: {ann_row[col]}")
            break
    pfam_vals = []
    for col in ["swissprot_Pfam", "Pfam", "pfam"]:
        if col in ann_row and pd.notna(ann_row[col]):
            pfam_vals = [d.strip() for d in str(ann_row[col]).split(",") if d.strip()]
            break
    if pfam_vals:
        pfam_text = ", ".join(pfam_vals[:3])
        if len(pfam_vals) > 3:
            pfam_text += f" (+{len(pfam_vals) - 3})"
        extras.append(f"Pfam: {pfam_text}")
    if extras:
        hover_text += "<br>" + " | ".join(extras[:3])

    return hover_text


def create_annotation_traces(df, x_range=None):
    """Create Plotly traces for annotation tracks with feature-specific styling."""
    if df is None or df.empty:
        return [], [], [], 0

    ftypes = sorted(df["ftype"].dropna().unique())
    ftype_to_y = {ft: i for i, ft in enumerate(ftypes)}

    traces = []
    x_min, x_max = sorted(x_range) if x_range else (None, None)

    for ftype in ftypes:
        color = FEATURE_COLORS.get(ftype, FEATURE_COLORS["default"])
        style = FEATURE_STYLES.get(ftype, FEATURE_STYLES["_default"])
        y_center = ftype_to_y[ftype]

        for _, ann in df[df["ftype"] == ftype].sort_values("start").iterrows():
            if x_min is not None and (ann["stop"] < x_min or ann["start"] > x_max):
                continue

            arrow = create_ann_arrow(
                int(ann["start"]), int(ann["stop"]), ann["strand"], y_center
            )
            hover_text = _ann_hover(ann)

            # Arrow body (filled polygon, hover disabled — avoids trace-name clutter)
            traces.append(
                go.Scatter(
                    x=arrow["x"],
                    y=arrow["y"],
                    mode="lines",
                    fill="toself",
                    fillcolor=color,
                    line=dict(color=color, width=style["line_width"]),
                    hoverinfo="skip",
                    showlegend=False,
                    opacity=style["opacity"],
                    name="",
                )
            )

            # Invisible centre point that carries the hover tooltip
            ann_center = (int(ann["start"]) + int(ann["stop"])) / 2
            traces.append(
                go.Scatter(
                    x=[ann_center],
                    y=[y_center],
                    mode="markers",
                    marker=dict(size=10, opacity=0.2, color=color),
                    hoverinfo="text",
                    text=[hover_text],
                    hoverlabel=dict(
                        namelength=-1,
                        bgcolor="rgba(255, 255, 255, 0.95)",
                        font=dict(size=12, family="Arial, sans-serif", color="#1a1a1a"),
                        bordercolor=color,
                    ),
                    showlegend=False,
                    name="",
                )
            )

    return traces, list(range(len(ftypes))), ftypes, len(ftypes)


# =============================================================================
# FIGURE BUILDING
# =============================================================================


def build_figure(
    tree_data,
    tip_order,
    tree_max_xy,
    df_intervals,
    bin_edges=None,
    colors=None,
    df_annotations=None,
    query_id=None,
    interval_lw=3.0,
    tree_lw=1.5,
    uirevision="base",
    y_range=None,
    x_range=None,
    x_limits=None,
    y_limits=None,
    leaf_distances=None,
    enum_only=True,
    color_by="dist",
    dist_range=(0.0, 0.5),
    scheme=COLORSCALE_DEFAULT,
    flip_colorscale=False,
):
    """Assemble the complete figure with tree, intervals, and optional annotations."""
    has_annot = df_annotations is not None and not df_annotations.empty

    # Determine colorbar position and length
    if has_annot:
        row1_bottom = ANNOTATION_ROW_HEIGHT + 0.05
        row1_top = 1.0
    else:
        row1_bottom, row1_top = 0.0, 1.0
    cb_y = (row1_bottom + row1_top) / 2.0
    cb_len = (row1_top - row1_bottom) * 0.85

    # Create subplot layout
    if has_annot:
        fig = make_subplots(
            rows=2,
            cols=2,
            shared_yaxes=True,
            shared_xaxes=True,
            column_widths=[PANEL_TREE_WIDTH, PANEL_INTERVAL_WIDTH],
            row_heights=[1 - ANNOTATION_ROW_HEIGHT, ANNOTATION_ROW_HEIGHT],
            horizontal_spacing=PANEL_SPACING,
            vertical_spacing=0.05,
        )
    else:
        fig = make_subplots(
            rows=1,
            cols=2,
            shared_yaxes=True,
            column_widths=[PANEL_TREE_WIDTH, PANEL_INTERVAL_WIDTH],
            horizontal_spacing=PANEL_SPACING,
        )

    # Colorbar in interval panel
    if enum_only:
        fig.add_trace(
            make_colorbar(
                bin_edges, colors, y_center=cb_y, length=cb_len, flip=flip_colorscale
            ),
            row=1,
            col=2,
        )
    else:
        # Compute continuous color range from data
        cmin, cmax = dist_range
        if not df_intervals.empty:
            cv = compute_color_values(df_intervals, color_by, dist_range)
            if cv is not None:
                cmin, cmax = cv["cmin"], cv["cmax"]
        fig.add_trace(
            make_continuous_colorbar(
                cmin,
                cmax,
                color_by,
                scheme,
                y_center=cb_y,
                length=cb_len,
                flip=flip_colorscale,
            ),
            row=1,
            col=2,
        )

    fig.update_layout(
        autosize=True,
        margin=FIG_MARGIN,
        plot_bgcolor=COLORS["plot_bg"],
        paper_bgcolor=COLORS["paper_bg"],
        hovermode="closest",
        uirevision=uirevision,
        font=dict(family=FIG_FONT, size=FIG_SIZE),
        dragmode="zoom",  # Enable drag-to-zoom by default
        showlegend=False,  # Cleaner look without legend
    )

    # Tree panel (col 1)
    fig.add_trace(
        go.Scattergl(
            x=tree_data["x"],
            y=tree_data["y"],
            hovertext=tree_data["text"],
            hovertemplate="%{hovertext}<extra></extra>",
            mode="lines",
            line=dict(color=COLORS["tree"], width=tree_lw),
            showlegend=False,
        ),
        row=1,
        col=1,
    )

    # Interval panel (col 2)
    if not df_intervals.empty:
        if enum_only:
            interval_traces = batch_by_color(
                df_intervals, bin_edges, colors, leaf_distances, flip=flip_colorscale
            )
        else:
            interval_traces, _, _ = batch_by_continuous_color(
                df_intervals,
                color_by,
                dist_range,
                scheme,
                leaf_distances,
                flip=flip_colorscale,
            )
        for color, data in interval_traces.items():
            fig.add_trace(
                go.Scattergl(
                    x=data["x"],
                    y=data["y"],
                    mode="lines",
                    line=dict(width=interval_lw, color=color),
                    text=data["text"],
                    hovertemplate="%{text}<extra></extra>",
                    showlegend=False,
                ),
                row=1,
                col=2,
            )

    # Annotation panel (row 2)
    ann_y_vals, ann_y_text, ann_n = [], [], 0
    if has_annot:
        ann_df = filter_annotations(df_annotations, query_id, x_range)
        if ann_df is not None and not ann_df.empty:
            ann_traces, ann_y_vals, ann_y_text, ann_n = create_annotation_traces(
                ann_df, x_range
            )
            for trace in ann_traces:
                fig.add_trace(trace, row=2, col=2)

    # ---- Axes ----
    tree_x_range = [
        -(tree_max_xy[0] * (TREE_X_MARGIN - 1)),
        tree_max_xy[0] * TREE_X_MARGIN,
    ]
    x_range = x_range if x_range is not None else x_limits

    # Tree X-axis (fixed)
    fig.update_xaxes(
        title="Branch length",
        title_font=dict(size=AXIS_TITLE_SIZE),
        range=tree_x_range,
        fixedrange=True,
        side="bottom",
        showline=True,
        showgrid=False,
        showticklabels=True,
        linewidth=AXIS_LINE_WIDTH,
        linecolor=COLORS["axis"],
        mirror=False,
        tickfont=dict(size=AXIS_TICK_SIZE),
        ticks="outside",
        nticks=5,
        row=1,
        col=1,
    )
    fig.update_yaxes(
        autorange=y_range is None,
        range=y_range if y_range else y_limits,
        showline=False,
        zeroline=False,
        showgrid=False,
        showticklabels=False,
        minallowed=y_limits[0],
        maxallowed=y_limits[1],
        row=1,
        col=1,
    )

    # Interval X-axis
    fig.update_xaxes(
        title="Position (bp)" if not has_annot else "",
        title_font=dict(size=AXIS_TITLE_SIZE),
        range=x_range,
        tickmode="auto",
        nticks=AXIS_X_NTICKS,
        showgrid=True,
        gridcolor=COLORS["grid"],
        gridwidth=AXIS_GRID_WIDTH,
        showline=True,
        linewidth=AXIS_LINE_WIDTH,
        linecolor=COLORS["axis"],
        mirror=False,
        tickfont=dict(size=AXIS_TICK_SIZE),
        showticklabels=not has_annot,
        ticks="outside" if not has_annot else "",
        minallowed=x_limits[0],
        maxallowed=x_limits[1],
        row=1,
        col=2,
    )
    tickvals, ticktext = compute_y_ticks(tip_order, y_range)
    fig.update_yaxes(
        autorange=y_range is None,
        range=y_range if y_range else y_limits,
        tickvals=tickvals,
        ticktext=ticktext,
        tickfont=dict(size=AXIS_TICK_SIZE),
        showgrid=False,
        zeroline=False,
        showline=True,
        linewidth=AXIS_LINE_WIDTH,
        linecolor=COLORS["axis"],
        mirror=False,
        minallowed=y_limits[0],
        maxallowed=y_limits[1],
        row=1,
        col=2,
    )

    # Annotation panel axes
    if has_annot:
        ann_y_range = [-0.5, max(ann_n - 0.5, 0.5)]

        fig.update_xaxes(
            title="Position (bp)",
            title_font=dict(size=AXIS_TITLE_SIZE),
            range=x_range,
            tickmode="auto",
            nticks=AXIS_X_NTICKS,
            showgrid=True,
            gridcolor=COLORS["grid"],
            gridwidth=AXIS_GRID_WIDTH,
            showline=True,
            linewidth=AXIS_LINE_WIDTH,
            linecolor=COLORS["axis"],
            mirror=False,
            tickfont=dict(size=AXIS_TICK_SIZE),
            showticklabels=True,
            ticks="outside",
            minallowed=x_limits[0],
            maxallowed=x_limits[1],
            row=2,
            col=2,
        )
        fig.update_yaxes(
            autorange=False,
            fixedrange=True,
            range=ann_y_range,
            tickvals=ann_y_vals,
            ticktext=list(ann_y_text),
            tickfont=dict(size=12, color="#555"),
            side="right",
            showline=False,
            zeroline=False,
            showgrid=True,
            gridcolor=COLORS["grid"],
            gridwidth=1,
            showticklabels=True,
            row=2,
            col=2,
        )
        # Row 2, col 1 (hidden, under tree)
        fig.update_yaxes(
            showline=False,
            showgrid=False,
            showticklabels=False,
            fixedrange=True,
            zeroline=False,
            row=2,
            col=1,
        )
        fig.update_xaxes(
            showline=False,
            showgrid=False,
            showticklabels=False,
            title="",
            fixedrange=True,
            row=2,
            col=1,
        )

    return fig


# =============================================================================
# HELPERS
# =============================================================================


def nearest_value(value, values):
    """Find nearest value in sorted list."""
    idx = np.searchsorted(values, value)
    if idx == 0:
        return values[0]
    if idx == len(values):
        return values[-1]
    left, right = values[idx - 1], values[idx]
    return left if (value - left) < (right - value) else right


def make_slider_marks(values):
    """Build slider marks. Interior marks use zero-width space so Dash
    renders a tick mark without visible text (empty string hides the tick)."""
    if not values:
        return {}
    return {
        values[0]: str(values[0]),
        values[-1]: str(values[-1]),
        **{v: "\u200b" for v in values[1:-1]},
    }


# =============================================================================
# UI BUILDING
# =============================================================================


def control_label(text):
    return html.Label(
        text,
        style={
            "fontWeight": "600",
            "marginRight": f"{UI_LABEL_MARGIN}px",
            "fontSize": UI_LABEL_SIZE,
            "color": UI_LABEL_COLOR,
            "whiteSpace": "nowrap",
        },
    )


def nav_button_style(position="middle"):
    radius = {"left": "4px 0 0 4px", "right": "0 4px 4px 0", "middle": "4px"}[position]
    style = {
        "padding": UI_NAV_PADDING,
        "fontSize": UI_NAV_SIZE,
        "cursor": "pointer",
        "border": f"{UI_NAV_BORDER}px solid {COLORS['button_border']}",
        "borderRadius": radius,
        "backgroundColor": COLORS["button_bg"],
        "color": UI_TOGGLE_COLOR_ACTIVE,
        "transition": "background-color 0.12s ease",
        "whiteSpace": "nowrap",
    }
    if position == "left":
        style["borderRight"] = "none"
    return style


def export_button_style():
    return {
        "padding": UI_EXPORT_PADDING,
        "fontSize": UI_EXPORT_FONT_SIZE,
        "cursor": "pointer",
        "border": f"1px solid {COLORS['button_border']}",
        "borderRadius": "4px",
        "backgroundColor": COLORS["button_bg"],
        "fontWeight": "600",
        "color": UI_TOGGLE_COLOR_ACTIVE,
        "transition": "background-color 0.12s ease",
        "whiteSpace": "nowrap",
    }


def export_input_style(width):
    return {
        "width": width,
        "fontSize": UI_EXPORT_FONT_SIZE,
        "padding": "3px 5px",
        "border": f"1px solid {COLORS['button_border']}",
        "borderRadius": "4px",
        "textAlign": "center",
        "color": UI_TOGGLE_COLOR_ACTIVE,
    }


def toggle_style(is_active, position="middle"):
    radius = {"left": "4px 0 0 4px", "right": "0 4px 4px 0", "middle": "0"}[position]
    style = {
        "padding": UI_TOGGLE_PADDING,
        "fontSize": UI_TOGGLE_SIZE,
        "cursor": "pointer",
        "border": f"{UI_TOGGLE_BORDER}px solid {COLORS['button_border']}",
        "borderRadius": radius,
        "backgroundColor": (
            COLORS["button_bg"] if is_active else COLORS["button_bg_inactive"]
        ),
        "fontWeight": "600" if is_active else "normal",
        "color": UI_TOGGLE_COLOR_ACTIVE if is_active else UI_TOGGLE_COLOR_INACTIVE,
        "transition": "background-color 0.12s ease, color 0.12s ease",
        "whiteSpace": "nowrap",
    }
    if position != "right":
        style["borderRight"] = "none"
    return style


def divider():
    return html.Div(
        style={
            "width": "1px",
            "alignSelf": "stretch",
            "backgroundColor": COLORS["button_border"],
            "opacity": "0.35",
            "margin": "0 4px",
        }
    )


def control_panel(children):
    return html.Div(
        children,
        style={
            "display": "flex",
            "alignItems": "center",
            "gap": UI_PANEL_GAP,
            "padding": UI_PANEL_PADDING,
            "backgroundColor": COLORS["panel_bg"],
            "borderRadius": UI_PANEL_RADIUS,
            "marginBottom": UI_PANEL_MARGIN_BOTTOM,
            "boxShadow": UI_PANEL_SHADOW,
            "flexWrap": "wrap",
        },
    )


def build_layout(
    seq_ids,
    dist_ths,
    has_pruned,
    has_strand=False,
    has_is_rc=False,
    has_qvalue=False,
    initial_prune=True,
    enum_only=True,
):
    dmin = dist_ths[0] if dist_ths else 0.0
    dmax = dist_ths[-1] if dist_ths else 0.0

    return html.Div(
        [
            # ── Stores ──────────────────────────────────────────────────────
            dcc.Store(id="y-range-store"),
            dcc.Store(id="x-range-store"),
            dcc.Store(id="tree-view-store", data="phylogeny"),
            dcc.Store(id="prune-store", data="pruned" if initial_prune else "full"),
            dcc.Store(id="strand-store", data="both"),
            dcc.Store(id="direction-store", data="both"),
            dcc.Store(id="mode-store", data="focus"),
            dcc.Store(id="color-by-store", data="dist"),
            dcc.Store(id="sig-metric-store", data="q" if has_qvalue else "p"),
            dcc.Store(id="fold-store", data="both"),
            dcc.Store(id="colorscale-flip-store", data=False),
            dcc.Store(id="interval-count-store", data=0),
            html.Div(id="keyboard-listener", tabIndex="0",
                      style={"position": "fixed", "top": 0, "left": 0,
                             "width": "100%", "height": "100%",
                             "zIndex": -1, "opacity": 0}),
            dcc.Store(id="keyboard-nav-store", data=0),
            dcc.Download(id="download-data"),
            dcc.Download(id="download-plot"),
            # ── Control panel ────────────────────────────────────────────────
            control_panel(
                [
                    # Query navigator: ◀ [dropdown] ▶
                    html.Div(
                        [
                            control_label("Query:"),
                            html.Button(
                                "◀",
                                id="prev-query-btn",
                                n_clicks=0,
                                style=nav_button_style("left"),
                            ),
                            dcc.Dropdown(
                                id="query-dropdown",
                                options=[{"label": q, "value": q} for q in seq_ids],
                                value=seq_ids[0] if seq_ids else None,
                                clearable=False,
                                style={"width": UI_QUERY_DROPDOWN_WIDTH},
                            ),
                            html.Button(
                                "▶",
                                id="next-query-btn",
                                n_clicks=0,
                                style=nav_button_style("right"),
                            ),
                        ],
                        style={"display": "flex", "alignItems": "center", "gap": 0},
                    ),
                    # Distance controls: slider + focus/overlap mode toggle (enum-only)
                    html.Div(
                        [
                            divider(),
                            html.Div(
                                [
                                    control_label("Distance:"),
                                    html.Div(
                                        dcc.Slider(
                                            id="dist-slider-focus",
                                            min=dmin,
                                            max=dmax,
                                            step=None,
                                            marks=make_slider_marks(dist_ths),
                                            value=dist_ths[0] if dist_ths else dmin,
                                            tooltip={
                                                "placement": "bottom",
                                                "always_visible": True,
                                            },
                                            updatemode="mouseup",
                                            included=False,
                                            className="focus-slider",
                                        ),
                                        id="slider-focus-wrap",
                                        style={
                                            "minWidth": UI_SLIDER_MIN_WIDTH,
                                            "flexShrink": 1,
                                        },
                                    ),
                                    html.Div(
                                        dcc.RangeSlider(
                                            id="dist-slider-overlap",
                                            min=dmin,
                                            max=dmax,
                                            step=None,
                                            marks=make_slider_marks(dist_ths),
                                            value=[dmin, dmax],
                                            tooltip={
                                                "placement": "bottom",
                                                "always_visible": False,
                                            },
                                            updatemode="mouseup",
                                            allowCross=False,
                                        ),
                                        id="slider-overlap-wrap",
                                        style={
                                            "minWidth": UI_SLIDER_MIN_WIDTH,
                                            "flexShrink": 1,
                                            "display": "none",
                                        },
                                    ),
                                    html.Div(
                                        [
                                            html.Button(
                                                "focus",
                                                id="mode-focus-btn",
                                                n_clicks=0,
                                                style=toggle_style(True, "left"),
                                            ),
                                            html.Button(
                                                "overlap",
                                                id="mode-overlap-btn",
                                                n_clicks=0,
                                                style=toggle_style(False, "right"),
                                            ),
                                        ],
                                        style={"display": "flex", "marginLeft": 8},
                                    ),
                                ],
                                style={
                                    "display": "flex",
                                    "alignItems": "center",
                                    "gap": 4,
                                },
                            ),
                        ],
                        style={
                            "display": "flex" if enum_only else "none",
                            "alignItems": "center",
                        },
                    ),
                    # Color-by toggle: dist / p-value / q-value (non-enum only)
                    html.Div(
                        [
                            divider(),
                            control_label("Color by:"),
                            html.Button(
                                "dist",
                                id="colorby-dist-btn",
                                n_clicks=0,
                                style=toggle_style(True, "left"),
                            ),
                            html.Button(
                                "p-value",
                                id="colorby-pval-btn",
                                n_clicks=0,
                                style=toggle_style(False, "middle"),
                            ),
                            html.Button(
                                "q-value",
                                id="colorby-qval-btn",
                                n_clicks=0,
                                style=toggle_style(False, "right"),
                            ),
                        ],
                        style={
                            "display": "flex" if not enum_only else "none",
                            "alignItems": "center",
                            "gap": 0,
                        },
                    ),
                    # Fold toggle: <1 / both / >1 (non-enum only)
                    html.Div(
                        [
                            divider(),
                            control_label("Fold:"),
                            html.Button(
                                "<1",
                                id="fold-lt-btn",
                                n_clicks=0,
                                style=toggle_style(False, "left"),
                            ),
                            html.Button(
                                "both",
                                id="fold-both-btn",
                                n_clicks=0,
                                style=toggle_style(True, "middle"),
                            ),
                            html.Button(
                                ">1",
                                id="fold-gt-btn",
                                n_clicks=0,
                                style=toggle_style(False, "right"),
                            ),
                        ],
                        style={
                            "display": "flex" if not enum_only else "none",
                            "alignItems": "center",
                            "gap": 0,
                        },
                    ),
                    # Significance filter slider (non-enum only)
                    html.Div(
                        [
                            divider(),
                            html.Span(
                                "q \u2264" if has_qvalue else "p \u2264",
                                id="sig-filter-label",
                                style={
                                    "fontWeight": "600",
                                    "marginRight": f"{UI_LABEL_MARGIN}px",
                                    "fontSize": UI_LABEL_SIZE,
                                    "color": UI_LABEL_COLOR,
                                    "whiteSpace": "nowrap",
                                },
                            ),
                            html.Div(
                                [
                                    html.Button(
                                        "p",
                                        id="sig-p-btn",
                                        n_clicks=0,
                                        style=toggle_style(not has_qvalue, "left"),
                                    ),
                                    html.Button(
                                        "q",
                                        id="sig-q-btn",
                                        n_clicks=0,
                                        style=toggle_style(has_qvalue, "right"),
                                    ),
                                ],
                                style={"display": "flex", "marginRight": 6},
                            ),
                            html.Div(
                                dcc.Slider(
                                    id="pval-slider",
                                    min=-10,
                                    max=0,
                                    step=None,
                                    marks={
                                        0: "1",
                                        -10: "10⁻¹⁰",
                                        **{
                                            v: "\u200b"
                                            for v in [-1, -1.301, -2, -3, -4, -5, -7]
                                        },
                                    },
                                    value=0,
                                    tooltip={
                                        "placement": "bottom",
                                        "always_visible": False,
                                    },
                                    updatemode="mouseup",
                                    included=True,
                                ),
                                style={"minWidth": 240, "flexShrink": 1},
                            ),
                        ],
                        style={
                            "display": "flex" if not enum_only else "none",
                            "alignItems": "center",
                            "gap": 4,
                        },
                    ),
                    divider(),
                    # Tree: phylogeny/cladogram + optional pruned/full
                    html.Div(
                        [
                            control_label("Tree:"),
                            html.Button(
                                "phylogeny",
                                id="phylogeny-btn",
                                n_clicks=0,
                                style=toggle_style(True, "left"),
                            ),
                            html.Button(
                                "cladogram",
                                id="cladogram-btn",
                                n_clicks=0,
                                style=toggle_style(False, "right"),
                            ),
                            (
                                html.Div(
                                    [
                                        html.Button(
                                            "pruned",
                                            id="pruned-btn",
                                            n_clicks=0,
                                            style=toggle_style(initial_prune, "left"),
                                        ),
                                        html.Button(
                                            "filtered",
                                            id="query-btn",
                                            n_clicks=0,
                                            style=toggle_style(False, "middle"),
                                        ),
                                        html.Button(
                                            "full",
                                            id="full-btn",
                                            n_clicks=0,
                                            style=toggle_style(
                                                not initial_prune, "right"
                                            ),
                                        ),
                                    ],
                                    style={"display": "flex", "marginLeft": 8},
                                )
                                if has_pruned
                                else html.Div()
                            ),
                        ],
                        style={"display": "flex", "alignItems": "center", "gap": 0},
                    ),
                    # Enum: filter by physical strand (IS_RC: fw / rc)
                    html.Div(
                        [
                            divider(),
                            control_label("Direction:"),
                            html.Button(
                                "fw",
                                id="direction-fw-btn",
                                n_clicks=0,
                                style=toggle_style(False, "left"),
                            ),
                            html.Button(
                                "both",
                                id="direction-both-btn",
                                n_clicks=0,
                                style=toggle_style(True, "middle"),
                            ),
                            html.Button(
                                "rc",
                                id="direction-rc-btn",
                                n_clicks=0,
                                style=toggle_style(False, "right"),
                            ),
                        ],
                        style={
                            "display": "flex" if enum_only and has_is_rc else "none",
                            "alignItems": "center",
                            "gap": 0,
                        },
                    ),
                    # Continuous: filter by closer/farther (STRAND + / -)
                    html.Div(
                        [
                            divider(),
                            control_label("Strand:"),
                            html.Button(
                                "+ (closer)",
                                id="strand-plus-btn",
                                n_clicks=0,
                                style=toggle_style(False, "left"),
                            ),
                            html.Button(
                                "both",
                                id="strand-both-btn",
                                n_clicks=0,
                                style=toggle_style(True, "middle"),
                            ),
                            html.Button(
                                "- (farther)",
                                id="strand-minus-btn",
                                n_clicks=0,
                                style=toggle_style(False, "right"),
                            ),
                        ],
                        style={
                            "display": "flex" if (not enum_only and has_strand) else "none",
                            "alignItems": "center",
                            "gap": 0,
                        },
                    ),
                    # Color scale + Export — pushed to the right with marginLeft: auto
                    html.Div(
                        [
                            control_label("Color:"),
                            dcc.Dropdown(
                                id="colorscheme-dropdown",
                                options=[
                                    {"label": cs.capitalize(), "value": cs}
                                    for cs in COLORSCALE_OPTIONS
                                ],
                                value=COLORSCALE_DEFAULT,
                                clearable=False,
                                style={"width": UI_COLORSCALE_DROPDOWN_WIDTH},
                            ),
                            html.Button(
                                "⇄",
                                id="flip-colorscale-btn",
                                n_clicks=0,
                                style={
                                    **toggle_style(False, "right"),
                                    "marginLeft": 4,
                                    "padding": "0 6px",
                                    "fontSize": "12px",
                                },
                            ),
                            divider(),
                            control_label("Export:"),
                            html.Button(
                                "Data",
                                id="export-btn",
                                n_clicks=0,
                                style=export_button_style(),
                            ),
                            html.Span(
                                "W",
                                style={
                                    "fontSize": UI_EXPORT_FONT_SIZE - 1,
                                    "marginLeft": 6,
                                    "marginRight": 2,
                                    "color": "#888",
                                },
                            ),
                            dcc.Input(
                                id="export-width",
                                type="number",
                                value=UI_EXPORT_W_DEFAULT,
                                min=100,
                                max=10000,
                                step=10,
                                style=export_input_style(UI_EXPORT_INPUT_WIDTH_W),
                            ),
                            html.Span(
                                "H",
                                style={
                                    "fontSize": UI_EXPORT_FONT_SIZE - 1,
                                    "margin": "0 2px 0 5px",
                                    "color": "#888",
                                },
                            ),
                            dcc.Input(
                                id="export-height",
                                type="number",
                                value=UI_EXPORT_H_DEFAULT,
                                min=100,
                                max=10000,
                                step=10,
                                style=export_input_style(UI_EXPORT_INPUT_WIDTH_H),
                            ),
                            html.Button(
                                "PDF",
                                id="export-pdf-btn",
                                n_clicks=0,
                                style={**export_button_style(), "marginLeft": 4},
                            ),
                        ],
                        style={
                            "display": "flex",
                            "alignItems": "center",
                            "gap": 4,
                            "marginLeft": "auto",
                        },
                    ),
                ]
            ),
            # ── Query title ──────────────────────────────────────────────────
            html.Div(
                id="query-title",
                style={
                    "textAlign": "center",
                    "fontWeight": "700",
                    "fontSize": UI_TITLE_SIZE,
                    "color": UI_TITLE_COLOR,
                    "marginBottom": 6,
                    "letterSpacing": "0.02em",
                },
            ),
            # ── Main graph ───────────────────────────────────────────────────
            html.Div(
                dcc.Graph(
                    id="graph",
                    style={"height": UI_GRAPH_HEIGHT},
                    config={
                        "doubleClick": "reset",
                        "displayModeBar": True,
                        "displaylogo": False,
                        "modeBarButtonsToRemove": ["lasso2d", "select2d"],
                        "scrollZoom": True,
                    },
                ),
                style={"overflowY": "auto", "height": UI_CONTAINER_HEIGHT},
            ),
        ],
        style={"fontFamily": FIG_FONT, "padding": "6px 10px"},
    )


# =============================================================================
# CALLBACK FACTORIES
# =============================================================================


def _make_toggle_callbacks(
    app: Dash,
    store_id: str,
    button_ids: list[str],
    values: list[str],
    default_value: str,
    reset_ranges: bool = True,
) -> None:
    """Register toggle button callbacks for a group of mutually exclusive buttons.

    Args:
        app: Dash application instance.
        store_id: ID of the dcc.Store to hold the toggle state.
        button_ids: List of button component IDs.
        values: List of values corresponding to each button.
        default_value: Value to use when no button triggered the callback.
        reset_ranges: If True, also reset y-range and x-range stores on toggle.
    """
    outputs = [Output(store_id, "data")]
    if reset_ranges:
        outputs.extend(
            [
                Output("y-range-store", "data", allow_duplicate=True),
                Output("x-range-store", "data", allow_duplicate=True),
            ]
        )

    @app.callback(
        *outputs,
        *[Input(bid, "n_clicks") for bid in button_ids],
        prevent_initial_call=True,
    )
    def _toggle_callback(*clicks):
        triggered = ctx.triggered_id
        for bid, val in zip(button_ids, values):
            if triggered == bid:
                if reset_ranges:
                    return val, None, None
                return val
        # Default: return default_value
        if reset_ranges:
            return default_value, None, None
        return default_value

    @app.callback(
        *[Output(bid, "style") for bid in button_ids], Input(store_id, "data")
    )
    def _style_callback(current_value):
        n = len(button_ids)
        positions = (
            ["left"] + ["middle"] * (n - 2) + ["right"] if n > 2 else ["left", "right"]
        )
        return [
            toggle_style(current_value == val, pos)
            for val, pos in zip(values, positions)
        ]


def _make_two_button_toggle(
    app: Dash,
    store_id: str,
    btn_a: str,
    btn_b: str,
    val_a: str,
    val_b: str,
    default_a: bool = True,
) -> None:
    """Register callbacks for a simple two-button toggle.

    Args:
        app: Dash application instance.
        store_id: ID of the dcc.Store to hold the toggle state.
        btn_a: ID of first button.
        btn_b: ID of second button.
        val_a: Value when first button is active.
        val_b: Value when second button is active.
        default_a: Whether first button is default active.
    """

    @app.callback(
        Output(store_id, "data"), Input(btn_a, "n_clicks"), Input(btn_b, "n_clicks")
    )
    def _set_value(a_clicks, b_clicks):
        return val_b if ctx.triggered_id == btn_b else val_a

    @app.callback(
        Output(btn_a, "style"),
        Output(btn_b, "style"),
        Input(store_id, "data"),
    )
    def _style(value):
        return toggle_style(value == val_a, "left"), toggle_style(
            value == val_b, "right"
        )


# =============================================================================
# APP FACTORY
# =============================================================================


def create_app(
    input_path: str,
    tree_path: str,
    query: Optional[str] = None,
    annotation_path: Optional[str] = None,
    enum_only: Optional[bool] = None,
) -> Dash:
    """Create and configure the Dash application."""
    df = load_gdiff_tsv(input_path)
    plot_mode = detect_plot_mode(df, enum_only=enum_only)
    enum_ui = is_enum_plot_mode(plot_mode)
    df_annot = load_annotations(annotation_path)

    if plot_mode == "legacy_enum":
        required = {
            "QUERY_ID",
            "REF_ID",
            "INTERVAL_START",
            "INTERVAL_END",
            "SEQ_LEN",
            "DIST_TH",
        }
    elif plot_mode == "enum":
        required = {
            "QUERY_ID",
            "REF_ID",
            "INTERVAL_START",
            "INTERVAL_END",
            "SEQ_LEN",
            "STRAND",
            "DIST",
            "MASK",
            "D_INTERVAL",
        }
    else:
        required = {
            "QUERY_ID",
            "REF_ID",
            "INTERVAL_START",
            "INTERVAL_END",
            "SEQ_LEN",
            "STRAND",
            "DIST",
            "PERCENTILE",
        }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {', '.join(sorted(missing))}")

    if not (df["INTERVAL_START"] != df["INTERVAL_END"]).any():
        raise ValueError("No intervals found (all INTERVAL_START == INTERVAL_END)")

    # Tree setup
    full_tree = load_tree(tree_path)
    leaf_names = {leaf.name for leaf in full_tree.iter_leaves()}

    if query and query not in leaf_names:
        raise ValueError(
            f"Query '{query}' not found. Available: {', '.join(sorted(leaf_names)[:10])}..."
        )

    retained = get_retained_leaves(df)
    has_pruned = len(retained) > 0
    pruned_tree = prune_tree(full_tree, retained) if has_pruned else None

    # Query-specific tree (only leaves with matches for current query)
    query_retained = get_query_retained_leaves(df, query) if query else set()
    has_query_pruned = len(query_retained) > 0
    query_tree = prune_tree(full_tree, query_retained) if has_query_pruned else None

    def tree_data_for(t, q):
        distances = compute_path_distances(t, q) if q else None
        leaf_dist = (
            {l.name: distances[l] for l in t.iter_leaves() if l in (distances or {})}
            if distances
            else None
        )
        phylo = compute_tree_layout(t, "phylogeny", distances)
        clado = compute_tree_layout(t, "cladogram", distances)
        return {
            "phylo_layout": phylo,
            "clado_layout": clado,
            "tip_order": phylo[1],
            "leaf_distances": leaf_dist,
        }

    full_data = tree_data_for(full_tree, query)
    pruned_data = tree_data_for(pruned_tree, query) if has_pruned else full_data
    query_data = tree_data_for(query_tree, query) if has_query_pruned else full_data

    full_df = add_tip_order(df, full_data["tip_order"])
    pruned_df = add_tip_order(df, pruned_data["tip_order"]) if has_pruned else full_df
    query_df = (
        add_tip_order(df, query_data["tip_order"]) if has_query_pruned else full_df
    )

    if enum_ui:
        dist_ths = tuple(get_distance_thresholds(df))
        if not dist_ths:
            raise ValueError("No distance thresholds in data.")
        global_dist_range = (D_EPS, min(0.5, D_UB - D_EPS))
    else:
        dist_ths = ()
        global_dist_range = compute_distance_range(df)

    seq_ids = get_sequence_identifiers(df)
    has_strand = has_informative_strand(df)
    has_is_rc = has_is_rc_column(df)
    has_qvalue = (
        "QVALUE" in df.columns
        and (df["QVALUE"].notna() & np.isfinite(df["QVALUE"])).any()
    )

    filter_full = make_cached_filter(full_df, enum_ui)
    filter_pruned = (
        make_cached_filter(pruned_df, enum_ui) if has_pruned else filter_full
    )
    filter_query = (
        make_cached_filter(query_df, enum_ui) if has_query_pruned else filter_full
    )

    app = Dash(__name__)

    # Store the full tree for dynamic query pruning
    app.full_tree = full_tree
    app.df = df
    app.enum_only = enum_ui
    app.plot_mode = plot_mode
    app.layout = build_layout(
        seq_ids,
        dist_ths,
        has_pruned,
        has_strand=has_strand,
        has_is_rc=has_is_rc,
        has_qvalue=has_qvalue,
        enum_only=enum_ui,
    )

    # ---- Navigation callbacks ----
    app.clientside_callback(
        """
        function(id) {
            document.addEventListener('keydown', function(e) {
                if (e.target.tagName === 'INPUT' || e.target.tagName === 'TEXTAREA'
                    || e.target.tagName === 'SELECT') return;
                var prev = document.getElementById('prev-query-btn');
                var next = document.getElementById('next-query-btn');
                if (e.key === 'ArrowLeft' && prev) { prev.click(); }
                else if (e.key === 'ArrowRight' && next) { next.click(); }
            });
            return window.dash_clientside.no_update;
        }
        """,
        Output("keyboard-nav-store", "data"),
        Input("keyboard-listener", "id"),
    )

    @app.callback(
        Output("query-dropdown", "value"),
        Input("prev-query-btn", "n_clicks"),
        Input("next-query-btn", "n_clicks"),
        State("query-dropdown", "value"),
        prevent_initial_call=True,
    )
    def navigate_query(prev, next_, current):
        if current is None or not seq_ids:
            return no_update
        idx = seq_ids.index(current) if current in seq_ids else 0
        if ctx.triggered_id == "prev-query-btn":
            idx = (idx - 1) % len(seq_ids)
        else:
            idx = (idx + 1) % len(seq_ids)
        return seq_ids[idx]

    @app.callback(
        Output("query-title", "children"),
        Input("query-dropdown", "value"),
        Input("interval-count-store", "data"),
    )
    def update_title(seq_id, n_intervals):
        if seq_id is None:
            return ""
        idx = seq_ids.index(seq_id) if seq_id in seq_ids else 0
        count_str = f"  [{n_intervals} intervals]" if n_intervals else ""
        return f"{seq_id}  ({idx + 1}/{len(seq_ids)}){count_str}"

    # ---- Range stores (clientside for zero-latency zoom/pan/scroll) ----
    app.clientside_callback(
        """
        function(relayoutData) {
            var no_update = window.dash_clientside.no_update;
            if (!relayoutData) return [no_update, no_update];

            var keys = Object.keys(relayoutData);

            // Double-click / autorange reset → clear both stores
            if (keys.some(function(k) { return k.endsWith('.autorange'); })) {
                return [null, null];
            }

            // Y range (shared axis is always yaxis)
            var y = no_update;
            if ('yaxis.range[0]' in relayoutData && 'yaxis.range[1]' in relayoutData) {
                y = [relayoutData['yaxis.range[0]'], relayoutData['yaxis.range[1]']];
            } else if ('yaxis.range' in relayoutData) {
                y = relayoutData['yaxis.range'];
            }

            // X range — check both possible axis positions (1-col vs 2-col layout)
            var xaxes = ['xaxis2', 'xaxis4'];
            for (var i = 0; i < xaxes.length; i++) {
                var key = xaxes[i];
                if ((key + '.range[0]') in relayoutData && (key + '.range[1]') in relayoutData) {
                    return [y, [relayoutData[key + '.range[0]'], relayoutData[key + '.range[1]']]];
                }
                if ((key + '.range') in relayoutData) {
                    return [y, relayoutData[key + '.range']];
                }
            }

            return [y, no_update];
        }
        """,
        Output("y-range-store", "data"),
        Output("x-range-store", "data"),
        Input("graph", "relayoutData"),
    )

    # ---- Tree view toggle ----
    _make_two_button_toggle(
        app,
        "tree-view-store",
        "phylogeny-btn",
        "cladogram-btn",
        "phylogeny",
        "cladogram",
    )

    # ---- Prune toggle ----
    if has_pruned:
        _make_toggle_callbacks(
            app,
            "prune-store",
            ["pruned-btn", "query-btn", "full-btn"],
            ["pruned", "filtered", "full"],
            "full",
        )

    # ---- Direction toggle (enum: fw/rc via IS_RC) ----
    _make_toggle_callbacks(
        app,
        "direction-store",
        ["direction-fw-btn", "direction-both-btn", "direction-rc-btn"],
        ["fw", "both", "rc"],
        "both",
        reset_ranges=False,
    )

    # ---- Strand toggle (continuous: closer/farther via STRAND) ----
    _make_toggle_callbacks(
        app,
        "strand-store",
        ["strand-plus-btn", "strand-both-btn", "strand-minus-btn"],
        ["+", "both", "-"],
        "both",
        reset_ranges=False,
    )

    # ---- Mode (focus/overlap) toggle ----
    _make_two_button_toggle(
        app, "mode-store", "mode-focus-btn", "mode-overlap-btn", "focus", "overlap"
    )

    @app.callback(
        Output("slider-focus-wrap", "style"),
        Output("slider-overlap-wrap", "style"),
        Input("mode-store", "data"),
    )
    def style_mode_sliders(mode):
        is_focus = mode == "focus"
        base = {"minWidth": UI_SLIDER_MIN_WIDTH, "flexShrink": 1}
        return base if is_focus else {**base, "display": "none"}, (
            base if not is_focus else {**base, "display": "none"}
        )

    # ---- Color-by toggle (non-enum mode) ----
    _make_toggle_callbacks(
        app,
        "color-by-store",
        ["colorby-dist-btn", "colorby-pval-btn", "colorby-qval-btn"],
        ["dist", "pval", "qval"],
        "dist",
        reset_ranges=False,
    )

    # ---- Significance metric toggle (p vs q) ----
    _make_toggle_callbacks(
        app,
        "sig-metric-store",
        ["sig-p-btn", "sig-q-btn"],
        ["p", "q"],
        "q" if has_qvalue else "p",
        reset_ranges=False,
    )

    @app.callback(
        Output("sig-filter-label", "children"),
        Input("sig-metric-store", "data"),
    )
    def update_sig_label(metric):
        return "q \u2264" if metric == "q" else "p \u2264"

    # ---- Fold toggle (<1 / both / >1) ----
    _make_toggle_callbacks(
        app,
        "fold-store",
        ["fold-lt-btn", "fold-both-btn", "fold-gt-btn"],
        ["<1", "both", ">1"],
        "both",
        reset_ranges=False,
    )

    # ---- Colorscale flip toggle ----
    @app.callback(
        Output("colorscale-flip-store", "data"),
        Input("flip-colorscale-btn", "n_clicks"),
        State("colorscale-flip-store", "data"),
    )
    def toggle_flip(n_clicks, current):
        if n_clicks and n_clicks > 0:
            return not current
        return current

    @app.callback(
        Output("flip-colorscale-btn", "style"),
        Input("colorscale-flip-store", "data"),
    )
    def style_flip(flip):
        return {
            **toggle_style(flip, "right"),
            "marginLeft": 4,
            "padding": "0 6px",
            "fontSize": "12px",
        }

    # ---- Shared helpers for callbacks ----
    @lru_cache(maxsize=CACHE_SIZE)
    def _filtered_tree_data(seq_id):
        """Cached per-query tree pruning for 'filtered' prune mode."""
        query_ret = get_query_retained_leaves(app.df, seq_id)
        if not query_ret:
            return None
        qt = prune_tree(app.full_tree, query_ret)
        qd = tree_data_for(qt, seq_id)
        qdf = add_tip_order(app.df, qd["tip_order"])
        return qd, make_cached_filter(qdf, enum_ui)

    def _resolve_prune(prune_mode, seq_id):
        """Return (cur_data, do_filter) for the given prune mode."""
        if prune_mode == "full":
            return full_data, filter_full
        elif prune_mode == "pruned":
            cd = pruned_data if has_pruned else full_data
            flt = filter_pruned if has_pruned else filter_full
            return cd, flt
        elif prune_mode == "filtered":
            result = _filtered_tree_data(seq_id)
            if result is not None:
                return result
            return full_data, filter_full
        return full_data, filter_full

    def _resolve_filter(
        do_filter,
        seq_id,
        tip_order,
        strand_val=None,
        direction_val=None,
        mode=None,
        focus_val=None,
        overlap_val=None,
        sig_log=None,
        sig_metric="p",
        fold_val=None,
    ):
        """Return filtered dataframe for current filter state."""
        if enum_ui:
            if mode == "focus":
                d = nearest_value(
                    focus_val if focus_val is not None else dist_ths[0], dist_ths
                )
                d_lo = d_hi = d
            else:
                if isinstance(overlap_val, (list, tuple)) and len(overlap_val) == 2:
                    d_lo = nearest_value(overlap_val[0], dist_ths)
                    d_hi = nearest_value(overlap_val[1], dist_ths)
                else:
                    d_lo, d_hi = dist_ths[0], dist_ths[-1]
            direction = (
                direction_val
                if has_is_rc and direction_val and direction_val != "both"
                else None
            )
            df_f = do_filter(seq_id, tuple(tip_order), d_hi, d_lo, direction)
            return df_f
        sig_col = "QVALUE" if sig_metric == "q" else "PERCENTILE"
        strand = strand_val if has_strand and strand_val and strand_val != "both" else None
        df_f = do_filter(
            seq_id,
            tuple(tip_order),
            strand,
            sig_log,
            sig_col,
            fold_val if fold_val and fold_val != "both" else None,
        )
        return df_f

    # ---- Main figure update ----
    @app.callback(
        Output("graph", "figure"),
        Output("interval-count-store", "data"),
        Input("query-dropdown", "value"),
        Input("dist-slider-focus", "value"),
        Input("dist-slider-overlap", "value"),
        Input("mode-store", "data"),
        Input("colorscheme-dropdown", "value"),
        Input("y-range-store", "data"),
        Input("x-range-store", "data"),
        Input("tree-view-store", "data"),
        Input("direction-store", "data"),
        Input("strand-store", "data"),
        Input("color-by-store", "data"),
        Input("sig-metric-store", "data"),
        Input("pval-slider", "value"),
        Input("fold-store", "data"),
        Input("colorscale-flip-store", "data"),
        Input("prune-store", "data") if has_pruned else State("prune-store", "data"),
    )
    def update_figure(
        seq_id,
        focus_val,
        overlap_val,
        mode,
        scheme,
        y_range,
        x_range,
        tree_mode,
        direction,
        strand,
        color_by,
        sig_metric,
        sig_log,
        fold_val,
        flip_colorscale,
        prune_mode,
    ):
        if seq_id is None:
            return go.Figure(), 0

        cur_data, do_filter = _resolve_prune(prune_mode, seq_id)
        tip_order = cur_data["tip_order"]
        leaf_dist = cur_data["leaf_distances"]
        strand_val = strand if (not enum_ui and has_strand and strand != "both") else None
        direction_val = direction if (enum_ui and has_is_rc) else None

        df_filt = _resolve_filter(
            do_filter, seq_id, tip_order,
            strand_val=strand_val,
            direction_val=direction_val,
            mode=mode, focus_val=focus_val, overlap_val=overlap_val,
            sig_log=sig_log, sig_metric=sig_metric, fold_val=fold_val,
        )
        bin_edges, colors = (
            get_binned_colors(dist_ths, scheme) if enum_ui else (None, None)
        )

        # Sequence length
        seq_len = get_seq_len(df_filt)
        if not seq_len or seq_len <= 0:
            if enum_ui:
                seq_len = get_seq_len(
                    do_filter(seq_id, tuple(tip_order), dist_ths[-1], dist_ths[0], None)
                )
            else:
                seq_len = get_seq_len(
                    do_filter(seq_id, tuple(tip_order), None)
                )
        if not seq_len or seq_len <= 0:
            return go.Figure(), 0

        n_tips = len(tip_order)

        # Enforce zoom constraints
        if y_range is not None:
            y_range = enforce_min_span(y_range, ZOOM_MIN_Y_SPAN, bounds=(0, n_tips - 1))
        if x_range is not None:
            x_range = enforce_min_span(
                x_range, min(ZOOM_MIN_X_SPAN, seq_len), bounds=(0, seq_len)
            )

        x_limits = [-AXIS_X_PAD, seq_len + AXIS_X_PAD]
        y_limits = [-AXIS_Y_PAD, n_tips - 1 + AXIS_Y_PAD]

        layout = (
            cur_data["clado_layout"]
            if tree_mode == "cladogram"
            else cur_data["phylo_layout"]
        )

        n_intervals = int(
            (df_filt["INTERVAL_START"] != df_filt["INTERVAL_END"]).sum()
        ) if not df_filt.empty else 0

        fig = build_figure(
            layout[0],
            tip_order,
            layout[2],
            df_filt,
            bin_edges=bin_edges,
            colors=colors,
            df_annotations=df_annot,
            query_id=seq_id,
            interval_lw=scaled_interval_width(y_range, n_tips),
            tree_lw=scaled_tree_width(y_range, n_tips),
            uirevision="query",
            y_range=y_range,
            x_range=x_range,
            x_limits=x_limits,
            y_limits=y_limits,
            leaf_distances=leaf_dist,
            enum_only=enum_ui,
            color_by=color_by,
            dist_range=global_dist_range,
            scheme=scheme,
            flip_colorscale=flip_colorscale,
        )
        return fig, n_intervals

    # ---- Data export ----
    @app.callback(
        Output("download-data", "data"),
        Input("export-btn", "n_clicks"),
        State("query-dropdown", "value"),
        State("dist-slider-focus", "value"),
        State("dist-slider-overlap", "value"),
        State("mode-store", "data"),
        State("y-range-store", "data"),
        State("x-range-store", "data"),
        State("prune-store", "data"),
        State("direction-store", "data"),
        State("strand-store", "data"),
        State("sig-metric-store", "data"),
        State("pval-slider", "value"),
        State("fold-store", "data"),
        prevent_initial_call=True,
    )
    def export_data(
        n,
        seq_id,
        focus_val,
        overlap_val,
        mode,
        y_rng,
        x_rng,
        prune_mode,
        direction,
        strand,
        sig_metric,
        sig_log,
        fold_val,
    ):
        if seq_id is None:
            return no_update

        cur_data, do_filter = _resolve_prune(prune_mode, seq_id)
        strand_val = strand if (not enum_ui and has_strand and strand != "both") else None
        direction_val = direction if (enum_ui and has_is_rc) else None
        df_exp = _resolve_filter(
            do_filter,
            seq_id,
            cur_data["tip_order"],
            strand_val=strand_val,
            direction_val=direction_val,
            mode=mode,
            focus_val=focus_val,
            overlap_val=overlap_val,
            sig_log=sig_log,
            sig_metric=sig_metric,
            fold_val=fold_val,
        )

        # Visible genomes
        if y_rng is None:
            y_min, y_max = 0, len(cur_data["tip_order"]) - 1
        else:
            y_min = max(0, int(np.floor(min(y_rng))))
            y_max = min(len(cur_data["tip_order"]) - 1, int(np.ceil(max(y_rng))))

        visible = cur_data["tip_order"][y_min : y_max + 1]
        df_exp = df_exp[df_exp["REF_ID"].isin(visible)].copy()

        # Visible x-range
        if x_rng is None:
            if enum_ui:
                seq_len = get_seq_len(df_exp) or get_seq_len(
                    do_filter(
                        seq_id,
                        tuple(cur_data["tip_order"]),
                        dist_ths[-1],
                        dist_ths[0],
                        None,
                    )
                )
            else:
                seq_len = get_seq_len(df_exp) or get_seq_len(
                    do_filter(seq_id, tuple(cur_data["tip_order"]), None)
                )
            x_min, x_max = 0, seq_len or 0
        else:
            x_min = int(np.floor(min(x_rng)))
            x_max = int(np.ceil(max(x_rng)))

        df_exp = df_exp[
            (df_exp["INTERVAL_END"] >= x_min) & (df_exp["INTERVAL_START"] <= x_max)
        ].copy()

        ref_to_y = {n: i for i, n in enumerate(cur_data["tip_order"])}
        df_exp["y_order"] = df_exp["REF_ID"].map(ref_to_y)
        df_exp = df_exp.sort_values(["y_order", "INTERVAL_START"]).drop(
            columns=["y_order", "y"]
        )

        tag = ""
        if enum_ui and has_is_rc and direction and direction != "both":
            tag = f"_{direction}"
        elif not enum_ui and has_strand and strand and strand != "both":
            tag = f"_{strand}"
        filename = f"{seq_id}{tag}_y{y_min}-{y_max}_x{x_min}-{x_max}.tsv"
        return dict(content=df_exp.to_csv(sep="\t", index=False), filename=filename)

    # ---- PDF export ----
    @app.callback(
        Output("download-plot", "data"),
        Input("export-pdf-btn", "n_clicks"),
        State("graph", "figure"),
        State("query-dropdown", "value"),
        State("export-width", "value"),
        State("export-height", "value"),
        prevent_initial_call=True,
    )
    def export_pdf(n, figure, seq_id, width, height):
        if figure is None:
            return no_update
        w, h = int(width) if width else 1600, int(height) if height else 900
        fig = go.Figure(figure)
        # Convert scattergl to scatter for PDF export
        new_traces = []
        for trace in fig.data:
            td = trace.to_plotly_json()
            if td.get("type") == "scattergl":
                td["type"] = "scatter"
            new_traces.append(td)
        export_fig = go.Figure(data=new_traces, layout=fig.layout)
        img_bytes = pio.to_image(export_fig, format="pdf", width=w, height=h, scale=2)
        name = seq_id or "plot"
        return dict(
            content=base64.b64encode(img_bytes).decode(),
            filename=f"{name}.pdf",
            base64=True,
        )

    return app


# =============================================================================
# ENTRY POINT
# =============================================================================


def parse_args():
    try:
        version = importlib.metadata.version("gdiff-plot")
    except importlib.metadata.PackageNotFoundError:
        version = "dev"

    p = argparse.ArgumentParser(
        prog="gdiff-plot",
        description="gdiff-plot — Interactive visualization of distance-based genome comparison results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--input",
        "-i",
        required=True,
        help="Path to TSV data file (QUERY_ID, REF_ID, INTERVAL_START, …)",
    )
    p.add_argument("--tree", "-t", required=True, help="Path to Newick tree file")
    p.add_argument(
        "--query",
        "-q",
        default=None,
        help="Default query leaf name (first sequence if omitted)",
    )
    p.add_argument(
        "--annotation",
        "-a",
        default=None,
        help="Annotation file: GFF3, GTF, or custom TSV",
    )
    p.add_argument(
        "--enum-only",
        action="store_true",
        default=False,
        help="Force enum UI (legacy DIST_TH or 15-column enum lite). "
        "Continuous/significance UI is auto-detected otherwise.",
    )
    p.add_argument("--port", "-p", type=int, default=8080, help="Port to serve on")
    p.add_argument(
        "--host",
        default="127.0.0.1",
        help="Host to bind to (use 0.0.0.0 to expose on LAN)",
    )
    p.add_argument(
        "--open", action="store_true", help="Open browser automatically after startup"
    )
    p.add_argument(
        "--debug", action="store_true", help="Enable Dash debug/hot-reload mode"
    )
    p.add_argument("--version", action="version", version=f"%(prog)s {version}")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    app = create_app(
        args.input,
        args.tree,
        query=args.query,
        annotation_path=args.annotation,
        enum_only=True if args.enum_only else None,
    )

    url = f"http://{args.host}:{args.port}"
    if args.open:
        import threading

        threading.Timer(1.2, lambda: webbrowser.open(url)).start()

    app.run(debug=args.debug, host=args.host, port=args.port)
