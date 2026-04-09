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


def get_distance_thresholds(df: pd.DataFrame) -> list[float]:
    """Return sorted unique DIST_TH values.

    Args:
        df: Input dataframe with DIST_TH column.

    Returns:
        Sorted list of unique distance thresholds.
    """
    return sorted(df["DIST_TH"].unique())


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
    strand: Optional[str] = None,
) -> pd.DataFrame:
    """Filter intervals by query, distance threshold, and strand."""
    df_q = df[df["QUERY_ID"] == seq_id].copy()

    if strand and strand != "both" and "STRAND" in df_q.columns:
        df_q = df_q[df_q["STRAND"] == strand]

    # Aggregate to unique intervals
    df_q = df_q.groupby(
        ["REF_ID", "y", "QUERY_ID", "INTERVAL_START", "INTERVAL_END", "SEQ_LEN"],
        as_index=False,
    )["DIST_TH"].min()

    df_q = df_q[df_q["DIST_TH"] <= dist_hi]
    if dist_lo is not None:
        df_q = df_q[df_q["DIST_TH"] >= dist_lo]

    return df_q.dropna(subset=["y"]).sort_values("y")


def filter_intervals_continuous(
    df: pd.DataFrame,
    seq_id: str,
    tip_order: list[str],
    strand: Optional[str] = None,
    pval_th: Optional[float] = None,
    sign: Optional[str] = None,
) -> pd.DataFrame:
    """Filter intervals by query, strand, p-value cutoff, and sign.

    strand: '+', '-', or 'both'. When 'both', pick strand with lowest DIST_CONTIG per interval.
    pval_th: if set, keep only rows with PERCENTILE <= pval_th.
    sign: "<" or ">" to filter on the SIGN column; None means no filter.
    """
    df_q = df[df["QUERY_ID"] == seq_id].copy()

    if pval_th is not None and "PERCENTILE" in df_q.columns:
        df_q = df_q[df_q["PERCENTILE"] <= pval_th]

    if sign is not None and "SIGN" in df_q.columns:
        df_q = df_q[df_q["SIGN"] == sign]

    if strand and strand != "both" and "STRAND" in df_q.columns:
        df_q = df_q[df_q["STRAND"] == strand]
    elif (
        strand == "both" and "STRAND" in df_q.columns and "DIST_CONTIG" in df_q.columns
    ):
        # Pick strand with lowest DIST_CONTIG for each (REF_ID, start, end)
        df_q = df_q.loc[
            df_q.groupby(["REF_ID", "INTERVAL_START", "INTERVAL_END"])[
                "DIST_CONTIG"
            ].idxmin()
        ].copy()
    elif strand == "both" and "STRAND" in df_q.columns and "DIST" in df_q.columns:
        # Fallback to DIST if DIST_CONTIG not available
        df_q = df_q.loc[
            df_q.groupby(["REF_ID", "INTERVAL_START", "INTERVAL_END"])["DIST"].idxmin()
        ].copy()

    return df_q.dropna(subset=["y"]).sort_values("y")


def make_cached_filter(df: pd.DataFrame, enum_only=True):
    if enum_only:

        @lru_cache(maxsize=CACHE_SIZE)
        def cached(seq_id, tip_order_tuple, dist_hi, dist_lo, strand):
            return filter_intervals(
                df, seq_id, list(tip_order_tuple), dist_hi, dist_lo, strand
            )

        return cached
    else:

        @lru_cache(maxsize=CACHE_SIZE)
        def cached_cont(seq_id, tip_order_tuple, strand, pval_th=None, sign=None):
            return filter_intervals_continuous(
                df, seq_id, list(tip_order_tuple), strand, pval_th, sign
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
        hi = center + min_span / 1.9999  # Slightly offset to handle rounding

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


def batch_by_color(df, bin_edges, colors, leaf_distances=None, flip=False):
    """Group intervals by color bin for efficient trace creation."""
    if df.empty:
        return {}
    mask = df["INTERVAL_START"].values != df["INTERVAL_END"].values
    if not mask.any():
        return {}

    # Reverse colors if flipped
    if flip:
        colors = list(reversed(colors))

    starts = df["INTERVAL_START"].values[mask]
    ends = df["INTERVAL_END"].values[mask]
    ys = df["y"].values[mask]
    dists = df["DIST_TH"].values[mask]
    refs = df["REF_ID"].values[mask]
    cidx = assign_color_indices(dists, bin_edges, len(colors))

    hovers = [
        f"{rid}<br>Pos: {s:,}-{e:,}<br>Dist: {d:.3f}"
        + (
            f"<br>Distance to query: {leaf_distances[rid]:.4f}"
            if leaf_distances and rid in leaf_distances
            else ""
        )
        for rid, s, e, d in zip(refs, starts, ends, dists)
    ]

    traces = {}
    for ci in np.unique(cidx)[::-1]:
        sel = cidx == ci
        n = sel.sum()
        x = np.empty(n * 3, dtype=object)
        x[0::3], x[1::3], x[2::3] = starts[sel], ends[sel], None

        y = np.empty(n * 3, dtype=object)
        y[0::3], y[1::3], y[2::3] = ys[sel], ys[sel], None

        sh = np.array(hovers, dtype=object)[sel]
        text = np.empty(n * 3, dtype=object)
        text[0::3], text[1::3], text[2::3] = sh, sh, None

        traces[colors[ci]] = {"x": x.tolist(), "y": y.tolist(), "text": text.tolist()}
    return traces


# Bins for p-value coloring: edges go from 1 down to 0; first bin (high p) is whitish
PVAL_BINS = (1.0, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001, 0.0)


def _pval_bin_colors(scheme):
    """Build colors for p-value bins: first bin whitish, rest from colorscale."""
    n = len(PVAL_BINS) - 1  # number of bins
    # Sample n-1 colors from the scheme for bins 1..n-1 (low p-values)
    fracs = [i / max(n - 2, 1) for i in range(n - 1)]
    sampled = px.colors.sample_colorscale(scheme, fracs)
    # First bin (p near 1) is whitish/light gray
    return ["rgba(230,230,230,0.6)"] + list(sampled)


def _pval_bin_index(pvals):
    """Assign p-values to bins. Bin 0 = [PVAL_BINS[1], PVAL_BINS[0]], etc."""
    # searchsorted on reversed edges (descending)
    edges = np.array(PVAL_BINS)  # descending
    idx = np.searchsorted(-edges, -pvals, side="left") - 1
    return np.clip(idx, 0, len(PVAL_BINS) - 2)


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

    # Select color column
    if color_by == "pval":
        raw_colors = pvals
        cmin, cmax = 0.0, 1.0
    else:  # "dist"
        raw_colors = dists
        cmin, cmax = dist_range

    return {
        "starts": starts,
        "ends": ends,
        "ys": ys,
        "refs": refs,
        "dists": dists,
        "pvals": pvals,
        "dist_contigs": dist_contigs,
        "dist_genomes": dist_genomes,
        "raw_colors": raw_colors,
        "cmin": cmin,
        "cmax": cmax,
    }


def batch_by_continuous_color(
    df, color_by, dist_range, scheme, leaf_distances=None, flip=False
):
    """Create traces with continuous or binned coloring for non-enum mode.

    For 'pval' color_by: uses fixed p-value bins with whitish first bin.
    For 'dist' color_by: 64-bin continuous discretization.
    """
    cv = compute_color_values(df, color_by, dist_range)
    if cv is None:
        return {}, 0.0, 1.0

    if color_by == "pval":
        # Binned p-value coloring
        bin_idx = _pval_bin_index(cv["pvals"])
        bin_colors = _pval_bin_colors(scheme)
        if flip:
            # Reverse: whitish goes to low p-values (high significance)
            bin_colors = list(reversed(bin_colors))
    else:
        # Continuous distance coloring
        n_bins = 64
        cmin, cmax = cv["cmin"], cv["cmax"]
        span = cmax - cmin if cmax > cmin else 1.0
        fracs = np.clip((cv["raw_colors"] - cmin) / span, 0.0, 1.0)
        bin_idx = np.clip((fracs * (n_bins - 1)).astype(int), 0, n_bins - 1)
        sample_fracs = [i / max(n_bins - 1, 1) for i in range(n_bins)]
        bin_colors = px.colors.sample_colorscale(scheme, sample_fracs)
        if flip:
            bin_colors = list(reversed(bin_colors))

    # Build hover text: dist + p-value + DIST_CONTIG + DIST_GENOME
    hovers = [
        f"{rid}<br>Pos: {s:,}-{e:,}<br>Dist: {d:.4f}<br>p-value: {p:.3e}"
        f"<br>Dist contig: {dc:.4f}<br>Dist genome: {dg:.4f}"
        + (
            f"<br>Tree dist: {leaf_distances[rid]:.4f}"
            if leaf_distances and rid in leaf_distances
            else ""
        )
        for rid, s, e, d, p, dc, dg in zip(
            cv["refs"],
            cv["starts"],
            cv["ends"],
            cv["dists"],
            cv["pvals"],
            cv["dist_contigs"],
            cv["dist_genomes"],
        )
    ]

    traces = {}
    for ci in np.unique(bin_idx)[::-1]:
        sel = bin_idx == ci
        n = sel.sum()
        x = np.empty(n * 3, dtype=object)
        x[0::3], x[1::3], x[2::3] = cv["starts"][sel], cv["ends"][sel], None
        y = np.empty(n * 3, dtype=object)
        y[0::3], y[1::3], y[2::3] = cv["ys"][sel], cv["ys"][sel], None
        sh = np.array(hovers, dtype=object)[sel]
        text = np.empty(n * 3, dtype=object)
        text[0::3], text[1::3], text[2::3] = sh, sh, None
        traces[bin_colors[ci]] = {
            "x": x.tolist(),
            "y": y.tolist(),
            "text": text.tolist(),
        }

    return traces, cv["cmin"], cv["cmax"]


def make_continuous_colorbar(
    cmin, cmax, color_by, scheme, y_center=0.5, length=None, flip=False
):
    """Create a colorbar trace for non-enum mode."""
    if color_by == "pval":
        # Binned p-value colorbar matching PVAL_BINS
        bins = PVAL_BINS
        n = len(bins) - 1
        colors = _pval_bin_colors(scheme)
        if flip:
            # Reverse colors: whitish goes to low p-values
            colors = list(reversed(colors))
        edges = np.linspace(0, 1, n + 1)
        scale = [[edges[i + j], colors[i]] for i in range(n) for j in (0, 1)]
        # Tick labels at bin boundaries
        tickvals = [edges[i] for i in range(n + 1)]
        ticktext = [f"{v:g}" for v in bins]
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
                    title="p-value",
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
                cmin=cmin,
                cmax=cmax,
                color=[cmin],
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


class ColorMapper:
    """Encapsulates colorscale mapping logic with optional flip support."""

    def __init__(self, scheme: str, flip: bool = False):
        self.scheme = scheme
        self.flip = flip
        self._binned_cache: dict[tuple, tuple] = {}
        self._colorscale_cache: dict[str, list] = {}

    def get_binned_colors(self, thresholds: tuple) -> tuple:
        """Return (bin_edges, colors) for discrete coloring with caching."""
        key = (thresholds, self.scheme, self.flip)
        if key not in self._binned_cache:
            n_bins = len(thresholds)
            frac = [i / max(n_bins - 1, 1) for i in range(n_bins)]
            colors = list(px.colors.sample_colorscale(self.scheme, frac))
            if self.flip:
                colors = list(reversed(colors))
            edges = (0.0,) + thresholds
            self._binned_cache[key] = (edges, tuple(colors))
        return self._binned_cache[key]

    def get_continuous_colors(self, n_bins: int = 64) -> list:
        """Return sampled colors for continuous coloring."""
        key = (self.scheme, n_bins, self.flip)
        if key not in self._colorscale_cache:
            fracs = [i / max(n_bins - 1, 1) for i in range(n_bins)]
            colors = px.colors.sample_colorscale(self.scheme, fracs)
            if self.flip:
                colors = list(reversed(colors))
            self._colorscale_cache[key] = colors
        return self._colorscale_cache[key]

    def get_plotly_colorscale(self) -> list:
        """Return full Plotly colorscale (for colorbars), optionally flipped."""
        if not self.flip:
            return self.scheme
        orig = px.colors.get_colorscale(self.scheme)
        return [[1 - pos, col] for pos, col in reversed(orig)]

    def reverse_colors(self, colors: list) -> list:
        """Reverse a color list if flip is enabled."""
        return list(reversed(colors)) if self.flip else colors


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
    seq_ids, dist_ths, has_pruned, has_strand=False, initial_prune=True, enum_only=True
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
            dcc.Store(id="mode-store", data="focus"),
            dcc.Store(id="color-by-store", data="dist"),
            dcc.Store(id="sign-store", data="both"),
            dcc.Store(id="colorscale-flip-store", data=False),
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
                    # Color-by toggle: dist / p-value (non-enum only)
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
                                style=toggle_style(False, "right"),
                            ),
                        ],
                        style={
                            "display": "flex" if not enum_only else "none",
                            "alignItems": "center",
                            "gap": 0,
                        },
                    ),
                    # Sign toggle: < / both / > (non-enum only)
                    html.Div(
                        [
                            divider(),
                            control_label("Sign:"),
                            html.Button(
                                "<",
                                id="sign-lt-btn",
                                n_clicks=0,
                                style=toggle_style(False, "left"),
                            ),
                            html.Button(
                                "both",
                                id="sign-both-btn",
                                n_clicks=0,
                                style=toggle_style(True, "middle"),
                            ),
                            html.Button(
                                ">",
                                id="sign-gt-btn",
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
                    # p-value filter slider (non-enum only)
                    html.Div(
                        [
                            divider(),
                            control_label("p ≤"),
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
                    # Strand toggle — only visible when strand data is present;
                    # leading divider keeps separation correct whether visible or not
                    html.Div(
                        [
                            divider(),
                            control_label("Strand:"),
                            html.Button(
                                "+ (fw)",
                                id="strand-fwd-btn",
                                n_clicks=0,
                                style=toggle_style(False, "left"),
                            ),
                            html.Button(
                                "± (both)",
                                id="strand-both-btn",
                                n_clicks=0,
                                style=toggle_style(True, "middle"),
                            ),
                            html.Button(
                                "− (rc)",
                                id="strand-rev-btn",
                                n_clicks=0,
                                style=toggle_style(False, "right"),
                            ),
                        ],
                        style={
                            "display": "flex" if has_strand else "none",
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


def extract_ranges(relayout):
    """Extract y and x ranges from relayoutData."""
    if not relayout:
        return no_update, no_update

    if any(k.endswith(".autorange") for k in relayout):
        return None, None

    # Y range
    if "yaxis.range[0]" in relayout and "yaxis.range[1]" in relayout:
        y = [relayout["yaxis.range[0]"], relayout["yaxis.range[1]"]]
    elif "yaxis.range" in relayout:
        y = relayout["yaxis.range"]
    else:
        y = no_update

    # X range (check both panel positions)
    for key in ("xaxis2", "xaxis4"):
        if f"{key}.range[0]" in relayout and f"{key}.range[1]" in relayout:
            return y, [relayout[f"{key}.range[0]"], relayout[f"{key}.range[1]"]]
        if f"{key}.range" in relayout:
            return y, relayout[f"{key}.range"]

    return y, no_update


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
    enum_only: bool = False,
) -> Dash:
    """Create and configure the Dash application."""
    df = pd.read_csv(input_path, sep="\t")
    df_annot = load_annotations(annotation_path)

    if enum_only:
        required = {
            "QUERY_ID",
            "REF_ID",
            "INTERVAL_START",
            "INTERVAL_END",
            "SEQ_LEN",
            "DIST_TH",
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

    if enum_only:
        dist_ths = tuple(get_distance_thresholds(df))
        if not dist_ths:
            raise ValueError("No distance thresholds in data.")
        # Adaptive distance range not needed for enum mode
        global_dist_range = (0.0, 0.5)
    else:
        dist_ths = ()  # Not used in continuous mode
        # Compute adaptive distance range from data
        d_min = float(df["DIST"].min()) if not df["DIST"].isna().all() else 0.0
        d_max = float(df["DIST"].max()) if not df["DIST"].isna().all() else 0.5
        if d_max <= d_min:
            d_max = d_min + 0.01
        global_dist_range = (d_min, d_max)

    seq_ids = get_sequence_identifiers(df)
    has_strand = "STRAND" in df.columns

    filter_full = make_cached_filter(full_df, enum_only)
    filter_pruned = (
        make_cached_filter(pruned_df, enum_only) if has_pruned else filter_full
    )
    filter_query = (
        make_cached_filter(query_df, enum_only) if has_query_pruned else filter_full
    )

    app = Dash(__name__)

    # Store the full tree for dynamic query pruning
    app.full_tree = full_tree
    app.df = df
    app.enum_only = enum_only
    app.layout = build_layout(
        seq_ids, dist_ths, has_pruned, has_strand=has_strand, enum_only=enum_only
    )

    # ---- Navigation callbacks ----
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

    @app.callback(Output("query-title", "children"), Input("query-dropdown", "value"))
    def update_title(seq_id):
        if seq_id is None:
            return ""
        idx = seq_ids.index(seq_id) if seq_id in seq_ids else 0
        return f"{seq_id}  ({idx + 1}/{len(seq_ids)})"

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

    # ---- Strand toggle ----
    _make_toggle_callbacks(
        app,
        "strand-store",
        ["strand-fwd-btn", "strand-both-btn", "strand-rev-btn"],
        ["+", "both", "-"],
        "both",
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
    _make_two_button_toggle(
        app, "color-by-store", "colorby-dist-btn", "colorby-pval-btn", "dist", "pval"
    )

    # ---- Sign toggle (< / both / >) ----
    _make_toggle_callbacks(
        app,
        "sign-store",
        ["sign-lt-btn", "sign-both-btn", "sign-gt-btn"],
        ["<", "both", ">"],
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

    # ---- Main figure update ----
    @app.callback(
        Output("graph", "figure"),
        Input("query-dropdown", "value"),
        Input("dist-slider-focus", "value"),
        Input("dist-slider-overlap", "value"),
        Input("mode-store", "data"),
        Input("colorscheme-dropdown", "value"),
        Input("y-range-store", "data"),
        Input("x-range-store", "data"),
        Input("tree-view-store", "data"),
        Input("strand-store", "data"),
        Input("color-by-store", "data"),
        Input("pval-slider", "value"),
        Input("sign-store", "data"),
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
        strand,
        color_by,
        pval_log,
        sign_val,
        flip_colorscale,
        prune_mode,
    ):
        if seq_id is None:
            return go.Figure()

        # Select tree data
        if prune_mode == "full":
            cur_data = full_data
            do_filter = filter_full
        elif prune_mode == "pruned":
            cur_data = pruned_data if has_pruned else full_data
            do_filter = filter_pruned if has_pruned else filter_full
        elif prune_mode == "filtered":
            # Dynamic filtered-specific pruning
            query_retained = get_query_retained_leaves(app.df, seq_id)
            if query_retained:
                query_tree = prune_tree(app.full_tree, query_retained)
                query_data_dynamic = tree_data_for(query_tree, seq_id)
                cur_data = query_data_dynamic
                query_df_dynamic = add_tip_order(
                    app.df, query_data_dynamic["tip_order"]
                )
                do_filter = make_cached_filter(query_df_dynamic, enum_only)
            else:
                cur_data = full_data
                do_filter = filter_full
        else:
            cur_data = full_data
            do_filter = filter_full

        tip_order = cur_data["tip_order"]
        leaf_dist = cur_data["leaf_distances"]
        strand_val = strand if has_strand and strand != "both" else None

        if enum_only:
            # Resolve distance thresholds
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

            df_filt = do_filter(seq_id, tuple(tip_order), d_hi, d_lo, strand_val)
            bin_edges, colors = get_binned_colors(dist_ths, scheme)
        else:
            # Convert log slider value to p-value cutoff (slider=0 → p≤1 → no filter)
            pval_th = 10.0**pval_log if pval_log is not None and pval_log < 0 else None
            sign_filter = sign_val if sign_val and sign_val != "both" else None
            df_filt = do_filter(
                seq_id, tuple(tip_order), strand_val, pval_th, sign_filter
            )
            bin_edges, colors = None, None

        # Sequence length
        seq_len = get_seq_len(df_filt)
        if not seq_len or seq_len <= 0:
            if enum_only:
                seq_len = get_seq_len(
                    do_filter(seq_id, tuple(tip_order), dist_ths[-1], dist_ths[0], None)
                )
            else:
                seq_len = get_seq_len(
                    do_filter(seq_id, tuple(tip_order), None, None, None)
                )
        if not seq_len or seq_len <= 0:
            return go.Figure()

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

        return build_figure(
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
            enum_only=enum_only,
            color_by=color_by,
            dist_range=global_dist_range,
            scheme=scheme,
            flip_colorscale=flip_colorscale,
        )

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
        State("strand-store", "data"),
        State("pval-slider", "value"),
        State("sign-store", "data"),
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
        strand,
        pval_log,
        sign_val,
    ):
        if seq_id is None:
            return no_update

        if prune_mode == "full":
            cur_data = full_data
            do_filter = filter_full
        elif prune_mode == "pruned":
            cur_data = pruned_data if has_pruned else full_data
            do_filter = filter_pruned if has_pruned else filter_full
        elif prune_mode == "filtered":
            query_retained = get_query_retained_leaves(app.df, seq_id)
            if query_retained:
                query_tree = prune_tree(app.full_tree, query_retained)
                query_data_dynamic = tree_data_for(query_tree, seq_id)
                cur_data = query_data_dynamic
                query_df_dynamic = add_tip_order(
                    app.df, query_data_dynamic["tip_order"]
                )
                do_filter = make_cached_filter(query_df_dynamic, enum_only)
            else:
                cur_data = full_data
                do_filter = filter_full
        else:
            cur_data = full_data
            do_filter = filter_full

        strand_val = strand if has_strand else None
        if enum_only:
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
            df_exp = do_filter(
                seq_id, tuple(cur_data["tip_order"]), d_hi, d_lo, strand_val
            )
        else:
            pval_th = 10.0**pval_log if pval_log is not None and pval_log < 0 else None
            sign_filter = sign_val if sign_val and sign_val != "both" else None
            df_exp = do_filter(
                seq_id, tuple(cur_data["tip_order"]), strand_val, pval_th, sign_filter
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
            if enum_only:
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
                    do_filter(seq_id, tuple(cur_data["tip_order"]), None, None, None)
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

        tag = f"_{strand}" if has_strand else ""
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
        help="Input is enum-only format (7-column TSV with DIST_TH instead of DIST/PERCENTILE)",
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
        enum_only=args.enum_only,
    )

    url = f"http://{args.host}:{args.port}"
    if args.open:
        import threading

        threading.Timer(1.2, lambda: webbrowser.open(url)).start()

    app.run(debug=args.debug, host=args.host, port=args.port)
