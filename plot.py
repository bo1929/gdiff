import argparse
import pandas as pd
import numpy as np
from ete3 import Tree
from dash import Dash, dcc, html, Input, Output, State, no_update, ctx
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots
import base64
from functools import lru_cache
from collections import defaultdict

# =============================================================================
# COLORS
# =============================================================================

COLORS = {
    "tree_line": "#2c3e50",
    "axis_line": "#34495e",
    "grid": "rgba(0,0,0,0.15)",
    "plot_bg": "#fbfbfb",
    "paper_bg": "#fafafa",
    "panel_bg": "#ecf0f1",
    "button_bg": "#ecf0f1",
    "button_bg_inactive": "#e8f4f8",
    "button_border": "#7f8c8d",
    "button_text": "#000",
}

# =============================================================================
# FIGURE
# =============================================================================

FIG_HEIGHT = 700
FIG_MARGIN = dict(l=10, r=100, t=40, b=50)
FIG_FONT_FAMILY = "Segoe UI, Arial, sans-serif"
FIG_FONT_SIZE = 13

# =============================================================================
# PANELS (subplots)
# =============================================================================

PANEL_TREE_WIDTH = 0.30
PANEL_INTERVAL_WIDTH = 0.70
PANEL_SPACING = 0.001

# =============================================================================
# TREE RENDERING
# =============================================================================

TREE_LINE_WIDTH_MIN = 0.5
TREE_LINE_WIDTH_MAX = 3.5
TREE_ROW_FILL = 0.20
TREE_X_MARGIN = 1.01
TREE_X_TITLE = "Branch length"

# =============================================================================
# INTERVAL RENDERING
# =============================================================================

INTERVAL_ROW_FILL = 0.65
INTERVAL_LINE_WIDTH_MIN = 0.01
INTERVAL_LINE_WIDTH_MAX = 25
INTERVAL_LINE_WIDTH_MAX_RATIO = 0.85
INTERVAL_X_TITLE = "Genomic Position (bp)"

# =============================================================================
# AXES
# =============================================================================

AXIS_LINE_WIDTH = 2.2
AXIS_TITLE_FONT_SIZE = 20
AXIS_TICK_FONT_SIZE = 16
AXIS_Y_MAX_TICKS = 30
AXIS_X_NTICKS = 12
AXIS_GRID_WIDTH = 2.5
AXIS_Y_PAD = 1.8  # extra space so top/bottom branches don't touch the axis lines
AXIS_X_PAD = 0.99

# =============================================================================
# COLORBAR & COLORSCALE
# =============================================================================

COLORBAR_LEN = 0.75
COLORBAR_X = 1.02
COLORBAR_THICKNESS = 24
COLORBAR_TITLE_FONT_SIZE = 20
COLORBAR_TICK_FONT_SIZE = 16
COLORSCALE_DEFAULT = "viridis"
COLORSCALE_OPTIONS = ["viridis", "plasma", "inferno", "magma", "cividis"]

# =============================================================================
# UI: CONTROL PANEL
# =============================================================================

UI_PANEL_PADDING = "6px 10px"
UI_PANEL_GAP = 10
UI_PANEL_BORDER_RADIUS = 4
UI_PANEL_MARGIN_BOTTOM = 24
UI_PANEL_BOX_SHADOW = "0 2px 4px rgba(0,0,0,0.1)"

# =============================================================================
# UI: LABELS & TITLE
# =============================================================================

UI_LABEL_FONT_SIZE = 16
UI_LABEL_MARGIN_RIGHT = 8
UI_TITLE_FONT_SIZE = 20
UI_TITLE_MARGIN_BOTTOM = 4

# =============================================================================
# UI: TOGGLE BUTTONS
# =============================================================================

UI_TOGGLE_PADDING = "6px 14px"
UI_TOGGLE_FONT_SIZE = 16
UI_TOGGLE_BORDER_WIDTH = 2
UI_TOGGLE_MARGIN_LEFT = 10

# =============================================================================
# UI: NAV BUTTONS (prev/next query)
# =============================================================================

UI_NAV_PADDING = "4px 8px"
UI_NAV_FONT_SIZE = 16
UI_NAV_BORDER_WIDTH = 1

# =============================================================================
# UI: CONTROLS (dropdowns, slider)
# =============================================================================

UI_QUERY_DROPDOWN_WIDTH = 200
UI_QUERY_MARGIN_RIGHT = 10
UI_SLIDER_MIN_WIDTH = 280
UI_SLIDER_GAP = 3
UI_COLORSCALE_DROPDOWN_WIDTH = 110
UI_COLORSCALE_GAP = 4
UI_COLORSCALE_MARGIN_LEFT = 10

# =============================================================================
# UI: GRAPH CONTAINER
# =============================================================================

UI_GRAPH_HEIGHT = "78vh"
UI_CONTAINER_HEIGHT = "80vh"
UI_APP_PADDING = 5

# =============================================================================
# ZOOM CONSTRAINTS & CACHE
# =============================================================================

ZOOM_MIN_Y_SPAN = 10
ZOOM_MIN_X_SPAN = 10000
LRU_CACHE_SIZE = 256

# =============================================================================
# ANNOTATION PANEL
# =============================================================================

PANEL_ANNOTATION_HEIGHT = 0.175  # fraction of figure height for annotation row
PANEL_TREE_INTERVAL_HEIGHT = 0.75  # fraction for tree+interval row
PANEL_VERTICAL_SPACING = 0.05  # gap between rows
PANEL_ANNOTATION_TRACK_HEIGHT = 0.20  # body half-height of gene arrows (track units ≤ 0.5)

# Arrow rendering
# Head occupies the last 25 % of each gene (body fills the first 75 %).
ARROW_HEAD_RATIO = 0.2  # fraction of gene length reserved for the arrowhead
ARROW_SIZE_MIN = 4  # minimum arrowhead size in bp (guards tiny features)

# Feature type colors — earthy/muted palette, no blue
FEATURE_COLORS = {
    "CDS": "#A0522D",  # Sienna
    "gene": "#4A7C4E",  # Forest green
    "rRNA": "#B44040",  # Brick red
    "tRNA": "#CC8800",  # Amber
    "misc_RNA": "#7B5EA7",  # Dusty purple
    "ncRNA": "#7B5EA7",  # Dusty purple
    "regulatory": "#C06020",  # Burnt orange
    "exon": "#7A8040",  # Olive
    "intron": "#909090",  # Gray
    "repeat": "#AA6688",  # Mauve
    "default": "#606060",  # Dark gray
}


# =============================================================================
# DATA HANDLING
# =============================================================================


def get_sequence_identifiers(df):
    return sorted(df["QUERY_ID"].unique())


def get_distance_thresholds(df):
    return sorted(df["DIST_TH"].unique())


def get_retained_leaves(df):
    df_with_matches = df[df["INTERVAL_START"] != df["INTERVAL_END"]]
    return set(df_with_matches["REF_ID"].unique())


def get_seq_len(df):
    return df["SEQ_LEN"].iloc[0] if not df.empty else None


def add_tip_order(df, tip_order):
    ref_to_y = {name: i for i, name in enumerate(tip_order)}
    df = df.copy()
    df["y"] = df["REF_ID"].map(ref_to_y)
    return df


def filter_data(df, seq_id, tip_order, dist_th_max=None, dist_th_min=None, strand=None):
    df_q = df[df["QUERY_ID"] == seq_id].copy()
    # Filter by strand before aggregation so we never mix + and - intervals
    if strand is not None and "STRAND" in df_q.columns:
        df_q = df_q[df_q["STRAND"] == strand]
    df_q = df_q.groupby(
        ["REF_ID", "y", "QUERY_ID", "INTERVAL_START", "INTERVAL_END", "SEQ_LEN"], as_index=False
    )["DIST_TH"].min()

    if dist_th_max is not None:
        df_q = df_q[df_q["DIST_TH"] <= dist_th_max]
    if dist_th_min is not None:
        df_q = df_q[df_q["DIST_TH"] >= dist_th_min]
    df_q = df_q.dropna(subset=["y"]).sort_values("y")
    return df_q


def make_cached_filter(df):
    @lru_cache(maxsize=LRU_CACHE_SIZE)
    def cached(seq_id, tip_order_tuple, dist_th_max, dist_th_min=None, strand=None):
        return filter_data(df, seq_id, list(tip_order_tuple), dist_th_max, dist_th_min, strand)

    return cached


# =============================================================================
# ANNOTATION HANDLING
# =============================================================================


def _parse_gff_attributes(attr_str):
    """Parse a GFF3 or GTF attribute string into a dict.

    GFF3 uses ``key=value`` pairs separated by ``;``.
    GTF uses ``key "value"`` pairs separated by ``;``.
    Handles both transparently.
    """
    attrs = {}
    if pd.isna(attr_str) or not str(attr_str).strip():
        return attrs
    for field in str(attr_str).split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            # GFF3: key=value
            key, _, val = field.partition("=")
            attrs[key.strip()] = val.strip()
        elif " " in field:
            # GTF: key "value"
            key, _, val = field.partition(" ")
            attrs[key.strip()] = val.strip().strip('"')
    return attrs


def _is_gff_format(path):
    """Heuristic: file is GFF/GTF if the first non-comment data line has >= 8
    tab-separated columns and column 4/5 look numeric (start/end)."""
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 8:
                try:
                    int(parts[3])
                    int(parts[4])
                    return True
                except ValueError:
                    pass
            return False
    return False


def _load_gff(path):
    """Load a GFF3 or GTF2.2 file and normalise it to the internal schema."""
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
            attr_str = cols[8] if len(cols) > 8 else ""
            try:
                start_i, end_i = int(start), int(end)
            except ValueError:
                continue
            attrs = _parse_gff_attributes(attr_str)
            # Derive locus_tag — prefer ID, then locus_tag, then Name, then gene_id
            locus_tag = (
                attrs.get("ID")
                or attrs.get("locus_tag")
                or attrs.get("Name")
                or attrs.get("gene_id")
                or f"{seqname}_{start_i}_{end_i}"
            )
            gene_name = attrs.get("gene") or attrs.get("gene_name") or attrs.get("Name")
            product = attrs.get("product") or attrs.get("description") or attrs.get("note")
            rows.append(
                {
                    "contig_id": seqname,
                    "locus_tag": locus_tag,
                    "ftype": feature if feature != "." else "misc",
                    "start": min(start_i, end_i),
                    "stop": max(start_i, end_i),
                    "strand": strand if strand in ("+", "-") else "+",
                    "gene_name": gene_name,
                    "product": product,
                    "ec_number": attrs.get("ec_number") or attrs.get("EC_number"),
                    "source": source,
                }
            )
    if not rows:
        return None
    return pd.DataFrame(rows)


def load_annotations(annotation_path):
    """Load annotations from GFF3, GTF, or custom TSV format.

    Detects format automatically.  Returns a DataFrame with at least
    ``contig_id, locus_tag, ftype, start, stop, strand`` or ``None``.
    """
    if annotation_path is None:
        return None

    try:
        # ── Try GFF/GTF first ────────────────────────────────────────────
        if _is_gff_format(annotation_path):
            return _load_gff(annotation_path)

        # ── Fallback: custom TSV ─────────────────────────────────────────
        df = pd.read_csv(annotation_path, sep="\t", na_values=["."], keep_default_na=True)

        # Drop a leading unnamed index column produced by tools like `cat -n`
        first_col = df.columns[0]
        if str(first_col) not in {"contig_id", "locus_tag", "ftype", "start", "stop", "strand"}:
            try:
                pd.to_numeric(df[first_col])
                df = df.drop(columns=[first_col])
            except (ValueError, TypeError):
                pass

        required_cols = {"contig_id", "locus_tag", "ftype", "start", "stop", "strand"}
        missing = required_cols - set(df.columns)
        if missing:
            raise ValueError(
                f"Annotation file missing required columns: {', '.join(sorted(missing))}"
            )

        df = df.copy()
        df["start"] = df[["start", "stop"]].min(axis=1)
        df["stop"] = df[["start", "stop"]].max(axis=1)

        # Coalesce product / EC / gene_name from multiple possible columns
        for target, candidates in [
            ("product", ["prokka_product", "swissprot_product", "product"]),
            ("ec_number", ["prokka_EC_number", "swissprot_EC_number", "ec_number"]),
            ("gene_name", ["prokka_gene", "swissprot_gene", "gene_name", "gene"]),
        ]:
            cols = [c for c in candidates if c in df.columns]
            if cols:
                df[target] = df[cols].apply(
                    lambda row: next((v for v in row if pd.notna(v)), None), axis=1
                )
            elif target not in df.columns:
                df[target] = None

        return df

    except Exception as e:
        print(f"Warning: Could not load annotation file: {e}")
        return None


def filter_annotations(df_annotations, query_id, x_range=None):
    """
    Filter annotations for a specific query and optional x-range.

    Args:
        df_annotations: Full annotation DataFrame
        query_id: Query sequence ID (contig_id)
        x_range: Optional [x_min, x_max] to filter by position

    Returns:
        Filtered annotation DataFrame
    """
    if df_annotations is None:
        return None

    # Filter by contig_id
    df_filtered = df_annotations[df_annotations["contig_id"] == query_id].copy()

    if df_filtered.empty:
        return None

    # Filter by x-range if provided
    if x_range is not None:
        x_min, x_max = sorted(x_range)
        df_filtered = df_filtered[
            (df_filtered["stop"] >= x_min) & (df_filtered["start"] <= x_max)
        ].copy()

    return df_filtered if not df_filtered.empty else None


# =============================================================================
# TREE & LAYOUT
# =============================================================================


def prune_tree(tree, retained_l):
    if retained_l is None or len(retained_l) == 0:
        return tree
    pruned_tree = tree.copy()
    leaves_to_keep = [leaf for leaf in pruned_tree.iter_leaves() if leaf.name in retained_l]
    if not leaves_to_keep:
        return tree
    pruned_tree.prune(leaves_to_keep, preserve_branch_length=True)
    return pruned_tree


def load_tree(newick_file):
    tree = Tree(newick_file, format=1)
    tree.ladderize()
    return tree


def compute_path_distances(tree, query_leaf_name):
    """Compute tree path distance from query leaf to all nodes."""
    query_leaf = None
    for leaf in tree.iter_leaves():
        if leaf.name == query_leaf_name:
            query_leaf = leaf
            break
    if query_leaf is None:
        return {}
    return {node: query_leaf.get_distance(node) for node in tree.traverse()}


def _node_hover(node, query_distances, count_offset=0):
    """Build hover text for a tree node."""
    branch_len = node.dist if node.dist is not None else 0.0
    if node.is_leaf():
        hover = f"{node.name}<br>Branch length: {branch_len:.4f}"
    else:
        name = node.name if node.name else "(internal)"
        size = len(node.get_leaves()) + count_offset
        hover = f"{name}<br>Subtree size: {size}<br>Branch length: {branch_len:.4f}"
    if query_distances and node in query_distances:
        hover += f"<br>Distance to query: {query_distances[node]:.4f}"
    return hover


def compute_tree_layout(tree, mode="phylogeny", query_distances=None):
    data = defaultdict(list)
    tip_order = [leaf.name for leaf in tree.iter_leaves()]
    y_pos = {name: i for i, name in enumerate(tip_order)}

    for leaf in tree.iter_leaves():
        leaf.add_feature("y", y_pos[leaf.name])
    for node in tree.traverse("postorder"):
        if not node.is_leaf():
            node.add_feature("y", np.mean([c.y for c in node.children]))

    # Assign x coordinates
    if mode == "cladogram":
        max_depth = 0
        for node in tree.traverse("postorder"):
            if node.is_leaf():
                node.add_feature("depth", 0)
            else:
                child_depths = [c.depth for c in node.children]
                node.add_feature("depth", max(child_depths) + 1 if child_depths else 1)
                max_depth = max(max_depth, node.depth)
        for node in tree.traverse("preorder"):
            node.add_feature("x", float(max_depth - node.depth))
        max_x = float(max_depth)
        count_offset = 1
    else:
        max_x = 0.0
        for node in tree.traverse("preorder"):
            if node.is_root():
                node.add_feature("x", 0)
            else:
                node.add_feature("x", node.up.x + (node.dist or 0))
            max_x = max(max_x, node.x)
        count_offset = 0

    max_y = len(tip_order) - 1

    # Horizontal lines
    for node in tree.traverse():
        if not node.is_root():
            data["x"].extend([node.up.x, node.x, None])
            data["y"].extend([node.y, node.y, None])
            hover = _node_hover(node, query_distances, count_offset)
            data["text"].extend([hover, hover, None])

    # Vertical lines
    for node in tree.traverse():
        if len(node.children) > 1:
            child_ys = [c.y for c in node.children]
            ix_min = int(np.argmin(child_ys))
            ix_max = int(np.argmax(child_ys))
            data["x"].extend([node.x, node.x, None])
            data["y"].extend([node.children[ix_min].y, node.children[ix_max].y, None])
            data["text"].extend(
                [
                    _node_hover(node.children[ix_min], query_distances, count_offset),
                    _node_hover(node.children[ix_max], query_distances, count_offset),
                    None,
                ]
            )

    return data, tip_order, (max_x, max_y)


# =============================================================================
# ZOOM & RANGE
# =============================================================================


def enforce_min_span(v_range, min_span, bounds=None):
    """Prevent zooming in beyond min_span; clamp to bounds if given."""
    if v_range is None:
        return None
    lo, hi = sorted(v_range)
    if hi - lo < min_span:
        center = (lo + hi) / 2
        lo = center - min_span / 2
        hi = center + min_span / 2
    if bounds is not None:
        b_lo, b_hi = bounds
        if hi - lo > b_hi - b_lo:
            lo, hi = b_lo, b_hi
        else:
            lo = max(b_lo, min(lo, b_hi - (hi - lo)))
            hi = lo + max(min_span, hi - lo)
            hi = min(hi, b_hi)
    lo, hi = int(np.floor(lo)), int(np.ceil(hi))
    return [lo, hi] if v_range[0] <= v_range[1] else [hi, lo]


def _visible_rows(y_range, total_leaves):
    if y_range is None:
        return total_leaves
    return max(1, abs(int(y_range[1] - y_range[0])))


def scaled_interval_width(y_range, total_leaves):
    visible = _visible_rows(y_range, total_leaves)
    row_h = FIG_HEIGHT / visible
    thickness = row_h * INTERVAL_ROW_FILL
    thickness = min(thickness, row_h * INTERVAL_LINE_WIDTH_MAX_RATIO)
    return min(max(INTERVAL_LINE_WIDTH_MIN, thickness), INTERVAL_LINE_WIDTH_MAX)


def scaled_tree_width(y_range, total_leaves):
    visible = _visible_rows(y_range, total_leaves)
    row_h = FIG_HEIGHT / visible
    thickness = row_h * TREE_ROW_FILL
    return min(max(TREE_LINE_WIDTH_MIN, thickness), TREE_LINE_WIDTH_MAX)


# =============================================================================
# COLOR MAPPING
# =============================================================================


@lru_cache(maxsize=LRU_CACHE_SIZE)
def get_binned_colors(dist_th_t, colorscheme):
    bin_edges = (0.0,) + dist_th_t
    n_bins = len(dist_th_t)
    colors = tuple(
        px.colors.sample_colorscale(
            colorscheme, [i / (n_bins - 1) for i in range(n_bins)] if n_bins > 1 else [0.5]
        )
    )
    return np.array(bin_edges), colors


def assign_color_indices(distances, bin_edges, n_colors):
    indices = np.searchsorted(bin_edges, distances, side="left") - 1
    return np.clip(indices, 0, n_colors - 1)


def _build_hover_texts(ref_ids, starts, ends, dist_ths, leaf_distances=None):
    """Vectorized hover text construction."""
    texts = [
        f"{rid}<br>Pos: {s:,}-{e:,}<br>Dist: {d:.3f}"
        for rid, s, e, d in zip(ref_ids, starts, ends, dist_ths)
    ]
    if leaf_distances:
        texts = [
            t + f"<br>Distance to query: {leaf_distances[rid]:.4f}" if rid in leaf_distances else t
            for t, rid in zip(texts, ref_ids)
        ]
    return texts


def batch_interval_by_color(df, bin_edges, colors, leaf_distances=None):
    if df.empty:
        return {}
    mask = df["INTERVAL_START"].values != df["INTERVAL_END"].values
    if not mask.any():
        return {}

    starts = df["INTERVAL_START"].values[mask]
    ends = df["INTERVAL_END"].values[mask]
    ys = df["y"].values[mask]
    dist_ths = df["DIST_TH"].values[mask]
    ref_ids = df["REF_ID"].values[mask]
    color_idx = assign_color_indices(dist_ths, bin_edges, len(colors))

    hovers = _build_hover_texts(ref_ids, starts, ends, dist_ths, leaf_distances)

    traces = {}
    for ci in np.unique(color_idx)[::-1]:
        sel = color_idx == ci
        n = sel.sum()
        # Build interleaved [start, end, None, ...] arrays
        x = np.empty(n * 3, dtype=object)
        x[0::3] = starts[sel]
        x[1::3] = ends[sel]
        x[2::3] = None
        y = np.empty(n * 3, dtype=object)
        y[0::3] = ys[sel]
        y[1::3] = ys[sel]
        y[2::3] = None
        sel_hovers = np.array(hovers, dtype=object)[sel]
        text = np.empty(n * 3, dtype=object)
        text[0::3] = sel_hovers
        text[1::3] = sel_hovers
        text[2::3] = None
        traces[colors[ci]] = {"x": x.tolist(), "y": y.tolist(), "text": text.tolist()}
    return traces


def make_colorbar(bin_edges, colors, colorbar_y=0.5, colorbar_len=None):
    n_bins = len(colors)
    norm_edges = np.linspace(0, 1, n_bins + 1)
    colorscale = [[norm_edges[i + j], colors[i]] for i in range(n_bins) for j in (0, 1)]

    tick_labels = [f"≤{bin_edges[i+1]:.3f}" for i in range(n_bins)]
    dummy_y = np.linspace(0, 1, n_bins)
    cb_len = colorbar_len if colorbar_len is not None else COLORBAR_LEN

    return go.Scatter(
        x=[None],
        y=[None],
        mode="markers",
        marker=dict(
            colorscale=colorscale,
            showscale=True,
            cmin=0,
            cmax=1,
            color=dummy_y,
            colorbar=dict(
                title="Distance",
                tickvals=[(norm_edges[i] + norm_edges[i + 1]) / 2 for i in range(n_bins)],
                ticktext=tick_labels,
                len=cb_len,
                x=COLORBAR_X,
                y=colorbar_y,
                yanchor="middle",
                thickness=COLORBAR_THICKNESS,
                title_font=dict(size=COLORBAR_TITLE_FONT_SIZE),
                tickfont=dict(size=COLORBAR_TICK_FONT_SIZE),
            ),
        ),
        hoverinfo="skip",
        showlegend=False,
    )


def compute_y_ticks(tip_order, y_range=None, max_ticks=AXIS_Y_MAX_TICKS):
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
            step = len(indices) // max_ticks
            indices = indices[::step]
    return indices, [tip_order[i] for i in indices]


# =============================================================================
# GENE VISUALIZATION (IGV.js CONVENTIONS)
# =============================================================================


def create_gene_arrow(start, end, strand, y_position, height=PANEL_ANNOTATION_TRACK_HEIGHT):
    """Arrow-shaped polygon for a gene at the given y_position (track index units).

    The body covers the first 75 % of the gene and the arrowhead the last 25 %,
    so the direction is obvious at every scale.  The head flares to 2.8 × the
    body half-height.
    """
    gene_length = max(end - start, 1)
    arrow_size = max(ARROW_SIZE_MIN, gene_length * ARROW_HEAD_RATIO)
    h = height * 1.5  # body half-height — noticeably thick
    H = height * 2.0  # arrowhead half-height — wider than body but not extreme

    if strand == "+":
        body_end = end - arrow_size
        if body_end <= start:
            # Gene too short for a body — just a triangle
            x = [start, end, start, start]
            y = [y_position - H, y_position, y_position + H, y_position - H]
        else:
            x = [start, body_end, body_end, end, body_end, body_end, start, start]
            y = [
                y_position - h,
                y_position - h,
                y_position - H,
                y_position,
                y_position + H,
                y_position + h,
                y_position + h,
                y_position - h,
            ]
    else:
        body_start = start + arrow_size
        if body_start >= end:
            x = [end, start, end, end]
            y = [y_position - H, y_position, y_position + H, y_position - H]
        else:
            x = [end, body_start, body_start, start, body_start, body_start, end, end]
            y = [
                y_position - h,
                y_position - h,
                y_position - H,
                y_position,
                y_position + H,
                y_position + h,
                y_position + h,
                y_position - h,
            ]

    return {"x": x, "y": y}


def build_gene_hover(gene_row):
    """
    Build hover text for a gene with comprehensive information.

    Args:
        gene_row: Series containing gene information

    Returns:
        Formatted hover string with gene information
    """
    # Calculate gene length
    length = int(gene_row["stop"]) - int(gene_row["start"])

    hover = f"<b>{gene_row['locus_tag']}</b><br>"
    hover += f"<b>Feature:</b> {gene_row['ftype']}<br>"
    hover += f"<b>Location:</b> {int(gene_row['start']):,} - {int(gene_row['stop']):,}<br>"
    hover += f"<b>Length:</b> {length:,} bp<br>"
    hover += f"<b>Strand:</b> {'Forward (+)' if gene_row['strand'] == '+' else 'Reverse (-)'}<br>"

    # Gene name
    if "gene_name" in gene_row and pd.notna(gene_row["gene_name"]):
        hover += f"<b>Gene:</b> {gene_row['gene_name']}<br>"

    # Product description
    if "product" in gene_row and pd.notna(gene_row["product"]):
        product = str(gene_row["product"])
        if len(product) > 100:
            product = product[:100] + "..."
        hover += f"<b>Product:</b> {product}<br>"

    # EC number
    if "ec_number" in gene_row and pd.notna(gene_row["ec_number"]):
        hover += f"<b>EC Number:</b> {gene_row['ec_number']}<br>"

    # Additional annotation sources
    for col, label in [
        ("swissprot_gene", "SwissProt Gene"),
        ("swissprot_KO", "KEGG Orthology"),
        ("swissprot_Pfam", "Pfam"),
        ("swissprot_CAZy", "CAZy"),
        ("swissprot_TIGRFAMs", "TIGRFAMs"),
    ]:
        if col in gene_row and pd.notna(gene_row[col]):
            val = str(gene_row[col])
            if len(val) > 50:
                val = val[:50] + "..."
            hover += f"<b>{label}:</b> {val}<br>"

    return hover


def create_annotation_traces(df_annotations, x_range=None):
    """
    Create Plotly traces for gene annotations, one horizontal track per ftype.

    Returns:
        traces        – list of go.Scatter trace objects
        y_tick_vals   – list of int y-positions for ftype labels
        y_tick_text   – list of ftype name strings
        n_tracks      – total number of ftype tracks
    """
    if df_annotations is None or df_annotations.empty:
        return [], [], [], 0

    # Assign a vertical track to each ftype (sorted alphabetically)
    ftypes = sorted(df_annotations["ftype"].dropna().unique())
    n_tracks = len(ftypes)
    ftype_to_y = {ft: i for i, ft in enumerate(ftypes)}

    traces = []
    x_min, x_max = sorted(x_range) if x_range is not None else (None, None)

    for ftype in ftypes:
        group = df_annotations[df_annotations["ftype"] == ftype].sort_values("start")
        color = FEATURE_COLORS.get(ftype, FEATURE_COLORS["default"])
        y_center = ftype_to_y[ftype]

        for _, gene in group.iterrows():
            if x_min is not None and (gene["stop"] < x_min or gene["start"] > x_max):
                continue

            arrow_data = create_gene_arrow(
                int(gene["start"]), int(gene["stop"]), gene["strand"], y_position=y_center
            )
            hover_text = build_gene_hover(gene)

            # Subtle edge: slightly darker shade of fill
            edge_color = color
            traces.append(
                go.Scatter(
                    x=arrow_data["x"],
                    y=arrow_data["y"],
                    mode="lines",
                    fill="toself",
                    fillcolor=color,
                    line=dict(color=edge_color, width=0.8),
                    hovertemplate=f"{hover_text}<extra></extra>",
                    showlegend=False,
                    name=str(gene.get("locus_tag", "")),
                    opacity=0.82,
                )
            )

    return traces, list(range(n_tracks)), ftypes, n_tracks


def build_figure(
    tree_data,
    tip_order,
    tree_max_xy,
    df_intervals,
    bin_edges,
    colors,
    df_annotations=None,
    query_id=None,
    interval_linewidth=3.0,
    tree_linewidth=1.5,
    uirevision="base",
    y_range=None,
    x_range=None,
    y_range_limits=None,
    x_range_limits=None,
    leaf_distances=None,
):
    """
    Build the complete figure with tree, interval, and optional annotation panels.

    Args:
        tree_data: Tree layout data
        tip_order: Ordered list of leaf names
        tree_max_xy: Maximum x,y coordinates for tree
        df_intervals: Interval data
        bin_edges: Color bin edges
        colors: Color list for bins
        df_annotations: Optional annotation DataFrame
        query_id: Current query sequence ID for filtering annotations
        interval_linewidth: Width of interval lines
        tree_linewidth: Width of tree lines
        uirevision: UI revision for preserving zoom state
        y_range: Current y-axis range
        x_range: Current x-axis range
        y_range_limits: Y-axis limits
        x_range_limits: X-axis limits
        leaf_distances: Distances from query to leaves

    Returns:
        Plotly figure object
    """
    # Determine if we have annotations to display
    has_annotations = df_annotations is not None and not df_annotations.empty

    if has_annotations:
        # 2x2 layout: tree + intervals (row 1), empty + annotations (row 2)
        # shared_yaxes: row 1 shares tree ↔ interval y-axes (row 2 shares
        #   empty ↔ annotation — harmless, overridden below).
        # shared_xaxes: col 2 shares interval ↔ annotation x-axes so genomic
        #   position zoom is synchronized.
        fig = make_subplots(
            rows=2,
            cols=2,
            shared_yaxes=True,
            shared_xaxes=True,
            column_widths=[PANEL_TREE_WIDTH, PANEL_INTERVAL_WIDTH],
            row_heights=[PANEL_TREE_INTERVAL_HEIGHT, PANEL_ANNOTATION_HEIGHT],
            horizontal_spacing=PANEL_SPACING,
            vertical_spacing=PANEL_VERTICAL_SPACING,
        )
    else:
        # 1x2 layout: tree + intervals only
        fig = make_subplots(
            rows=1,
            cols=2,
            shared_yaxes=True,
            column_widths=[PANEL_TREE_WIDTH, PANEL_INTERVAL_WIDTH],
            horizontal_spacing=PANEL_SPACING,
        )

    # Compute colorbar vertical center so it sits in the middle of the interval
    # row regardless of whether the annotation panel is present.
    # In paper coords: with annotations, row 1 top ≈ 1.0, row 1 bottom ≈
    #   PANEL_ANNOTATION_HEIGHT + PANEL_VERTICAL_SPACING.  Without annotations
    #   the single row spans [0, 1].
    if has_annotations:
        row1_bottom = PANEL_ANNOTATION_HEIGHT + PANEL_VERTICAL_SPACING
        row1_top = 1.0
    else:
        row1_bottom = 0.0
        row1_top = 1.0
    cb_y = (row1_bottom + row1_top) / 2.0
    cb_len = (row1_top - row1_bottom) * 0.85  # slightly shorter than the row

    # Add colorbar (in interval panel)
    fig.add_trace(
        make_colorbar(bin_edges, colors, colorbar_y=cb_y, colorbar_len=cb_len), row=1, col=2
    )

    fig.update_layout(
        autosize=True,
        margin=FIG_MARGIN,
        plot_bgcolor=COLORS["plot_bg"],
        paper_bgcolor=COLORS["paper_bg"],
        hovermode="closest",
        uirevision=uirevision,
        font=dict(family=FIG_FONT_FAMILY, size=FIG_FONT_SIZE),
    )

    # Tree panel (left)
    fig.add_trace(
        go.Scattergl(
            x=tree_data["x"],
            y=tree_data["y"],
            text=tree_data["text"],
            hovertemplate="%{text}<extra></extra>",
            mode="lines",
            line=dict(color=COLORS["tree_line"], width=tree_linewidth),
            showlegend=False,
        ),
        row=1,
        col=1,
    )

    # Interval panel (right)
    if not df_intervals.empty:
        traces_by_color = batch_interval_by_color(df_intervals, bin_edges, colors, leaf_distances)
        for color, data in traces_by_color.items():
            fig.add_trace(
                go.Scattergl(
                    x=data["x"],
                    y=data["y"],
                    mode="lines",
                    line=dict(width=interval_linewidth, color=color),
                    text=data["text"],
                    hovertemplate="%{text}<extra></extra>",
                    showlegend=False,
                ),
                row=1,
                col=2,
            )

    # Annotation panel (row 2, col 2)
    ann_y_tick_vals, ann_y_tick_text, ann_n_tracks = [], [], 0
    if has_annotations:
        query_annotations = filter_annotations(df_annotations, query_id, x_range)
        if query_annotations is not None and not query_annotations.empty:
            ann_traces, ann_y_tick_vals, ann_y_tick_text, ann_n_tracks = create_annotation_traces(
                query_annotations, x_range
            )
            for trace in ann_traces:
                fig.add_trace(trace, row=2, col=2)

    # ── Tree panel (row 1, col 1) ──────────────────────────────────────────────
    tree_range_min = -(tree_max_xy[0] * (TREE_X_MARGIN - 1))
    tree_range_max = tree_max_xy[0] * TREE_X_MARGIN

    # Resolve x_range before configuring axes
    if x_range is None:
        x_range = x_range_limits

    fig.update_xaxes(
        title=TREE_X_TITLE,
        title_font=dict(size=AXIS_TITLE_FONT_SIZE),
        range=[tree_range_min, tree_range_max],
        fixedrange=True,
        side="bottom",
        showline=True,
        showgrid=False,
        showticklabels=True,
        linewidth=AXIS_LINE_WIDTH,
        linecolor=COLORS["axis_line"],
        mirror=False,
        tickfont=dict(size=AXIS_TICK_FONT_SIZE),
        ticks="outside",
        nticks=5,
        row=1,
        col=1,
    )
    fig.update_yaxes(
        autorange=y_range is None,
        range=y_range if y_range else y_range_limits,
        showline=False,
        zeroline=False,
        showgrid=False,
        showticklabels=False,
        minallowed=y_range_limits[0],
        maxallowed=y_range_limits[1],
        row=1,
        col=1,
    )

    # ── Interval panel (row 1, col 2) ─────────────────────────────────────────
    tickvals, ticktext = compute_y_ticks(tip_order, y_range)

    # When annotations are present the shared x-axis renders tick labels at the
    # very bottom of column 2 (below row 2).  Suppress them on row 1 here to
    # avoid duplication; the title and ticks are set on row 2 further below.
    fig.update_xaxes(
        title=INTERVAL_X_TITLE if not has_annotations else "",
        title_font=dict(size=AXIS_TITLE_FONT_SIZE),
        range=x_range,
        tickmode="auto",
        nticks=AXIS_X_NTICKS,
        showgrid=True,
        gridcolor=COLORS["grid"],
        gridwidth=AXIS_GRID_WIDTH,
        showline=True,
        linewidth=AXIS_LINE_WIDTH,
        linecolor=COLORS["axis_line"],
        mirror=False,
        tickfont=dict(size=AXIS_TICK_FONT_SIZE),
        showticklabels=not has_annotations,
        ticks="outside" if not has_annotations else "",
        minallowed=x_range_limits[0],
        maxallowed=x_range_limits[1],
        row=1,
        col=2,
    )
    fig.update_yaxes(
        autorange=y_range is None,
        range=y_range if y_range else y_range_limits,
        tickvals=tickvals,
        ticktext=ticktext,
        tickfont=dict(size=AXIS_TICK_FONT_SIZE),
        showgrid=False,
        zeroline=False,
        showline=True,
        linewidth=AXIS_LINE_WIDTH,
        linecolor=COLORS["axis_line"],
        mirror=False,
        minallowed=y_range_limits[0],
        maxallowed=y_range_limits[1],
        row=1,
        col=2,
    )

    # ── Annotation panel (row 2) ─────────────────────────────────────────────
    if has_annotations:
        ann_y_range = [-0.5, max(ann_n_tracks - 0.5, 0.5)]

        # X-axis for annotation panel (row 2, col 2):
        # This is the ONLY place the genomic-position ticks and title appear.
        # (Row 1, col 2 has its ticks/title suppressed above.)
        fig.update_xaxes(
            title=INTERVAL_X_TITLE,
            title_font=dict(size=AXIS_TITLE_FONT_SIZE),
            range=x_range,
            tickmode="auto",
            nticks=AXIS_X_NTICKS,
            showgrid=True,
            gridcolor=COLORS["grid"],
            gridwidth=AXIS_GRID_WIDTH,
            showline=True,
            linewidth=AXIS_LINE_WIDTH,
            linecolor=COLORS["axis_line"],
            mirror=False,
            tickfont=dict(size=AXIS_TICK_FONT_SIZE),
            showticklabels=True,
            ticks="outside",
            minallowed=x_range_limits[0],
            maxallowed=x_range_limits[1],
            row=2,
            col=2,
        )

        # Y-axis: ftype track labels on the right, no vertical zoom
        fig.update_yaxes(
            autorange=False,
            fixedrange=True,
            range=ann_y_range,
            tickvals=ann_y_tick_vals,
            ticktext=list(ann_y_tick_text),
            tickfont=dict(size=10, color="#555"),
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

        # Row 2, col 1 (below tree) — completely hidden
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


def nearest_value(value, valid_values):
    idx = np.searchsorted(valid_values, value)
    if idx == 0:
        return valid_values[0]
    if idx == len(valid_values):
        return valid_values[-1]
    left = valid_values[idx - 1]
    right = valid_values[idx]
    return left if (value - left) < (right - value) else right


def create_slider_marks(dist_th_l):
    if not dist_th_l:
        return {}
    n = len(dist_th_l)
    marks = {dist_th_l[0]: str(dist_th_l[0]), dist_th_l[-1]: str(dist_th_l[-1])}
    for i in range(1, n - 1):
        marks[dist_th_l[i]] = "\u200b"
    return marks


# =============================================================================
# DASH APP
# =============================================================================


def control_label(text):
    """Create a standardized control label."""
    return html.Label(
        text,
        style={
            "fontWeight": "bold",
            "marginRight": f"{UI_LABEL_MARGIN_RIGHT}px",
            "fontSize": UI_LABEL_FONT_SIZE,
        },
    )


def get_toggle_button_style(is_active, position="middle"):
    """
    Get style dict for a toggle button.

    Args:
        is_active: Whether this button is the active one
        position: "left", "right", or "middle"
    """
    border_radius = {"left": "4px 0 0 4px", "right": "0 4px 4px 0", "middle": "0"}[position]

    return {
        "padding": UI_TOGGLE_PADDING,
        "fontSize": UI_TOGGLE_FONT_SIZE,
        "cursor": "pointer",
        "border": f"{UI_TOGGLE_BORDER_WIDTH}px solid {COLORS['button_border']}",
        "borderRight": "none" if position != "right" else None,
        "borderRadius": border_radius,
        "backgroundColor": COLORS["button_bg"] if is_active else COLORS["button_bg_inactive"],
        "fontWeight": "bold" if is_active else "normal",
        "color": COLORS["button_text"],
    }


def control_panel(children):
    return html.Div(
        children,
        style={
            "display": "flex",
            "alignItems": "center",
            "gap": UI_PANEL_GAP,
            "padding": UI_PANEL_PADDING,
            "backgroundColor": COLORS["panel_bg"],
            "borderRadius": UI_PANEL_BORDER_RADIUS,
            "marginBottom": UI_PANEL_MARGIN_BOTTOM,
            "boxShadow": UI_PANEL_BOX_SHADOW,
            "flexWrap": "wrap",
        },
    )


def build_layout(seq_id_l, dist_th_l, has_pruned_tree, has_strand=False, initial_prune=True):
    """
    Build the application layout.

    Args:
        seq_id_l: List of sequence IDs
        dist_th_l: List of distance thresholds
        has_pruned_tree: Whether a pruned tree is available
        initial_prune: Initial state for prune toggle (True = start with pruned view)
    """
    dmin = dist_th_l[0] if dist_th_l else 0.0
    dmax = dist_th_l[-1] if dist_th_l else 0.0
    return html.Div(
        [
            dcc.Store(id="y-range-store", data=None),
            dcc.Store(id="x-range-store", data=None),
            dcc.Store(id="tree-view-store", data="phylogeny"),
            dcc.Store(id="prune-store", data=initial_prune if has_pruned_tree else False),
            dcc.Store(id="strand-store", data="+"),
            dcc.Store(id="mode-store", data="focus"),
            dcc.Download(id="download-data"),
            dcc.Download(id="download-plot"),
            control_panel(
                [
                    # ── Query navigator ───────────────────────────────────────
                    html.Div(
                        [
                            control_label("Query:"),
                            html.Button(
                                "◀",
                                id="prev-query-btn",
                                n_clicks=0,
                                style={
                                    "padding": UI_NAV_PADDING,
                                    "fontSize": UI_NAV_FONT_SIZE,
                                    "cursor": "pointer",
                                    "border": f"{UI_NAV_BORDER_WIDTH}px solid {COLORS['button_border']}",
                                    "borderRadius": "4px 0 0 4px",
                                    "backgroundColor": COLORS["button_bg"],
                                },
                            ),
                            dcc.Dropdown(
                                id="query-dropdown",
                                options=[{"label": q, "value": q} for q in seq_id_l],
                                value=seq_id_l[0] if seq_id_l else None,
                                style={"width": UI_QUERY_DROPDOWN_WIDTH},
                            ),
                            html.Button(
                                "▶",
                                id="next-query-btn",
                                n_clicks=0,
                                style={
                                    "padding": UI_NAV_PADDING,
                                    "fontSize": UI_NAV_FONT_SIZE,
                                    "cursor": "pointer",
                                    "border": f"{UI_NAV_BORDER_WIDTH}px solid {COLORS['button_border']}",
                                    "borderRadius": "0 4px 4px 0",
                                    "backgroundColor": COLORS["button_bg"],
                                },
                            ),
                        ],
                        style={"display": "flex", "alignItems": "center", "gap": 0},
                    ),
                    # ── Vertical divider ─────────────────────────────────────
                    html.Div(
                        style={
                            "width": "1px",
                            "alignSelf": "stretch",
                            "backgroundColor": COLORS["button_border"],
                            "opacity": "0.35",
                            "margin": "0 4px",
                        }
                    ),
                    # ── Distance slider + focus/overlap mode ──────────────────
                    html.Div(
                        [
                            control_label("Distance:"),
                            # Focus mode: single Slider (select exactly one threshold)
                            html.Div(
                                dcc.Slider(
                                    id="dist-slider-focus",
                                    min=dmin,
                                    max=dmax,
                                    step=None,
                                    marks=create_slider_marks(dist_th_l),
                                    value=dist_th_l[0] if dist_th_l else dmin,
                                    tooltip={"placement": "bottom", "always_visible": False},
                                    updatemode="mouseup",
                                ),
                                id="dist-slider-focus-wrapper",
                                style={"minWidth": UI_SLIDER_MIN_WIDTH, "flexShrink": 1},
                            ),
                            # Overlap mode: RangeSlider (filter by range)
                            html.Div(
                                dcc.RangeSlider(
                                    id="dist-slider-overlap",
                                    min=dmin,
                                    max=dmax,
                                    step=None,
                                    marks=create_slider_marks(dist_th_l),
                                    value=[dmin, dmax],
                                    tooltip={"placement": "bottom", "always_visible": False},
                                    updatemode="mouseup",
                                    allowCross=False,
                                ),
                                id="dist-slider-overlap-wrapper",
                                style={
                                    "minWidth": UI_SLIDER_MIN_WIDTH,
                                    "flexShrink": 1,
                                    "display": "none",
                                },
                            ),
                            # focus | overlap toggle
                            html.Div(
                                [
                                    html.Button(
                                        "focus",
                                        id="mode-focus-btn",
                                        n_clicks=0,
                                        style=get_toggle_button_style(True, "left"),
                                    ),
                                    html.Button(
                                        "overlap",
                                        id="mode-overlap-btn",
                                        n_clicks=0,
                                        style=get_toggle_button_style(False, "right"),
                                    ),
                                ],
                                style={"display": "flex", "marginLeft": 8},
                            ),
                        ],
                        style={"display": "flex", "alignItems": "center", "gap": UI_SLIDER_GAP},
                    ),
                    # ── Vertical divider ─────────────────────────────────────
                    html.Div(
                        style={
                            "width": "1px",
                            "alignSelf": "stretch",
                            "backgroundColor": COLORS["button_border"],
                            "opacity": "0.35",
                            "margin": "0 4px",
                        }
                    ),
                    # ── Tree controls: view mode + scope ─────────────────────
                    # Phylogeny/Cladogram and Pruned/Full are both tree-display
                    # settings, so they share a single "Tree:" label.
                    html.Div(
                        [
                            control_label("Tree:"),
                            # view mode
                            html.Button(
                                "phylogeny",
                                id="phylogeny-btn",
                                n_clicks=0,
                                style=get_toggle_button_style(True, "left"),
                            ),
                            html.Button(
                                "cladogram",
                                id="cladogram-btn",
                                n_clicks=0,
                                style=get_toggle_button_style(False, "right"),
                            ),
                            # scope (pruned / full) — only when pruned tree exists
                            (
                                html.Div(
                                    [
                                        html.Button(
                                            "pruned",
                                            id="pruned-btn",
                                            n_clicks=0,
                                            style=get_toggle_button_style(initial_prune, "left"),
                                        ),
                                        html.Button(
                                            "full",
                                            id="full-tree-btn",
                                            n_clicks=0,
                                            style=get_toggle_button_style(
                                                not initial_prune, "right"
                                            ),
                                        ),
                                    ],
                                    style={
                                        "display": "flex" if has_pruned_tree else "none",
                                        "marginLeft": 6,
                                    },
                                )
                                if has_pruned_tree
                                else html.Div()
                            ),
                        ],
                        style={"display": "flex", "alignItems": "center", "gap": 0},
                    ),
                    # ── Strand toggle (only when STRAND column present) ────────
                    html.Div(
                        [
                            # thin divider before strand group
                            html.Div(
                                style={
                                    "width": "1px",
                                    "alignSelf": "stretch",
                                    "backgroundColor": COLORS["button_border"],
                                    "opacity": "0.35",
                                    "marginRight": 8,
                                }
                            ),
                            control_label("Strand:"),
                            html.Button(
                                "+ (fw)",
                                id="strand-fwd-btn",
                                n_clicks=0,
                                title="Forward strand",
                                style=get_toggle_button_style(True, "left"),
                            ),
                            html.Button(
                                "− (rc)",
                                id="strand-rev-btn",
                                n_clicks=0,
                                title="Reverse complement strand",
                                style=get_toggle_button_style(False, "right"),
                            ),
                        ],
                        style={
                            "display": "flex" if has_strand else "none",
                            "alignItems": "center",
                            "gap": 0,
                        },
                    ),
                    # ── Vertical divider ─────────────────────────────────────
                    html.Div(
                        style={
                            "width": "1px",
                            "alignSelf": "stretch",
                            "backgroundColor": COLORS["button_border"],
                            "opacity": "0.35",
                            "margin": "0 4px",
                        }
                    ),
                    # ── Color scheme ──────────────────────────────────────────
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
                                style={"width": UI_COLORSCALE_DROPDOWN_WIDTH},
                                clearable=False,
                            ),
                        ],
                        style={"display": "flex", "alignItems": "center", "gap": UI_COLORSCALE_GAP},
                    ),
                    # ── Export (pushed to far right) ──────────────────────────
                    html.Div(
                        [
                            # vertical divider
                            html.Div(
                                style={
                                    "width": "1px",
                                    "alignSelf": "stretch",
                                    "backgroundColor": COLORS["button_border"],
                                    "opacity": "0.35",
                                    "marginRight": 8,
                                }
                            ),
                            control_label("Export:"),
                            # Data TSV
                            html.Button(
                                "Data",
                                id="export-btn",
                                n_clicks=0,
                                style={
                                    "padding": "5px 10px",
                                    "fontSize": 12,
                                    "cursor": "pointer",
                                    "border": f"2px solid {COLORS['button_border']}",
                                    "borderRadius": 4,
                                    "backgroundColor": COLORS["button_bg"],
                                    "fontWeight": "bold",
                                    "color": COLORS["button_text"],
                                    "marginRight": 8,
                                },
                            ),
                            # PDF export with W / H dimensions
                            html.Span(
                                "W", style={"fontSize": 11, "marginRight": 2, "color": "#555"}
                            ),
                            dcc.Input(
                                id="export-width",
                                type="number",
                                value=1600,
                                min=100,
                                max=10000,
                                step=10,
                                style={
                                    "width": 60,
                                    "fontSize": 12,
                                    "padding": "3px 4px",
                                    "border": f"2px solid {COLORS['button_border']}",
                                    "borderRadius": 4,
                                    "textAlign": "center",
                                },
                            ),
                            html.Span(
                                "H",
                                style={"fontSize": 11, "margin": "0 2px 0 5px", "color": "#555"},
                            ),
                            dcc.Input(
                                id="export-height",
                                type="number",
                                value=900,
                                min=100,
                                max=10000,
                                step=10,
                                style={
                                    "width": 52,
                                    "fontSize": 12,
                                    "padding": "3px 4px",
                                    "border": f"2px solid {COLORS['button_border']}",
                                    "borderRadius": 4,
                                    "textAlign": "center",
                                },
                            ),
                            html.Button(
                                "PDF",
                                id="export-pdf-btn",
                                n_clicks=0,
                                style={
                                    "padding": "5px 10px",
                                    "fontSize": 12,
                                    "cursor": "pointer",
                                    "border": f"2px solid {COLORS['button_border']}",
                                    "borderRadius": 4,
                                    "backgroundColor": COLORS["button_bg"],
                                    "fontWeight": "bold",
                                    "color": COLORS["button_text"],
                                    "marginLeft": 5,
                                },
                            ),
                        ],
                        style={
                            "display": "flex",
                            "alignItems": "center",
                            "gap": 3,
                            "marginLeft": "auto",
                        },
                    ),
                ]
            ),
            html.Div(
                id="query-title",
                style={
                    "textAlign": "center",
                    "fontWeight": "bold",
                    "fontSize": UI_TITLE_FONT_SIZE,
                    "marginBottom": UI_TITLE_MARGIN_BOTTOM,
                },
            ),
            html.Div(
                dcc.Graph(
                    id="graph",
                    style={"height": UI_GRAPH_HEIGHT},
                    config={
                        "doubleClick": "reset",
                        "displayModeBar": True,
                        "displaylogo": False,
                        "modeBarButtonsToRemove": ["lasso2d", "select2d"],
                        "scrollZoom": False,
                    },
                ),
                style={"overflowY": "auto", "height": UI_CONTAINER_HEIGHT},
            ),
        ],
        style={"fontFamily": FIG_FONT_FAMILY, "padding": UI_APP_PADDING},
    )


def extract_y_range(relayout):
    """Extract y-axis range from relayoutData. Returns no_update if y was not changed."""
    if not relayout:
        return no_update
    if any(key.endswith(".autorange") for key in relayout):
        return None
    if "yaxis.range[0]" in relayout and "yaxis.range[1]" in relayout:
        return [relayout["yaxis.range[0]"], relayout["yaxis.range[1]"]]
    if "yaxis.range" in relayout:
        return relayout["yaxis.range"]
    return no_update


def extract_x_range(relayout):
    """Extract x-axis range from relayoutData. Returns no_update if x was not changed."""
    if not relayout:
        return no_update
    if any(key.endswith(".autorange") for key in relayout):
        return None

    # Check interval panel (xaxis2) first
    if "xaxis2.range[0]" in relayout and "xaxis2.range[1]" in relayout:
        return [relayout["xaxis2.range[0]"], relayout["xaxis2.range[1]"]]
    if "xaxis2.range" in relayout:
        return relayout["xaxis2.range"]

    # Also check annotation panel (xaxis4) - they should be synchronized
    if "xaxis4.range[0]" in relayout and "xaxis4.range[1]" in relayout:
        return [relayout["xaxis4.range[0]"], relayout["xaxis4.range[1]"]]
    if "xaxis4.range" in relayout:
        return relayout["xaxis4.range"]

    return no_update


# =============================================================================
# APP FACTORY
# =============================================================================


def create_app(input_path, tree_path, query=None, annotation_path=None):
    """
    Create the Dash application with support for toggling between full and pruned trees.

    Args:
        input_path: Path to TSV data file
        tree_path: Path to Newick tree file
        query: Optional query leaf name for distance calculations
        annotation_path: Optional path to annotation TSV file
    """
    df = pd.read_csv(input_path, sep="\t")

    # Load annotations if provided
    df_annotations = load_annotations(annotation_path)

    required_cols = {"QUERY_ID", "REF_ID", "INTERVAL_START", "INTERVAL_END", "SEQ_LEN", "DIST_TH"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Input file missing required columns: {', '.join(sorted(missing))}")

    if not (df["INTERVAL_START"] != df["INTERVAL_END"]).any():
        raise ValueError("Input file contains no intervals (all INTERVAL_START == INTERVAL_END).")

    full_tree = load_tree(tree_path)
    leaf_names = {leaf.name for leaf in full_tree.iter_leaves()}

    if query is not None and query not in leaf_names:
        raise ValueError(
            f"Query '{query}' not found in tree. "
            f"Available leaves ({len(leaf_names)}): {', '.join(sorted(leaf_names)[:10])}..."
        )

    retained_l = get_retained_leaves(df)
    has_pruned_tree = len(retained_l) > 0  # and len(retained_l) < len(leaf_names)
    pruned_tree = prune_tree(full_tree, retained_l) if has_pruned_tree else None
    if pruned_tree is not None:
        pruned_tree.ladderize()

    def compute_tree_data(tree, query):
        """Compute layouts, tip order, and distances for a tree."""
        query_distances = None
        leaf_distances = None
        if query:
            query_distances = compute_path_distances(tree, query)
            leaf_distances = {
                leaf.name: query_distances[leaf]
                for leaf in tree.iter_leaves()
                if leaf in query_distances
            }

        phylo_layout = compute_tree_layout(tree, "phylogeny", query_distances)
        clado_layout = compute_tree_layout(tree, "cladogram", query_distances)
        _, tip_order, _ = phylo_layout

        return {
            "phylo_layout": phylo_layout,
            "clado_layout": clado_layout,
            "tip_order": tip_order,
            "leaf_distances": leaf_distances,
        }

    full_tree_data = compute_tree_data(full_tree, query)
    pruned_tree_data = compute_tree_data(pruned_tree, query) if has_pruned_tree else None

    full_df = add_tip_order(df, full_tree_data["tip_order"])
    pruned_df = add_tip_order(df, pruned_tree_data["tip_order"]) if has_pruned_tree else full_df

    dist_th_l = tuple(get_distance_thresholds(df))
    if not dist_th_l:
        raise ValueError("No distance thresholds found in data.")

    seq_id_l = get_sequence_identifiers(df)
    has_strand = "STRAND" in df.columns

    obtain_filtered_full = make_cached_filter(full_df)
    obtain_filtered_pruned = (
        make_cached_filter(pruned_df) if has_pruned_tree else obtain_filtered_full
    )

    app = Dash(__name__)
    app.layout = build_layout(seq_id_l, dist_th_l, has_pruned_tree, has_strand=has_strand)

    @app.callback(
        Output("query-dropdown", "value"),
        Input("prev-query-btn", "n_clicks"),
        Input("next-query-btn", "n_clicks"),
        State("query-dropdown", "value"),
        prevent_initial_call=True,
    )
    def navigate_query(prev_clicks, next_clicks, current):
        """Cycle through queries with prev/next buttons."""
        if current is None or not seq_id_l:
            return no_update
        idx = seq_id_l.index(current) if current in seq_id_l else 0
        if ctx.triggered_id == "prev-query-btn":
            idx = (idx - 1) % len(seq_id_l)
        elif ctx.triggered_id == "next-query-btn":
            idx = (idx + 1) % len(seq_id_l)
        return seq_id_l[idx]

    @app.callback(Output("query-title", "children"), Input("query-dropdown", "value"))
    def update_query_title(seq_id):
        if seq_id is None:
            return ""
        idx = seq_id_l.index(seq_id) if seq_id in seq_id_l else 0
        return f"{seq_id}  ({idx + 1}/{len(seq_id_l)})"

    @app.callback(
        Output("y-range-store", "data"),
        Output("x-range-store", "data"),
        Input("graph", "relayoutData"),
        prevent_initial_call=False,
    )
    def update_view_stores(relayout):
        return extract_y_range(relayout), extract_x_range(relayout)

    @app.callback(
        Output("tree-view-store", "data"),
        Input("phylogeny-btn", "n_clicks"),
        Input("cladogram-btn", "n_clicks"),
        prevent_initial_call=False,
    )
    def update_tree_view(phylogeny_clicks, cladogram_clicks):
        """Handle tree view toggle buttons."""
        if ctx.triggered_id == "phylogeny-btn":
            return "phylogeny"
        elif ctx.triggered_id == "cladogram-btn":
            return "cladogram"
        return "phylogeny"

    @app.callback(
        Output("phylogeny-btn", "style"),
        Output("cladogram-btn", "style"),
        Input("tree-view-store", "data"),
    )
    def update_tree_view_button_styles(tree_view_mode):
        """Update tree view button styles (phylogeny/cladogram)."""
        return (
            get_toggle_button_style(tree_view_mode == "phylogeny", "left"),
            get_toggle_button_style(tree_view_mode == "cladogram", "right"),
        )

    # Add prune toggle callbacks only if pruned tree exists
    if has_pruned_tree:

        @app.callback(
            Output("prune-store", "data"),
            Output("y-range-store", "data", allow_duplicate=True),
            Output("x-range-store", "data", allow_duplicate=True),
            Input("pruned-btn", "n_clicks"),
            Input("full-tree-btn", "n_clicks"),
            prevent_initial_call=True,
        )
        def toggle_prune(pruned_clicks, full_clicks):
            """Handle prune/full tree toggle and reset view."""
            if ctx.triggered_id == "pruned-btn":
                return True, None, None
            elif ctx.triggered_id == "full-tree-btn":
                return False, None, None
            return no_update, no_update, no_update

        @app.callback(
            Output("pruned-btn", "style"),
            Output("full-tree-btn", "style"),
            Input("prune-store", "data"),
        )
        def update_prune_button_styles(is_pruned):
            """Update prune button styles."""
            return (
                get_toggle_button_style(is_pruned, "left"),
                get_toggle_button_style(not is_pruned, "right"),
            )

    @app.callback(
        Output("strand-store", "data"),
        Output("y-range-store", "data", allow_duplicate=True),
        Output("x-range-store", "data", allow_duplicate=True),
        Input("strand-fwd-btn", "n_clicks"),
        Input("strand-rev-btn", "n_clicks"),
        prevent_initial_call=True,
    )
    def toggle_strand(fwd_clicks, rev_clicks):
        """Toggle between forward (+) and reverse-complement (−) strand; reset view."""
        if ctx.triggered_id == "strand-fwd-btn":
            return "+", None, None
        elif ctx.triggered_id == "strand-rev-btn":
            return "-", None, None
        return no_update, no_update, no_update

    @app.callback(
        Output("strand-fwd-btn", "style"),
        Output("strand-rev-btn", "style"),
        Input("strand-store", "data"),
    )
    def update_strand_button_styles(strand):
        """Update strand button styles to reflect active strand."""
        return (
            get_toggle_button_style(strand == "+", "left"),
            get_toggle_button_style(strand == "-", "right"),
        )

    @app.callback(
        Output("mode-store", "data"),
        Input("mode-focus-btn", "n_clicks"),
        Input("mode-overlap-btn", "n_clicks"),
        prevent_initial_call=False,
    )
    def update_mode(focus_clicks, overlap_clicks):
        """Toggle between focus (single threshold) and overlap (range) mode."""
        if ctx.triggered_id == "mode-overlap-btn":
            return "overlap"
        return "focus"  # default

    @app.callback(
        Output("mode-focus-btn", "style"),
        Output("mode-overlap-btn", "style"),
        Output("dist-slider-focus-wrapper", "style"),
        Output("dist-slider-overlap-wrapper", "style"),
        Input("mode-store", "data"),
    )
    def update_mode_ui(mode):
        """Sync button styles and slider visibility with current mode."""
        is_focus = mode == "focus"
        slider_w = {"minWidth": UI_SLIDER_MIN_WIDTH, "flexShrink": 1}
        hidden = {**slider_w, "display": "none"}
        return (
            get_toggle_button_style(is_focus, "left"),
            get_toggle_button_style(not is_focus, "right"),
            slider_w if is_focus else hidden,
            slider_w if not is_focus else hidden,
        )

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
        Input("prune-store", "data") if has_pruned_tree else State("prune-store", "data"),
    )
    def update_figure(
        seq_id,
        focus_val,
        overlap_val,
        mode,
        colorscheme,
        y_range,
        x_range,
        tree_view_mode,
        strand,
        is_pruned,
    ):
        """Update the figure based on all current settings."""
        if seq_id is None:
            return go.Figure()

        # Select the appropriate tree data based on prune state
        if has_pruned_tree and not is_pruned:
            current_data = full_tree_data
            obtain_filtered = obtain_filtered_full
        else:
            current_data = pruned_tree_data if has_pruned_tree else full_tree_data
            obtain_filtered = obtain_filtered_pruned if has_pruned_tree else obtain_filtered_full

        tip_order = current_data["tip_order"]
        leaf_distances = current_data["leaf_distances"]
        phylo_layout = current_data["phylo_layout"]
        clado_layout = current_data["clado_layout"]

        # Resolve dist_lo / dist_hi from the active slider
        if mode == "focus":
            # Single threshold — show only intervals at that exact distance
            dist_val = nearest_value(
                focus_val if focus_val is not None else dist_th_l[0], dist_th_l
            )
            dist_lo = dist_val
            dist_hi = dist_val
        else:
            # Range — overlap all thresholds within [lo, hi]
            if isinstance(overlap_val, (list, tuple)) and len(overlap_val) == 2:
                dist_lo = nearest_value(overlap_val[0], dist_th_l)
                dist_hi = nearest_value(overlap_val[1], dist_th_l)
            else:
                dist_lo = dist_th_l[0]
                dist_hi = dist_th_l[-1]

        active_strand = strand if has_strand else None
        df_filtered = obtain_filtered(seq_id, tuple(tip_order), dist_hi, dist_lo, active_strand)
        bin_edges, colors = get_binned_colors(dist_th_l, colorscheme)

        # Calculate bounds and line width
        n_tips = len(tip_order)
        seq_len = get_seq_len(df_filtered)
        if seq_len is None or seq_len <= 0:
            seq_len = get_seq_len(
                obtain_filtered(seq_id, tuple(tip_order), dist_th_l[-1], dist_th_l[0], None)
            )
        if seq_len is None or seq_len <= 0:
            return go.Figure()

        if y_range is not None:
            y_range = enforce_min_span(y_range, ZOOM_MIN_Y_SPAN, bounds=(0, n_tips - 1))
        if x_range is not None:
            x_range = enforce_min_span(x_range, min(ZOOM_MIN_X_SPAN, seq_len), bounds=(0, seq_len))

        x_range_limits = [-AXIS_X_PAD, seq_len + AXIS_X_PAD]
        y_range_limits = [-AXIS_Y_PAD, n_tips - 1 + AXIS_Y_PAD]
        interval_lw = scaled_interval_width(y_range, n_tips)
        tree_lw = scaled_tree_width(y_range, n_tips)

        if tree_view_mode == "cladogram":
            tree_data, _, layout_xy = clado_layout
        else:
            tree_data, _, layout_xy = phylo_layout

        return build_figure(
            tree_data,
            tip_order,
            layout_xy,
            df_filtered,
            bin_edges,
            colors,
            df_annotations=df_annotations,
            query_id=seq_id,
            interval_linewidth=interval_lw,
            tree_linewidth=tree_lw,
            uirevision="query",
            y_range=y_range,
            x_range=x_range,
            x_range_limits=x_range_limits,
            y_range_limits=y_range_limits,
            leaf_distances=leaf_distances,
        )

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
        prevent_initial_call=True,
    )
    def export_visible_data(
        n_clicks, seq_id, focus_val, overlap_val, mode, y_range, x_range, is_pruned, strand
    ):
        """Export currently visible data (genomes and intervals in view)."""
        if seq_id is None:
            return no_update

        if has_pruned_tree and not is_pruned:
            current_data = full_tree_data
            obtain_filtered = obtain_filtered_full
        else:
            current_data = pruned_tree_data if has_pruned_tree else full_tree_data
            obtain_filtered = obtain_filtered_pruned if has_pruned_tree else obtain_filtered_full

        tip_order = current_data["tip_order"]

        # Resolve thresholds the same way as update_figure
        if mode == "focus":
            dist_val = nearest_value(
                focus_val if focus_val is not None else dist_th_l[0], dist_th_l
            )
            dist_lo = dist_hi = dist_val
        else:
            if isinstance(overlap_val, (list, tuple)) and len(overlap_val) == 2:
                dist_lo = nearest_value(overlap_val[0], dist_th_l)
                dist_hi = nearest_value(overlap_val[1], dist_th_l)
            else:
                dist_lo = dist_th_l[0]
                dist_hi = dist_th_l[-1]

        active_strand = strand if has_strand else None
        df_filtered = obtain_filtered(seq_id, tuple(tip_order), dist_hi, dist_lo, active_strand)

        # Determine visible genomes (y-axis)
        if y_range is None:
            y_min, y_max = 0, len(tip_order) - 1
            visible_indices = list(range(len(tip_order)))
        else:
            y_min, y_max = sorted(y_range)
            y_min = max(0, int(np.floor(y_min)))
            y_max = min(len(tip_order) - 1, int(np.ceil(y_max)))
            visible_indices = list(range(y_min, y_max + 1))

        visible_genomes = [tip_order[i] for i in visible_indices]
        df_export = df_filtered[df_filtered["REF_ID"].isin(visible_genomes)].copy()

        # Determine visible x-range
        if x_range is None:
            seq_len = get_seq_len(df_filtered)
            if seq_len is None:
                seq_len = get_seq_len(
                    obtain_filtered(seq_id, tuple(tip_order), dist_th_l[-1], dist_th_l[0], None)
                )
            x_min, x_max = 0, seq_len if seq_len else 0
        else:
            x_min, x_max = sorted(x_range)
            x_min = int(np.floor(x_min))
            x_max = int(np.ceil(x_max))

        df_export = df_export[
            (df_export["INTERVAL_END"] >= x_min) & (df_export["INTERVAL_START"] <= x_max)
        ].copy()

        df_export["y_order"] = df_export["REF_ID"].map(
            {name: i for i, name in enumerate(tip_order)}
        )
        df_export = df_export.sort_values(["y_order", "INTERVAL_START"])
        df_export = df_export.drop(columns=["y_order", "y"])

        tsv_content = df_export.to_csv(sep="\t", index=False)
        strand_tag = f"_{strand}" if has_strand else ""
        filename = f"{seq_id}{strand_tag}_y{y_min}-{y_max}_x{x_min}-{x_max}.tsv"
        return dict(content=tsv_content, filename=filename)

    @app.callback(
        Output("download-plot", "data"),
        Input("export-pdf-btn", "n_clicks"),
        State("graph", "figure"),
        State("query-dropdown", "value"),
        State("export-width", "value"),
        State("export-height", "value"),
        prevent_initial_call=True,
    )
    def export_plot(pdf_clicks, figure, seq_id, width, height):
        """Export the current plot view as PDF.

        Scattergl (WebGL) traces render blank in static image export,
        so we convert them to regular Scatter before calling pio.to_image.
        """
        if figure is None:
            return no_update
        w = int(width) if width else 1600
        h = int(height) if height else 900
        fig = go.Figure(figure)
        # Scattergl (WebGL) traces render blank in kaleido's headless exporter.
        # Rebuild the figure with Scattergl replaced by plain Scatter.
        new_traces = []
        for trace in fig.data:
            td = trace.to_plotly_json()
            if td.get("type") == "scattergl":
                td["type"] = "scatter"
            new_traces.append(td)
        export_fig = go.Figure(data=new_traces, layout=fig.layout)
        img_bytes = pio.to_image(export_fig, format="pdf", width=w, height=h, scale=2)
        encoded = base64.b64encode(img_bytes).decode()
        name = seq_id if seq_id else "plot"
        return dict(content=encoded, filename=f"{name}.pdf", base64=True)

    return app


def parse_args():
    p = argparse.ArgumentParser(description=None)
    p.add_argument("--input", "-i", required=True, help="Path to TSV data")
    p.add_argument("--tree", "-t", required=True, help="Path to Newick tree")
    p.add_argument("--query", "-q", default=None, help="Query leaf name for distance calculations")
    p.add_argument(
        "--annotation", "-a", default=None, help="Path to annotation TSV file (optional)"
    )
    p.add_argument("--port", "-p", type=int, default=8080, help="Port (default: 8080)")
    p.add_argument("--host", default="127.0.0.1", help="Host address (default: 127.0.0.1)")
    p.add_argument("--debug", action="store_true", default=False, help="Debug mode")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    app = create_app(args.input, args.tree, query=args.query, annotation_path=args.annotation)
    app.run(debug=args.debug, host=args.host, port=args.port)
