import argparse
import pandas as pd
import numpy as np
from ete3 import Tree
from dash import Dash, dcc, html, Input, Output, State, no_update, ctx
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
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
FIG_FONT_SIZE = 12

# =============================================================================
# PANELS (subplots)
# =============================================================================

PANEL_TREE_WIDTH = 0.30
PANEL_INTERVAL_WIDTH = 0.70
PANEL_SPACING = 0.001
PANEL_Y_DOMAIN = [0.01, 1.0]

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
AXIS_TITLE_FONT_SIZE = 16
AXIS_TICK_FONT_SIZE = 13
AXIS_Y_MAX_TICKS = 30
AXIS_X_NTICKS = 12
AXIS_GRID_WIDTH = 2.5
AXIS_Y_PAD = 0.99
AXIS_X_PAD = 0.99

# =============================================================================
# COLORBAR & COLORSCALE
# =============================================================================

COLORBAR_LEN = 0.75
COLORBAR_X = 1.02
COLORBAR_THICKNESS = 24
COLORBAR_TITLE_FONT_SIZE = 15
COLORBAR_TICK_FONT_SIZE = 13
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

UI_LABEL_FONT_SIZE = 14
UI_LABEL_MARGIN_RIGHT = 6
UI_TITLE_FONT_SIZE = 16
UI_TITLE_MARGIN_BOTTOM = 4

# =============================================================================
# UI: TOGGLE BUTTONS
# =============================================================================

UI_TOGGLE_PADDING = "6px 14px"
UI_TOGGLE_FONT_SIZE = 13
UI_TOGGLE_BORDER_WIDTH = 2
UI_TOGGLE_MARGIN_LEFT = 10

# =============================================================================
# UI: NAV BUTTONS (prev/next query)
# =============================================================================

UI_NAV_PADDING = "4px 8px"
UI_NAV_FONT_SIZE = 14
UI_NAV_BORDER_WIDTH = 1

# =============================================================================
# UI: CONTROLS (dropdowns, slider)
# =============================================================================

UI_QUERY_DROPDOWN_WIDTH = 200
UI_QUERY_MARGIN_RIGHT = 10
UI_SLIDER_MIN_WIDTH = 300
UI_SLIDER_GAP = 3
UI_COLORSCALE_DROPDOWN_WIDTH = 140
UI_COLORSCALE_GAP = 5
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


def filter_data(df, seq_id, tip_order, dist_th=None):
    df_q = df[df["QUERY_ID"] == seq_id].copy()
    df_q = df_q.groupby(
        ["REF_ID", "y", "QUERY_ID", "INTERVAL_START", "INTERVAL_END", "SEQ_LEN"], as_index=False
    )["DIST_TH"].min()

    if dist_th is not None:
        df_q = df_q[df_q["DIST_TH"] <= dist_th]
    df_q = df_q.dropna(subset=["y"]).sort_values("y")
    return df_q


def make_cached_filter(df):
    @lru_cache(maxsize=LRU_CACHE_SIZE)
    def cached(seq_id, tip_order_tuple, dist_th):
        return filter_data(df, seq_id, list(tip_order_tuple), dist_th)

    return cached


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
    if node.is_leaf():
        hover = node.name
    else:
        size = len(node.get_leaves()) + count_offset
        hover = f"{node.name}<br>Subtree size: {size}<br>Branch length: {node.dist:.3f}"
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


def make_colorbar(bin_edges, colors):
    n_bins = len(colors)
    norm_edges = np.linspace(0, 1, n_bins + 1)
    colorscale = [[norm_edges[i + j], colors[i]] for i in range(n_bins) for j in (0, 1)]

    tick_labels = [f"≤{bin_edges[i+1]:.3f}" for i in range(n_bins)]
    dummy_y = np.linspace(0, 1, n_bins)

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
                len=COLORBAR_LEN,
                x=COLORBAR_X,
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


def build_figure(
    tree_data,
    tip_order,
    tree_max_xy,
    df_intervals,
    bin_edges,
    colors,
    interval_linewidth=3.0,
    tree_linewidth=1.5,
    uirevision="base",
    y_range=None,
    x_range=None,
    y_range_limits=None,
    x_range_limits=None,
    leaf_distances=None,
):
    fig = make_subplots(
        rows=1,
        cols=2,
        shared_yaxes=True,
        column_widths=[PANEL_TREE_WIDTH, PANEL_INTERVAL_WIDTH],
        horizontal_spacing=PANEL_SPACING,
    )
    fig.add_trace(make_colorbar(bin_edges, colors), row=1, col=2)
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

    # Interval panel
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

    # Left panel axes
    tree_range_min = -(tree_max_xy[0] * (TREE_X_MARGIN - 1))
    tree_range_max = tree_max_xy[0] * TREE_X_MARGIN
    fig.update_xaxes(
        title=TREE_X_TITLE,
        title_font=dict(size=AXIS_TITLE_FONT_SIZE),
        range=[tree_range_min, tree_range_max],
        fixedrange=True,
        side="bottom",
        showline=True,
        showgrid=False,
        linewidth=AXIS_LINE_WIDTH,
        linecolor=COLORS["axis_line"],
        mirror=True,
        tickfont=dict(size=AXIS_TICK_FONT_SIZE),
        ticks="outside",
        anchor="free",
        position=0.0,
        row=1,
        col=1,
    )
    fig.update_yaxes(
        autorange=y_range is None,
        showline=False,
        zeroline=False,
        domain=PANEL_Y_DOMAIN,
        showgrid=False,
        showticklabels=False,
        minallowed=y_range_limits[0],
        maxallowed=y_range_limits[1],
        row=1,
        col=1,
        range=y_range if y_range else y_range_limits,
    )

    # Right panel axes
    tickvals, ticktext = compute_y_ticks(tip_order, y_range)
    if x_range is None:
        x_range = x_range_limits

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
        mirror=True,
        tickfont=dict(size=AXIS_TICK_FONT_SIZE),
        minallowed=x_range_limits[0],
        maxallowed=x_range_limits[1],
        row=1,
        col=2,
    )
    fig.update_yaxes(
        autorange=y_range is None,
        tickvals=tickvals,
        ticktext=ticktext,
        tickfont=dict(size=AXIS_TICK_FONT_SIZE),
        showgrid=False,
        zeroline=False,
        showline=True,
        linewidth=AXIS_LINE_WIDTH,
        linecolor=COLORS["axis_line"],
        mirror=False,
        domain=PANEL_Y_DOMAIN,
        minallowed=y_range_limits[0],
        maxallowed=y_range_limits[1],
        row=1,
        col=2,
        range=y_range if y_range else y_range_limits,
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


def build_layout(seq_id_l, dist_th_l, has_pruned_tree, initial_prune=True):
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
            dcc.Download(id="download-data"),  # Download component for export
            control_panel(
                [
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
                        style={
                            "display": "flex",
                            "alignItems": "center",
                            "marginRight": UI_QUERY_MARGIN_RIGHT,
                        },
                    ),
                    html.Div(
                        [
                            control_label("Distance threshold:"),
                            html.Div(
                                dcc.Slider(
                                    id="dist-slider",
                                    min=dmin,
                                    max=dmax,
                                    step=None,
                                    marks=create_slider_marks(dist_th_l),
                                    value=dmax,
                                    tooltip={"placement": "bottom", "always_visible": False},
                                    updatemode="mouseup",
                                ),
                                style={"minWidth": UI_SLIDER_MIN_WIDTH, "flexShrink": 1},
                            ),
                        ],
                        style={"display": "flex", "alignItems": "center", "gap": UI_SLIDER_GAP},
                    ),
                    html.Div(
                        [
                            control_label("Color Scheme:"),
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
                        style={
                            "display": "flex",
                            "alignItems": "center",
                            "gap": UI_COLORSCALE_GAP,
                            "marginLeft": UI_COLORSCALE_MARGIN_LEFT,
                        },
                    ),
                    # Tree view toggle buttons (phylogeny/cladogram)
                    html.Div(
                        [
                            html.Button(
                                "Phylogeny",
                                id="phylogeny-btn",
                                n_clicks=0,
                                style=get_toggle_button_style(True, "left"),
                            ),
                            html.Button(
                                "Cladogram",
                                id="cladogram-btn",
                                n_clicks=0,
                                style=get_toggle_button_style(False, "right"),
                            ),
                        ],
                        style={"display": "flex", "gap": 0},
                    ),
                    # Prune/Full tree toggle (only if pruned tree available)
                    (
                        html.Div(
                            [
                                html.Button(
                                    "Pruned",
                                    id="pruned-btn",
                                    n_clicks=0,
                                    style=get_toggle_button_style(initial_prune, "left"),
                                ),
                                html.Button(
                                    "Full Tree",
                                    id="full-tree-btn",
                                    n_clicks=0,
                                    style=get_toggle_button_style(not initial_prune, "right"),
                                ),
                            ],
                            style={
                                "display": "flex" if has_pruned_tree else "none",
                                "gap": 0,
                                "marginLeft": UI_TOGGLE_MARGIN_LEFT,
                            },
                        )
                        if has_pruned_tree
                        else html.Div()
                    ),
                    # Export button (rightmost)
                    html.Button(
                        "Export View",
                        id="export-btn",
                        n_clicks=0,
                        style={
                            "padding": "10px 20px",
                            "fontSize": 14,
                            "cursor": "pointer",
                            "border": f"2px solid {COLORS['button_border']}",
                            "borderRadius": "4px",
                            "backgroundColor": COLORS["button_bg"],
                            "fontWeight": "bold",
                            "color": COLORS["button_text"],
                            "marginLeft": "auto",  # Push to right
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
    if "xaxis2.range[0]" in relayout and "xaxis2.range[1]" in relayout:
        return [relayout["xaxis2.range[0]"], relayout["xaxis2.range[1]"]]
    if "xaxis2.range" in relayout:
        return relayout["xaxis2.range"]
    return no_update


# =============================================================================
# APP FACTORY
# =============================================================================


def create_app(input_path, tree_path, query=None):
    """
    Create the Dash application with support for toggling between full and pruned trees.

    Args:
        input_path: Path to TSV data file
        tree_path: Path to Newick tree file
        query: Optional query leaf name for distance calculations
    """
    df = pd.read_csv(input_path, sep="\t")

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

    obtain_filtered_full = make_cached_filter(full_df)
    obtain_filtered_pruned = (
        make_cached_filter(pruned_df) if has_pruned_tree else obtain_filtered_full
    )

    app = Dash(__name__)
    app.layout = build_layout(seq_id_l, dist_th_l, has_pruned_tree)

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
        Output("graph", "figure"),
        Input("query-dropdown", "value"),
        Input("dist-slider", "value"),
        Input("colorscheme-dropdown", "value"),
        Input("y-range-store", "data"),
        Input("x-range-store", "data"),
        Input("tree-view-store", "data"),
        Input("prune-store", "data") if has_pruned_tree else State("prune-store", "data"),
    )
    def update_figure(
        seq_id, slider_value, colorscheme, y_range, x_range, tree_view_mode, is_pruned
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

        # Filter data
        dist_th = nearest_value(slider_value, dist_th_l)
        df_filtered = obtain_filtered(seq_id, tuple(tip_order), dist_th)
        bin_edges, colors = get_binned_colors(dist_th_l, colorscheme)

        # Calculate bounds and line width
        n_tips = len(tip_order)
        seq_len = get_seq_len(df_filtered)
        assert seq_len is not None and seq_len > 0

        if y_range is not None:
            y_range = enforce_min_span(y_range, ZOOM_MIN_Y_SPAN, bounds=(0, n_tips - 1))
        if x_range is not None:
            x_range = enforce_min_span(x_range, min(ZOOM_MIN_X_SPAN, seq_len), bounds=(0, seq_len))

        x_range_limits = [-AXIS_X_PAD, seq_len + AXIS_X_PAD]
        y_range_limits = [-AXIS_Y_PAD, n_tips - 1 + AXIS_Y_PAD]
        interval_lw = scaled_interval_width(y_range, n_tips)
        tree_lw = scaled_tree_width(y_range, n_tips)

        # Select tree layout (phylogeny vs cladogram)
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
        State("dist-slider", "value"),
        State("y-range-store", "data"),
        State("x-range-store", "data"),
        State("prune-store", "data"),
        prevent_initial_call=True,
    )
    def export_visible_data(n_clicks, seq_id, slider_value, y_range, x_range, is_pruned):
        """Export currently visible data (genomes and intervals in view)."""
        if seq_id is None:
            return no_update

        # Get the appropriate tree data based on prune state
        if has_pruned_tree and not is_pruned:
            current_data = full_tree_data
            obtain_filtered = obtain_filtered_full
        else:
            current_data = pruned_tree_data if has_pruned_tree else full_tree_data
            obtain_filtered = obtain_filtered_pruned if has_pruned_tree else obtain_filtered_full

        tip_order = current_data["tip_order"]

        # Filter data
        dist_th = nearest_value(slider_value, dist_th_l)
        df_filtered = obtain_filtered(seq_id, tuple(tip_order), dist_th)

        # Determine visible genomes (y-axis)
        if y_range is None:
            # Full view - all genomes visible
            visible_indices = list(range(len(tip_order)))
        else:
            # Get visible range
            y_min, y_max = sorted(y_range)
            y_min = max(0, int(np.floor(y_min)))
            y_max = min(len(tip_order) - 1, int(np.ceil(y_max)))
            visible_indices = list(range(y_min, y_max + 1))

        visible_genomes = [tip_order[i] for i in visible_indices]

        # Filter to visible genomes
        df_export = df_filtered[df_filtered["REF_ID"].isin(visible_genomes)].copy()

        # Determine visible x-range (genomic positions)
        if x_range is None:
            # Full view - use sequence length
            seq_len = get_seq_len(df_filtered)
            x_min, x_max = 0, seq_len
        else:
            x_min, x_max = sorted(x_range)
            x_min = int(np.floor(x_min))
            x_max = int(np.ceil(x_max))

        # Filter to intervals that overlap with visible x-range
        df_export = df_export[
            (df_export["INTERVAL_END"] >= x_min) & (df_export["INTERVAL_START"] <= x_max)
        ].copy()

        # Sort by genome order and position
        df_export["y_order"] = df_export["REF_ID"].map(
            {name: i for i, name in enumerate(tip_order)}
        )
        df_export = df_export.sort_values(["y_order", "INTERVAL_START"])
        df_export = df_export.drop(columns=["y_order", "y"])

        # Generate TSV content
        tsv_content = df_export.to_csv(sep="\t", index=False)

        # Create filename with query and ranges
        filename = f"{seq_id}_y{y_min}-{y_max}_x{x_min}-{x_max}.tsv"

        return dict(content=tsv_content, filename=filename)

    return app


def parse_args():
    p = argparse.ArgumentParser(description=None)
    p.add_argument("--input", "-i", required=True, help="Path to TSV data")
    p.add_argument("--tree", "-t", required=True, help="Path to Newick tree")
    p.add_argument("--query", "-q", default=None, help="Query leaf name for distance calculations")
    p.add_argument("--port", "-p", type=int, default=8080, help="Port (default: 8080)")
    p.add_argument("--host", default="127.0.0.1", help="Host address (default: 127.0.0.1)")
    p.add_argument("--debug", action="store_true", default=False, help="Debug mode")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    app = create_app(args.input, args.tree, query=args.query)
    app.run(debug=args.debug, host=args.host, port=args.port)
