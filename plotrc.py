"""Configuration constants for the plotting application."""

__all__ = [
    # Layout
    "FIG_HEIGHT",
    "FIG_MARGIN",
    "FIG_FONT",
    "FIG_SIZE",
    "PANEL_TREE_WIDTH",
    "PANEL_INTERVAL_WIDTH",
    "PANEL_SPACING",
    # Axes
    "AXIS_LINE_WIDTH",
    "AXIS_TITLE_SIZE",
    "AXIS_TICK_SIZE",
    "AXIS_X_NTICKS",
    "AXIS_GRID_WIDTH",
    "AXIS_Y_PAD",
    "AXIS_X_PAD",
    "AXIS_Y_MAX_TICKS",
    # Tree
    "TREE_LINE_MIN",
    "TREE_LINE_MAX",
    "TREE_ROW_FILL",
    "TREE_X_MARGIN",
    # Intervals
    "INTERVAL_ROW_FILL",
    "INTERVAL_LINE_MIN",
    "INTERVAL_LINE_MAX",
    "INTERVAL_LINE_MAX_RATIO",
    # Colors
    "COLORS",
    "COLORBAR_LEN",
    "COLORBAR_X",
    "COLORBAR_THICKNESS",
    "COLORBAR_TITLE_SIZE",
    "COLORBAR_TICK_SIZE",
    "COLORSCALE_DEFAULT",
    "COLORSCALE_OPTIONS",
    "FEATURE_COLORS",
    "FEATURE_STYLES",
    # Zoom & Cache
    "ZOOM_MIN_Y_SPAN",
    "ZOOM_MIN_X_SPAN",
    "CACHE_SIZE",
    # UI Panel
    "UI_PANEL_PADDING",
    "UI_PANEL_GAP",
    "UI_PANEL_RADIUS",
    "UI_PANEL_MARGIN_BOTTOM",
    "UI_PANEL_SHADOW",
    # UI Labels
    "UI_LABEL_SIZE",
    "UI_LABEL_MARGIN",
    "UI_LABEL_COLOR",
    "UI_TITLE_SIZE",
    "UI_TITLE_COLOR",
    # UI Toggles
    "UI_TOGGLE_PADDING",
    "UI_TOGGLE_SIZE",
    "UI_TOGGLE_BORDER",
    "UI_TOGGLE_COLOR_ACTIVE",
    "UI_TOGGLE_COLOR_INACTIVE",
    # UI Navigation
    "UI_NAV_PADDING",
    "UI_NAV_SIZE",
    "UI_NAV_BORDER",
    # UI Export
    "UI_EXPORT_PADDING",
    "UI_EXPORT_FONT_SIZE",
    "UI_EXPORT_INPUT_WIDTH_W",
    "UI_EXPORT_INPUT_WIDTH_H",
    "UI_EXPORT_W_DEFAULT",
    "UI_EXPORT_H_DEFAULT",
    # UI Controls
    "UI_QUERY_DROPDOWN_WIDTH",
    "UI_SLIDER_MIN_WIDTH",
    "UI_COLORSCALE_DROPDOWN_WIDTH",
    "UI_GRAPH_HEIGHT",
    "UI_CONTAINER_HEIGHT",
    # Annotation
    "ANNOTATION_ROW_HEIGHT",
    "ANNOTATION_TRACK_HEIGHT",
    "ARROW_HEAD_RATIO",
    "ARROW_SIZE_MIN",
    "ARROW_HEIGHT_MAX",
    "ARROW_BODY_HEIGHT_RATIO",
    "ARROW_STEM_WIDTH",
]

# =============================================================================
# LAYOUT
# =============================================================================

FIG_HEIGHT = 700
FIG_MARGIN = dict(l=10, r=100, t=40, b=50)
FIG_FONT = "Segoe UI, Arial, sans-serif"
FIG_SIZE = 13

PANEL_TREE_WIDTH = 0.30
PANEL_INTERVAL_WIDTH = 0.70
PANEL_SPACING = 0.001

# =============================================================================
# AXES
# =============================================================================

AXIS_LINE_WIDTH = 2.2
AXIS_TITLE_SIZE = 20
AXIS_TICK_SIZE = 16
AXIS_X_NTICKS = 12
AXIS_GRID_WIDTH = 2.5
AXIS_Y_PAD = 1.8
AXIS_X_PAD = 0.99
AXIS_Y_MAX_TICKS = 30

# =============================================================================
# TREE VISUALIZATION
# =============================================================================

TREE_LINE_MIN = 0.5
TREE_LINE_MAX = 2.5
TREE_ROW_FILL = 0.20
TREE_X_MARGIN = 1.01

# =============================================================================
# INTERVAL VISUALIZATION
# =============================================================================

INTERVAL_ROW_FILL = 0.65
INTERVAL_LINE_MIN = 0.01
INTERVAL_LINE_MAX = 25
INTERVAL_LINE_MAX_RATIO = 0.85

# =============================================================================
# COLOR SCHEME
# =============================================================================

COLORS = {
    "tree": "#2c3e50",
    "axis": "#34495e",
    "grid": "rgba(0,0,0,0.15)",
    "plot_bg": "#fbfbfb",
    "paper_bg": "#fafafa",
    "panel_bg": "#ecf0f1",
    "button_bg": "#ecf0f1",
    "button_bg_inactive": "#e8f4f8",
    "button_border": "#7f8c8d",
}

COLORBAR_LEN = 0.75
COLORBAR_X = 1.02
COLORBAR_THICKNESS = 24
COLORBAR_TITLE_SIZE = 20
COLORBAR_TICK_SIZE = 16
COLORSCALE_DEFAULT = "viridis"
COLORSCALE_OPTIONS = ["viridis", "plasma", "inferno", "magma", "cividis"]

# =============================================================================
# ZOOM & CACHING
# =============================================================================

ZOOM_MIN_Y_SPAN = 10
ZOOM_MIN_X_SPAN = 10000
CACHE_SIZE = 256

# =============================================================================
# UI STYLING - CONTROL PANEL
# =============================================================================

UI_PANEL_PADDING = "7px 12px"
UI_PANEL_GAP = 10
UI_PANEL_RADIUS = 6
UI_PANEL_MARGIN_BOTTOM = 20
UI_PANEL_SHADOW = "0 1px 4px rgba(0,0,0,0.12)"

# Labels
UI_LABEL_SIZE = 14
UI_LABEL_MARGIN = 6
UI_LABEL_COLOR = "#34495e"
UI_TITLE_SIZE = 18
UI_TITLE_COLOR = "#2c3e50"

# Toggle buttons (phylo/clado, focus/overlap, strand, prune)
UI_TOGGLE_PADDING = "5px 12px"
UI_TOGGLE_SIZE = 14
UI_TOGGLE_BORDER = 1
UI_TOGGLE_COLOR_ACTIVE = "#2c3e50"
UI_TOGGLE_COLOR_INACTIVE = "#6c7a89"

# Navigation buttons (◀ ▶)
UI_NAV_PADDING = "4px 8px"
UI_NAV_SIZE = 14
UI_NAV_BORDER = 1

# Export buttons and inputs
UI_EXPORT_PADDING = "5px 11px"
UI_EXPORT_FONT_SIZE = 13
UI_EXPORT_INPUT_WIDTH_W = 62
UI_EXPORT_INPUT_WIDTH_H = 54
UI_EXPORT_W_DEFAULT = 1600
UI_EXPORT_H_DEFAULT = 900

# Dropdowns and sliders
UI_QUERY_DROPDOWN_WIDTH = 210
UI_SLIDER_MIN_WIDTH = 270
UI_COLORSCALE_DROPDOWN_WIDTH = 115

# Graph sizing
UI_GRAPH_HEIGHT = "78vh"
UI_CONTAINER_HEIGHT = "80vh"

# =============================================================================
# ANNOTATION PANEL
# =============================================================================

ANNOTATION_ROW_HEIGHT = 0.175
ANNOTATION_TRACK_HEIGHT = 0.25
ARROW_HEAD_RATIO = 0.125  # smaller head for better proportions
ARROW_SIZE_MIN = 6
ARROW_HEIGHT_MAX = 0.40  # maximum arrow half-height in track units (0.5 max)
ARROW_BODY_HEIGHT_RATIO = 0.6  # body height relative to track height
ARROW_STEM_WIDTH = 1.5  # width of arrow stem lines

# =============================================================================
# BIOLOGICAL FEATURE COLORS (IGV / UCSC Genome Browser conventions)
# =============================================================================

FEATURE_COLORS = {
    # Protein-coding
    "CDS": "#8B4513",  # saddle brown
    "gene": "#A0522D",  # sienna
    "ann": "#4A7C4E",  # muted green (legacy key)
    # RNA genes
    "rRNA": "#0066CC",  # bright blue
    "tRNA": "#9370DB",  # medium purple
    "ncRNA": "#4682B4",  # steel blue
    "misc_RNA": "#708090",  # slate gray
    "snRNA": "#BA55D3",  # medium orchid
    "snoRNA": "#DDA0DD",  # plum
    # Regulatory
    "regulatory": "#228B22",  # forest green
    "promoter": "#32CD32",  # lime green
    "enhancer": "#90EE90",  # light green
    "terminator": "#FF8C00",  # dark orange
    # Structural
    "exon": "#FF69B4",  # hot pink
    "intron": "#C0C0C0",  # silver
    "UTR": "#FFA07A",  # light salmon
    "five_prime_UTR": "#FFB6C1",  # light pink
    "three_prime_UTR": "#FFA07A",  # light salmon
    # Repetitive / mobile elements
    "repeat": "#DC143C",  # crimson
    "transposon": "#8B008B",  # dark magenta
    "mobile_element": "#9932CC",  # dark orchid
    # Other
    "misc_feature": "#808080",
    "source": "#2F4F4F",
    "region": "#696969",
    "gap": "#000000",
    "default": "#696969",
}

# Per-feature line width and opacity used in create_annotation_traces
FEATURE_STYLES = {
    "CDS": {"line_width": 2.0, "opacity": 0.80},
    "gene": {"line_width": 2.0, "opacity": 0.80},
    "rRNA": {"line_width": 1.8, "opacity": 0.75},
    "tRNA": {"line_width": 1.6, "opacity": 0.75},
    "ncRNA": {"line_width": 1.4, "opacity": 0.75},
    "regulatory": {"line_width": 1.5, "opacity": 0.75},
    "promoter": {"line_width": 1.5, "opacity": 0.75},
    "exon": {"line_width": 1.8, "opacity": 0.75},
    "intron": {"line_width": 1.2, "opacity": 0.75},
    "repeat": {"line_width": 1.6, "opacity": 0.70},
    "transposon": {"line_width": 1.6, "opacity": 0.75},
    # fallback
    "_default": {"line_width": 1.5, "opacity": 0.75},
}
