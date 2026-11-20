"""Static plots: exact-coordinate interval tracks, empirical CDFs, pair matrices.

Notation follows the Gdiff paper: D-hat is a distance estimate, pi the tail
level, and theta_pi / theta_{1-pi} the lower/upper distance thresholds.
"""

import base64
import heapq
import math
from html import escape
from io import BytesIO
from pathlib import Path
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import FuncFormatter, MaxNLocator
import numpy as np

from .data import display_names

GRAY = "#9aa4ae"
TITLES = {"map": "Outlier intervals", "roll": "Local distance profile",
          "dist": "Genome-wide distances", "samples": "Window distance samples"}


def label(key, names):
    descriptions = {"low": "lower tail", "high": "upper tail"}
    parts = [names[key[0]], *[descriptions.get(s, s) for s in key[1:] if s]]
    return textwrap.fill(" · ".join(parts), 65)


def position_scale(extent):
    end = max(abs(extent[0]), abs(extent[1]))
    return (1e6, "Mb") if end >= 1e6 else (1e3, "kb") if end >= 1e3 else (1, "bp")


def style(ax, grid="y"):
    ax.set_facecolor("white")
    ax.grid(axis=grid, color="#e8edf1", lw=0.6, zorder=0)
    ax.tick_params(length=3, color="#aab4bd", labelsize=9)
    for spine in ax.spines.values():
        spine.set_color("#c7d0d8")
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    try:
        ax.ticklabel_format(style="plain", useOffset=False, axis="both")
    except AttributeError:
        pass  # log or categorical axes keep their own formatter


def thin_style(ax, text):
    ax.set_yticks([])
    ax.set_ylabel(text, fontsize=9, rotation=0, ha="right", va="center", labelpad=14)
    for spine in ax.spines.values():
        spine.set_color("#c7d0d8")


def pack_genes(features, limit=8):
    active, free, lanes, packed = [], [], 0, []
    for feature in sorted(features, key=lambda f: (f.start, f.end)):
        while active and active[0][0] < feature.start:
            _, lane = heapq.heappop(active)
            heapq.heappush(free, lane)
        if free:
            lane = heapq.heappop(free)
        else:
            lane = lanes
            lanes += 1
        if lanes > limit:
            return None
        heapq.heappush(active, (feature.end, lane))
        packed.append((feature, lane))
    return packed


def track_label(key, kind):
    """Minimal y label: plain tail names for outlier intervals, D-hat for windows."""
    if kind == "roll":
        return r"$\hat{D}$"
    return SIDE_LABELS.get(key[1] if len(key) > 1 else "", r"$\hat{D}$")


def draw_order(segments):
    """Highest distance first, so lower distances are drawn on top where they overlap."""
    return sorted(segments, key=lambda s: (not math.isfinite(s[2]),
                                           -s[2] if math.isfinite(s[2]) else 0.0))


def positional_figure(report, index=0):
    data = report.data
    seqdata = data["sequences"][index]
    start, end = seqdata["extent"]
    scale, unit = position_scale((start, end))
    tracks = list(seqdata["tracks"].items())
    overlays = list(seqdata["overlay"].items())
    genes = seqdata.get("genes") or []
    names = display_names([key[0] for key, _ in tracks] + [key[0] for key, _ in overlays])

    finite = [d for _, segments in tracks + overlays for _, _, d in segments if math.isfinite(d)]
    # Fixed scale: the estimator is bounded to (0, 0.5), so colors stay
    # comparable across figures and the legend never looks truncated.
    high = 0.5 if finite else 0.05
    norm, cmap = Normalize(0, high), plt.get_cmap("viridis")

    rows, heads, drawn = [], set(), set()
    for key, segments in tracks:
        if rows and rows[-1][1][0] != key[0]:
            rows.append(("spacer", key, None, None))
        if not rows or rows[-1][0] == "spacer":
            heads.add(len(rows))
        rows.append(("track", key, segments, report.kind))
        for okey, osegments in overlays:
            if okey[0] == key[0] and okey[0] not in drawn:
                rows.append(("track", okey, osegments, data["overlay_kind"]))
                drawn.add(okey[0])
    if genes:
        if rows:
            rows.append(("spacer", None, None, None))
        heads.add(len(rows))
        rows.append(("genes", None, None, None))

    heights = [{"spacer": 0.25, "genes": 0.6}.get(kind, 0.55) for kind, _, _, _ in rows]
    height = max(2.6, 1.35 * sum(heights))
    fig, axes = plt.subplots(len(rows), 1, figsize=(11, height), sharex=True, squeeze=False,
                             gridspec_kw={"height_ratios": heights})
    fig.subplots_adjust(left=0.11, right=0.93, bottom=0.6 / height, top=1 - 0.55 / height,
                        hspace=0.35)
    axes = axes[:, 0]
    formatter = FuncFormatter(lambda x, pos: f"{x / scale:,.0f}")
    for index, (ax, (kind, key, segments, source)) in enumerate(zip(axes, rows)):
        ax.set_xlim(start - 1, end)
        ax.xaxis.set_major_formatter(formatter)
        ax.xaxis.set_major_locator(MaxNLocator(7))
        if kind == "spacer":
            ax.axis("off")
            continue
        if kind == "genes":
            packed = pack_genes(genes)
            lanes = 1
            if packed is not None:
                lanes = max(lane for _, lane in packed) + 1
                for feature, lane in packed:
                    ax.plot([feature.start - 1, feature.end], [-lane, -lane],
                            color="#5f7280", lw=4, solid_capstyle="butt")
                if len(genes) <= 24:
                    for feature, lane in packed:
                        if feature.end - feature.start + 1 >= 0.02 * (end - start):
                            ax.text((feature.start - 1 + feature.end) / 2, -lane + 0.35,
                                    feature.name, fontsize=6.5, ha="center", va="bottom",
                                    color="#465767")
            else:
                for feature in genes:
                    ax.plot([feature.start - 1, feature.end], [0, 0], color="#5f7280", lw=4,
                            alpha=0.55, solid_capstyle="butt")
            ax.set_ylim(-lanes + 0.2, 0.9 if packed is not None and len(genes) <= 24 else 0.6)
            thin_style(ax, "genes")
            continue
        ordered = draw_order(segments)
        spans = [(a - 1, b - (a - 1)) for a, b, _ in ordered]
        colors = [cmap(norm(d)) if math.isfinite(d) else GRAY for _, _, d in ordered]
        # White edges keep the boundaries of overlapping intervals visible.
        ax.broken_barh(spans, (-0.45, 0.9), facecolors=colors, edgecolor="white", linewidth=0.7)
        ax.set_ylim(-0.7, 0.7)
        thin_style(ax, track_label(key, source))
        if index in heads:
            ax.set_title(label(key, names), loc="left", fontsize=10.5, pad=8)
    fig.colorbar(ScalarMappable(norm, cmap), ax=list(axes), pad=0.015,
                 label=r"distance $\hat{D}$")
    fig.suptitle(report.heading or TITLES[report.kind], x=0.11, y=1 - 0.08 / height, ha="left",
                 va="top", fontsize=14)
    fig.supxlabel(f"Position on {seqdata['seq']} ({unit})", y=0.08 / height, fontsize=10)
    return fig


def samples_figure(report):
    data = report.data
    ecdfs = data["ecdfs"]
    names = display_names(name for key in ecdfs for name in key[:2])
    pairs = {}
    for key in ecdfs:
        pairs.setdefault((*sorted(key[:2]), key[2]), []).append(key)
    if not pairs:
        fig, ax = plt.subplots(figsize=(8, 3), layout="constrained")
        ax.text(0.5, 0.5, "No distance estimates", transform=ax.transAxes, ha="center",
                color="#637382")
        ax.set_yticks([])
        style(ax, "x")
        fig.suptitle(report.heading or TITLES["samples"], x=0.04, ha="left", fontsize=15)
        return fig
    fig, axes = plt.subplots(len(pairs), 1, figsize=(10, max(3.2, 2.6 * len(pairs))),
                             sharex=True, squeeze=False, layout="constrained")
    for ax, (pair, keys) in zip(axes[:, 0], pairs.items()):
        style(ax)
        a, b, config = pair
        ax.set_title(textwrap.fill(f"{names[a]}  ↔  {names[b]}"
                                   + (f" · {config}" if config != "gdiff" else ""), 75),
                     loc="left", fontsize=10.5, pad=10)
        for key in sorted(keys):
            values = ecdfs[key]
            color = "#2875a7" if key[0] == a else "#ca7939"
            description = f"{names[key[0]]} → {names[key[1]]}"
            if values.size:
                ax.plot(np.r_[values[0], values], np.r_[0, (np.arange(values.size) + 1) / values.size],
                        drawstyle="steps-post", color=color, lw=1.6, label=description)
            else:
                ax.plot([], [], color=color, label=description)
        if not any(ecdfs[key].size for key in keys):
            ax.text(0.5, 0.5, "No distance estimates", transform=ax.transAxes,
                    ha="center", color="#637382")
        ax.axhline(0.5, color="#b9c3cc", ls="--", lw=0.7)
        present = [ecdfs[key] for key in keys if ecdfs[key].size]
        all_values = np.concatenate(present) if present else np.array([0.0, 1.0])
        ax.set_xlim(0, all_values.max() * 1.025)
        ax.set_ylim(0, 1.025)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
        ax.set_ylabel("Cumulative fraction", fontsize=9)
        ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=8)
    fig.suptitle(report.heading or TITLES["samples"], x=0.04, ha="left", fontsize=15)
    fig.supxlabel(r"window distance $\hat{D}$", fontsize=10)
    return fig


def _matrix_panel(ax, fig, data, names, value_label):
    matrix = data["matrix"]
    n_rows, n_cols = matrix.shape
    cmap = plt.get_cmap("YlGnBu").with_extremes(bad="#edf0f3")
    finite = matrix[np.isfinite(matrix)]
    high = max(0.01, math.ceil((float(finite.max()) if finite.size else 0.01) / 0.05) * 0.05)
    im = ax.imshow(np.ma.masked_invalid(matrix), cmap=cmap, vmin=0, vmax=high,
                   aspect="auto", interpolation="nearest")
    font = max(5, min(10, 12 - max(n_rows, n_cols) // 8))
    ax.set_xticks(range(n_cols), [textwrap.fill(names[n], 24) for n in data["x_names"]],
                  rotation=45 if n_cols > 5 else 0, ha="right" if n_cols > 5 else "center",
                  fontsize=font)
    ax.set_yticks(range(n_rows), [textwrap.fill(names[n], 24) for n in data["y_names"]],
                  fontsize=font)
    ax.tick_params(length=0, pad=8)
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlabel("Genome" if data["symmetric"] else "Genome B", labelpad=12)
    ax.set_ylabel("Genome" if data["symmetric"] else "Genome A", labelpad=12)
    fig.colorbar(im, ax=ax, orientation="horizontal", label=value_label, shrink=0.8, pad=0.1)
    if matrix.size <= 144:
        for (i, j), value in np.ndenumerate(matrix):
            ax.text(j, i, f"{value:.3g}" if np.isfinite(value) else "—", ha="center",
                    va="center", color="white" if np.isfinite(value) and value / high > 0.55
                    else "#334452", fontsize=10)


def distance_figure(report):
    data = report.data
    names = display_names(data["names"])
    value_label = r"genome-wide $\hat{D}$"
    if data["matrix"] is None:
        closest = data["closest"]
        fig, ax = plt.subplots(figsize=(9, max(3.0, 0.32 * len(closest) + 1.6)),
                               layout="constrained")
        if closest:
            labels = [textwrap.fill(f"{a} ↔ {b}", 40) for _, a, b in closest]
            ax.barh(range(len(closest)), [d for d, _, _ in closest], color="#2875a7", height=0.65)
            ax.set_yticks(range(len(closest)), labels, fontsize=7)
            ax.invert_yaxis()
            style(ax, "x")
        else:
            ax.text(0.5, 0.5, "No finite estimates", transform=ax.transAxes, ha="center")
            ax.set_yticks([])
        ax.set(xlabel=r"$\hat{D}$", title="Closest reported comparisons")
    else:
        matrix = data["matrix"]
        fig, ax = plt.subplots(figsize=(max(6, min(12, matrix.shape[1] * 0.55 + 3)),
                                        max(3.6, min(12, matrix.shape[0] * 0.45 + 2))),
                               layout="constrained")
        _matrix_panel(ax, fig, data, names, value_label)
    fig.suptitle(report.heading or TITLES["dist"], x=0.04, ha="left", fontsize=15)
    return fig


def figure(report, index=0):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "text.color": "#263746", "axes.labelcolor": "#465767",
                         "xtick.color": "#637382", "ytick.color": "#637382",
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.fonttype": "none"})
    if report.kind in {"map", "roll"}:
        return positional_figure(report, index)
    return samples_figure(report) if report.kind == "samples" else distance_figure(report)


def _encode(fig):
    buffer = BytesIO()
    fig.savefig(buffer, format="png", dpi=160, facecolor="white")
    return base64.b64encode(buffer.getvalue()).decode("ascii")


VIRIDIS = ["#440154", "#482878", "#3e4989", "#31688e", "#26828e",
           "#1f9e89", "#35b779", "#6ece58", "#b5de2b", "#fde725"]
SIDE_LABELS = {"low": "lower tail", "high": "upper tail",
               "forward": "forward windows", "reverse complement": "reverse windows"}


def track_symbol(key, kind, overlayed):
    symbol = SIDE_LABELS.get(key[1] if len(key) > 1 else "", r"$\hat{D}$")
    return symbol + (f" ({kind})" if overlayed else "")


def viewer_payload(report):
    """JSON for the interactive HTML viewer: exact segments, no rendering decisions."""
    data = report.data
    sequences = []
    for sd in data["sequences"]:
        short = display_names([key[0] for key in sd["tracks"]] + [key[0] for key in sd["overlay"]])
        tracks = [{"key": list(key), "ref": short[key[0]], "kind": report.kind,
                   "segs": [[a, b, d if math.isfinite(d) else None] for a, b, d in segments]}
                  for key, segments in sd["tracks"].items()]
        for key, segments in sd["overlay"].items():
            tracks.append({"key": list(key), "ref": short[key[0]], "kind": data["overlay_kind"],
                           "segs": [[a, b, d if math.isfinite(d) else None] for a, b, d in segments]})
        sequences.append({"seq": sd["seq"], "extent": [sd["extent"][0], sd["extent"][1]],
                          "tracks": tracks,
                          "genes": [[f.start, f.end, f.name] for f in sd["genes"]]})
    return {"kind": report.kind, "vmax": 0.5,  # estimator bound; keeps colors comparable
            "gray": GRAY, "viridis": VIRIDIS, "side_labels": SIDE_LABELS, "sequences": sequences}


VIEWER_JS = r"""
const DATA = JSON.parse(document.getElementById('gvizz-data').textContent);
const $ = id => document.getElementById(id);
const L = 170, R = 980, ROW = 16, GAP = 14;
const state = {i: 0, refs: new Set(), sides: new Set(), genes: true, zoom: null};
const cur = () => DATA.sequences[state.i];
const hex = h => [1, 3, 5].map(i => parseInt(h.slice(i, i + 2), 16));
const STOPS = DATA.viridis.map(hex);
function ramp(t) {
  t = Math.min(1, Math.max(0, t));
  const x = t * (STOPS.length - 1), i = Math.min(STOPS.length - 2, Math.floor(x)), f = x - i;
  const c = STOPS[i].map((v, j) => Math.round(v + f * (STOPS[i + 1][j] - v)));
  return `rgb(${c[0]},${c[1]},${c[2]})`;
}
const fmt = v => {
  const span = cur().extent[1];
  return span >= 1e6 ? (v / 1e6).toLocaleString(undefined, {maximumFractionDigits: 2}) + ' Mb'
       : span >= 1e3 ? (v / 1e3).toLocaleString(undefined, {maximumFractionDigits: 1}) + ' kb'
       : Math.round(v).toLocaleString() + ' bp';
};
const zoomNow = () => state.zoom || [cur().extent[0] - 1, cur().extent[1]];
function setZoom(a, b) {
  const e0 = cur().extent[0] - 1, e1 = cur().extent[1];
  const span = Math.min(Math.max(b - a, 500), e1 - e0);
  a = Math.min(Math.max(a, e0), e1 - span);
  state.zoom = (a <= e0 && a + span >= e1) ? null : [a, a + span];
  render();
}
function resetFilters() {
  const s = cur();
  state.refs = new Set(s.tracks.map(t => t.key[0]));
  state.sides = new Set(s.tracks.map(t => t.key[1]).filter(Boolean));
  state.genes = true; state.zoom = null;
}
function controls() {
  const s = cur();
  const refs = [...new Set(s.tracks.map(t => t.key[0]))];
  const refLabel = Object.fromEntries(s.tracks.map(t => [t.key[0], t.ref]));
  const sides = [...new Set(s.tracks.map(t => t.key[1]).filter(Boolean))];
  $('filters').innerHTML =
    refs.map(r => `<label><input type="checkbox" data-ref="${r}" ${state.refs.has(r) ? 'checked' : ''}> ${refLabel[r]}</label>`).join('') +
    sides.map(sd => `<label><input type="checkbox" data-side="${sd}" ${state.sides.has(sd) ? 'checked' : ''}> ${DATA.side_labels[sd] || sd}</label>`).join('') +
    (s.genes.length ? `<label><input type="checkbox" data-genes="1" ${state.genes ? 'checked' : ''}> genes</label>` : '') +
    `<span class="hint">scroll to zoom · drag to pan · double-click to reset</span>`;
  $('filters').onchange = e => {
    const t = e.target;
    if (t.dataset.ref) t.checked ? state.refs.add(t.dataset.ref) : state.refs.delete(t.dataset.ref);
    if (t.dataset.side) t.checked ? state.sides.add(t.dataset.side) : state.sides.delete(t.dataset.side);
    if (t.dataset.genes) state.genes = t.checked;
    render();
  };
}
function render() {
  const s = cur();
  const [z0, z1] = zoomNow();
  const rows = [];
  const byRef = new Map();
  for (const t of s.tracks) {
    if (!byRef.has(t.key[0])) byRef.set(t.key[0], []);
    byRef.get(t.key[0]).push(t);
  }
  let firstGroup = true;
  for (const [ref, ts] of byRef) {
    const visible = ts.filter(t => state.refs.has(ref) && (!t.key[1] || state.sides.has(t.key[1])));
    if (!visible.length) continue;
    if (!firstGroup) rows.push({gap: true});
    firstGroup = false;
    visible.forEach((t, i) => rows.push({track: t, ref: i === 0}));
  }
  if (state.genes && s.genes.length) { if (rows.length) rows.push({gap: true}); rows.push({genes: true}); }
  const H = 46 + rows.reduce((a, r) => a + (r.gap ? GAP : ROW), 0) + 34;
  const X = v => L + (v - z0) / (z1 - z0) * (R - L);
  let y = 40, body = '';
  for (const r of rows) {
    if (r.gap) { y += GAP; continue; }
    if (r.genes) {
      body += `<text x="8" y="${y + 11}" class="lab">genes</text>`;
      for (const [a, b, name] of s.genes) {
        if (b < z0 || a - 1 > z1) continue;
        const x0 = X(Math.max(a - 1, z0)), x1 = X(Math.min(b, z1));
        body += `<rect x="${x0}" y="${y + 3}" width="${Math.max(1, x1 - x0)}" height="8" fill="#5f7280" stroke="#ffffff" stroke-width="0.7"/>`;
        if (x1 - x0 > 40) body += `<text x="${(x0 + x1) / 2}" y="${y}" class="gene" text-anchor="middle">${name}</text>`;
      }
      y += ROW; continue;
    }
    const t = r.track;
    body += `<text x="8" y="${y + 11}" class="lab">${r.ref ? t.ref : ''}</text>`;
    body += `<text x="128" y="${y + 11}" class="sym">${DATA.side_labels[t.key[1]] || 'D̂'}${t.kind !== DATA.kind ? ' ⁺' : ''}</text>`;
    for (const [a, b, d] of t.segs) {
      if (b < z0 || a - 1 > z1) continue;
      const x0 = X(Math.max(a - 1, z0)), x1 = X(Math.min(b, z1));
      body += `<rect x="${x0}" y="${y + 3}" width="${Math.max(0.6, x1 - x0)}" height="9" stroke="#ffffff" stroke-width="0.7" fill="${d === null ? DATA.gray : ramp(d / DATA.vmax)}"/>`;
    }
    y += ROW;
  }
  const ticks = [...Array(6).keys()].map(i => z0 + i * (z1 - z0) / 5);
  body += `<line x1="${L}" y1="${y + 8}" x2="${R}" y2="${y + 8}" stroke="#c7d0d8"/>`;
  body += ticks.map(v => `<line x1="${X(v)}" y1="${y + 8}" x2="${X(v)}" y2="${y + 13}" stroke="#aab4bd"/>` +
                         `<text x="${X(v)}" y="${y + 26}" class="tick" text-anchor="middle">${fmt(v)}</text>`).join('');
  const svg = $('plot');
  svg.setAttribute('viewBox', `0 0 1000 ${H}`);
  svg.innerHTML = `<defs><linearGradient id="cb">${DATA.viridis.map((c, i) =>
    `<stop offset="${i / (DATA.viridis.length - 1)}" stop-color="${c}"/>`).join('')}</linearGradient></defs>` +
    `<rect x="${R - 190}" y="12" width="150" height="8" fill="url(#cb)" stroke="#c7d0d8" stroke-width="0.5"/>` +
    `<text x="${R - 196}" y="19" class="tick" text-anchor="end">0</text>` +
    `<text x="${R - 34}" y="19" class="tick">D̂ ${DATA.vmax}</text>` +
    body;
  $('seqname').textContent = s.seq;
  $('counter').textContent = `${state.i + 1}/${DATA.sequences.length}`;
  $('prev').disabled = state.i === 0;
  $('next').disabled = state.i === DATA.sequences.length - 1;
  $('zoomlabel').textContent = state.zoom ? `zoom ${fmt(z0)}–${fmt(z1)}` : '';
}
const svgEl = $('plot');
const pxOf = e => {
  const rect = svgEl.getBoundingClientRect();
  return (e.clientX - rect.left) * (1000 / rect.width);
};
svgEl.addEventListener('wheel', e => {
  e.preventDefault();
  const [z0, z1] = zoomNow();
  const c = (pxOf(e) - L) / (R - L) * (z1 - z0) + z0;
  const f = e.deltaY > 0 ? 1.25 : 0.8;
  setZoom(c - (c - z0) * f, c + (z1 - c) * f);
}, {passive: false});
let pan = null;
svgEl.addEventListener('mousedown', e => {
  const [z0, z1] = zoomNow();
  pan = {x: pxOf(e), z0, z1};
  svgEl.style.cursor = 'grabbing';
});
window.addEventListener('mousemove', e => {
  if (!pan) return;
  const d = (pan.x - pxOf(e)) / (R - L) * (pan.z1 - pan.z0);
  const e0 = cur().extent[0] - 1, e1 = cur().extent[1];
  let a = pan.z0 + d, b = pan.z1 + d;
  if (a < e0) { b += e0 - a; a = e0; }
  if (b > e1) { a -= b - e1; b = e1; }
  state.zoom = [a, b];
  render();
});
window.addEventListener('mouseup', () => { pan = null; svgEl.style.cursor = ''; });
svgEl.addEventListener('dblclick', () => { state.zoom = null; render(); });
$('zin').onclick = () => { const [z0, z1] = zoomNow(), c = (z0 + z1) / 2, h = (z1 - z0) * 0.4; setZoom(c - h, c + h); };
$('zout').onclick = () => { const [z0, z1] = zoomNow(), c = (z0 + z1) / 2, h = (z1 - z0) * 0.625; setZoom(c - h, c + h); };
$('zreset').onclick = () => { state.zoom = null; render(); };
$('prev').onclick = () => { if (state.i) { state.i--; resetFilters(); controls(); render(); } };
$('next').onclick = () => { if (state.i < DATA.sequences.length - 1) { state.i++; resetFilters(); controls(); render(); } };
resetFilters(); controls(); render();
"""


def write_report(report, table, output):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.suffix.lower() != ".html":
        fig = figure(report)
        try:
            fig.savefig(output, dpi=160, facecolor="white")
        finally:
            plt.close(fig)
        return
    notes = "".join(f"<li>{escape(note)}</li>" for note in report.notes)
    provenance = escape("\n".join(table.provenance))
    heading = escape(report.heading or TITLES[report.kind])
    details = (f'<details><summary>Plot details and input provenance</summary><ul>{notes}</ul>'
               f'<pre>{provenance or "No provenance comments in input."}</pre></details>')
    head = f'''<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>gvizz · {heading}</title>
<style>
body {{font:15px/1.6 system-ui,sans-serif;color:#263746;background:#f7f9fb;max-width:1180px;margin:36px auto;padding:0 24px}}
header {{margin:0 0 14px}} h1 {{font-size:25px;font-weight:600;letter-spacing:-.4px;margin:5px 0}}
.kind {{font-size:11px;font-weight:700;letter-spacing:1.5px;color:#2875a7}}
.meta {{font-size:13px;color:#687987;overflow-wrap:anywhere;margin:0}}
.toolbar {{display:flex;align-items:center;gap:12px;margin:0 0 8px;font-size:14px}}
.toolbar button {{font-size:15px;line-height:1;color:#2875a7;padding:2px 10px;border:1px solid #dfe6ec;border-radius:6px;background:white;cursor:pointer}}
.toolbar button:disabled {{color:#c3ccd4;cursor:default}}
#seqname {{font-weight:600}} #counter, #zoomlabel {{color:#687987;font-size:12px}}
#filters {{display:flex;flex-wrap:wrap;gap:6px 16px;align-items:center;margin:0 0 6px;font-size:12.5px;color:#465767}}
#filters label {{display:inline-flex;gap:5px;align-items:center}}
#filters .hint {{color:#98a4ae;font-size:11.5px}}
svg {{display:block;width:100%;background:white;border:1px solid #e1e7ec;border-radius:10px;user-select:none}}
.lab {{font:11px system-ui,sans-serif;fill:#465767}} .sym {{font:10px system-ui,sans-serif;fill:#637382}}
.tick {{font:9px system-ui,sans-serif;fill:#637382}} .gene {{font:8px system-ui,sans-serif;fill:#465767}}
figure {{margin:0 0 20px;padding:14px;background:white;border:1px solid #e1e7ec;border-radius:10px}}
img {{display:block;width:100%;height:auto}} details {{font-size:13px;color:#526575;border-top:1px solid #dde4ea;padding-top:14px;margin-top:18px}}
summary {{cursor:pointer;font-weight:600}} li {{margin:6px 0}} pre {{white-space:pre-wrap;overflow-wrap:anywhere;background:white;padding:14px}}
</style>
<header><div class="kind">GVIZZ / {escape(report.kind.upper())}</div>
<h1>{heading}</h1><p class="meta">{escape(table.path.name)} · {escape(report.title)}</p></header>
'''
    if report.kind in {"map", "roll"}:
        import json
        payload = json.dumps(viewer_payload(report), separators=(",", ":")).replace("</", "<\\/")
        output.write_text(head + f'''<div class="toolbar"><button id="prev">◀</button><span id="seqname"></span>
<span id="counter"></span><button id="next">▶</button>
<button id="zin" title="zoom in">+</button><button id="zout" title="zoom out">−</button><button id="zreset" title="reset zoom">reset</button>
<span id="zoomlabel"></span></div>
<div id="filters"></div>
<svg id="plot" role="img" aria-label="interval tracks"></svg>
{details}
<script id="gvizz-data" type="application/json">{payload}</script>
<script>{VIEWER_JS}</script></html>''', encoding="utf-8")
        return
    fig = figure(report)
    try:
        encoded = _encode(fig)
    finally:
        plt.close(fig)
    output.write_text(head + f'<figure><img alt="{heading}" src="data:image/png;base64,{encoded}"></figure>'
                             f'{details}</html>', encoding="utf-8")
