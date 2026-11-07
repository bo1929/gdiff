"""Static plots with aligned reference panels and bounded interval tracks."""

import base64
import heapq
import math
from collections import Counter
from html import escape
from io import BytesIO
from pathlib import Path
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import MaxNLocator
import numpy as np

BLUE, ORANGE, GRAY = "#2875a7", "#ca7939", "#9aa4ae"
TITLES = {"map": "Mapped intervals", "roll": "Local sequence distances",
          "dist": "Genome distance comparison", "samples": "Sampled distance distributions"}


def short(name):
    for suffix in (".gz", ".fasta", ".fastq", ".fna", ".fa", "_contigs"):
        name = name.removesuffix(suffix)
    return name


def display_names(names):
    """Remove common suffixes only when the remaining names are unambiguous."""
    names = set(names)
    counts = Counter(short(name) for name in names)
    return {name: short(name) if counts[short(name)] == 1 else name for name in names}


def label(key, names):
    descriptions = {"low": "closer", "high": "farther"}
    parts = [names[key[0]], *[descriptions.get(s, s) for s in key[1:] if s]]
    return textwrap.fill(" · ".join(parts), 65)


def position_scale(edges):
    end = max(abs(edges[0]), abs(edges[-1]))
    return (1e6, "Mb") if end >= 1e6 else (1e3, "kb") if end >= 1e3 else (1, "bp")


def style(ax, grid="y"):
    ax.set_facecolor("white")
    ax.grid(axis=grid, color="#e8edf1", lw=0.6, zorder=0)
    ax.tick_params(length=3, color="#aab4bd", labelsize=9)
    for spine in ax.spines.values():
        spine.set_color("#c7d0d8")
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.ticklabel_format(style="plain", useOffset=False, axis="both")


def pack_intervals(intervals, limit=12):
    """Greedy lane packing for small tracks; crowded sets use density instead."""
    if intervals is None:
        return None
    active, free, packed = [], [], []
    lanes = 0
    for start, end, value in sorted(intervals, key=lambda item: (item[0], item[1])):
        while active and active[0][0] < start:
            _, lane = heapq.heappop(active)
            heapq.heappush(free, lane)
        if free:
            lane = heapq.heappop(free)
        else:
            lane = lanes
            lanes += 1
        if lanes > limit:
            return None
        heapq.heappush(active, (end, lane))
        packed.append((start, end, value, lane))
    return packed


def positional_figure(report):
    items = list(report.data["series"].items())
    names = display_names(key[0] for key, _ in items)
    height = max(3.4, 2.1 * len(items))
    fig, axes = plt.subplots(len(items), 1, figsize=(11, height), sharex=True, squeeze=False)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.6 / height,
                        top=1 - 0.85 / height, hspace=0.7)
    axes = axes[:, 0]
    edges = report.data["edges"]
    scale, unit = position_scale(edges)
    maxima = [float(np.nanmax(v["maximum"])) for _, v in items
              if "maximum" in v and np.isfinite(v["maximum"]).any()]
    upper = max(0.05, math.ceil(max(maxima, default=0.05) * 1.1 / 0.05) * 0.05)
    for ax, (key, values) in zip(axes, items):
        style(ax, "x" if report.kind == "map" else "y")
        rows, missing = values["rows"], values["missing_rows"]
        heading = label(key, names)
        noun = "interval" if report.kind == "map" else "window"
        details = f"{rows:,} {noun}" + ("s" if rows != 1 else "")
        if report.kind == "map" and values["intervals"] is not None:
            distances = [d for _, _, d in values["intervals"] if np.isfinite(d)]
            if distances:
                details += f" · d {min(distances):.3g}–{max(distances):.3g}"
        if missing:
            details += f" · {missing:,} without distance"
        ax.set_title(heading, loc="left", fontsize=10.5, pad=25)
        ax.text(0, 1.035, details, transform=ax.transAxes, fontsize=8, color="#637382")
        if report.kind == "roll":
            observed = values["observed"]
            x = values["positions"] / scale
            ax.plot(x[observed], values["mean"][observed], color=BLUE, lw=1.3, marker=".", ms=3)
            ax.fill_between(x[observed], values["minimum"][observed], values["maximum"][observed],
                            color=BLUE, alpha=0.14, linewidth=0)
            missing_bins = observed & (values["missing"] > 0)
            ax.plot(x[missing_bins], np.full(missing_bins.sum(), 0.025), "|", color=GRAY,
                    ms=9, transform=ax.get_xaxis_transform(), clip_on=False)
            if not np.isfinite(values["mean"]).any():
                ax.text(0.5, 0.5, "No distance estimates", transform=ax.transAxes,
                        ha="center", color="#637382")
            ax.set_ylim(0, upper)
            ax.set_ylabel("Distance", fontsize=9)
        else:
            color = BLUE if len(key) < 2 or key[1] == "low" else ORANGE if key[1] == "high" else GRAY
            packed = pack_intervals(values["intervals"])
            if packed is not None:
                segments = [[((a - 1) / scale, -lane), (b / scale, -lane)] for a, b, _, lane in packed]
                colors = [color if np.isfinite(d) else GRAY for _, _, d, _ in packed]
                ax.add_collection(LineCollection(segments, colors=colors, linewidths=3.5, capstyle="butt"))
                endpoints = [position / scale for a, b, _, _ in packed for position in (a - 1, b)]
                lane_positions = [-lane for _, _, _, lane in packed for _ in range(2)]
                ax.plot(endpoints, lane_positions, "|", color="#52616e", ms=6, mew=0.7)
                lanes = max((p[3] for p in packed), default=0) + 1
                ax.set_ylim(-lanes + 0.2, 0.9)
                ax.set_yticks([])
                ax.spines["left"].set_visible(False)
            else:
                ax.stairs(values["overlap"], edges / scale, color=color, lw=1.2)
                ax.stairs(values["overlap"], edges / scale, color=color, fill=True, alpha=0.14)
                ax.set_ylim(bottom=0)
                ax.set_ylabel("Mean overlap", fontsize=9)
                ax.text(0.01, 0.92, "Density summary", transform=ax.transAxes,
                        fontsize=8, color="#637382", va="top")
        ax.set_xlim(edges[0] / scale, edges[-1] / scale)
    fig.suptitle(TITLES[report.kind], x=0.08, y=1 - 0.1 / height, ha="left", va="top", fontsize=15)
    fig.supxlabel(f"Position on {report.title} ({unit})", y=0.08 / height, fontsize=10)
    return fig


def samples_figure(report):
    data = report.data
    names = display_names(name for key in data["histograms"] for name in key[:2])
    pairs = {}
    for key, counts in data["histograms"].items():
        pair = (*sorted(key[:2]), key[2])
        pairs.setdefault(pair, []).append((key, counts))
    fig, axes = plt.subplots(len(pairs), 1, figsize=(10, max(3.5, 2.8 * len(pairs))),
                             sharex=True, squeeze=False, layout="constrained")
    for ax, (pair, directions) in zip(axes[:, 0], pairs.items()):
        style(ax)
        a, b, config = pair
        title = f"{names[a]}  ↔  {names[b]}" + (f" · {config}" if config != "gdiff" else "")
        ax.set_title(textwrap.fill(title, 75), loc="left", fontsize=10.5, pad=10)
        for key, counts in directions:
            n, missing = int(counts.sum()), data["missing"][key]
            color = BLUE if key[0] == a else ORANGE
            description = f"{names[key[0]]} → {names[key[1]]}\n{n:,} estimates · {missing:,} missing"
            values = np.r_[0, 0, np.cumsum(counts) / n] if n else np.zeros(len(counts) + 2)
            if n:
                ax.step(np.r_[0, data["edges"]], values, where="post", color=color,
                        lw=1.6, label=description)
            else:
                ax.plot([], [], color=color, label=description)
        if not any(counts.any() for _, counts in directions):
            ax.text(0.5, 0.5, "No distance estimates", transform=ax.transAxes,
                    ha="center", color="#637382")
        ax.axhline(0.5, color="#b9c3cc", ls="--", lw=0.7)
        ax.set_xlim(0, data["edges"][-1] * 1.025)
        ax.set_ylim(0, 1.025)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
        ax.set_ylabel("Cumulative fraction", fontsize=9)
        ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=8)
    fig.suptitle(TITLES["samples"], x=0.04, ha="left", fontsize=15)
    fig.supxlabel("Window distance", fontsize=10)
    return fig


def distance_figure(report):
    data = report.data
    names = display_names(data["names"])
    if data["matrix"] is None:
        fig, (ax, ranked) = plt.subplots(1, 2, figsize=(12, 6), layout="constrained")
        style(ax)
        ax.stairs(data["histogram"], data["edges"], fill=True, color=BLUE, alpha=0.65)
        ax.set(xlabel="Distance", ylabel="Pair rows", title="All reported distances")
        style(ranked, "x")
        closest = data["closest"]
        if closest:
            labels = [textwrap.fill(f"{names[a]} ↔ {names[b]}", 38) for _, a, b in closest]
            ranked.barh(range(len(closest)), [d for d, _, _ in closest], color=BLUE, height=0.65)
            ranked.set_yticks(range(len(closest)), labels, fontsize=7)
            ranked.invert_yaxis()
        else:
            ranked.text(0.5, 0.5, "No finite estimates", transform=ranked.transAxes, ha="center")
            ranked.set_yticks([])
        ranked.set(xlabel="Distance", title="Closest reported comparisons")
    else:
        matrix = data["matrix"]
        n_rows, n_cols = matrix.shape
        fig, ax = plt.subplots(figsize=(max(6, min(12, n_cols * 0.55 + 3)),
                                        max(3.6, min(12, n_rows * 0.45 + 2))), layout="constrained")
        cmap = plt.get_cmap("YlGnBu").with_extremes(bad="#edf0f3")
        finite = matrix[np.isfinite(matrix)]
        maximum = float(finite.max()) if finite.size else 0.01
        high = max(0.01, math.ceil(maximum / 0.05) * 0.05)
        im = ax.imshow(np.ma.masked_invalid(matrix), cmap=cmap, vmin=0, vmax=high,
                       aspect="auto", interpolation="nearest")
        font = max(5, min(10, 12 - max(n_rows, n_cols) // 8))
        ax.set_xticks(range(n_cols), [textwrap.fill(names[n], 24) for n in data["x_names"]],
                      rotation=45 if n_cols > 5 else 0, ha="right" if n_cols > 5 else "center", fontsize=font)
        ax.set_yticks(range(n_rows), [textwrap.fill(names[n], 24) for n in data["y_names"]], fontsize=font)
        ax.tick_params(length=0, pad=8)
        ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=1.5)
        ax.tick_params(which="minor", length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xlabel("Genome" if data["symmetric"] else "Genome B", labelpad=12)
        ax.set_ylabel("Genome" if data["symmetric"] else "Genome A", labelpad=12)
        fig.colorbar(im, ax=ax, orientation="horizontal", label="Genome-wide distance", shrink=0.8, pad=0.1)
        if matrix.size <= 144:
            for (i, j), value in np.ndenumerate(matrix):
                ax.text(j, i, f"{value:.3g}" if np.isfinite(value) else "—", ha="center", va="center",
                        color="white" if np.isfinite(value) and value / high > 0.55 else "#334452", fontsize=10)
    fig.suptitle(TITLES["dist"], x=0.04, ha="left", fontsize=15)
    return fig


def figure(report):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "text.color": "#263746", "axes.labelcolor": "#465767",
                         "xtick.color": "#637382", "ytick.color": "#637382",
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.fonttype": "none"})
    if report.kind in {"map", "roll"}:
        return positional_figure(report)
    return samples_figure(report) if report.kind == "samples" else distance_figure(report)


def write_report(report, table, output):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig = figure(report)
    try:
        if output.suffix.lower() != ".html":
            fig.savefig(output, dpi=160, facecolor="white")
            return
        buffer = BytesIO()
        fig.savefig(buffer, format="png", dpi=160, facecolor="white")
        encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
        notes = "".join(f"<li>{escape(note)}</li>" for note in report.notes)
        provenance = escape("\n".join(table.provenance))
        context = f" · {escape(report.title)}" if report.kind in {"map", "roll"} else ""
        output.write_text(f'''<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>gviz · {escape(TITLES[report.kind])}</title>
<style>
body {{font:15px/1.6 system-ui,sans-serif;color:#263746;background:#f7f9fb;max-width:1180px;margin:36px auto;padding:0 24px}}
header {{margin:0 0 22px}} h1 {{font-size:25px;font-weight:600;letter-spacing:-.4px;margin:5px 0}}
.kind {{font-size:11px;font-weight:700;letter-spacing:1.5px;color:#2875a7}}
.meta {{font-size:13px;color:#687987;overflow-wrap:anywhere;margin:0}}
figure {{margin:0 0 20px;padding:18px;background:white;border:1px solid #e1e7ec;border-radius:10px}}
img {{display:block;width:100%;height:auto}} details {{font-size:13px;color:#526575;border-top:1px solid #dde4ea;padding-top:14px}}
summary {{cursor:pointer;font-weight:600}} li {{margin:6px 0}} pre {{white-space:pre-wrap;overflow-wrap:anywhere;background:white;padding:14px}}
</style>
<header><div class="kind">GVIZ / {escape(report.kind.upper())}</div>
<h1>{escape(TITLES[report.kind])}</h1><p class="meta">{escape(table.path.name)}{context}<br>
{report.rows:,} selected rows · {report.total:,} input rows</p></header>
<figure><img alt="{escape(TITLES[report.kind], quote=True)}" src="data:image/png;base64,{encoded}"></figure>
<details><summary>Plot details and input provenance</summary><ul>{notes}</ul>
<pre>{provenance or 'No provenance comments in input.'}</pre></details></html>''', encoding="utf-8")
    finally:
        plt.close(fig)
