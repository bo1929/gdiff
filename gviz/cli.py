"""Command line for gviz.py: inputs, output path, sequence, annotations. Nothing else."""

import argparse
from pathlib import Path
import sys

from .data import Selection, Table, read_annotations


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Static, offline plots of gdiff map, roll, dist, and sampled distances.")
    parser.add_argument("input", type=Path, nargs="+",
                        help="gdiff TSV or TSV.gz; an optional second positional file adds its "
                             "tracks to the same position axis (map over roll, or roll over map)")
    parser.add_argument("-o", "--output", type=Path,
                        help="HTML (default, next to the input), PNG, or SVG")
    parser.add_argument("--seq", help="exact sequence name for positional plots; without it the "
                                      "HTML report traverses every sequence with arrow links")
    parser.add_argument("--annotations", type=Path,
                        help="gene track: TSV with seq, start, end (plus optional strand and name "
                             "columns) or a GFF3 file; drawn under positional plots")
    args = parser.parse_args(argv)
    if len(args.input) > 2:
        parser.error("at most two input tables (a primary and one overlay)")
    output = args.output or args.input[0].with_name(args.input[0].name + ".html")
    if output.suffix.lower() not in {".html", ".png", ".svg"}:
        parser.error("output must end in .html, .png, or .svg")
    if output.resolve() == args.input[0].resolve():
        parser.error("output must differ from input")
    try:
        from .aggregate import summarize
        from .render import write_report
        table = Table(args.input[0])
        overlay = Table(args.input[1]) if len(args.input) == 2 else None
        if table.kind in {"map", "roll"} and not args.seq and output.suffix.lower() != ".html":
            parser.error("PNG/SVG output needs --seq; the HTML report traverses every "
                         "sequence with arrow links")
        if overlay is not None and {table.kind, overlay.kind} - {"map", "roll"}:
            raise ValueError("overlay inputs must be map or roll outputs")
        if table.kind == "dist" and args.seq:
            raise ValueError("--seq applies to positional or sampled-window output")
        selection = Selection(args.seq)
        annotations = read_annotations(args.annotations) if args.annotations else None
        report = summarize(table, selection, overlay=overlay, annotations=annotations)
        write_report(report, table, output)
    except ImportError as error:
        parser.exit(2, f"gviz: {error}; install dependencies with python -m pip install -r gviz/requirements.txt\n")
    except (ValueError, OSError, EOFError) as error:
        parser.exit(2, f"gviz: {error}\n")
    print(f"{output}: {table.kind}, {report.rows:,}/{report.total:,} rows; {output.stat().st_size / 1024:.0f} KiB")
    for note in report.notes:
        print(note)
    return 0


if __name__ == "__main__":
    sys.exit(main())
