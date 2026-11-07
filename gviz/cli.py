"""Command line for gviz.py."""

import argparse
from pathlib import Path
import sys
import webbrowser

from .data import Selection, Table


def region(value):
    try:
        start, end = map(int, value.split(":"))
        if 1 <= start <= end:
            return start, end
    except ValueError:
        pass
    raise argparse.ArgumentTypeError("region must be START:END with 1 <= START <= END")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Static, offline plots of gdiff map, roll, dist, and sampled distances.")
    parser.add_argument("input", type=Path, help="gdiff TSV or TSV.gz")
    parser.add_argument("-o", "--output", type=Path, help="HTML (default), PNG, or SVG")
    parser.add_argument("--seq", help="exact sequence name; positional plots default to the longest sequence")
    parser.add_argument("--region", type=region, help="1-based inclusive START:END; requires --seq")
    parser.add_argument("--reference", action="append", default=[], metavar="NAME",
                        help="exact reference name (repeatable); dist/samples select pairs involving NAME")
    parser.add_argument("--bins", type=int, default=1000, help="maximum positional bins [1000], range 1–10000")
    parser.add_argument("--open", action="store_true", help="open the output in a browser")
    args = parser.parse_args(argv)
    if not 1 <= args.bins <= 10000:
        parser.error("--bins must be between 1 and 10000")
    if args.region and not args.seq:
        parser.error("--region requires --seq")
    output = args.output or args.input.with_name(args.input.name + ".html")
    if output.suffix.lower() not in {".html", ".png", ".svg"}:
        parser.error("output must end in .html, .png, or .svg")
    if output.resolve() == args.input.resolve():
        parser.error("output must differ from input")
    try:
        from .aggregate import summarize
        from .render import write_report
        table = Table(args.input)
        if table.kind == "dist" and args.seq:
            raise ValueError("--seq and --region apply to positional or sampled-window output")
        selection = Selection(args.seq, args.region, tuple(args.reference))
        report = summarize(table, selection, args.bins)
        write_report(report, table, output)
    except ImportError as error:
        parser.exit(2, f"gviz: {error}; install dependencies with python -m pip install -r gviz/requirements.txt\n")
    except (ValueError, OSError, EOFError) as error:
        parser.exit(2, f"gviz: {error}\n")
    print(f"{output}: {table.kind}, {report.rows:,}/{report.total:,} rows; {output.stat().st_size / 1024:.0f} KiB")
    for note in report.notes:
        print(note)
    if args.open:
        webbrowser.open(output.resolve().as_uri())
    return 0


if __name__ == "__main__":
    sys.exit(main())
