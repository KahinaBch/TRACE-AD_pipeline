#!/usr/bin/env python3
"""run_report.txt: parameters, tool versions, sanity check, FLARE statistics,
and how each target variant was resolved."""
import argparse
from datetime import date


def grep(path, keys):
    """Lines of a log that start (after whitespace) with one of `keys`."""
    return [l.strip() for l in open(path) if l.strip().startswith(keys)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--protein", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--sanity", required=True)
    ap.add_argument("--notes", required=True)
    ap.add_argument("--beagle-log", required=True)
    ap.add_argument("--flare-log", required=True)
    ap.add_argument("--pipeline-version", default="unknown")
    ap.add_argument("--out", default="run_report.txt")
    args = ap.parse_args()

    sections = [
        f"TRACE-AD local ancestry pipeline -- {args.cohort} / {args.protein} (chr{args.chrom})",
        f"date: {date.today()}   pipeline commit: {args.pipeline_version}",
        f"input: {args.input}",
        f"reference panel: {args.panel}",
        "",
        "NOTE: the global ancestry file is FLARE's average over chromosome "
        f"{args.chrom} only, not a genome-wide estimate.",
        "",
        "== tool versions",
        *grep(args.beagle_log, ("beagle.", "Program", "Copyright")),
        *grep(args.flare_log, ("Program",)),
        "",
        "== sanity check",
        open(args.sanity).read().rstrip(),
        "",
        "== FLARE run",
        *grep(args.flare_log, ("reference samples", "target samples", "markers", "Wallclock")),
        "",
        "== target variants (local ancestry file)",
        open(args.notes).read().rstrip(),
    ]
    open(args.out, "w").write("\n".join(sections) + "\n")


if __name__ == "__main__":
    main()
