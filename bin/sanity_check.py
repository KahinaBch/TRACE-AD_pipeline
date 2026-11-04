#!/usr/bin/env python3
"""Pre-flight checks on the cohort VCF, run before any processing.

Checks, on the protein's chromosome only:
  1. the chromosome is present (named 'N' or 'chrN')
  2. the file carries genotypes (FORMAT/GT) and has no duplicate sample IDs
  3. genome build is GRCh38, from two independent signals:
       - header: the chromosome's ##contig length (differs between GRCh37/38)
       - data: positions/alleles agree with the GRCh38 reference panel
         (~99% of shared positions on GRCh38 data, ~35% by chance on GRCh37)
     The data signal decides; a GRCh37-looking header alone is only a warning.
  4. no strand flips / REF-ALT swaps (they would be silently dropped later --
     how A4 lost rs7412); if found, the pipeline stops and suggests fix_strands
  5. each target variant: present in the cohort and in the panel, alleles as
     listed in assets/targets.csv

Always exits 0 and writes the outcome to the report (STATUS: PASS/FAIL), so the
report is kept even on failure; the workflow reads STATUS (and BUILD /
RECOVERABLE_FLIPS, to suggest liftover or fix_strands) and stops the run.
"""
import argparse
import csv
import sys

import pysam

from strand_utils import CATEGORIES, RECOVERABLE, classify, compare, find_contig

MIN_CONSISTENT = 0.90  # share of shared positions with consistent alleles expected on GRCh38
CONTIG_LENGTH = {      # chromosomes with a published phased panel
    "2":  {242193529: "GRCh38", 243199373: "GRCh37"},
    "6":  {170805979: "GRCh38", 171115067: "GRCh37"},
    "19": {58617616: "GRCh38", 59128983: "GRCh37"},
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--targets", required=True)
    ap.add_argument("--protein", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--ignore-flips", action="store_true")
    ap.add_argument("--out", default="sanity_report.txt")
    args = ap.parse_args()

    targets = [t for t in csv.DictReader(open(args.targets)) if t["protein"].upper() == args.protein.upper()]
    lines, fails, warns = [], [], []
    build = "GRCh38"

    header = pysam.VariantFile(args.vcf).header
    samples = list(header.samples)
    lines.append(f"input VCF: {args.vcf}")
    lines.append(f"samples: {len(samples)}")

    contig = find_contig(args.vcf, args.chrom)
    if contig is None:
        fails.append(f"chromosome {args.chrom} not found (looked for '{args.chrom}' and 'chr{args.chrom}')")
    if "GT" not in header.formats:
        fails.append("no FORMAT/GT field: the VCF carries no genotypes")
    dupes = len(samples) - len(set(samples))
    if dupes:
        fails.append(f"{dupes} duplicate sample IDs")

    if contig is not None:
        # --- genome build, header signal
        length = header.contigs[contig].length if contig in header.contigs else None
        header_build = CONTIG_LENGTH.get(args.chrom, {}).get(length)
        lines.append(f"chromosome {args.chrom} named '{contig}' in input")
        lines.append(f"build from header (contig length {length}): {header_build or 'not determinable'}")

        # --- compare with the GRCh38 reference panel
        panel_contig = find_contig(args.panel, args.chrom)
        counts, _, cohort_at, panel_at = compare(args.vcf, contig, args.panel, panel_contig,
                                                 [int(t["pos"]) for t in targets])
        shared = sum(counts[c] for c in CATEGORIES)
        cohort_sites = counts["cohort_sites"]
        lines.append(f"biallelic SNP positions in input: {cohort_sites}, shared with the GRCh38 "
                     f"reference panel: {shared} ({100 * shared / cohort_sites if cohort_sites else 0:.1f}%)")
        for cat in CATEGORIES:
            pct = 100 * counts[cat] / shared if shared else 0
            lines.append(f"  {cat:22s} {counts[cat]:9d}  ({pct:.2f}%)")

        # --- genome build, data signal (decides)
        consistent = shared - counts["allele_mismatch"]
        if shared == 0 or consistent / shared < MIN_CONSISTENT:
            build = "GRCh37" if header_build == "GRCh37" else "not GRCh38 (most likely GRCh37)"
            fails.append(f"input is {build}: only {100 * consistent / shared if shared else 0:.1f}% of positions "
                         f"shared with the GRCh38 panel have consistent alleles (GRCh38 data: >{100 * MIN_CONSISTENT:.0f}%)")
        elif header_build == "GRCh37":
            warns.append("header contig length is GRCh37's, but positions/alleles match GRCh38: "
                         "treating the input as GRCh38 (header probably not updated after liftover)")
        lines.append(f"BUILD: {build}")

        # --- strand flips (only meaningful once the build is right)
        recoverable = sum(counts[c] for c in RECOVERABLE) if build == "GRCh38" else 0
        lines.append(f"RECOVERABLE_FLIPS: {recoverable}")
        if recoverable:
            msg = (f"{recoverable} strand flips / REF-ALT swaps found; these sites would be dropped "
                   "when matching to the reference panel")
            (warns if args.ignore_flips else fails).append(msg)

        # --- target variants
        lines.append("\ntarget variants:")
        for t in targets:
            pos, label = int(t["pos"]), t["label"]
            exp = (t["ref"], t["alt"])
            if pos not in panel_at:
                status = "ABSENT from reference panel -> will be reported as NA"
                warns.append(f"{label} not in reference panel")
            elif classify(*exp, *panel_at[pos]) != "match":
                status = f"panel has {panel_at[pos][0]}/{panel_at[pos][1]}, targets.csv says {exp[0]}/{exp[1]}"
                fails.append(f"{label}: REF/ALT in targets.csv disagree with the reference panel")
            elif pos not in cohort_at:
                status = "ABSENT from input (not genotyped/imputed, or filtered upstream) -> NA"
                warns.append(f"{label} not in input")
            else:
                cat = classify(*exp, *cohort_at[pos])
                status = "OK" if cat == "match" else f"{cat} in input ({cohort_at[pos][0]}/{cohort_at[pos][1]})"
            lines.append(f"  {label:10s} {t['rsid']:12s} {args.chrom}:{pos} {exp[0]}>{exp[1]}  {status}")

    lines.append("")
    lines += [f"WARNING: {w}" for w in warns]
    lines += [f"FAIL: {f}" for f in fails]
    lines.append(f"STATUS: {'FAIL' if fails else 'PASS'}")
    report = "\n".join(lines) + "\n"
    open(args.out, "w").write(report)
    sys.stdout.write(report)


if __name__ == "__main__":
    main()
