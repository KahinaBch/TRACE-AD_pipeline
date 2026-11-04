#!/usr/bin/env python3
"""Correct REF/ALT swaps and strand flips against the reference panel, on one
chromosome -- same logic as sanity_check/fix_strand_flips.py.

  simple_swap           REF/ALT set to the panel's; GT allele codes 0<->1 swapped
  strand_flip           REF/ALT set to the panel's; GT unchanged
  strand_flip_and_swap  REF/ALT set to the panel's; GT allele codes 0<->1 swapped

Everything else (matches, ambiguous A/T-C/G sites, true mismatches, indels,
multiallelics) passes through unchanged. At swapped sites only GT is kept:
other per-genotype fields (DS, GP...) and INFO (AF...) would otherwise still
describe the old ALT allele. The pipeline itself only uses GT.

Output: <out>.vcf.gz (that chromosome only) + .tbi index.
"""
import argparse

import pysam

from strand_utils import RECOVERABLE, compare, find_contig

SWAP_GT = str.maketrans("01", "10")


def fix_line(line, new_ref, new_alt, swap):
    f = line.split("\t")
    f[3], f[4] = new_ref, new_alt
    if swap:
        gt_idx = f[8].split(":").index("GT")
        f[7], f[8] = ".", "GT"
        f[9:] = [s.split(":")[gt_idx].translate(SWAP_GT) for s in f[9:]]
    return "\t".join(f)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--out", required=True, help="output path, ending in .vcf.gz")
    args = ap.parse_args()

    contig = find_contig(args.vcf, args.chrom)
    if contig is None:
        raise SystemExit(f"ERROR: chromosome {args.chrom} not found in {args.vcf}")
    _, mismatches, _, _ = compare(args.vcf, contig, args.panel, find_contig(args.panel, args.chrom))
    fixes = {pos: (c_ref, c_alt, p_ref, p_alt, cat != "strand_flip")
             for pos, (c_ref, c_alt, p_ref, p_alt, cat) in mismatches.items() if cat in RECOVERABLE}
    print(f"correcting {len(fixes)} sites on chromosome {args.chrom}")

    tbx = pysam.TabixFile(args.vcf)
    n_done = 0
    with pysam.BGZFile(args.out, "wb") as out:
        for h in tbx.header:
            out.write((h + "\n").encode())
        for line in tbx.fetch(contig):
            pos_ref_alt = line.split("\t", 5)
            fix = fixes.get(int(pos_ref_alt[1]))
            # only rewrite the exact record that was classified (not a second record at that position)
            if fix and (pos_ref_alt[3].upper(), pos_ref_alt[4].upper()) == fix[:2]:
                line = fix_line(line, fix[2], fix[3], fix[4])
                n_done += 1
            out.write((line + "\n").encode())
    pysam.tabix_index(args.out, preset="vcf", force=True)
    print(f"corrected {n_done} sites -> {args.out}")


if __name__ == "__main__":
    main()
