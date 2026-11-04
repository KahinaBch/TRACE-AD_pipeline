"""Shared by sanity_check.py and fix_strands.py.

Reads biallelic SNP sites from an indexed VCF and classifies every position
shared with the reference panel -- same categories and logic as
sanity_check/check_strand_flips.py:

  match                 cohort REF/ALT == panel REF/ALT
  simple_swap           same alleles, REF and ALT reversed
  strand_flip           complement of the panel's alleles
  strand_flip_and_swap  complement, reversed
  ambiguous_palindromic A/T or C/G site: swap and flip can't be told apart
  allele_mismatch       a genuinely different variant
"""
import pysam

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
CATEGORIES = ("match", "simple_swap", "strand_flip", "strand_flip_and_swap",
              "ambiguous_palindromic", "allele_mismatch")
RECOVERABLE = ("simple_swap", "strand_flip", "strand_flip_and_swap")


def find_contig(vcf, chrom):
    """Name this VCF uses for chromosome `chrom` ('6' or 'chr6'), or None."""
    contigs = set(pysam.TabixFile(vcf).contigs)
    for name in (str(chrom), f"chr{chrom}"):
        if name in contigs:
            return name
    return None


def biallelic_snps(vcf, contig):
    """Yield (pos, ref, alt) for biallelic SNPs; only the first 5 columns are parsed."""
    for line in pysam.TabixFile(vcf).fetch(contig):
        f = line.split("\t", 5)
        ref, alt = f[3].upper(), f[4].upper()
        if ref in COMPLEMENT and alt in COMPLEMENT:
            yield int(f[1]), ref, alt


def classify(ref_r, alt_r, ref_c, alt_c):
    """Category of a cohort site (ref_c/alt_c) relative to the panel (ref_r/alt_r)."""
    if (ref_c, alt_c) == (ref_r, alt_r):
        return "match"
    is_swap = (ref_c, alt_c) == (alt_r, ref_r)
    cr, ca = COMPLEMENT[ref_r], COMPLEMENT[alt_r]
    is_flip = (ref_c, alt_c) == (cr, ca)
    is_flip_swap = (ref_c, alt_c) == (ca, cr)
    if {ref_r, alt_r} in ({"A", "T"}, {"C", "G"}) and (is_swap or is_flip):
        return "ambiguous_palindromic"
    if is_swap:
        return "simple_swap"
    if is_flip:
        return "strand_flip"
    if is_flip_swap:
        return "strand_flip_and_swap"
    return "allele_mismatch"


def _unique_positions(sites):
    """Drop the second record at an already-seen position (split multiallelics)."""
    last = None
    for pos, ref, alt in sites:
        if pos != last:
            last = pos
            yield pos, ref, alt


def compare(sample_vcf, sample_contig, panel_vcf, panel_contig, keep_positions=()):
    """Sorted-merge join of cohort and panel sites.

    Returns:
      counts      {category: n} over all shared positions, plus 'cohort_sites' (all cohort
                  biallelic SNP positions, shared or not)
      mismatches  {pos: (cohort_ref, cohort_alt, panel_ref, panel_alt, category)} for non-matches
      sample_kept {pos: (ref, alt)} cohort sites found at keep_positions
      panel_kept  {pos: (ref, alt)} panel sites found at keep_positions
    """
    keep = set(keep_positions)
    sample_kept, panel_kept = {}, {}
    n_cohort = [0]

    def tap(sites, kept, counter=None):
        for pos, ref, alt in sites:
            if counter is not None:
                counter[0] += 1
            if pos in keep:
                kept[pos] = (ref, alt)
            yield pos, ref, alt

    s_iter = tap(_unique_positions(biallelic_snps(sample_vcf, sample_contig)), sample_kept, n_cohort)
    r_iter = tap(_unique_positions(biallelic_snps(panel_vcf, panel_contig)), panel_kept)
    counts = dict.fromkeys(CATEGORIES, 0)
    mismatches = {}
    s, r = next(s_iter, None), next(r_iter, None)
    while s is not None and r is not None:
        if s[0] < r[0]:
            s = next(s_iter, None)
        elif s[0] > r[0]:
            r = next(r_iter, None)
        else:
            cat = classify(r[1], r[2], s[1], s[2])
            counts[cat] += 1
            if cat != "match":
                mismatches[s[0]] = (s[1], s[2], r[1], r[2], cat)
            s, r = next(s_iter, None), next(r_iter, None)
    for _ in s_iter:  # drain both streams so every keep_position gets captured
        pass
    for _ in r_iter:
        pass
    counts["cohort_sites"] = n_cohort[0]
    return counts, mismatches, sample_kept, panel_kept
