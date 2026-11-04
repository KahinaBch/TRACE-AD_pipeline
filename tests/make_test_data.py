#!/usr/bin/env python3
"""Generate the synthetic dataset used by `-profile test` (tests/data/).

Users don't need to run this: the generated files are part of the repository.
It only documents how they were made (and regenerates them identically: fixed seed).

  - 1 Mb of chr19 around APOE (44.4-45.4 Mb, GRCh38), ~1500 SNPs plus the two APOE
    target variants of assets/targets.csv (rs429358 T>C, rs7412 C>T)
  - reference panel: 3 ancestries (AFR, EUR, EAS) x 100 individuals, phased,
    chromosome named '19' -- like the real panels
  - cohort: 150 admixed individuals whose haplotypes switch between ancestries
    every few hundred kb, unphased, chromosome named 'chr19' -- like a TOPMed VCF
Allele frequencies differ between ancestries (Balding-Nichols, F=0.15) so that
FLARE has ancestry signal to work with.

Requires pysam:  python3 tests/make_test_data.py
"""
import random
from pathlib import Path

import pysam

random.seed(2026)
OUT = Path(__file__).parent / "data"
CHR19_GRCH38_LENGTH = 58617616
START, END = 44_400_000, 45_400_000
ANCESTRIES = ["AFR", "EUR", "EAS"]
N_REF_PER_ANC, N_COHORT = 100, 150
COHORT_ANC_WEIGHTS = [0.3, 0.6, 0.1]   # AFR, EUR, EAS share of tract starts
SWITCH_RATE = 3e-6                      # ancestry switches per bp (~3 per Mb)
BASES = "ACGT"
TARGETS = {  # pos: (ref, alt, ALT frequency per ancestry AFR/EUR/EAS)
    44908684: ("T", "C", (0.25, 0.15, 0.10)),  # rs429358 (APOE e4)
    44908822: ("C", "T", (0.10, 0.08, 0.06)),  # rs7412   (APOE e2)
}


def balding_nichols(p, f=0.15):
    a, b = p * (1 - f) / f, (1 - p) * (1 - f) / f
    return min(max(random.betavariate(a, b), 0.02), 0.98)


positions = sorted(set(random.sample(range(START, END), 1500)) - set(TARGETS))
sites = []
for pos in positions:
    ref = random.choice(BASES)
    alt = random.choice([b for b in BASES if b != ref])
    p = random.uniform(0.05, 0.95)
    sites.append((pos, ref, alt, [balding_nichols(p) for _ in ANCESTRIES]))
for pos, (ref, alt, freqs) in TARGETS.items():
    sites.append((pos, ref, alt, list(freqs)))
sites.sort()


def draw(freq):
    return int(random.random() < freq)


def mosaic_haplotype():
    """Ancestry index at each site: switches along the chromosome."""
    anc, last_pos, path = random.choices(range(3), COHORT_ANC_WEIGHTS)[0], START, []
    for pos, *_ in sites:
        if random.random() < SWITCH_RATE * (pos - last_pos):
            anc = random.choices(range(3), COHORT_ANC_WEIGHTS)[0]
        path.append(anc)
        last_pos = pos
    return path


def write_vcf(path, contig, samples, genotypes, phased):
    sep = "|" if phased else "/"
    lines = ["##fileformat=VCFv4.2",
             f"##contig=<ID={contig},length={CHR19_GRCH38_LENGTH}>",
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]
    for i, (pos, ref, alt, _) in enumerate(sites):
        calls = []
        for h0, h1 in genotypes:
            a, b = h0[i], h1[i]
            if not phased:
                a, b = sorted((a, b))
            calls.append(f"{a}{sep}{b}")
        lines.append(f"{contig}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(calls))
    with pysam.BGZFile(str(path), "wb") as out:
        out.write(("\n".join(lines) + "\n").encode())
    pysam.tabix_index(str(path), preset="vcf", force=True)


OUT.mkdir(exist_ok=True)

ref_samples, ref_gts = [], []
for a, anc in enumerate(ANCESTRIES):
    for k in range(1, N_REF_PER_ANC + 1):
        ref_samples.append(f"REF_{anc}_{k:03d}")
        ref_gts.append(tuple([draw(s[3][a]) for s in sites] for _ in range(2)))
write_vcf(OUT / "ref_19_phased.vcf.gz", "19", ref_samples, ref_gts, phased=True)
with open(OUT / "test_sample_map.txt", "w") as fh:
    fh.writelines(f"{s}\t{s.split('_')[1]}\n" for s in ref_samples)

cohort_samples, cohort_gts = [f"TEST_{k:03d}" for k in range(1, N_COHORT + 1)], []
for _ in cohort_samples:
    haps = []
    for _ in range(2):
        path = mosaic_haplotype()
        haps.append([draw(s[3][anc]) for s, anc in zip(sites, path)])
    cohort_gts.append(tuple(haps))
write_vcf(OUT / "test_cohort.vcf.gz", "chr19", cohort_samples, cohort_gts, phased=False)

print(f"wrote {OUT}: {len(sites)} sites, {len(ref_samples)} reference and {len(cohort_samples)} cohort individuals")
