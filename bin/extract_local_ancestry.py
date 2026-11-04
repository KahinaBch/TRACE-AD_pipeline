#!/usr/bin/env python3
"""Wide local-ancestry table at the protein's target variants.

One row per participant; for each target variant four columns:
  <label>_hap0_allele  <label>_hap0_anc  <label>_hap1_allele  <label>_hap1_anc

Alleles come from the phased cohort VCF (Beagle output), ancestry from FLARE's
AN1/AN2. FLARE only calls ancestry at markers above its own min-maf/min-mac
filter, so a rare target (e.g. TREM2 R47H) has no AN1/AN2 of its own: its
ancestry is then taken from the nearest FLARE marker within --max-distance bp
(local-ancestry tracts span tens of kb to Mb). The marker used and its
distance are written to the notes file. A target missing from the phased
data (not genotyped, or not in the reference panel) gets NA.
"""
import argparse
import csv

import pysam


def ancestry_codes(vcf):
    for line in str(vcf.header).splitlines():
        if line.startswith("##ANCESTRY=<"):
            body = line[len("##ANCESTRY=<"):].rstrip(">")
            return {int(code): name for name, code in (tok.split("=") for tok in body.split(","))}
    raise SystemExit("ERROR: no ##ANCESTRY header line in the FLARE output")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phased", required=True, help="phased cohort VCF (indexed)")
    ap.add_argument("--anc", required=True, help="FLARE .anc.vcf.gz")
    ap.add_argument("--targets", required=True)
    ap.add_argument("--protein", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--max-distance", type=int, default=20000)
    ap.add_argument("--out", required=True)
    ap.add_argument("--notes", default="targets_notes.tsv")
    args = ap.parse_args()

    pysam.tabix_index(args.anc, preset="vcf", force=True)
    phased = pysam.VariantFile(args.phased)
    anc = pysam.VariantFile(args.anc)
    codes = ancestry_codes(anc)
    samples = list(anc.header.samples)
    targets = [t for t in csv.DictReader(open(args.targets)) if t["protein"].upper() == args.protein.upper()]

    columns = {s: {} for s in samples}
    header, notes = ["participant_id"], []
    for t in targets:
        label, pos = t["label"], int(t["pos"])
        cols = [f"{label}_hap0_allele", f"{label}_hap0_anc", f"{label}_hap1_allele", f"{label}_hap1_anc"]
        header += cols

        gt_rec = next((r for r in phased.fetch(args.chrom, pos - 1, pos)
                       if r.pos == pos and r.ref == t["ref"] and r.alts == (t["alt"],)), None)
        window = anc.fetch(args.chrom, max(0, pos - 1 - args.max_distance), pos + args.max_distance)
        anc_rec = min(window, key=lambda r: abs(r.pos - pos), default=None)

        if gt_rec is None or anc_rec is None:
            reason = "not in phased data" if gt_rec is None else f"no FLARE marker within {args.max_distance} bp"
            notes.append([label, t["rsid"], pos, "NA", "", "", reason])
            for s in samples:
                columns[s].update(dict.fromkeys(cols, "NA"))
            continue

        n_carrier_haps = 0
        for s in samples:
            alleles = [gt_rec.alleles[i] if i is not None else "NA" for i in gt_rec.samples[s]["GT"]]
            ancs = [codes[anc_rec.samples[s]["AN1"]], codes[anc_rec.samples[s]["AN2"]]]
            n_carrier_haps += alleles.count(t["alt"])
            columns[s].update(zip(cols, [alleles[0], ancs[0], alleles[1], ancs[1]]))
        how = "own FLARE call" if anc_rec.pos == pos else "nearest FLARE marker"
        notes.append([label, t["rsid"], pos, "OK", anc_rec.pos, anc_rec.pos - pos,
                      f"{how}; {n_carrier_haps} haplotypes carry {t['alt']}"])

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for s in samples:
            w.writerow([s] + [columns[s][c] for c in header[1:]])
    with open(args.notes, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["label", "rsid", "pos", "status", "ancestry_marker_pos", "distance_bp", "detail"])
        w.writerows(notes)


if __name__ == "__main__":
    main()
