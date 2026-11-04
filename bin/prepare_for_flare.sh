#!/bin/bash
set -euo pipefail

## FINAL VERIFICATION (query + reference pair, ready for FLARE)
#
# Verifies the phased query (from phasing_flare.sh) and its matching
# pre-intersected reference subset (from pre_phasing_flare.sh) are
# actually FLARE-ready: phased, complete, same marker count, correct
# contig naming. Does not transform the files further -- by construction
# (pre_phasing_flare.sh intersected them, phasing_flare.sh phased the
# query against that exact reference subset) they should already match;
# this step exists to confirm that, not assume it.
#
# No file paths are hardcoded here; all inputs are arguments.

usage() {
  echo "Usage: bash prepare_for_flare.sh -gt phased_query.vcf.gz -ref ref_com.vcf.gz -chr N"
  echo "  -gt    phased query VCF (phasing_flare.sh's output)"
  echo "  -ref   pre-intersected, already-phased reference VCF (pre_phasing_flare.sh's *_ref_com.vcf.gz)"
  echo "  -chr   chromosome, for the contig-naming check (e.g. 19)"
  exit 1
}

GT=""; REF=""; CHR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -gt) [[ -z "${2:-}" || "${2:-}" == -* ]] && { echo "ERROR: -gt requires an argument"; usage; }; GT="$2"; shift 2 ;;
    -ref) [[ -z "${2:-}" || "${2:-}" == -* ]] && { echo "ERROR: -ref requires an argument"; usage; }; REF="$2"; shift 2 ;;
    -chr) [[ -z "${2:-}" || "${2:-}" == -* ]] && { echo "ERROR: -chr requires an argument"; usage; }; CHR="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown argument: $1"; usage ;;
  esac
done

[[ -z "$GT" || -z "$REF" || -z "$CHR" ]] && { echo "ERROR: -gt, -ref and -chr are required"; usage; }
for f in "$GT" "$REF"; do
  [[ -f "$f" ]] || { echo "ERROR: file not found: $f"; exit 1; }
done

echo "[$(date)] checking gt= ($GT)"
N_GT=$(bcftools index -n "$GT")
UNPHASED_GT=$(bcftools query -f '[%GT\n]' "$GT" | grep -vc '|' || true)
MISSING_GT=$(bcftools query -f '[%GT\n]' "$GT" | grep -c '\.' || true)
echo "  variants: $N_GT   unphased calls: $UNPHASED_GT   missing calls: $MISSING_GT"

echo "[$(date)] checking ref= ($REF)"
N_REF=$(bcftools index -n "$REF")
UNPHASED_REF=$(bcftools query -f '[%GT\n]' "$REF" | grep -vc '|' || true)
MISSING_REF=$(bcftools query -f '[%GT\n]' "$REF" | grep -c '\.' || true)
echo "  variants: $N_REF   unphased calls: $UNPHASED_REF   missing calls: $MISSING_REF"

CONTIG_GT=$(bcftools view -h "$GT" | grep -m1 "^##contig=<ID=${CHR}[,>]" || echo "NOT FOUND")
CONTIG_REF=$(bcftools view -h "$REF" | grep -m1 "^##contig=<ID=${CHR}[,>]" || echo "NOT FOUND")
echo "  gt= contig: $CONTIG_GT"
echo "  ref= contig: $CONTIG_REF"

FAIL=0
[[ "$UNPHASED_GT" == "0" ]] || { echo "FAIL: gt= has unphased genotypes"; FAIL=1; }
[[ "$MISSING_GT" == "0" ]] || { echo "FAIL: gt= has missing genotypes"; FAIL=1; }
[[ "$UNPHASED_REF" == "0" ]] || { echo "FAIL: ref= has unphased genotypes"; FAIL=1; }
[[ "$MISSING_REF" == "0" ]] || { echo "FAIL: ref= has missing genotypes"; FAIL=1; }
[[ "$N_GT" == "$N_REF" ]] || { echo "FAIL: marker counts differ (gt=$N_GT vs ref=$N_REF)"; FAIL=1; }
[[ "$CONTIG_GT" != "NOT FOUND" && "$CONTIG_REF" != "NOT FOUND" ]] || { echo "FAIL: contig ${CHR} missing from gt= or ref= header"; FAIL=1; }

if [[ "$FAIL" == "1" ]]; then
  echo "[$(date)] NOT ready for FLARE -- see FAIL lines above"
  exit 1
fi
echo "[$(date)] ready for FLARE: gt=$GT  ref=$REF  ($N_GT markers)"
