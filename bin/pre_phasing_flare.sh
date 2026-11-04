#!/bin/bash
set -euo pipefail

## PRE-PHASING (query prep for FLARE)
#
# Subsets a chromosome from a raw query sample VCF, optionally renames its
# contigs to match the reference, and intersects its markers against the
# phased reference panel -- mirroring this project's existing
# pre_phasing.sh (used for the RFMix pipeline), just targeting the
# already-phased FLARE reference instead of the raw RFMix one. Does NOT
# phase anything -- that's phasing_flare.sh, next.
#
# No file paths are hardcoded here; all inputs are arguments.

usage() {
  echo "Usage: bash pre_phasing_flare.sh -sample sample.vcf.gz -chr N -ref ref_phased.vcf.gz -out prefix [-contigs contigs.txt] [-threads 8]"
  echo "  -sample    raw (or whole-genome) query VCF/VCF.gz"
  echo "  -chr       chromosome to extract (e.g. 19)"
  echo "  -ref       phased reference VCF to intersect against (e.g. reference_phasing/ref_19_phased.vcf.gz)"
  echo "  -out       output prefix; writes <prefix>_query_com.vcf.gz and <prefix>_ref_com.vcf.gz"
  echo "  -contigs   optional chr-rename map (old<TAB>new), applied to -sample before intersecting"
  echo "  -threads   default 8"
  exit 1
}

SAMPLE=""; CHR=""; REF=""; OUT=""; CONTIGS=""; THREADS=8

while [[ $# -gt 0 ]]; do
  case "$1" in
    -sample) [[ -z "${2:-}" || "${2:-}" == -* ]] && { echo "ERROR: -sample requires an argument"; usage; }; SAMPLE="$2"; shift 2 ;;
    -chr) [[ -z "${2:-}" || "${2:-}" == -* ]] && { echo "ERROR: -chr requires an argument"; usage; }; CHR="$2"; shift 2 ;;
    -ref) [[ -z "${2:-}" || "${2:-}" == -* ]] && { echo "ERROR: -ref requires an argument"; usage; }; REF="$2"; shift 2 ;;
    -out) [[ -z "${2:-}" || "${2:-}" == -* ]] && { echo "ERROR: -out requires an argument"; usage; }; OUT="$2"; shift 2 ;;
    -contigs) [[ -z "${2:-}" || "${2:-}" == -* ]] && { echo "ERROR: -contigs requires an argument"; usage; }; CONTIGS="$2"; shift 2 ;;
    -threads) THREADS="${2:-}"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown argument: $1"; usage ;;
  esac
done

[[ -z "$SAMPLE" || -z "$CHR" || -z "$REF" || -z "$OUT" ]] && { echo "ERROR: -sample, -chr, -ref and -out are required"; usage; }
[[ -f "$SAMPLE" ]] || { echo "ERROR: sample file not found: $SAMPLE"; exit 1; }
[[ -f "$REF" ]] || { echo "ERROR: reference file not found: $REF"; exit 1; }
[[ -n "$CONTIGS" && ! -f "$CONTIGS" ]] && { echo "ERROR: contigs file not found: $CONTIGS"; exit 1; }

echo "[$(date)] indexing sample (force, to be safe across machines -- same as this project's pre_phasing.sh)"
bcftools index -f --threads "$THREADS" "$SAMPLE"

echo "[$(date)] subsetting chr${CHR} from sample"
WORK="$SAMPLE"
if [[ -n "$CONTIGS" ]]; then
  echo "[$(date)] renaming sample contigs per $CONTIGS"
  bcftools annotate --rename-chrs "$CONTIGS" --threads "$THREADS" -Oz -o "${OUT}_renamed.vcf.gz" "$SAMPLE"
  bcftools index -ft --threads "$THREADS" "${OUT}_renamed.vcf.gz"
  WORK="${OUT}_renamed.vcf.gz"
fi
bcftools view -r "$CHR" --threads "$THREADS" -Oz -o "${OUT}_chr${CHR}.vcf.gz" "$WORK"
bcftools index -ft --threads "$THREADS" "${OUT}_chr${CHR}.vcf.gz"
N_SAMPLE=$(bcftools index -n "${OUT}_chr${CHR}.vcf.gz")
echo "[$(date)] sample chr${CHR} variants: $N_SAMPLE"

echo "[$(date)] intersecting with reference ($REF)"
bcftools isec -n=2 -w1 --threads "$THREADS" -Oz -o "${OUT}_query_com.vcf.gz" "${OUT}_chr${CHR}.vcf.gz" "$REF"
bcftools isec -n=2 -w2 --threads "$THREADS" -Oz -o "${OUT}_ref_com.vcf.gz" "${OUT}_chr${CHR}.vcf.gz" "$REF"
bcftools index -ft --threads "$THREADS" "${OUT}_query_com.vcf.gz"
bcftools index -ft --threads "$THREADS" "${OUT}_ref_com.vcf.gz"
N_COM=$(bcftools index -n "${OUT}_query_com.vcf.gz")
N_REF_COM=$(bcftools index -n "${OUT}_ref_com.vcf.gz")
[[ "$N_COM" == "$N_REF_COM" ]] || { echo "ERROR: intersected counts differ ($N_COM vs $N_REF_COM)"; exit 1; }
echo "[$(date)] shared markers: $N_COM (of $N_SAMPLE sample / $(bcftools index -n "$REF") reference)"
echo "[$(date)] pre-phasing done: ${OUT}_query_com.vcf.gz, ${OUT}_ref_com.vcf.gz"
