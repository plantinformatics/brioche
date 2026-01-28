#!/bin/bash
set -euo pipefail

#
# Convert a VCF to a Brioche input table: ID, Sequence, Target.bp, Target.base
#
# ARGUMENTS:
#   --vcfin FILE.vcf[.gz]      Path to input VCF (plain .vcf or .vcf.gz). Records after the '#CHROM  POS' header are used.
#                              Columns used: 3=ID, 4=REF, 5=ALT (first ALT taken if multiallelic).
#   --refgenome ref.fa[.gz]    Reference FASTA corresponding to the VCF; will be copied & uncompressed in outdir and indexed.
#   --fragmentsize N           (Deprecated convenience) Half-window for both sides; equivalent to --left N --right N.
#   --left N                   Left flank size (bp) for extraction around POS.
#   --right N                  Right flank size (bp) for extraction around POS.
#   --chrommapping file.tsv    (Optional) Two-column TSV mapping VCF chromosome -> reference chromosome:
#                                vcfchromname<TAB>genomechromname
#                                e.g.,   1<TAB>Chr1
#   --outdir /path/to/outdir   Output directory for intermediates and results.
#
# BEHAVIOUR:
#   - Extracts CHROM, POS, ID; optionally remaps CHROM via --chrommapping.
#   - Builds a window using LEFT/RIGHT around POS, extracts sequence only (no FASTA headers) in input order.
#   - Replaces the central base (LEFT+1) with 'D'.
#   - Emits:
#       Target.bp   = "[REF/ALT]" (ALT = first allele if multiallelic)
#       Target.base = LEFT+1  (if truncated at contig edge, sequence may be shorter)
#   - Final table: ${OUTPUTDIR}/results/Brioche_inputfile.tsv
#
# Example:
#   bash run_convert_vcf_to_brioche_input.sh \
#     --vcfin data.vcf.gz \
#     --refgenome refgenome.fa.gz \
#     --left 100 --right 200 \
#     --chrommapping chrom_mapping.tsv \
#     --outdir ./brioche_out
#

# -------------------
# Defaults & helpers
# -------------------
VCFIN=""
REFGENOME=""
FRAGMENTSIZE=""     # optional legacy knob
LEFT=""
RIGHT=""
CHROMCHROMMAPPING=""
OUTPUTDIR=""

print_usage() {
  echo "Usage: $0 --vcfin FILE.vcf[.gz] --refgenome ref.fa[.gz] [--fragmentsize N | --left L --right R] [--chrommapping map.tsv] --outdir OUTDIR"
}

# -------------------
# Parse args
# -------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcfin)          VCFIN="$2"; shift 2;;
    --refgenome)      REFGENOME="$2"; shift 2;;
    --fragmentsize)   FRAGMENTSIZE="$2"; shift 2;;
    --left)           LEFT="$2"; shift 2;;
    --right)          RIGHT="$2"; shift 2;;
    --chrommapping)   CHROMCHROMMAPPING="$2"; shift 2;;
    --outdir)         OUTPUTDIR="$2"; shift 2;;
    -h|--help)        print_usage; exit 0;;
    *) echo "Unknown arg: $1"; print_usage; exit 1;;
  esac
done

if [[ -z "$VCFIN" || -z "$REFGENOME" || -z "$OUTPUTDIR" ]]; then
  echo "Missing required arguments."; print_usage; exit 1
fi

# Back-compat: --fragmentsize sets both LEFT and RIGHT if not provided
if [[ -n "${FRAGMENTSIZE:-}" ]]; then
  if [[ -z "${LEFT:-}"  ]]; then LEFT="$FRAGMENTSIZE"; fi
  if [[ -z "${RIGHT:-}" ]]; then RIGHT="$FRAGMENTSIZE"; fi
fi

# Defaults if neither was provided
if [[ -z "${LEFT:-}" ]];  then LEFT="150";  fi
if [[ -z "${RIGHT:-}" ]]; then RIGHT="150"; fi

# Validate that LEFT/RIGHT are non-negative integers
if ! [[ "$LEFT" =~ ^[0-9]+$ ]] || ! [[ "$RIGHT" =~ ^[0-9]+$ ]]; then
  echo "ERROR: --left and --right must be non-negative integers." >&2
  exit 1
fi
if [[ "$LEFT" -eq 0 || "$RIGHT" -eq 0 ]]; then
  echo "[INFO] Using zero-length flank on one side (LEFT=$LEFT, RIGHT=$RIGHT)."
fi

if [[ -z "${CHROMCHROMMAPPING:-}" ]]; then
  echo "No chrom mapping file was given; assuming VCF chrom names match the reference."
fi

# -------------------
# Conda env (unchanged)
# -------------------
module load Miniconda3 || true

eval "$(conda shell.bash hook)"
conda_path=$(conda info --base)
conda config --set channel_priority flexible >/dev/null 2>&1 || true

ENV_NAME="brioche-vcf"
if conda env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
  echo "Conda env '$ENV_NAME' exists."
else
  echo "Conda env '$ENV_NAME' not found; creating..."
  mkdir -p "$conda_path"/envs/brioche-vcf
  env_path="$conda_path"/envs/brioche-vcf
  conda env create --file brioche-vcf.yaml --prefix "$env_path" --yes
  echo "Conda environment created at: $env_path"
fi
conda activate "$ENV_NAME"

# -------------------
# Dirs
# -------------------
mkdir -p "${OUTPUTDIR}/referencegenome" "${OUTPUTDIR}/temp" "${OUTPUTDIR}/results"

# -------------------
# Reference: copy & ensure uncompressed
# -------------------
OUTREF_DIR="${OUTPUTDIR}/referencegenome"
LOCAL_REF="${OUTREF_DIR}/reference.fa"

echo "[INFO] Copying reference to working directory and ensuring uncompressed FASTA..."
case "$REFGENOME" in
  *.fa.gz|*.fna.gz|*.fasta.gz|*.fa.bgz|*.fna.bgz|*.fasta.bgz)
    gzip -cd -- "$REFGENOME" > "$LOCAL_REF"
    ;;
  *)
    # Already uncompressed: copy as-is
    cp -f -- "$REFGENOME" "$LOCAL_REF"
    ;;
esac

# Index & lengths
echo "[INFO] Indexing reference..."
samtools faidx "$LOCAL_REF"
cut -f1,2 "${LOCAL_REF}.fai" > "${OUTREF_DIR}/ref.genome"

# -------------------
# Extract CHROM, POS, ID from VCF (skip header)
# -------------------
VCF_BODY="${OUTPUTDIR}/temp/chrom_pos_marker.tsv"
if [[ "$VCFIN" = *.gz ]]; then
  zcat -- "$VCFIN" | awk '/^#CHROM[[:space:]]+POS/{seen=1; next} seen' | cut -f1,2,3 > "$VCF_BODY"
elif [[ "$VCFIN" = *.vcf ]]; then
  awk '/^#CHROM[[:space:]]+POS/{seen=1; next} seen' "$VCFIN" | cut -f1,2,3 > "$VCF_BODY"
else
  echo "Error: please provide a .vcf or .vcf.gz file" >&2; exit 1
fi
# Columns: CHROM \t POS \t ID

# -------------------
# Apply chromosome-name mapping if provided
# -------------------
if [[ -n "${CHROMCHROMMAPPING:-}" ]]; then
  echo "[INFO] Applying chromosome mapping from: ${CHROMCHROMMAPPING}"
  awk -v OFS='\t' '
    NR==FNR {
      # Skip optional header
      if (FNR==1 && ($1=="vcfchromname" || $1=="VCFChrom" || $1=="vcf")) next
      map[$1]=$2; next
    }
    {
      if ($1 in map) $1=map[$1];
      print
    }' "${CHROMCHROMMAPPING}" "$VCF_BODY" > "${VCF_BODY}.mapped"
  mv "${VCF_BODY}.mapped" "$VCF_BODY"
fi

echo "[INFO] Wrote: $VCF_BODY  (mapped if mapping provided)"

# -------------------
# Make 0-based point BED, clamp start = 0
# -------------------
POINT_BED="${OUTPUTDIR}/temp/markers.point.bed"
awk -v OFS='\t' '{
  s=$2-1; if (s<0) s=0;
  print $1, s, $2, $3
}' "$VCF_BODY" > "$POINT_BED"

# -------------------
# Expand windows with bedtools slop using LEFT/RIGHT
# -------------------
SLopBED="${OUTPUTDIR}/temp/markers.L${LEFT}_R${RIGHT}.bed"
bedtools slop \
  -l "$LEFT" \
  -r "$RIGHT" \
  -g "${OUTREF_DIR}/ref.genome" \
  -i "$POINT_BED" \
  > "$SLopBED"

# -------------------
# Extract sequences
# -------------------
SEQ_TXT="${OUTPUTDIR}/temp/seq_L${LEFT}_R${RIGHT}.txt"
bedtools getfasta \
  -fi "$LOCAL_REF" \
  -bed "$SLopBED" \
  -name -tab \
| cut -f2 > "$SEQ_TXT"
echo "[INFO] Wrote sequences to: $SEQ_TXT"

# -------------------
# Pull ID, REF, ALT (first ALT)
# -------------------
IDS_TXT="${OUTPUTDIR}/temp/ids.txt"
REF_TXT="${OUTPUTDIR}/temp/ref.txt"
ALT_TXT="${OUTPUTDIR}/temp/alt.txt"

if [[ "$VCFIN" = *.gz ]]; then
  zcat -- "$VCFIN" | awk '/^#CHROM[[:space:]]+POS/{seen=1; next} seen' | cut -f3 > "$IDS_TXT"
  zcat -- "$VCFIN" | awk '/^#CHROM[[:space:]]+POS/{seen=1; next} seen' | cut -f4 > "$REF_TXT"
  zcat -- "$VCFIN" | awk '/^#CHROM[[:space:]]+POS/{seen=1; next} seen' | cut -f5 > "$ALT_TXT"
else
  awk '/^#CHROM[[:space:]]+POS/{seen=1; next} seen' "$VCFIN" | cut -f3 > "$IDS_TXT"
  awk '/^#CHROM[[:space:]]+POS/{seen=1; next} seen' "$VCFIN" | cut -f4 > "$REF_TXT"
  awk '/^#CHROM[[:space:]]+POS/{seen=1; next} seen' "$VCFIN" | cut -f5 > "$ALT_TXT"
fi

TARGET_BP="${OUTPUTDIR}/temp/target_bp.txt"
paste "$REF_TXT" "$ALT_TXT" \
| awk -F'\t' '{split($2,a,","); printf "[%s/%s]\n",$1,a[1]}' \
> "$TARGET_BP"

# -------------------
# Replace center with D (center = LEFT+1)
# -------------------
CENTER=$(( LEFT + 1 ))
SEQ_D_TXT="${OUTPUTDIR}/temp/seq_D.txt"
awk -v C="$CENTER" '{
  if (length($0) >= C) {
    printf "%sD%s\n", substr($0,1,C-1), substr($0,C+1)
  } else {
    # Edge-truncated sequence: leave as-is
    print $0
  }
}' "$SEQ_TXT" > "$SEQ_D_TXT"

TARGET_BASE="${OUTPUTDIR}/temp/target_base.txt"
awk -v C="$CENTER" '{print C}' "$SEQ_D_TXT" > "$TARGET_BASE"

# -------------------
# Assemble output
# -------------------
OUT_TSV="${OUTPUTDIR}/results/Brioche_inputfile.tsv"
echo -e "ID\tSequence\tTarget.bp\tTarget.base" > "$OUT_TSV"
paste "$IDS_TXT" "$SEQ_D_TXT" "$TARGET_BP" "$TARGET_BASE" >> "$OUT_TSV"

echo "[OK] Wrote: $OUT_TSV"

# -------------------
# Cleanup local reference copy (to save space)
# -------------------
echo "[CLEANUP] Removing working reference copy and index..."
rm -f -- "$LOCAL_REF" "${LOCAL_REF}.fai"

echo "[DONE]"
