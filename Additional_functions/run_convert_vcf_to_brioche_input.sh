#!/bin/bash
set -euo pipefail

#
# Convert a VCF to a Brioche input table: ID, Sequence, Target.bp, Target.base
#
# ARGUMENTS:
#   --vcfin FILE.vcf[.gz]      Path to input VCF (plain .vcf or .vcf.gz). Records after the '#CHROM  POS' header are used.
#                              Columns used: 3=ID, 4=REF, 5=ALT (first ALT taken if multiallelic).
#   --refgenome ref.fa         Reference FASTA corresponding to the VCF; will create/use ref.fa.fai index.
#   --fragmentsize N           Half-window size (bp) for extraction around POS. Total length = 2*N + 1. Default: 150.
#   --chrommapping file.tsv    (Optional) Two-column TSV mapping VCF chromosome name -> reference chromosome name:
#                                vcfchromname<TAB>genomechromname
#                                e.g.,   1<TAB>Chr1
#   --outdir /path/to/outdir   Output directory for intermediates and results.
#
# BEHAVIOUR:
#   - Extracts CHROM, POS, ID; optionally remaps CHROM via --chrommapping.
#   - Builds a ±FRAGMENTSIZE window around POS, extracts sequence only (no FASTA headers) in input order.
#   - Replaces the central base (FRAGMENTSIZE+1) with 'D'.
#   - Emits:
#       Target.bp   = "[REF/ALT]" (ALT = first allele if multiallelic)
#       Target.base = FRAGMENTSIZE+1 (will extract equal length on either side of marker, in exception case where sequence can't extend fragmentsize location is adjusted automatically and a shorter sequence is returned)
#   - Final table: ${OUTPUTDIR}/results/Brioche_inputfile.tsv
#
# Run example
#   bash run_convert_vcf_to_brioche_input.sh \
#     --vcfin mergedImputed.vcf.gz \
#     --refgenome genome1.fa \
#     --fragmentsize 150 \
#     --chrommapping chrom_mapping.tsv \
#     --outdir ./brioche_out



VCFIN=""
REFGENOME=""
FRAGMENTSIZE="150"
CHROMCHROMMAPPING=""
OUTPUTDIR=""

module load Miniconda3

print_usage() {
  echo "Usage: $0 --vcfin FILE.vcf[.gz] --refgenome ref.fa --fragmentsize N --chrommapping chrompairings.tsv --outdir /path/to/outdir"
}

#Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcfin)        VCFIN="$2"; shift 2;;
    --refgenome)    REFGENOME="$2"; shift 2;;
    --fragmentsize) FRAGMENTSIZE="$2"; shift 2;;
    --chrommapping) CHROMCHROMMAPPING="$2"; shift 2;;
    --outdir)       OUTPUTDIR="$2"; shift 2;;
    -h|--help)      print_usage; exit 0;;
    *) echo "Unknown arg: $1"; print_usage; exit 1;;
  esac
done


if [[ -z "$VCFIN" || -z "$REFGENOME" || -z "$OUTPUTDIR" ]]; then
  echo "Missing required arguments."; print_usage; exit 1
fi

if [[ -z "${CHROMCHROMMAPPING:-}" ]]; then
  echo "No chrom mapping file was given; assuming VCF chrom names match the reference."
fi

### In case first time running, a new conda environment for the R anchoring scripts is required. env yaml in same directory as this file
eval "$(conda shell.bash hook)"
conda_path=$(conda info --base)
conda_base=$(conda info --base)
conda config --set channel_priority flexible

ENV_NAME="brioche-vcf"

if conda env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
  echo "Conda env '$ENV_NAME' exists proceeding with anchoring"
else
  echo "Conda env '$ENV_NAME' not found will begin download of the environment before beginning anchoring"

  mkdir -p "$conda_path"/envs/brioche-vcf

  env_path="$conda_path"/envs/brioche-vcf

  conda env create --file brioche-vcf.yaml --prefix $env_path --yes

  echo "Conda environment created at: $env_path"
  echo "Conda environment brioche-vcf fully installed"

fi

conda activate brioche-vcf

# dirs
mkdir -p "${OUTPUTDIR}/referencegenome" "${OUTPUTDIR}/temp" "${OUTPUTDIR}/results"

# prepare reference genome index + genome size
samtools faidx "${REFGENOME}"                      # creates ${REFGENOME}.fai next to the FASTA
cut -f1,2 "${REFGENOME}.fai" > "${OUTPUTDIR}/referencegenome/ref.genome"

# extract CHROM, POS, ID from VCF and skip header rows 
if [[ "$VCFIN" = *.gz ]]; then
  zcat -- "$VCFIN" | awk '/^#CHROM[[:space:]]+POS/{seen=1; next} seen' | cut -f1,2,3 > "${OUTPUTDIR}/temp/chrom_pos_marker.tsv"
elif [[ "$VCFIN" = *.vcf ]]; then
  cat  -- "$VCFIN" | awk '/^#CHROM[[:space:]]+POS/{seen=1; next} seen' | cut -f1,2,3 > "${OUTPUTDIR}/temp/chrom_pos_marker.tsv"
else
  echo "Error: please provide a .vcf or .vcf.gz file" >&2; exit 1
fi
# chrom_pos_marker.tsv columns: CHROM \t POS \t ID

# Apply chromosome-name mapping if provided 

if [[ -n "${CHROMCHROMMAPPING:-}" ]]; then
  echo "Applying chromosome mapping from: ${CHROMCHROMMAPPING}"
  awk -v OFS='\t' '
    NR==FNR {
      # Skip an optional header that literally says "vcfchromname  genomechromname"
      if (FNR==1 && ($1=="vcfchromname" || $1=="VCFChrom" || $1=="vcf" )) next
      map[$1]=$2; next
    }
    {
      if ($1 in map) $1=map[$1];  # replace CHROM if a mapping exists; leave untouched otherwise
      print
    }' "${CHROMCHROMMAPPING}" "${OUTPUTDIR}/temp/chrom_pos_marker.tsv" \
    > "${OUTPUTDIR}/temp/chrom_pos_marker.mapped.tsv"

  mv "${OUTPUTDIR}/temp/chrom_pos_marker.mapped.tsv" "${OUTPUTDIR}/temp/chrom_pos_marker.tsv"
fi

echo "Wrote: ${OUTPUTDIR}/temp/chrom_pos_marker.tsv  (mapped if mapping provided)"

awk -v OFS='\t' '{print $1, $2-1, $2, $3}' \
  "${OUTPUTDIR}/temp/chrom_pos_marker.tsv" > "${OUTPUTDIR}/temp/markers.point.bed"

bedtools slop \
  -b "${FRAGMENTSIZE}" \
  -g "${OUTPUTDIR}/referencegenome/ref.genome" \
  -i "${OUTPUTDIR}/temp/markers.point.bed" \
  > "${OUTPUTDIR}/temp/markers.pm${FRAGMENTSIZE}.bed"


bedtools getfasta \
  -fi "${REFGENOME}" \
  -bed "${OUTPUTDIR}/temp/markers.pm${FRAGMENTSIZE}.bed" \
  -name -tab \
| cut -f2 > "${OUTPUTDIR}/temp/seq_pm${FRAGMENTSIZE}bp.txt"

echo "Wrote sequences to: ${OUTPUTDIR}/temp/seq_pm${FRAGMENTSIZE}bp.txt"


if [[ "$VCFIN" = *.gz ]]; then
  zcat -- "$VCFIN" | awk '/^#CHROM[[:space:]]+POS/{seen=1; next} seen' | cut -f3 > "${OUTPUTDIR}/temp/ids.txt"
  zcat -- "$VCFIN" | awk '/^#CHROM[[:space:]]+POS/{seen=1; next} seen' | cut -f4 > "${OUTPUTDIR}/temp/ref.txt"
  zcat -- "$VCFIN" | awk '/^#CHROM[[:space:]]+POS/{seen=1; next} seen' | cut -f5 > "${OUTPUTDIR}/temp/alt.txt"
elif [[ "$VCFIN" = *.vcf ]]; then
  cat  -- "$VCFIN" | awk '/^#CHROM[[:space:]]+POS/{seen=1; next} seen' | cut -f3 > "${OUTPUTDIR}/temp/ids.txt"
  cat  -- "$VCFIN" | awk '/^#CHROM[[:space:]]+POS/{seen=1; next} seen' | cut -f4 > "${OUTPUTDIR}/temp/ref.txt"
  cat  -- "$VCFIN" | awk '/^#CHROM[[:space:]]+POS/{seen=1; next} seen' | cut -f5 > "${OUTPUTDIR}/temp/alt.txt"
else
  echo "Error: please provide a .vcf or .vcf.gz file" >&2; exit 1
fi

paste "${OUTPUTDIR}/temp/ref.txt" "${OUTPUTDIR}/temp/alt.txt" \
| awk -F'\t' '{split($2,a,","); printf "[%s/%s]\n",$1,a[1]}' \
> "${OUTPUTDIR}/temp/target_bp.txt"

CENTER=$(( FRAGMENTSIZE + 1 ))

awk -v C="$CENTER" '{
  if (length($0) >= C) {
    printf "%sD%s\n", substr($0,1,C-1), substr($0,C+1)
  } else {
    print $0  
  }
}' "${OUTPUTDIR}/temp/seq_pm${FRAGMENTSIZE}bp.txt" > "${OUTPUTDIR}/temp/seq_D.txt"

awk -v C="$CENTER" '{print C}' "${OUTPUTDIR}/temp/seq_D.txt" > "${OUTPUTDIR}/temp/target_base.txt"

echo -e "ID\tSequence\tTarget.bp\tTarget.base" > "${OUTPUTDIR}/results/Brioche_inputfile.tsv"
paste ${OUTPUTDIR}/temp/ids.txt ${OUTPUTDIR}/temp/seq_D.txt  ${OUTPUTDIR}/temp/target_bp.txt ${OUTPUTDIR}/temp/target_base.txt >> "${OUTPUTDIR}/results/Brioche_inputfile.tsv"