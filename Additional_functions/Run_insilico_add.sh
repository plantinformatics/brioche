#!/bin/bash

# Multithreaded (SMP) job: must run on one node
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="example_insilicorunning"

# The project ID which this job should run under:
#SBATCH --account="default"

# The partition in which to submit the job:
#SBATCH --partition="batch"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# The total amount of memory in megabytes in the job:
#SBATCH --mem=60GB

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL

# Use this email address:
#SBATCH --mail-user=user123@mail.mail

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=72:0:00


############### Brioche settings ################################

#Build target table for the new brioche format
# this file points to the input markers file for brioche input. 
target=/filepath/briochetargetdesign.tsv

# Path of working project directory where nextflow main and modules is kept 
projectpath=/filepath/brioche/

# Shortname of targets file and its function e.g., is it the targets of a SNPchip, WGS etc
chipname="probename"

# list of the absolute path to genomes to apply insilico as a single column txt file to e.g., /path/genome1.fasta  (can be gz or not) 
genomelist="/filepath/genomes_insilico.txt"

# Path to main.nf in Brioche folder (points to main.nf located in the landing folder of the brioche git 
nextflowpath="/filepath/brioche/main.nf"

# Location of parameter file (TODO: MAYBE: have this as a list so that different genomes can have different optimised params?) 
params="/filepath/brioche/params.config"

# email address for updates (will give an update when each Nextflow run completes/if they fail)
emailaddress="user123@mail.mail"

# whether to do orientation iterative mapping (anchor each marker relative to their first unique reference (required for detection of tri+ allelic states)
doorientation="Yes"

# File path of starting orientation can be prefilled with known orientations or set as unknown. Two column tsv with header ID, Orientation  with nrows=nmarkers in targets file and Orientation values of "Unknown"|"plus"|"minus"
orientationfile="/filepath/brioche/Orientation_file.tsv"

# Whether to do biallelic mapping exclusively. Will force brioche to orient everything to the targets Target.base so that only biallelic states will be found and any trialleles are assumed to be orientation errors in the reference genome and coerced to same orientation as Target.base
dobiallelic="No"


######## Output master folder #######
outputmasterfolder="/filepath/brioche/resultsinsilicobrioche"


##################### Anchoring settings #############

# Absolute path to the location of the raw genotypes file you want anchored with the brioche results. Note, make sure you update the below boolean variables for what this file is. e.g., is it a Dart genotypes file, a vcf file, or the SNPchipgenotypes file 
Genotypesfile="/filepath/Rawgenotypes.tsv"


# Absolute path to anchor Anchoring script (R) and scondary bash script to add reference col
ANCHOR_SCRIPT="/filepath/brioche/Additional_functions/Anchoring_script.R"
ADD_REF_COL_SCRIPT="/filepath/brioche/Additional_functions/add_reference_sample_from_header.sh"

# Booleans for first anchoring stage (R script flags)
# Use lowercase true/false strings
isvcf=false # set to true if the genotypes file is a vcf file
issnp=true # set to true if the genotypes file is in SNPchip genotypes format (e.g., col1=Name [marker name] col2=REF, col3=ALT, cols4-N=names of samples, rows 1-n=markers
isdart=false # set to true if the genotypes file is a DArT RawGenotypes file (Note must be in 1 row format, to allow for faster analysis files are read in as chuncks so the script can't have markers which bleed over multiple rows)
isbiallelic=true # set to true if brioche was run in biallelic mode.

####################### Final Brioche settings ################


# Final (true) reference genome for last Brioche run. This is what you want to have everything mapped to in the end
finalgenomepath="/filepathgenome/"
finalgenome="genome_1.fasta"
finalgenomename="genome_1"
finalgenomeACC="GCF_0000001234.1"






##### Runscript ##


module load Nextflow
module load Miniconda3
module load R
module load git


require() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH." >&2; exit 1; }; }
require nextflow
require Rscript
[[ -s "$genomelist" ]] || { echo "ERROR: genomelist not found or empty: $genomelist" >&2; exit 1; }
[[ -s "$target" ]]     || { echo "ERROR: target TSV not found: $target" >&2; exit 1; }
[[ -s "$params" ]]     || { echo "ERROR: params file not found: $params" >&2; exit 1; }
[[ -s "$ANCHOR_SCRIPT" ]] || { echo "ERROR: Anchoring script not found: $ANCHOR_SCRIPT" >&2; exit 1; }
[[ -r "$Genotypesfile" ]] || { echo "ERROR: Genotypes file not readable: $Genotypesfile" >&2; exit 1; }
[[ -r "$finalgenomepath/$finalgenome" ]] || { echo "ERROR: Final genome FASTA not readable: $finalgenomepath/$finalgenome" >&2; exit 1; }



 # remove extensions in files function (save calling it several times) 
# strip known FASTA/FASTQ extensions (optionally .gz)
strip_ext() {
  # lower-case the extension portion and strip
  local x="$1"
  x="${x%.gz}"
  echo "$x" | sed -E 's/\.(fasta|fa|fna|fastq|fq)$//I'
}

########################### Create folders ###############

mkdir -p "${outputmasterfolder}"/{intermediate_brioche,anchoring,final/brioche,final/anchored_results}


  while IFS= read -r row || [[ -n "${row:-}" ]]; do
    # skip blanks and comments
    [[ -z "${row// }" ]] && continue
    [[ "${row:0:1}" == "#" ]] && continue

    fullgenomename="$(basename "$row")"
    genomename="$(strip_ext "$fullgenomename")"

    outdir="${outputmasterfolder}/intermediate_brioche/${genomename}_insilico"
    mkdir -p "$outdir"
    mkdir "$outdir"/brioche-results
    cp "$params" "$outdir"/brioche-results

    echo "[INFO] Brioche: ${genomename}"
    nextflow run "$nextflowpath" -profile slurm \
      --emailaddress "$emailaddress" \
      --genomefasta "$row" \
      --genomename "$genomename" \
      --probename "$chipname" \
      --targetdesign "$target" \
      --paramfile "$params" \
      -c "$params" \
      --resultsdir "${outdir}/" \
      --workdir "${outdir}" \
      --markercharacter "D" \
      --coverage 70 \
      --pident 90 \
      --maximumhits 10 \
      --mode prod \
      --keepduplicates false

  done < "$genomelist"



############################
# 2) Collect mapping CSVs IN GENOME ORDER (from $genomelist)
############################
echo "[INFO] Collecting mapping CSVs in genome-list order…"

declare -a anchorfiles        # ordered CSV paths for anchoring
declare -a anchor_genomes     # matching genomename per CSV
idx=0

while IFS= read -r row || [[ -n "${row:-}" ]]; do
  # skip blanks and comments
  [[ -z "${row// }" ]] && continue
  [[ "${row:0:1}" == "#" ]] && continue

  fullgenomename="$(basename "$row")"
  genomename="$(strip_ext "$fullgenomename")"

  # the CSV produced by Brioche for this genome
  csv_src="$(find "${outputmasterfolder}/intermediate_brioche/${genomename}_insilico" \
                   -type f -name '*Brioche_all_markers1to1stagingforvcf.csv' -print -quit)"

  if [[ -z "$csv_src" ]]; then
    echo "[WARN] No mapping CSV found for ${genomename}; skipping."
    continue
  fi

  ((idx+=1))
  ord=$(printf '%04d' "$idx")
  dst="${outputmasterfolder}/anchoring/${ord}_${genomename}_Brioche_all_markers1to1stagingforvcf.csv"

  cp -f "$csv_src" "$dst"
  anchorfiles+=("$dst")
  anchor_genomes+=("$genomename")
done < "$genomelist"

if (( ${#anchorfiles[@]} == 0 )); then
  echo "ERROR: No mapping CSVs were copied into anchoring/; aborting." >&2
  exit 1
fi

echo "[INFO] Prepared ${#anchorfiles[@]} mapping files (order preserved)."


#######################
# 3) Iterative anchoring (preserve order)
#######################
prev_vcf="$Genotypesfile"


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

for ((i=0; i<${#anchorfiles[@]}; i++)); do
  targetfile="${anchorfiles[$i]}"
  genomename="${anchor_genomes[$i]}"
  ord=$(printf '%04d' $((i+1)))
  out_vcf="${outputmasterfolder}/anchoring/Building_insilico_${ord}_${genomename}.vcf"

  printf '[INFO] Anchoring (%d/%d): %s\n' "$((i+1))" "${#anchorfiles[@]}" "$genomename"

  if (( i == 0 )); then
    # First round: use the user-provided genotypes (table/VCF/DArT per your flags)
    Rscript "$ANCHOR_SCRIPT" \
      --Rawgenotypes "$prev_vcf" \
      --Briochemappings "$targetfile" \
      --IsVCFraws "$isvcf" \
      --IsSNPchip "$issnp" \
      --isDArTfile "$isdart" \
      --Outputfilename "$out_vcf" \
      --genomename "$genomename"
  else
    # Subsequent rounds: always feed previous VCF
    Rscript "$ANCHOR_SCRIPT" \
      --Rawgenotypes "$prev_vcf" \
      --Briochemappings "$targetfile" \
      --IsVCFraws true \
      --IsSNPchip false \
      --isDArTfile false \
      --Outputfilename "$out_vcf" \
      --genomename "$genomename"
  fi

  # Add reference-sample column and chain to next round
  bash "$ADD_REF_COL_SCRIPT" "$out_vcf" > "${out_vcf%.vcf}_coladd.vcf"
  prev_vcf="${out_vcf%.vcf}_coladd.vcf"
done

last_iter_vcf="$prev_vcf"
echo "[INFO] Iterative anchoring complete. Last VCF: $last_iter_vcf"


############################
# 4) Brioche on final reference genome
############################
final_outdir="${outputmasterfolder}/final/brioche/${finalgenomename}_insilico"
mkdir -p "$final_outdir"

echo "[INFO] Brioche on final reference: $finalgenomename"
    nextflow run "$nextflowpath" -profile slurm \
      --emailaddress "$emailaddress" \
      --genomefasta "${finalgenomepath}/${finalgenome}" \
      --genomename "$finalgenomename" \
      --probename "$chipname" \
      --targetdesign "$target" \
      --paramfile "$params" \
      -c "$params" \
      --resultsdir "${final_outdir}/" \
      --workdir "${final_outdir}" \
      --markercharacter "D" \
      --coverage 70 \
      --pident 90 \
      --maximumhits 10 \
      --mode prod \
      --keepduplicates false

mapfile -t finalmappings < <(
  find "${final_outdir}/brioche-results" -type f \
       -name '*Brioche_all_markers1to1stagingforvcf.csv' -print | sort
)

if (( ${#finalmappings[@]} == 0 )); then
  echo "ERROR: No final mapping CSV produced for final reference." >&2
  exit 1
fi

final_csv="${finalmappings[0]}"
final_vcf="${outputmasterfolder}/final/anchored_results/Building_final_mapped_against_${finalgenomename}.vcf"

echo "[INFO] Final anchoring onto ${finalgenomename}"
Rscript "$ANCHOR_SCRIPT" \
  --Rawgenotypes "$last_iter_vcf" \
  --Briochemappings "$final_csv" \
  --IsVCFraws true \
  --IsSNPchip "$issnp" \
  --isDArTfile "$isdart" \
  --Outputfilename "$final_vcf" \
  --genomename "$finalgenomename" \
  --Accession "$finalgenomeACC"

echo "[DONE] Final VCF: $final_vcf"