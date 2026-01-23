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


############# End Brioche settings ################

######## Output master folder #####################

outputmasterfolder="/filepath/brioche/resultsinsilicobrioche"

#############################################################

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

##################End anchoring settings #############

####################### Final Brioche settings ################


# Final (true) reference genome for last Brioche run. This is what you want to have everything mapped to in the end
finalgenomepath="/filepathgenome/"
finalgenome="genome_1.fasta"
finalgenomename="genome_1"
finalgenomeACC="GCF_0000001234.1"



################### End final Brioche settings ###############

################### Modules ##################

module load Nextflow
module load Miniconda3
module load R
module load git


################### End Modules #################


##### Runscript ##




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



#################### Step 1 brioche mapping per reference genome #######################

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
      --maximumhits 10 \
      --mode prod \
      --keepduplicates false

  done < "$genomelist"




############################
# 2 Collect mapping CSVs IN GENOME ORDER (from $genomelist)
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
  # The priors tsv produced by brioche
  tsvpriors_src="$(find "${outputmasterfolder}/intermediate_brioche/${genomename}_insilico" \
                   -type f -name '*priors_informed_strictmapping.tsv' -print -quit)"
  # The dups tsv produced by brioche
  tsvdups_src="$(find "${outputmasterfolder}/intermediate_brioche/${genomename}_insilico" \
                   -type f -name '*marker_localdups_NULLS_counts.tsv' -print -quit)"


  if [[ -z "$csv_src" ]]; then
    echo "[WARN] No mapping CSV found for ${genomename}; skipping."
    continue
  fi

  ((idx+=1))
  ord=$(printf '%04d' "$idx")
  dst="${outputmasterfolder}/anchoring/${ord}_${genomename}_Brioche_all_markers1to1stagingforvcf.csv"
  dst2="${outputmasterfolder}/anchoring/${ord}_${genomename}_priors_informed_strictmapping.tsv"
  dst3="${outputmasterfolder}/anchoring/${ord}_${genomename}_marker_localdups_NULLS_counts.tsv"

  cp -f "$csv_src" "$dst"
  cp -f "$tsvpriors_src" "$dst2"
  cp -f "$tsvdups_src" "$dst3"
  
  anchorfiles+=("$dst")
  priorsfiles+=("$dst2")
  dupsfiles+=("$dst3")
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
  priorsfile="${priorsfiles[$i]}"
  dupfile="${dupsfiles[$i]}"
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
      --genomename "$genomename" \
      --mappriors "$priorsfile" 

    outvcfbasename="$(basename "$out_vcf")"
    outvcffulldata="$out_vcf"
    keep2=$(bcftools query -l "$out_vcf" | head -n 2 | paste -sd, -)

    if [[ -n "${keep2}" ]]; then
      bcftools view -s "$keep2" -Ov -o "${outputmasterfolder}/anchoring/${outvcfbasename}_subset.vcf" "$out_vcf"
    else
      # no samples? just copy
      cp "$out_vcf" "${outputmasterfolder}/anchoring/${outvcfbasename}_subset.vcf"
    fi

    bash "$ADD_REF_COL_SCRIPT" "${outputmasterfolder}/anchoring/${outvcfbasename}_subset.vcf" "$dupfile" > "${out_vcf%.vcf}_coladd.vcf"
    prev_vcf="${out_vcf%.vcf}_coladd.vcf"

  else
    # Subsequent rounds: always feed previous VCF
    Rscript "$ANCHOR_SCRIPT" \
      --Rawgenotypes "$prev_vcf" \
      --Briochemappings "$targetfile" \
      --IsVCFraws true \
      --IsSNPchip false \
      --isDArTfile false \
      --Outputfilename "$out_vcf" \
      --genomename "$genomename" \
      --mappriors "$priorsfile" \
      --droplist ""

    bash "$ADD_REF_COL_SCRIPT" "$out_vcf" "$dupfile" > "${out_vcf%.vcf}_coladd.vcf"
    prev_vcf="${out_vcf%.vcf}_coladd.vcf"

  fi

  # Add reference-sample column and chain to next round
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
      --coverage 80 \
      --pident 95 \
      --maximumhits 10 \
      --mode prod \
      --keepduplicates false

mapfile -t finalmappings < <(
  find "${final_outdir}/brioche-results" -type f \
       -name '*Brioche_all_markers1to1stagingforvcf.csv' -print | sort
)

mapfile -t finalmappings2 < <(
  find "${final_outdir}/brioche-results" -type f \
       -name '*priors_informed_strictmapping.tsv' -print | sort
)

# Pretty sure below is redundant 
mapfile -t finalmappings3 < <(
  find "${final_outdir}/brioche-results" -type f \
       -name '*marker_localdups_NULLS_counts.tsv' -print | sort
)


if (( ${#finalmappings[@]} == 0 )); then
  echo "ERROR: No final mapping CSV produced for final reference." >&2
  exit 1
fi

final_csv="${finalmappings[0]}"
final_tsvpriors="${finalmappings2[0]}"
final_tsvdups="${finalmappings3[0]}" # Pretty sure this is redundant now as we aren't extracting final gt ref
final_vcfrefs="${outputmasterfolder}/final/anchored_results/Building_final_mapped_against_${finalgenomename}_refs.vcf"
final_vcfseqs="${outputmasterfolder}/final/anchored_results/Building_final_mapped_against_${finalgenomename}_seqs.vcf"
final_vcf="${outputmasterfolder}/final/anchored_results/Building_final_mapped_against_${finalgenomename}.vcf"

echo "[INFO] Final anchoring onto ${finalgenomename}"
Rscript "$ANCHOR_SCRIPT" \
  --Rawgenotypes "$last_iter_vcf" \
  --Briochemappings "$final_csv" \
  --IsVCFraws true \
  --IsSNPchip "$issnp" \
  --isDArTfile "$isdart" \
  --Outputfilename "$final_vcfrefs" \
  --genomename "$finalgenomename" \
  --Accession "$finalgenomeACC" \
  --mappriors "$final_tsvpriors" \
  --droplist ""


Rscript "$ANCHOR_SCRIPT" \
  --Rawgenotypes "$outvcffulldata" \
  --Briochemappings "$final_csv" \
  --IsVCFraws true \
  --IsSNPchip "$issnp" \
  --isDArTfile "$isdart" \
  --Outputfilename "$final_vcfseqs" \
  --genomename "$finalgenomename" \
  --Accession "$finalgenomeACC" \
  --mappriors "$final_tsvpriors" \
  --droplist ""



bcftools query -l "$final_vcfseqs" | sort -u > vcf_seq_samplenames.samples.txt
bcftools query -l "$final_vcfrefs" | sort -u > vcf_ref_samplenames.samples.txt
comm -12 vcf_seq_samplenames.samples.txt vcf_ref_samplenames.samples.txt > dup.samples.txt           # overlap
comm -23 vcf_ref_samplenames.samples.txt dup.samples.txt > vcf_ref_samplenames.keep.samples.txt      # unique-to-vcf2

# Keep samples from vcf_ref_samplenames.keep.samples.txt while dropping any other in vcf_ref_samplenames.samples.txt
bcftools view -S vcf_ref_samplenames.keep.samples.txt  -Oz -o vcf2.nodupSamples.vcf.gz "$final_vcfrefs"

bcftools sort -O v -o vcf2.nodupSamples_srt.vcf.gz vcf2.nodupSamples.vcf.gz
bcftools sort -O v -o vcfseqs_srt.vcf "$final_vcfseqs"
bcftools index -c vcf2.nodupSamples_srt.vcf.gz

bgzip vcfseqs_srt.vcf
bcftools index -c vcfseqs_srt.vcf.gz

# Merge samples across identical sites
bcftools merge -m none -Ov -o "$final_vcf".gz vcfseqs_srt.vcf.gz vcf2.nodupSamples_srt.vcf.gz

gunzip "$final_vcf".gz

# Remove excess files 

rm vcf_seq_samplenames.samples.txt vcf_ref_samplenames.samples.txt dup.samples.txt vcf_ref_samplenames.keep.samples.txt vcf2.nodupSamples.vcf.gz vcf2.nodupSamples_srt.vcf.gz vcfseqs_srt.vcf.gz vcfseqs_srt.vcf.gz.csi vcf2.nodupSamples_srt.vcf.gz.csi

echo "[DONE] Final VCF: $final_vcf"

echo "Preparing BCF compatible versions of vcf file for $final_vcf"

final_vcfbasename=$(basename "$final_vcf")

# --------------------------------
# STep 1: If user provided a chrom:numeric mapping file use it
# -------------------------------
used_user_map=0
if [[ -n "${ChromChrommappingfile:-}" && -s "$ChromChrommappingfile" ]]; then
  used_user_map=1

  awk -v MAP="$ChromChrommappingfile" '
    BEGIN{
      FS=OFS="\t"
      while ((getline line < MAP) > 0) {
        gsub(/\r$/, "", line)                # tolerate CRLF
        if (line ~ /^[[:space:]]*$/) continue
        split(line, f, "\t")
        m[f[1]] = f[2]
      }
      close(MAP)
    }

    # contig header: replace ID if mapped (works with "," or ">" so will work with any brioche anchoring script output)
    /^##contig=<ID=/{
      line=$0
      sub(/^##contig=<ID=/, "", line)
      split(line, a, /[,>]/)
      id=a[1]
      if (id in m) {
        sub("##contig=<ID=" id ",", "##contig=<ID=" m[id] ",", $0)
        sub("##contig=<ID=" id ">", "##contig=<ID=" m[id] ">", $0)
      }
      print
      next
    }

    /^#/ { print; next }

    { if ($1 in m) $1=m[$1]; print }
  ' "$final_vcf" > "${outputmasterfolder}/final/anchored_results/${final_vcfbasename}_numericchroms.vcf"
fi


# -------------------------
# Steps 2 and 3: Only run if user did not provide ChromChrommappingfile
# ----------------------
if [[ $used_user_map -eq 0 ]]; then

  # We will generate a mapping file locally
  ChromChrommappingfile="chrommatchingtxt.tsv"


  # STEP 2: pattern matching

  step2_ok=1

  awk -F'\t' '
    !/^#/ { if (!seen[$1]++) print $1 }
  ' "$final_vcf" |
  awk '
    BEGIN { OFS="\t" }
    function die(msg) { print msg > "/dev/stderr"; exit 2 }

    {
      orig=$0
      o=tolower(orig)

      # Normalize to underscore-separated tokens
      s=o
      gsub(/[^0-9a-z]+/, "_", s)
      sub(/^_+/, "", s); sub(/_+$/, "", s)

      # Unknown/unplaced/NA-like => map to Un (do NOT trigger fallback)
      # This is checked early so chrUn_random etc becomes Un.
      if (s ~ /^(chromosome|chrom|chr)?_?(un|unk|unknown|na|nan|unplaced|unlocali[sz]ed|random)(_|$)/) {
        print orig, 3, 0, "", 0
        next
      }

      # Accession-like contigs => force fallback (e.g., NW_0123.1 -> nw_0123_1 after normalization)
      if (s ~ /^[a-z]{2}_[0-9]+(_[0-9]+)?$/) {
        die("Unparseable accession-like contig (fallback to faidx): " orig)
      }

      # Strip common leading chromosome-ish prefixes (repeat to handle nesting)
      for (i=0; i<5; i++) {
        if (s ~ /^(chromosome|chrom|chr|linkage_group|linkagegroup|linkage|lg|group)_?/) {
          sub(/^(chromosome|chrom|chr|linkage_group|linkagegroup|linkage|lg|group)_?/, "", s)
          sub(/^_+/, "", s)
        } else break
      }

      # Sex chromosomes
      if (s=="x" || s=="y" || s=="w" || s=="z") {
        sexrank=(s=="x"?1:(s=="y"?2:(s=="w"?3:4)))
        print orig, 1, 0, "", sexrank
        next
      }

      # Mito
      if (s=="mt" || s=="m" || s=="mitochondrion" || s=="mitochondrial" || s=="mitochondria") {
        print orig, 2, 0, "", 0
        next
      }

      # Autosome-like: 0–999 + optional suffix letters ONLY
      if (s ~ /^0*[0-9]{1,3}[a-z]*$/) {
        tmp=s
        sub(/^0+/, "", tmp); if (tmp=="") tmp="0"

        num=tmp
        sub(/[^0-9].*$/, "", num)   # leading digits
        suf=tmp
        sub(/^[0-9]+/, "", suf)     # trailing letters (may be empty)

        print orig, 0, num+0, suf, 0
        next
      }

      # Anything else => abandon Step 2
      die("Unparseable contig (fallback to faidx): " orig)
    }
  ' |
  LC_ALL=C sort -t$'\t' -k2,2n -k3,3n -k4,4 -k5,5n -k1,1 |
  awk -F'\t' '
    BEGIN{OFS="\t"; i=0}
    {
      orig=$1; cat=$2
      if (cat==3) { print orig, "Un" }
      else { i++; print orig, i }
    }
  ' > "$ChromChrommappingfile" || step2_ok=0

  # If step2 produced nothing, also treat as fail
  if [[ ! -s "$ChromChrommappingfile" ]]; then
    step2_ok=0
  fi
  set +o pipefail


  # STEP 3: fallback to FASTA order

  if [[ $step2_ok -eq 0 ]]; then
    ref="${finalgenomepath}/${finalgenome}"

    if [[ ! -f "${ref}.fai" ]]; then
      samtools faidx "$ref"
    fi

    # Extract contig names in fasta order and assign numeric 1..N.
    # Unknown-like contigs are mapped to "Un" (not incrementing i).
    awk '
      BEGIN{OFS="\t"; i=0}

      function norm(x,   t) {
        t=tolower(x)
        gsub(/[^0-9a-z]+/, "_", t)
        sub(/^_+/, "", t); sub(/_+$/, "", t)
        return t
      }

      {
        name=$1
        t=norm(name)

        # Same unknown rule (broad)
        if (t ~ /^(chr|chrom|chromosome)?_?(un|unk|unknown|na|nan|unplaced|unlocali[sz]ed|random)(_|$)/) {
          print name, "Un"
        } else {
          i++
          print name, i
        }
      }
    ' "${ref}.fai" > "$ChromChrommappingfile"
  fi

#Apply new chrom names to vcf file 

  awk -v MAP="$ChromChrommappingfile" '
    BEGIN{
      FS=OFS="\t"
      while ((getline line < MAP) > 0) {
        gsub(/\r$/, "", line)
        if (line ~ /^[[:space:]]*$/) continue
        split(line, f, "\t")
        m[f[1]] = f[2]
      }
      close(MAP)
    }

    /^##contig=<ID=/{
      line=$0
      sub(/^##contig=<ID=/, "", line)
      split(line, a, /[,>]/)
      id=a[1]
      if (id in m) {
        sub("##contig=<ID=" id ",", "##contig=<ID=" m[id] ",", $0)
        sub("##contig=<ID=" id ">", "##contig=<ID=" m[id] ">", $0)
      }
      print
      next
    }

    /^#/ { print; next }

    { if ($1 in m) $1=m[$1]; print }
  ' "$final_vcf" > "${outputmasterfolder}/final/anchored_results/${final_vcfbasename}_numericchroms.vcf"

fi


bcftools sort -O v -o "${outputmasterfolder}/final/anchored_results/${final_vcfbasename}_sorted.vcf" "${outputmasterfolder}/final/anchored_results/${final_vcfbasename}_numericchroms.vcf"


  awk '
    /^#/ { print; next }
    /MAPSTATUS=Failed_to_map_uniquely/ { print }
  ' "$final_vcfrefs" > "${outputmasterfolder}/final/anchored_results/${final_vcfbasename}_failedtomap.vcf"
  cut -f 3 "${outputmasterfolder}/final/anchored_results/${final_vcfbasename}_failedtomap.vcf" > "${outputmasterfolder}/final/anchored_results/unmapped_markers.tsv"

