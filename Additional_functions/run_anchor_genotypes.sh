#!/bin/bash

# Multithreaded (SMP) job: must run on one node
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="example_anchoring"

# The project ID which this job should run under:
#SBATCH --account="default"

# The partition in which to submit the job:
#SBATCH --partition="batch"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

# The total amount of memory in megabytes in the job:
#SBATCH --mem=60GB

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL

# Use this email address:
#SBATCH --mail-user=user123@mail.mail

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=4:0:00


################
# To run this script, 
# 1. Update the SBATCH commands with your information 
# 2. Set do next stages to no to skip downstream filtering or yes to update chromosome names, sort, and also output a mapped markers only vcf output alongside all markers output
# 3. skip down to Rscript Anchoring_script.R and update the variables for the required files 
# 4.  if donextstages is set to "yes" fill in the name of the final vcf output (The name to output from the R anchoring script in setting --Outputfilename) and the path/name of a reference genome, and if one is present a 
## tsv file with two columns, first column the name of markers in the reference genome used to map against and second column, the numeric numbers from 1-N for chromosomes + one final line for chrUnk\tUn for unknown chromosome 
### or leave ChromChrommappingfile blank for attempted automatic chrom number assigning (relies on 1, pattern matching, 2, reference genome fasta header order)
# NOTE: next stages assumes that the output files are in the working directory where this is run
###############

###############Variables to set ###################

donextstages="no" # This setting determines whether to stop after anchoring or to also create a new vcf file which removes all unmapped markers and converts the chromosome numbers to numeric values and sorts the VCF

finalgenomepath="/filepathgenome/" # filepath to reference genome
finalgenome="genome_1.fasta" # reference genome fasta file
ChromChrommappingfile="" # Chromosome name to numeric comparison file
Rawgenotypesfile="Rawgenotypesfile.tsv" # Raw genotypes file (either SNP, VCF, or DArT format)
briochemappingfile="Briocheresultsdatasetx_all_markers1to1stagingforvcf.csv" # brioche1to1 mappings path + file from brioche-results
Brioche_priors="Briocheresultsdatasetx_priors_informed_strictmapping.tsv" # brioche priorsfile path +file from brioche-results

VCFraws="false" # If the input genotypes file is a VCF set to "true" if not set to "false"
SNPchip="true" # If the input genotypes file is in Name REF ALT SAMPLE1, 0,1,2 raw genotypes format set to "true" if not set to "false"
dartfile="false" # If the input genotypes file is a DArT 1 row format result report set to "true" if not set to "false". Note this setting has not been tested as extensively as the other two
final_vcf="Reanchored_vcf_datasetx_against_refgenomeY.vcf" # What to save the output file/path as
refgenomename="Species X whole genome v1.1" # The name of the reference genome used to remap markers (taken from the NCBI title of the genome)
ACCcode="GCFXXX123" # The NCBI Accession code of the reference genome used


######################End variables to set ###################


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


date

Rscript Anchoring_script.R --Rawgenotypes "$Rawgenotypesfile" --Briochemappings "$briochemappingfile" --mappriors "$Brioche_priors" --IsVCFraws "$VCFraws" --IsSNPchip "$SNPchip" --isDArTfile "$dartfile" --Outputfilename "$final_vcf" --genomename "$refgenomename" --Accession "$ACCcode" --droplist ""

date


# If updating chroms, sorting of vcf and removal of unmapped markers is wanted
# The following is nearly identical copy to the script in run insilico for this 
# Update the variables below to run

if [[ $donextstages == "yes" ]]; then


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
  ' "$final_vcf" > "${final_vcfbasename}_numericchroms.vcf"
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
  ' "$final_vcf" > "${final_vcfbasename}_numericchroms.vcf"

fi


bcftools sort -O v -o "${final_vcfbasename}_sorted.vcf" "${final_vcfbasename}_numericchroms.vcf"


  awk '
    /^#/ { print; next }
    /MAPSTATUS=Failed_to_map_uniquely/ { print }
  ' "$final_vcf" > "${final_vcfbasename}_failedtomap.vcf"
  cut -f 3 ${final_vcfbasename}_failedtomap.vcf > unmapped_markers.tsv
  awk '
    /^#/ { print; next }
    !/MAPSTATUS=Failed_to_map_uniquely/ { print }
  ' "${final_vcfbasename}_sorted.vcf" > "${final_vcfbasename}_sorted_uniquelymappedmarkers.vcf"


fi




