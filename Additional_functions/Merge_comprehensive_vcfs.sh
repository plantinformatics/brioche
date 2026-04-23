#!/usr/bin/env bash
set -euo pipefail


# Simple bash script for the merging of vcf/vcf.gz files. 
# Script will merge Brioche outputs which contain chrUnk markers without the risk of markers being lost due to pos=0 information. 
# Script will detect the number and order of markers in each vcf and confirm they are compatible for merging, attempt to sort if required and conform REF/ALT snp calls if required due to insilico results. 
# Script will output a vcf with nMarkers=nMarkersFile1=nMarkersFile2 and nSamples= nSamplesFile1 +nSamplesFile2
# The user should update all settings in the madatory user inputs section and confirm that 


############################################

# Mandatory user inputs (edit these, or I'll make it a command with argparse, not sure yet)

############################################

# First vcf/vcf.gz file to merge. If one of the files is a vcf generated from insilico genotyping set it to File1 as incongruities will be settled with the insilico results. If neither are, can set the files in whatever order
File1="File1.vcf.gz"

 # Second vcf/vcf.gz file to merge. If one of the files is a vcf generated from insilico genotyping do not set that files as File 2. 
File2="File2.vcf"

# Output file name. Files output will be a bgzip vcf file and it's index e.g., merged_vcf_refs_comprehensive.vcf.gz merged_vcf_refs_comprehensive.vcf.gz.csi
Outputfilename="merged_vcf_refs_comprehensive"


# If chrUnk REF/ALT states differ between File1 and File2:
#   Yes = automatically transform File2 to match File1 using the R script
#   No  = stop and exit for manual review

update_REFALT_state="Yes"

###########################################
# End mandatory user inputs
###########################################

#######################################

# Optional user input (You can update the file location but it comes in the same directory as this file so is fine unless things have been moved)

#######################################

# Path to the chrUnk REF/ALT harmonisation script
transform_chrunk_script="transform_chrUnk_base_to_insilico.R"

###########################################
# End optional user inputs
###########################################


#########################################

# Load required modules/programs 

#########################################


# Update this for whatever you are running so that conda is active 
# e.g., use a module or change to conda activate base if conda is installed locally for you or any other method so that the conda env created by brioche can be built/used

module load Miniconda3



#########################################
# end load required modules/programs 
#########################################


###################################################################################################################################################################

################################################ NO FURTHER EDITS REQUIRED BY USER PAST THIS POINT ################################################################

###################################################################################################################################################################


#########################################

 ########## START RUN_SCRIPT ##########

#########################################

# check if conda env is installed, install/activate as necessary 
eval "$(conda shell.bash hook)"
conda_path="$(conda info --base)"
conda config --set channel_priority flexible

ENV_NAME="brioche-vcf"

if conda env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
  echo "Conda env '$ENV_NAME' exists; proceeding with anchoring."
else
  echo "Conda env '$ENV_NAME' not found; creating environment before beginning anchoring..."
  conda env create --file "brioche-vcf.yaml" --name "$ENV_NAME" --yes
  echo "Conda environment '$ENV_NAME' fully installed."
fi

conda activate "$ENV_NAME"
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
echo "Active env: ${CONDA_DEFAULT_ENV}"

# confirm basic programs are there
for exe in bcftools bgzip awk grep tail cut diff mktemp zcat wc basename sort join head; do
  command -v "$exe" >/dev/null 2>&1 || { echo "ERROR: missing '$exe' in PATH"; exit 127; }
done


# Create temp working directory

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT


# Compare and transform vcf file in case of anchoring mismatch from insilico into chrUnk



############################################
reconcile_chrunk_refalt() {
  local insilico_vcf="$1"
  local base_vcf="$2"

  local transform_out="$tmpdir/$(basename "${base_vcf%.gz}").transformed.vcf"
  local f1_chrunk="$tmpdir/$(basename "${insilico_vcf%.gz}").chrunk.tsv"
  local f2_chrunk="$tmpdir/$(basename "${base_vcf%.gz}").chrunk.tsv"
  local f1_ids="$tmpdir/$(basename "${insilico_vcf%.gz}").chrunk.ids"
  local f2_ids="$tmpdir/$(basename "${base_vcf%.gz}").chrunk.ids"
  local joined="$tmpdir/chrunk_join.tsv"
  local mismatches="$tmpdir/chrunk_refalt_mismatches.tsv"

  echo "Checking chrUnk REF/ALT compatibility between:"
  echo "  File1 (insilico): $insilico_vcf"
  echo "  File2 (base):     $base_vcf"

  # Extract chrUnk marker ID + REF + ALT from each VCF
  awk 'BEGIN{OFS="\t"} !/^#/ && $1 ~ /chrUnk/ { print $3, toupper($4), toupper($5) }' "$insilico_vcf" | sort -k1,1 > "$f1_chrunk"
  awk 'BEGIN{OFS="\t"} !/^#/ && $1 ~ /chrUnk/ { print $3, toupper($4), toupper($5) }' "$base_vcf"     | sort -k1,1 > "$f2_chrunk"

  # Confirm the chrUnk marker sets are identical
  cut -f1 "$f1_chrunk" > "$f1_ids"
  cut -f1 "$f2_chrunk" > "$f2_ids"

  if ! diff -q "$f1_ids" "$f2_ids" >/dev/null; then
    echo "ERROR: chrUnk marker IDs differ between the two numbered VCFs."
    echo "Please confirm both VCFs contain the same marker set before merge."
    exit 2
  fi

  # Join on marker ID and find REF/ALT mismatches
  join -t $'\t' -j 1 "$f1_chrunk" "$f2_chrunk" > "$joined"
  awk 'BEGIN{OFS="\t"} $2 != $4 || $3 != $5 { print $1, $2, $3, $4, $5 }' "$joined" > "$mismatches"

  local mismatch_count
  mismatch_count=$(wc -l < "$mismatches")
  mismatch_count="${mismatch_count//[[:space:]]/}"

  if [[ "$mismatch_count" -eq 0 ]]; then
    echo "No chrUnk REF/ALT discrepancies detected."
    return 0
  fi

  echo "Detected $mismatch_count chrUnk REF/ALT discrepancies."
  echo "Examples (ID  insilico_REF  insilico_ALT  base_REF  base_ALT):"
  head -n 10 "$mismatches"

  case "${update_REFALT_state,,}" in
    yes)
      echo "One or more unanchored markers has had their REF/ALT state anchored by a different reference genome."
      echo "update_REFALT_state=Yes, so File2 will now be transformed to match File1."
      echo "Please confirm that File1 is the insilico results file."

      Rscript "$transform_chrunk_script" \
        --insilico_vcf "$insilico_vcf" \
        --base_vcf "$base_vcf" \
        --Outputfilename "$transform_out"

      mv -f "$transform_out" "$base_vcf"

      # If you want to keep the sidecar reports in tmpdir, leave these in place.
      # If not needed, you can remove these mv lines entirely.
      [[ -f "${transform_out}.audit" ]]   && mv -f "${transform_out}.audit"   "${base_vcf}.audit"
      [[ -f "${transform_out}.summary" ]] && mv -f "${transform_out}.summary" "${base_vcf}.summary"

      echo "chrUnk REF/ALT transformation complete: $base_vcf"
      ;;
    no)
      echo "update_REFALT_state=No. Exiting without merge."
      echo "Please either set update_REFALT_state=Yes, or manually review the VCFs before merging."
      exit 2
      ;;
    *)
      echo "ERROR: update_REFALT_state must be set to Yes or No (current value: '$update_REFALT_state')."
      exit 2
      ;;
  esac
}

##################




# Create uncompressed versions of files in necessary
File1_plain="$tmpdir/file1.vcf"
File2_plain="$tmpdir/file2.vcf"

echo "Preparing uncompressed working copy for File1..."
if [[ "$File1" =~ \.gz$ ]]; then
  # Try bgzip -dc first; fall back to zcat if needed
  if bgzip -dc "$File1" > "$File1_plain" 2>/dev/null; then
    :
  else
    zcat "$File1" > "$File1_plain"
  fi
else
  cp -f "$File1" "$File1_plain"
fi

echo "Preparing uncompressed working copy for File2..."
if [[ "$File2" =~ \.gz$ ]]; then
  if bgzip -dc "$File2" > "$File2_plain" 2>/dev/null; then
    :
  else
    zcat "$File2" > "$File2_plain"
  fi
else
  cp -f "$File2" "$File2_plain"
fi


File1basename="$(basename "${File1%.gz}")"
File2basename="$(basename "${File2%.gz}")"

#  Checking compatibility of files
# 1. Same number of markers?

file1length="$(grep -vc '^#' "$File1_plain")"
file2length="$(grep -vc '^#' "$File2_plain")"

echo "Record counts: File1=$file1length, File2=$file2length"

if [[ "$file1length" -ne "$file2length" ]]; then
  echo "ERROR: The two VCF files provided have different numbers of markers."
  echo "Are you sure they are both anchored with the same Brioche mappings file?"
  echo " Terminating merge as these files may not be compatible"
  exit 2
fi


# 2. Compare chrUnk order last 50 marker IDs
final50file1="$tmpdir/final50_file1.txt"
final50file2="$tmpdir/final50_file2.txt"

grep -v '^#' "$File1_plain" | tail -n 50 | cut -f3 > "$final50file1"
grep -v '^#' "$File2_plain" | tail -n 50 | cut -f3 > "$final50file2"

if diff -q "$final50file1" "$final50file2" >/dev/null; then
  echo "Final 50 markers of both files are in the same order."
  echo "Now naming chrUnk sequentially as files are identically ordered."

  ##########################################
  # Sequential numbering for both VCFs
  ##########################################
  File1_numbered="$tmpdir/${File1basename}_numbered.vcf"
  File2_numbered="$tmpdir/${File2basename}_numbered.vcf"

  awk '
  BEGIN { OFS="\t"; unk=1 }
  /^##/ { print; next }
  /^#CHROM/ { print; next }
  {
      if ($1 ~ /chrUnk/) {
          $2 = unk++
      }
      print
  }
  ' "$File1_plain" > "$File1_numbered"

  awk '
  BEGIN { OFS="\t"; unk=1 }
  /^##/ { print; next }
  /^#CHROM/ { print; next }
  {
      if ($1 ~ /chrUnk/) {
          $2 = unk++
      }
      print
  }
  ' "$File2_plain" > "$File2_numbered"


  # Confirm no REF/ALT discrepencies
  reconcile_chrunk_refalt "$File1_numbered" "$File2_numbered"



  File1_numbered_gz="$tmpdir/${File1basename}_numbered.vcf.gz"
  File2_numbered_gz="$tmpdir/${File2basename}_numbered.vcf.gz"

  bgzip -c "$File1_numbered" > "$File1_numbered_gz"
  bgzip -c "$File2_numbered" > "$File2_numbered_gz"

  bcftools index -c "$File1_numbered_gz"
  bcftools index -c "$File2_numbered_gz"

  merged_gz="$tmpdir/merge.vcf.gz"
  bcftools merge -Oz -o "$merged_gz" "$File1_numbered_gz" "$File2_numbered_gz"


  merged_plain="$tmpdir/merge.vcf"
  bgzip -dc "$merged_gz" > "$merged_plain"

  merged_chrunk_plain="$tmpdir/merge_chrunk.vcf"

  awk '
  BEGIN { OFS="\t" }
  /^##/ { print; next }
  /^#CHROM/ { print; next }
  {
      if ($1 ~ /chrUnk/) {
          $2 = 0
      }
      print
  }
  ' "$merged_plain" > "$merged_chrunk_plain"

  ##########################################
  # Output name handling
  ##########################################
  out="$Outputfilename"
  if [[ "$out" != *.vcf.gz ]]; then
    out="${out}.vcf.gz"
  fi

  bgzip -c "$merged_chrunk_plain" > "$out"
  bcftools index -c "$out" || true

  echo "Done. Output: $out"

else
  echo "Final 50 markers are NOT the same between VCFs."
  echo "Assuming one or both files are not sorted; assigning chrUnk POS based on marker name mapping from File1."

  # 1.Build marker map from File1 chrUnk markers: markerID -> sequential index

  mapfile="$tmpdir/chrUnk_marker_map.tsv"

  awk '
  BEGIN{ OFS="\t"; i=0 }
  /^#/ { next }
  $1 ~ /chrUnk/ {
    i++
    print $3, i
  }
  ' "$File1_plain" > "$mapfile"

  # 2.Apply marker map to File1
  File1_numbered="$tmpdir/${File1basename}_numbered.vcf"

  awk -v MAP="$mapfile" '
  BEGIN {
    OFS="\t"
    while ((getline < MAP) > 0) {
      map[$1] = $2
    }
    close(MAP)
  }
  /^##/ { print; next }
  /^#CHROM/ { print; next }
  {
    if ($1 ~ /chrUnk/) {
      if (!($3 in map)) {
        print "ERROR: chrUnk marker has no mapping in File1: " $3 > "/dev/stderr"
        exit 3
      }
      $2 = map[$3]
    }
    print
  }
  ' "$File1_plain" > "$File1_numbered"


  #3. Apply marker map to File2 (same mapping)

  File2_numbered="$tmpdir/${File2basename}_numbered.vcf"

  awk -v MAP="$mapfile" '
  BEGIN {
    OFS="\t"
    while ((getline < MAP) > 0) {
      map[$1] = $2
    }
    close(MAP)
  }
  /^##/ { print; next }
  /^#CHROM/ { print; next }
  {
    if ($1 ~ /chrUnk/) {
      if (!($3 in map)) {
        print "ERROR: chrUnk marker has no mapping in File2: " $3 > "/dev/stderr"
        exit 3
      }
      $2 = map[$3]
    }
    print
  }
  ' "$File2_plain" > "$File2_numbered"


  # Confirm no REF/ALT discrepencies
  reconcile_chrunk_refalt "$File1_numbered" "$File2_numbered"



  File1_numbered_gz="$tmpdir/${File1basename}_numbered.vcf.gz"
  File2_numbered_gz="$tmpdir/${File2basename}_numbered.vcf.gz"

  bgzip -c "$File1_numbered" > "$File1_numbered_gz"
  bgzip -c "$File2_numbered" > "$File2_numbered_gz"

  File1_sorted_gz="$tmpdir/${File1basename}_numbered_srt.vcf.gz"
  File2_sorted_gz="$tmpdir/${File2basename}_numbered_srt.vcf.gz"

  bcftools sort -Oz -o "$File1_sorted_gz" "$File1_numbered_gz"
  bcftools sort -Oz -o "$File2_sorted_gz" "$File2_numbered_gz"

  bcftools index -c "$File1_sorted_gz"
  bcftools index -c "$File2_sorted_gz"

  finalmarkerfile1="$(bgzip -dc "$File1_sorted_gz" | grep -v '^#' | tail -n 1 | cut -f3)"
  finalmarkerfile2="$(bgzip -dc "$File2_sorted_gz" | grep -v '^#' | tail -n 1 | cut -f3)"
  finalposfile1="$(bgzip -dc "$File1_sorted_gz" | grep -v '^#' | tail -n 1 | cut -f2)"
  finalposfile2="$(bgzip -dc "$File2_sorted_gz" | grep -v '^#' | tail -n 1 | cut -f2)"

  if [[ "$finalmarkerfile1" == "$finalmarkerfile2" && "$finalposfile1" == "$finalposfile2" ]]; then
    echo "Files have been sorted and are compatible for merging."
  else
    echo "ERROR: After sorting, the final row marker and position are not identical between VCFs."
    echo "File1 last: marker=$finalmarkerfile1 pos=$finalposfile1"
    echo "File2 last: marker=$finalmarkerfile2 pos=$finalposfile2"
    echo "Please confirm the VCFs have identical markers and came from the same brioche mapping file!"
    echo "Terminating merge."
    exit 2
  fi



  
  # 5. Merge finals
  
  merged_gz="$tmpdir/merge.vcf.gz"
  bcftools merge -Oz -o "$merged_gz" "$File1_sorted_gz" "$File2_sorted_gz"

  #5.  Decompress merged, set chrUnk POS back to 0
  
  merged_plain="$tmpdir/merge.vcf"
  bgzip -dc "$merged_gz" > "$merged_plain"

  merged_chrunk_plain="$tmpdir/merge_chrunk.vcf"

  awk '
  BEGIN { OFS="\t" }
  /^##/ { print; next }
  /^#CHROM/ { print; next }
  {
      if ($1 ~ /chrUnk/) {
          $2 = 0
      }
      print
  }
  ' "$merged_plain" > "$merged_chrunk_plain"

  
  #6. Output name handling
  
  out="$Outputfilename"
  if [[ "$out" != *.vcf.gz ]]; then
    out="${out}.vcf.gz"
  fi

  bgzip -c "$merged_chrunk_plain" > "$out"
  bcftools index -c "$out" || true

  echo "Done. Output: $out"
fi
