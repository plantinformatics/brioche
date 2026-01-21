#!/bin/bash
set -euo pipefail
# Usage:
#   add_reference_sample_from_header.sh input.vcf [marker_meta.tsv] > output.vcf
#   zcat input.vcf.gz | add_reference_sample_from_header.sh - marker_meta.tsv > output.vcf
#
# Behavior
# - Sample name from last '##reference=' (basename; strips .fa/.fna/.fasta[.gz]); falls back to "REFERENCE".
# - Appends that sample to the #CHROM header unless it already exists.
# - Rows with FILTER=LowQual (or containing LowQual) OR CHROM==chrUnk:
#       reference sample: GT=./., NU=1, DU='.' (if DU in FORMAT), all other FORMAT keys '.'
# - Otherwise (good rows):
#       reference sample: GT=0/0, numeric keys=0 (incl. NU by default), non-numeric '.'
#       If a TSV is provided (qaccver,copy_number,consensusbase,keep,nulcall):
#         * When TSV.nulcall == "Yes" for that ID, set NU=1 (overrides default 0).
#         * If DU present in FORMAT: set DU to (copy_number-1) when copy_number >=1; else '.'
#       If DU is NOT present in FORMAT, we append ':DU' to FORMAT and add ':.' to ALL existing
#       sample fields on the line so the row stays valid; then we set DU for the reference sample.
#       Populating for the missing Genotype formats is probably excessive given this is an internal script though so will potentially strip down later.
#

infile="${1:-/dev/stdin}"
tsv="${2:-/dev/null}"

awk -v HAVE_TSV="$( [ -s "${tsv}" ] && echo 1 || echo 0 )" '
  BEGIN{
    FS = "\t"; OFS = "\t"
    add_sample = 1
    refln = ""; newsample = ""

    # Tags that default to numeric 0 on good rows
    delete numeric
    numeric["NU"]=1; numeric["DP"]=1; numeric["GQ"]=1
    numeric["AD"]=1; numeric["PL"]=1; numeric["MIN_DP"]=1; numeric["MQ"]=1

    # maps from TSV
    delete nul_by_id
    delete copy_by_id
  }

  ## TSV
  HAVE_TSV && FNR==NR {
    # Expect header at first line; flexible column order
    if (FNR==1) {
      for (i=1; i<=NF; i++){
        k=tolower($i)
        if (k=="qaccver") qi=i
        else if (k=="copy_number") ci=i
        else if (k=="nulcall") ni=i
      }
      next
    }
    id = (qi ? $qi : "")
    if (id=="") next
    # store nulcall flag (Yes/True/1 => 1; else 0)
    raw = (ni ? $ni : "")
    gsub(/\r/,"",raw)
    l = tolower(raw)
    nul_by_id[id] = (l=="yes" || l=="true" || l=="1") ? 1 : 0

    # store copy_number (as integer if possible)
    cn = (ci ? $ci : "")
    gsub(/\r/,"",cn)
    if (cn ~ /^[0-9]+([.][0-9]+)?$/) {
      # integerize; only accept >=1 for DU (we subtract 1 later)
      val = int(cn+0)
      copy_by_id[id] = val
    } else {
      copy_by_id[id] = ""  # unknown / NA
    }
    next
  }

  ## ---------- Pass 2: VCF ----------
  # pass meta; remember last ##reference=
  /^##reference=/ { refln = $0; print; next }
  /^##/           { print; next }

  # header line: derive sample name, append if missing
  /^#CHROM/ {
    sub(/\r$/, "", $0)
    name = refln
    sub(/^##reference=/, "", name)
    gsub(/^file:\/\//, "", name)
    sub(/.*\//, "", name)
    sub(/\.(fa|fna|fasta)(\.gz)?$/, "", name)
    if (name == "") name = "REFERENCE"
    newsample = name

    n = split($0, H, "\t")
    dup = 0
    for (i=10; i<=n; i++) if (H[i] == newsample) { dup=1; break }

    if (dup) { add_sample=0; print $0 }
    else     { add_sample=1; print $0, newsample }
    next
  }

  # data rows
  {
    sub(/\r$/, "", $0)
    if (add_sample == 0) { print $0; next }

    # bad if LowQual OR chrUnk
    lowq = ($7 == "LowQual" || $7 ~ /(^|;)LowQual($|;)/) ? 1 : 0
    unk  = ($1 == "chrUnk") ? 1 : 0
    bad  = (lowq || unk) ? 1 : 0

    id = $3

    fmt = $9
    if (fmt == "." || fmt == "") {
      # no FORMAT: just emit GT
      print $0, (bad ? "./." : "0/0")
      next
    }

    # Ensure DU tag exists if we have a TSV and wish to record copy_number
    has_DU = (fmt ~ /(^|:)DU($|:)/) ? 1 : 0
    if (!has_DU) {
      # Extend FORMAT and all existing sample columns with a placeholder for DU
      $9 = (fmt == "." ? "DU" : fmt ":DU")
      for (i=10; i<=NF; i++) {
        if ($i == "." || $i == "") $i = "."
        else $i = $i ":."
      }
      has_DU = 1
      fmt = $9
    }

    nf = split(fmt, F, ":")

    # Build the reference sample value
    sampleval = ""
    for (i=1; i<=nf; i++) {
      tag = F[i]
      v = "."
      if (bad) {
        if      (tag == "GT") v = "./."
        else if (tag == "NU") v = "1"
        else if (tag == "DU") v = "."
        else                  v = "."
      } else {
        if      (tag == "GT")    v = "0/0"
        else if (tag == "NU") {
          # default 0; override to 1 if TSV says nulcall Yes
          v = "0"
          if (id in nul_by_id && nul_by_id[id] == 1) v = "1"
        }
        else if (tag == "DU") {
          # default "."; if copy_number >=1 then DU = copy_number - 1
          if (id in copy_by_id) {
            cn = copy_by_id[id]
            if (cn ~ /^[0-9]+$/ && cn+0 >= 1) v = (cn-1)
            else v = "."
          } else v = "."
        }
        else if (tag in numeric) v = "0"
        else                     v = "."
      }
      sampleval = (i==1 ? v : sampleval ":" v)
    }

    print $0, sampleval
  }
' "${tsv}" "${infile}"
