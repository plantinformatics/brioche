#!/bin/bash
set -euo pipefail
# Usage
#  add_reference_sample_from_header.sh input.vcf > output.vcf
#  zcat input.vcf.gz | add_reference_sample_from_header.sh - > output.vcf
#
# Behavior
# - Sample name from last '##reference=' (basename; strips .fa/.fna/.fasta[.gz]).
# - Appends that sample to the #CHROM header unless it already exists.
# - Rows with FILTER=LowQual (or containing LowQual) **OR CHROM==chrUnk**:
#       GT=./., NC=1, all other FORMAT keys '.'
# - Otherwise:
#       GT=0/0, NU=0 for numeric keys, others '.'
# - Does not edit INFO tags like MAF/AC/AN (intermediate only).

infile="${1:-/dev/stdin}"

awk '
  BEGIN{
    FS = "\t"; OFS = "\t"
    add_sample = 1
    refln = ""; newsample = ""

    # Tags that default to numeric 0 on good rows
    numeric["NU"]=1; numeric["DP"]=1; numeric["GQ"]=1
    numeric["AD"]=1; numeric["PL"]=1; numeric["MIN_DP"]=1; numeric["MQ"]=1
  }

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

    fmt = $9
    if (fmt == "." || fmt == "") {
      # no FORMAT: just emit GT
      print $0, (bad ? "./." : "0/0")
      next
    }

    nf = split(fmt, F, ":")

    if (bad) {
      # bad rows: GT=./., NU=1, others "."
      sampleval = ""
      for (i=1; i<=nf; i++) {
        tag = F[i]
        if      (tag == "GT") v = "./."
        else if (tag == "NU") v = "1"
        else                  v = "."
        sampleval = (i==1 ? v : sampleval ":" v)
      }
      print $0, sampleval
      next
    }

    # good rows: GT=0/0, numeric=0 (incl. NU), others "."
    sampleval = ""
    for (i=1; i<=nf; i++) {
      tag = F[i]
      if      (tag == "GT")         v = "0/0"
      else if (tag in numeric)      v = "0"
      else                          v = "."
      sampleval = (i==1 ? v : sampleval ":" v)
    }
    print $0, sampleval
  }
' "${infile}"
