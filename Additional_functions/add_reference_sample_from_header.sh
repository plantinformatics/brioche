#!/bin/bash
set -euo pipefail
# Usage
#  add_reference_sample_from_header.sh input.vcf > output.vcf
#   zcat input.vcf.gz | add_reference_sample_from_header.sh - > output.vcf
#
# Behavior
# Sample name from last '##reference=' (basename; strips .fa/.fna/.fasta[.gz]).
# Appends that sample to the #CHROM header.
# for rows where FILTER has 'LowQual' -> GT=./., NU=1, all other FORMAT keys '.'
#   (so with FORMAT=GT:NU you get exactly './.:0' with the NU=0 indicating an NA vs a Not called drop out [see Null allele info in anchoring script.).
# - Otherwise -> GT=0/0, NU='0', others '.'.
# - If the sample already exists in the sample col, output will be unchanged and no new columns added.
# script does not update the MAF or AC etc values in the vcf as it is intened to be an intermediate file! 

infile="${1:-/dev/stdin}"

awk '
  BEGIN{
    FS = "\t"; OFS = "\t"
    add_sample = 1
    refln = ""; newsample = ""

    # Tags that should default to numeric 0 on non-LowQual rows
    numeric["NU"]=1; numeric["DP"]=1; numeric["GQ"]=1
    numeric["AD"]=1; numeric["PL"]=1; numeric["MIN_DP"]=1; numeric["MQ"]=1
  }

  # find last ##reference= and pass all meta lines
  /^##reference=/ { refln = $0; print; next }
  /^##/           { print; next }

  # Header line: derive sample name, check duplicates in other samplenames, append if new
  /^#CHROM/ {
    # normalize \r if present
    sub(/\r$/, "", $0)

    name = refln
    sub(/^##reference=/, "", name)
    gsub(/^file:\/\//, "", name)        # strip file://
    sub(/.*\//, "", name)               # basename
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

  # Variant rows
  {
    
    sub(/\r$/, "", $0)

    if (add_sample == 0) { print $0; next }

    # Detect LowQual in FILTER (col 7, semicolon-delimited)
    lowq = ($7 == "LowQual" || $7 ~ /(^|;)LowQual($|;)/) ? 0 : 1

    fmt = $9
    if (fmt == "." || fmt == "") {
      print $0, (lowq ? "./." : "0/0")
      next
    }

    nf = split(fmt, F, ":")

    if (lowq) {
      sampleval = ""
      for (i=1; i<=nf; i++) {
        tag = F[i]
        if (tag == "GT")      v = "./."
        else if (tag == "NU") v = "0"
        else                  v = "."
        sampleval = (i==1 ? v : sampleval ":" v)
      }
      print $0, sampleval
      next
    }

    # Non-LowQual path: GT=0/0, numeric=0, others '.'
    sampleval = ""
    for (i=1; i<=nf; i++) {
      tag = F[i]
      if (tag == "GT")            v = "0/0"
      else if (tag in numeric)    v = "0"
      else                        v = "."
      sampleval = (i==1 ? v : sampleval ":" v)
    }
    print $0, sampleval
  }
' "${infile}"
