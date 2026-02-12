/*
 * Process 1: Create a blast database with makeblastdb
*/
process BUILD_TARGET {
    tag "${file(targetdesign).baseName}"
    executor 'local'
    input:
        val targetdesign
    output:
        tuple path("target.table"),
              path("target.fasta")
    script:
        mincols = 4
        """
        #!/usr/bin/env bash
        set -euo pipefail  # Enable strict error handling
        # Check if the file exists
        if [[ ! -e "${targetdesign}" ]]; then
            echo "File ${targetdesign} not found"
            exit 1
        fi

        # Detect delimiter (comma or tab)
        format=""
        if [[ ! -z \$(head -n1 "${targetdesign}" | awk -F, '{if(NF>1) print \$NF}') ]]; then
            format=","
        elif [[ ! -z \$(head -n1 "${targetdesign}" | awk -F'\t' '{if(NF>1) print \$NF }') ]]; then
            format="\t"
        else
            echo "${targetdesign} is not a CSV or tab-separated file."
            exit 1
        fi

        # Check for minimum number of columns
        ncols=\$(head -n1 "${targetdesign}" | awk -F"\${format}" '{print NF}')
        if [[ "\${ncols}" -lt "${mincols}" ]]; then
            echo "${targetdesign} does not have the correct number of columns, minimum expected is ${mincols}."
            exit 1
        fi

        # Check if the second column contains a valid sequence
        seq1=\$(head "${targetdesign}" | tail -n1 | awk -F"\${format}" '{print \$2}')
        nucleotides="ACTGUN${params.markercharacter}"
        if [[ ! "\${seq1}" =~ ^[\${nucleotides}]+\$ ]]; then  # Use ^ and \$ for exact match
            echo "Error! Second column in the target file is supposed to be the probe/marker sequence"
            exit 1
        fi

        # Check for duplicate IDs
        if [[ ! -z \$(awk -F"\${format}" 'BEGIN{print } {print \$1}' "${targetdesign}" | sort | uniq -d) ]]; then
            echo "ID column in ${targetdesign} cannot have duplicates"
            exit 1
        fi

        echo "Build target and fasta files"

        # Create target table
        if [[ "\${ncols}" -ge 5 ]]; then
            awk -F"\${format}" 'BEGIN{print "ID\tTarget.Chr\tTarget.bp\tTarget.base"} !/Target/{print \$1"\t"\$5"\t"\$3"\t"\$4}' "${targetdesign}" > target.table
        else
            awk -F"\${format}" 'BEGIN{print "ID\tTarget.Chr\tTarget.bp\tTarget.base"} !/Target/{print \$1"\t-\t"\$3"\t"\$4}' "${targetdesign}" > target.table
        fi

        # Create probe fasta file
        grep -v "Sequence" "${targetdesign}" | awk -F"\${format}" '{print ">"\$1"\\n"\$2}' > target.fasta
        """
}

process SPLIT_TARGET_FASTA {
  tag { file(fasta).baseName }
  executor 'local'
  input:
    path fasta
    val  chunk_size
  output:
    path "split/*.fa"
  script:
  """
  set -euo pipefail
  mkdir -p split

  # sanity check
  if [[ \$(grep -c '^>' "${fasta}") -eq 0 ]]; then
    echo "No sequences found in ${fasta}" >&2
    exit 1
  fi

  awk -v n=${chunk_size} '
    BEGIN{ c=0; f="" }
    /^>/{
      c++
      if ((c-1)%n==0) {
        if (f!="") close(f)
        f=sprintf("split/q_%06d.fa", int((c-1)/n)+1)
      }
    }
    { print >> f }
    END{
      if (f!="") close(f)
    }
  ' "${fasta}"
  """
}

process PREP_REFGENOME {
  input:
    path ref_in

  output:
    path "ref_full.fa"

  """
  if [[ "${ref_in}" == *.gz ]] || gzip -t "${ref_in}" &>/dev/null; then
    gzip -dc "${ref_in}" > ref_full.fa
  else
    ln -sf "${ref_in}" ref_full.fa
  fi
  """
}

process BUILD_BLASTDb {
  tag   "${file(fasta).baseName.replaceAll('-split','')}"
  label 'large'
  input:
    path fasta
    val  ref_origin
    val  genomename
  output:
    val  "${blastdatabase}"

  script:
  def fastadir       = new File(ref_origin as String).getParent()
  def blastsubdir    = new File("${fastadir}/BLAST")
  def dbstem         = file(fasta).getName().replaceAll("-split.*","")
  def refgenomename  = new File(ref_origin as String).getName()
  def genomeTag      = (genomename as String).replaceAll(/[^A-Za-z0-9._-]/,'_')
  def writable = new File(fastadir).canWrite() && ( !blastsubdir.exists() || blastsubdir.canWrite() )
  def runWorkDir = workflow.workDir.toString()
  def blastdbBase = writable ? blastsubdir.toString()
                             : "${runWorkDir}/BLAST/${genomeTag}"

  if( dbstem != refgenomename ) {
    if( writable )
      blastdbBase += "/brioche-blastdb${genomeTag}"
    // else: fallback path already segregated by ${genomeTag} keep as-is
  }

  blastdatabase = "${blastdbBase}/${dbstem}"

  """
  set -euo pipefail

  mkdir -p "\$(dirname "${blastdatabase}")"

  # Build only if missing (cover classic and multi-volume layouts)
  if [[ ! -e "${blastdatabase}.ndb" && ! -e "${blastdatabase}.nin" && ! -e "${blastdatabase}.00.nin" ]]; then
    echo "Creating BLAST DB for ${fasta} -> ${blastdatabase}"
    makeblastdb -in "${fasta}" -out "${blastdatabase}" -dbtype nucl -title "${genomename}"
  else
    echo "BLAST DB already exists at ${blastdatabase}"
  fi
  """
}


process SPLIT_REFGENOME {
  tag "${file(refgenome).baseName}"
  input:
    path refgenome
    val  resultsdir
    val  chromstoexclude
    val  buildblastdbonly
  output:
    path "*-split.fasta"

  // Pre-compute a few safe strings in Groovy
  script:
  def excl       = (chromstoexclude ?: '').replace(',', '|')   // "a,b,c" -> "a|b|c"
  def restrict   = params.restrict2chrom ?: ''
  def buildonly  = (buildblastdbonly instanceof Boolean ? buildblastdbonly : "${buildblastdbonly}".toBoolean()).toString()
  def path2briocheR = params.pathtobriocher
  def genomefasta   = "${file(refgenome).getName()}-split.fasta"

  """
  set -euo pipefail

  mkdir -p "${resultsdir}"

  # Use the staged input directly (no self-symlink)
  REF="${refgenome}"

  # Dictionary next to REF using bash param expansion (note the backslashes)
  dict="\${REF%.*}.dict"
  if [[ ! -e "\${dict}" ]]; then
    picard -Xms4g -Xmx10g CreateSequenceDictionary -R "\${REF}" -O "\${dict}"
  else
    echo "Dictionary file found"
  fi

  # Optional exclusion pattern from Groovy
  exclude_array="DC"
  if [[ -n "${excl}" ]]; then
    exclude_array="${excl}"
  fi

  # List contigs and decide whether to split
  chromosomes="\$(grep '@SQ' "\${dict}" | cut -f2 | sed -e 's/SN://g' | grep -v -P "^\${exclude_array}" || true)"
  split="\$(echo "\${chromosomes}" | awk -v maxchrom=30 '{print (NF <= maxchrom && NF>0) ? "yes" : "no"}')"

  if [[ "\${split}" == "yes" && "${buildonly}" == "false" ]]; then
    if [[ -n "\${chromosomes}" ]]; then
      for chr in \${chromosomes}; do
        samtools faidx "\${REF}" "\$chr" > "\${chr}-split.fasta"
      done
    else
      echo "No chromosomes found in fasta file"
      exit 1
    fi
  elif [[ -n "${restrict}" ]]; then
    samtools faidx "\${REF}" "${restrict}" > "${restrict}-split.fasta"
  else
    # No split: expose a single *-split.fasta pointing at REF
    ln -sf "\${REF}" "${genomefasta}"
  fi

  # Install local briocheR tarball if needed
  Rscript -e "if(!require('briocheR')){
    install.packages('${path2briocheR}', repos=NULL);
    print('BriocheR installed')
  }"
  """
}



/*
 * Step2: Process to run blast of the probes provided against reference genome
 */
process RUN_BLAST {
  tag "${file(probefasta).baseName}"
  label 'large'

  input:
    path probefasta
    val  blastdb
    val  blastoutformat

  // Single BLAST output file per run; name decided in the script
  output:
    path "*.blast.txt"

  script:
  """
  dbbase=\$(basename "${blastdb}")
  qbase=\$(basename "${probefasta}")
  outprefix="\${dbbase}_\${qbase}"

  blastn -db ${blastdb} -query ${probefasta} \\
         -out "\${outprefix}.blast.txt" \\
         -num_threads $task.cpus \\
         -outfmt ${blastoutformat} \\
         -max_target_seqs 20 \\
         -max_hsps 20 \\
         -evalue ${params.evalue} ${params.otherblastoptions} \\
         -dust ${params.dust}
  """
}

/*
* Step 3 Process  to pre-process the results from blastn
*/
process PROCESS_BLAST_RESULTS {
  tag "${blastout.baseName}"
  label 'large'
  input:
    path blastout
    val  targetdesign
    val  minmismatches
    val  extendablebps
    path genomefasta
    val  blastoutformat
    val  istarget3primeend
    val  markerchar
  output:
    path "${rds}"
    path "${snpmapped}"
  script:
  snpmapped = "${blastout.baseName}_samtools_coordinates_mapped.tsv"
  rds       = "SNPcoordinates/${blastout.baseName}_blast_results.RDS"
  """
  echo ${blastout}
  snpcoords="\$(pwd)/SNPcoordinates/"
  R --version
  Rscript -e "briocheR::processBlastResults(
                blast.file='${blastout}',
                min.length=${minmismatches},
                path.2save.coords = '\${snpcoords}',
                extendable.site.bps = ${extendablebps},
                targets='${targetdesign}',
                outfmt=strsplit(${blastoutformat},' ')[[1]][-1],
                istarget3primeend='${istarget3primeend}',
                marker.char='${markerchar}')"
  mv \${snpcoords}/processed_blast_results.RDS ${rds}

  samtools faidx ${genomefasta} -r "\${snpcoords}/samtools_marker_positions.txt" | awk '
    BEGIN {FS = "\\n"}
    /^>/ {
      if (seq != "") {
        print id "\\t" seq
      }
      id = substr(\$1, 2)
      seq = ""
    }
    !/^>/ {
      seq = seq \$1
    }
    END {
      print id "\\t" seq
    }
  ' > "${snpmapped}"
  """
}


/*
Step4: Perform filtering and generate annotations
*/
process ADD_MARKERINFO {
    tag "${blast_results.baseName}"
    input:
        path blast_results
        path coordinates
        val  probename
        val  genomename
        val  keepduplicates
        val  coveragefilter
        val  pidentfilter
        val  istarget3primeend
    /*
     * Outputs per BLAST chunk to optimise for ultra large file sizes:
     */
    output:
        tuple path(allmappings),
              path(filteredmappings),
              path(mappings),
              path(filtered_mappings_tsv),
              path(rawsallmappings)
    script:
    def random_factor        = new Random().nextInt()
    def raw_prefix           = "${blast_results.baseName}_${random_factor}"
    allmappings              = "${raw_prefix}_all_mappings.csv"
    mappings                 = "${raw_prefix}_mapping.tsv"
    rawsallmappings          = "${raw_prefix}_complete_unfiltered_blast_results.csv"
    filteredmappings         = "${raw_prefix}_filtered_mappings.csv"
    filtered_mappings_tsv    = "${raw_prefix}_filtered_mappings.tsv"

    """
    # 1) Run R addSNPdetails
    R --version
    Rscript -e "briocheR::addSNPdetails(
           blast.path='${blast_results}',
           reference.bases='${coordinates}',
           probe.name='${probename}',
           genome.name='${genomename}',
           output.path = '\$(pwd)',
           pident.filter=${pidentfilter},
           coverage.filter=${coveragefilter},
           is3prime=${istarget3primeend})"
    mv "${probename}_with_${genomename}_all_mappings.csv" ${allmappings}
    mv "${probename}_with_${genomename}_mapping.tsv"      ${mappings}
    mv "${probename}_with_${genomename}_complete_unfiltered_blast_results.csv" ${rawsallmappings}

    # 2 CSV: filter all_mappings -> filteredmappings

    header_csv=\$(head -n1 ${allmappings})
    snpcol_csv=\$(echo "\$header_csv" | awk -F, '{for(i=1;i<=NF;i++) if(\$i=="qaccver")  print i}')
    hybcol_csv=\$(echo "\$header_csv" | awk -F, '{for(i=1;i<=NF;i++) if(\$i=="Hybridized") print i}')
    covcol_csv=\$(echo "\$header_csv" | awk -F, '{for(i=1;i<=NF;i++) if(\$i=="Coverage") print i}')
    pidcol_csv=\$(echo "\$header_csv" | awk -F, '{for(i=1;i<=NF;i++) if(\$i=="pident") print i}')

    if [ -z "\$snpcol_csv" ]; then
        echo "ERROR: qaccver column not found in ${allmappings}" 1>&2
        exit 1
    fi
    if [ -z "\$covcol_csv" ] || [ -z "\$pidcol_csv" ]; then
        echo "ERROR: Coverage and/or pident column not found in ${allmappings}" 1>&2
        exit 1
    fi

    #  Sort allmappings by qaccver then Coverage, and pident
    tail -n +2 "${allmappings}" \
      | sort -t, -k\${snpcol_csv},\${snpcol_csv} -k\${covcol_csv},\${covcol_csv}nr -k\${pidcol_csv},\${pidcol_csv}nr \
      > "${allmappings}.sorted"

    {
      echo "\$header_csv"
      cat "${allmappings}.sorted"
    } > "${allmappings}.tmp"

    mv "${allmappings}.tmp" "${allmappings}"
    rm -f "${allmappings}.sorted"

    # count hits per qaccver(not including hybridisation)
    awk -F, -v col="\$snpcol_csv" -v hybcol="\$hybcol_csv" -v mh=${params.maximumhits} '
        NR>1 {
            snp = \$col
            hybrid = \$hybcol
            gsub("^[[:space:]]+|[[:space:]]+\$", "", snp)
            gsub("^[[:space:]]+|[[:space:]]+\$", "", hybrid)
            if (snp != "" && hybrid != "No") cnt[snp]++
        }
        END {
            for (k in cnt)
                if (cnt[k] > mh)
                    print k
        }
    ' "${allmappings}" > drop_snps_allmaps.txt

    echo "Alternate_SNP_ID,\$header_csv" > ${filteredmappings}

    awk -F, -v OFS=, \
        -v col_snp="\$snpcol_csv" \
        -v col_hyb="\$hybcol_csv" \
        -v dropfile="drop_snps_allmaps.txt" '
        BEGIN{
            while ((getline line < dropfile) > 0) {
                gsub("^[[:space:]]+|[[:space:]]+\$", "", line)
                if (line != "") drop[line] = 1
            }
        }
        NR==1 { next }  # header already written
        {
            # Trim whitespace on all fields
            for (i=1; i<=NF; i++)
                gsub("^[[:space:]]+|[[:space:]]+\$", "", \$i)

            snp = \$col_snp
            if (snp == "") next

            # SNP count over maxhits thresh for qaccver count
            if (snp in drop) next

            # Hybridised filter
            hyb = (col_hyb > 0 ? \$col_hyb : "")
            gsub("^[[:space:]]+|[[:space:]]+\$", "", hyb)
            if (col_hyb > 0 && hyb == "No") next

            # Per-SNP hit index for Alternate_SNP_ID
            idx[snp]++
            alt = snp "." idx[snp]

            # Write filtered row with Alternate_SNP_ID
            print alt, \$0
        }
    ' "${allmappings}" >> ${filteredmappings}

    rm -f drop_snps_allmaps.txt

    # Add Alternate_SNP_ID to full all_mappings CSV (so all csvs have the same columns regardless of analysis stage filtering stopped [maybe redundant but helps in case of returned empty datasets?])
    awk -F, -v OFS=, -v col_snp="\$snpcol_csv" '
      NR==1 {
        for (i=1; i<=NF; i++)
            gsub("^[[:space:]]+|[[:space:]]+\$", "", \$i)
        print "Alternate_SNP_ID", \$0
        next
      }
      {
        for (i=1; i<=NF; i++)
            gsub("^[[:space:]]+|[[:space:]]+\$", "", \$i)
        snp = \$col_snp
        if (snp == "") next
        idx[snp]++
        alt = snp "." idx[snp]
        print alt, \$0
      }
    ' "${allmappings}" > "${allmappings}.tmp"

    mv "${allmappings}.tmp" "${allmappings}"

    # Add Alternate_SNP_ID to unfiltered BLAST CSV (as above)

    header_raw=\$(head -n1 ${rawsallmappings})
    snpcol_raw=\$(echo "\$header_raw" | awk -F, '{for(i=1;i<=NF;i++) if(\$i=="qaccver") print i}')
    covcol_raw=\$(echo "\$header_raw" | awk -F, '{for(i=1;i<=NF;i++) if(\$i=="Coverage") print i}')
    pidcol_raw=\$(echo "\$header_raw" | awk -F, '{for(i=1;i<=NF;i++) if(\$i=="pident") print i}')

    if [ -z "\$snpcol_raw" ]; then
        echo "ERROR: qaccver column not found in ${rawsallmappings}" 1>&2
        exit 1
    fi
    if [ -z "\$covcol_raw" ] || [ -z "\$pidcol_raw" ]; then
        echo "ERROR: Coverage and/or pident column not found in ${rawsallmappings}" 1>&2
        exit 1
    fi


    tail -n +2 "${rawsallmappings}" \
      | sort -t, -k\${snpcol_raw},\${snpcol_raw} -k\${covcol_raw},\${covcol_raw}nr -k\${pidcol_raw},\${pidcol_raw}nr \
      > "${rawsallmappings}.sorted"

    {
      echo "\$header_raw"
      cat "${rawsallmappings}.sorted"
    } > "${rawsallmappings}.tmp1"

    mv "${rawsallmappings}.tmp1" "${rawsallmappings}"
    rm -f "${rawsallmappings}.sorted"

    awk -F, -v OFS=, -v col_snp="\$snpcol_raw" '
      NR==1 {
        for (i=1; i<=NF; i++)
            gsub("^[[:space:]]+|[[:space:]]+\$", "", \$i)
        print "Alternate_SNP_ID", \$0
        next
      }
      {
        for (i=1; i<=NF; i++)
            gsub("^[[:space:]]+|[[:space:]]+\$", "", \$i)
        snp = \$col_snp
        if (snp == "") next
        idx[snp]++
        alt = snp "." idx[snp]
        print alt, \$0
      }
    ' "${rawsallmappings}" > "${rawsallmappings}.tmp2"

    mv "${rawsallmappings}.tmp2" "${rawsallmappings}"

    # TSV stage

    header_tsv=\$(head -n1 ${mappings})
    snpcol_tsv=\$(echo "\$header_tsv" | awk '{for(i=1;i<=NF;i++) if(\$i=="SNP_ID") print i}')
    hybcol_tsv=\$(echo "\$header_tsv" | awk '{for(i=1;i<=NF;i++) if(\$i=="Hybridised") print i}')
    covcol_tsv=\$(echo "\$header_tsv" | awk '{for(i=1;i<=NF;i++) if(\$i=="Coverage") print i}')
    pidcol_tsv=\$(echo "\$header_tsv" | awk '{for(i=1;i<=NF;i++) if(\$i=="pident") print i}')

    if [ -z "\$snpcol_tsv" ] || [ -z "\$hybcol_tsv" ] || [ -z "\$covcol_tsv" ] || [ -z "\$pidcol_tsv" ]; then
        echo "ERROR: required columns (SNP_ID/Hybridised/Coverage/pident) not found in ${mappings}" 1>&2
        exit 1
    fi


    echo "\$header_tsv" | awk -v col_snp="\$snpcol_tsv" 'BEGIN{FS="[[:space:]]+"; OFS=" "}{
        first=1
        for (i=1; i<=NF; i++) {
            if (first) { first=0 } else { printf OFS }
            printf "%s", \$i
            if (i == col_snp) printf OFS "Alternate_SNP_ID"
        }
        printf ORS
    }' > ${filtered_mappings_tsv}

    tail -n +2 "${mappings}" \
      | sort -k\${snpcol_tsv},\${snpcol_tsv} -k\${covcol_tsv},\${covcol_tsv}nr -k\${pidcol_tsv},\${pidcol_tsv}nr \
      | awk -v FS="[[:space:]]+" -v OFS=" " \
            -v col_snp="\$snpcol_tsv" \
            -v col_hyb="\$hybcol_tsv" \
            -v mh=${params.maximumhits} '
        {
            snp = \$col_snp
            gsub("^[[:space:]]+|[[:space:]]+\$", "", snp)
            if (snp == "") next

            hyb = \$col_hyb
            gsub("^[[:space:]]+|[[:space:]]+\$", "", hyb)
            if (hyb == "No") next

            # limit to mh hits per SNP_ID (after Hybridised filter)
            idx[snp]++
            if (idx[snp] > mh) next

            alt = snp "." idx[snp]

            # Print row with Alternate_SNP_ID inserted after SNP_ID
            first=1
            for (i=1; i<=NF; i++) {
                if (first) { first=0 } else { printf OFS }
                if (i == col_snp) {
                    printf "%s", \$i        # SNP_ID
                    printf OFS "%s", alt    # Alternate_SNP_ID
                } else {
                    printf "%s", \$i
                }
            }
            printf ORS
        }
    ' >> ${filtered_mappings_tsv}
    """
}

/*
Step5: Process to merge individual mapping results into a single file
*/
process MERGE_MAPPINGS {
    publishDir "${params.resultsdirectory}/", mode: 'copy'
    label 'merge'
    tag 'Merge_all_CSVs'

    input:
        path allmappings_in          // list of *_all_mappings.csv
        path filteredmappings_in     // list of *_filtered_mappings.csv
        path rawsallmappings_in      // list of *_complete_unfiltered_blast_results.csv
        path resultsdir
        val  probename
        val  genomename

    output:
        tuple path(filteredmappings_out),
              path(allmappings_out),
              path(allmappingsunfiltered_out)

    script:
    // Pick one file from each group to grab the header
    def oneFilt = filteredmappings_in[0]
    def oneAll  = allmappings_in[0]
    def oneRaw  = rawsallmappings_in[0]

    // String lists of filenames for the for-loops
    def filtList = filteredmappings_in.collect { it.getName() }.join(' ')
    def allList  = allmappings_in.collect { it.getName() }.join(' ')
    def rawList  = rawsallmappings_in.collect { it.getName() }.join(' ')

    filteredmappings_out       = "${probename}_with_${genomename}_filtered_mappings.csv"
    allmappings_out            = "${probename}_with_${genomename}_all_mappings.csv"
    allmappingsunfiltered_out  = "${probename}_with_${genomename}_blastmappings_unfiltered.csv"

    """
    header_filtered=\$(head -n1 "${oneFilt.getName()}")
    header_all=\$(head -n1 "${oneAll.getName()}")
    header_raw=\$(head -n1 "${oneRaw.getName()}")

    echo "##brioche parameters" > header_params.txt
    echo "##evalue            = ${params.evalue} #The number of expected hits of similar quality (score) that could be found just by chance." >> header_params.txt
    echo "##minimumlength     = ${params.minlength} #minimum length of HSP to be considered fully hybridized" >> header_params.txt
    echo "##extendablebps     =${params.extendablebps}  #number of matching base pairs from the 3 prime end for a probe to be considered as extendable" >> header_params.txt
    echo "##maximumgaps       = ${params.maxgaps} #Maximum number of gaps allowed per HSP" >> header_params.txt
    echo "##coverage          = ${params.coverage} #Option to filter any hits with coverage less than the provided percentage" >> header_params.txt
    echo "##pident            = ${params.pident}  #Option to filter any hits with pident less than the provided percentage" >> header_params.txt
    echo "##maximumhits       = ${params.maximumhits} #Filter probes with hits more than maximumhits" >> header_params.txt

    # Merge filtered mappings (maxhits + Hybridized)
    {
      cat header_params.txt
      echo "\$header_filtered"
      for f in ${filtList}; do
        [ -s "\$f" ] || continue
        tail -n +2 "\$f"
      done
    } > "${filteredmappings_out}"

    # merge pident/coverage-filtered mappings
    {
      cat header_params.txt
      echo "\$header_all"
      for f in ${allList}; do
        [ -s "\$f" ] || continue
        tail -n +2 "\$f"
      done
    } > "${allmappings_out}"

    # Merge fully unfiltered BLAST mappings ( may need to optimise further depending on size limitations of sponge)
    {
      cat header_params.txt
      echo "\$header_raw"
      for f in ${rawList}; do
        [ -s "\$f" ] || continue
        tail -n +2 "\$f"
      done
    } > "${allmappingsunfiltered_out}"

    if command -v sponge >/dev/null 2>&1; then
      LC_ALL=C tr -d '\\t\\r' < "${filteredmappings_out}"      | sponge "${filteredmappings_out}"
      LC_ALL=C tr -d '\\t\\r' < "${allmappings_out}"           | sponge "${allmappings_out}"
      LC_ALL=C tr -d '\\t\\r' < "${allmappingsunfiltered_out}" | sponge "${allmappingsunfiltered_out}"
    fi
    """
}



/*
Step5: Process to merge individual mapping results into a single file
*/
process MERGE_MIN_MAPPINGS {
    publishDir "${params.resultsdirectory}/", mode: 'copy'
    label 'merge'
    tag 'Merge_all_TSVs'

    input:
        path filtered_mappings_tsv_in   // *_filtered_mappings.tsv from ADD_MARKERINFO
        path resultsdir
        val  probename
        val  genomename

    output:
        path(outfile)

    script:
    // Use one file to grab the TSV header
    def oneFiltTsv = filtered_mappings_tsv_in[0]
    outfile        = "${probename}_with_${genomename}_mappings.tsv"

    """
    header_tsv=\$(head -n1 "${oneFiltTsv}")

    echo "##brioche parameters" > header_params.txt
    echo "##evalue            = ${params.evalue} #The number of expected hits of similar quality (score) that could be found just by chance." >> header_params.txt
    echo "##minimumlength     = ${params.minlength} #minimum length of HSP to be considered fully hybridized" >> header_params.txt
    echo "##extendablebps     =${params.extendablebps}  #number of matching base pairs from the 3 prime end for a probe to be considered as extendable" >> header_params.txt
    echo "##maximumgaps       = ${params.maxgaps} #Maximum number of gaps allowed per HSP" >> header_params.txt
    echo "##coverage          = ${params.coverage} #Option to filter any hits with coverage less than the provided percentage" >> header_params.txt
    echo "##pident            = ${params.pident}  #Option to filter any hits with pident less than the provided percentage" >> header_params.txt
    echo "##maximumhits       = ${params.maximumhits} #Filter probes with hits more than maximumhits" >> header_params.txt

    {
      cat header_params.txt
      echo "\$header_tsv"
      for f in *_filtered_mappings.tsv; do
        [ -s "\$f" ] || continue
        # Skip each file's header row
        tail -n +2 "\$f"
      done
    } | awk 'BEGIN{
               FS="[[:space:]]+";
               OFS="\\t"
             }
             /^##/ { print; next }
             {
               # Convert whitespace-delimited header/body to tab-delimited
               n = split(\$0, f, /[[:space:]]+/);
               for (i=1; i<=n; i++) {
                 if (i>1) printf OFS;
                 printf "%s", f[i];
               }
               printf ORS;
             }' > "${outfile}"
    """
}


/*
Step6: Perform filtering and generate annotations
*/
process INTERMEDIATE_FILTERING {
    tag "Apply_intermediate_filters"
    label 'large'
    publishDir "${params.resultsdirectory}/", mode: 'copy', pattern: '*_intermediate_filtering_*'
    input:
        path filteredmappings
        path outfile
        val resultsdir
        val probename
        val genomename
        val chrominfo
        val chrompath
        val targetinf
        val usemapmarkers
        val similarmarkersmaps
        val useldedge
        val ldmapp
        val usegeneticmap
        val geneticmap
        val dupdist
        val keepdups
    output:
        tuple path("${filteredmappretzelcsv}"), path("${intermediatefilteredmap}"), path("${dupmapinter}")
    script:
    filteredmappretzelcsv="${probename}_with_${genomename}_intermediate_filtering_hits.csv"
    intermediatefilteredmap="${probename}_with_${genomename}_intermediate_filtering_mappings.tsv"
    dupmapinter="${probename}_with_${genomename}_marker_localduplications_counts.tsv"
    """
    R --version
        Rscript -e "briocheR::DoIntermediatefiltering(
            probe.name='${probename}',
            genome.name='${genomename}',
            output.path='${resultsdir}',
            dotargetfilt='${chrominfo}',
            chrom.comp='${chrompath}',
            target.data='${targetinf}',
            domapmarkers='${usemapmarkers}',
            similarity.maps='${similarmarkersmaps}',
            doldedge='${useldedge}',
            ldmapp='${ldmapp}',
            dogeneticmap='${usegeneticmap}',
            geneticmap.file='${geneticmap}',
            blast.hits='${filteredmappings}',
            mappings.file='${outfile}',
            dup.dist='${dupdist}',
            keeplocalduppos='${keepdups}')"
    #mv "name of intermediate mapping blasts" filteredmappretzelcsv 
    #mv "name of intermediate mapping maps" intermediatefilteredmap
    """
}

/*
Step7: Perform advanced filtering and generate annotations
*/
process ADVANCED_FILTERING {
    tag "Apply_strict_filters"
    label 'large'
    publishDir "${params.resultsdirectory}/", mode: 'copy'

    input:
        path filteredmappretzelcsv
        path intermediatefilteredmap
        path dupmapinter
        val  resultsdir
        val  probename
        val  genomename
        val  chrominfo
        val  chrompath
        val  targetinf
        val  usemapmarkers
        val  similarmarkersmaps
        val  useldedge
        val  ldmapp
        val  usegeneticmap
        val  geneticmap

    output:
        tuple path("${strictmappedcsv}"),
              path("${strictmappedctsv}"),
              path("${pretzelfile}"),
              path("${priors_informatives_strict}")

    script:
    strictmappedcsv   = "${probename}_with_${genomename}_strict_filtering_hits.csv"
    strictmappedctsv  = "${probename}_with_${genomename}_strict_filtering_mappings.tsv"
    pretzelfile       = "${probename}_with_${genomename}_mappings-pretzel-alignment.xlsx"
    priors_informatives_strict  = "${probename}_with_${genomename}_priors_informed_strictmapping.tsv"
    """
    R --version
      Rscript -e "briocheR::DoStrictfiltering(
            probe.name='${probename}',
            genome.name='${genomename}',
            output.path='${resultsdir}',
            dotargetfilt='${chrominfo}',
            chrom.comp='${chrompath}',
            target.data='${targetinf}',
            domapmarkers='${usemapmarkers}',
            similarity.maps='${similarmarkersmaps}',
            doldedge='${useldedge}',
            ldmapp='${ldmapp}',
            blast.hits='${filteredmappretzelcsv}',
            mappings.file='${intermediatefilteredmap}',
            dogeneticmap='${usegeneticmap}',
            geneticmap.file='${geneticmap}',
            dupmapinter.file='${dupmapinter}')"
    #mv "name of intermediate mapping blasts" filteredmappretzelcsv 
    #mv "name of intermediate mapping maps" intermediatefilteredmap
    #Generate pretzel files
    echo  ${pretzelfile}
    #Format of alignment file
    #Name	Chromosome	Start	End	#Comments
    #Rows which start with # such as this one are filtered out and not loaded into the database.
    #Columns with a column header which starts with #, e.g. #Comments in column 5, are filtered out and not loaded into the database.
    #The given column headings are the minimum for an Alignment dataset;  other columns may be added and will be loaded into the database if they have a column heading.
    #There values will be stored in Feature.values.columnName, and will be displayed in the brushed/selected Features table in the right panel of Pretzel.
    displayName="ChangeMe-Genome-${genomename}-${probename}-Markers"
    parent="ChangeMe"
    commonName="ChangeMe"
    Description="${probename.replaceAll('_',' ')} markers mapped to ${genomename} using brioche"
    datasetname=\$(echo "Alignment|${genomename}_${probename}"|awk '{print substr(\$0,1,30)}')
    Rscript -e "
    data<-read.delim('${intermediatefilteredmap}',header=T,sep='\\\\t',comment.char = '#',stringsAsFactors=F,strip.white=T);
    cols=colnames(data)[c(1,4,5,5,6,7)]
    cols2=colnames(data)
    print(cols)
    cols2=cols2[!cols2 %in% cols]
    data=cbind(data[,cols],data[,cols2])
    #rearrange columns
    colnames(data)[1:6]<-c('Name','Chromosome','Start','End','REF','ALT');
    # Identify duplicate columns
    duplicate_cols <- colnames(data)[duplicated(lapply(data, function(x) as.character(x)))]
    # Keep only unique columns
    keep.cols<-!colnames(data) %in% duplicate_cols
    keep.cols[colnames(data) %in% 'End']=T
    data <- data[, keep.cols]
    if(!require('openxlsx'))
    {
        install.packages('openxlsx',repos = 'https://cloud.r-project.org');
    }
    data <-list('\${datasetname}'=data,Metadata = data.frame(
            Field = c('displayName','parentName', 'species', 'Description','citation','Data source','Uploaded by','Contact','Licencing')));
            data[['Metadata']][['\${datasetname}']] <-c('\${displayName}','\${commonName}','\${parent}','\${Description}','Change me','Change me','Change me','Change me','Change me');
    openxlsx::write.xlsx(data, file='${pretzelfile}')
    "
    """
}

process COLLECT_SUMSTATS {
  tag "Collect_results"
  label 'large'
  publishDir "${params.resultsdirectory}/", mode: 'copy'

  input:
    path strictmappedcsv
    path strictmappedtsv
    path dupmapinter
    path priors_informatives_strict
    val  resultsdir
    val  probename
    val  genomename
    val  targetdesign
    val  brioche_version
    val  coverage
    val  pident
    val  otherblastoptions
    val  used_targetchrom
    val  used_sharedmap
    val  used_ldedge
    val  used_geneticmap
    val  brioche_repo_url
    val  brioche_branch      
  output:
      tuple path("${summaryfiltering}"),
            path("${duphitsfile}")
  script:
  summaryfiltering = "${probename}_with_${genomename}_summary_filtering.csv"
  duphitsfile = "${probename}_with_${genomename}_marker_localdups_NULLS_counts.tsv"
  """
  R --version
    Rscript -e "briocheR::Dosumstats(
      probe.name='${probename}',
      genome.name='${genomename}',
      markers_priorsinformed.hits='${priors_informatives_strict}',
      duplicationmap.hits='${dupmapinter}',
      output.path='${resultsdir}',
      blast.hits='${strictmappedcsv}',
      mappings.file='${strictmappedtsv}',
      target.file='${targetdesign}',
      metadata=list(
        brioche_version='${brioche_version}',
        brioche_repo_url='${brioche_repo_url}',
        brioche_commit_url='${brioche_branch}',
        coverage=as.numeric(${coverage}),
        pident=as.numeric(${pident}),
        otherblastoptions='${otherblastoptions.replace("'", "\\'")}',
        usetargetchrom='${used_targetchrom}',
        usesharedmarkersmap='${used_sharedmap}',
        useldedgemap='${used_ldedge}',
        usegeneticmap='${used_geneticmap}',
        run_date=Sys.time()
      )
    )"
  """
}
/*
Step9: make metadata 
*/
process MAKE_RUN_META {
  tag "make_run_meta"
  label 'small'
  publishDir "${params.resultsdirectory ?: 'results'}", mode: 'copy'

  input:
  val run_info

  output:
  path "run_meta.json"

  // Use base64 to avoid any heredoc edge cases (quotes/newlines)
  script:
  """
  META_B64='${groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(run_info)).bytes.encodeBase64().toString()}'
  printf "%s" "\$META_B64" | base64 -d > run_meta.json
  """
}

/*
Step10: Build summary 
*/
process BUILD_SUMMARY_CORE {
  tag "summary_core"
  label 'small'
  publishDir "${params.resultsdirectory ?: 'results'}/summary", mode: 'copy'

  input:
  path META_JSON
  path FIGS_DIR
  path SUMMARY_DONE

  output:
  path "summary_core.html"
  path "summary_core.Rmd"

  // Precompute outdir to avoid nested quotes in the script GString
  def outdir_fallback = params.resultsdirectory ?: projectDir + '/results'

  script:
  """
  set -euo pipefail

  mkdir -p figures

  # Copy any images from the figs directory (may be empty; that's fine)
  if [ -d "${FIGS_DIR}" ]; then
      shopt -s nullglob dotglob
      files=( "${FIGS_DIR}"/* )
      if [ \${#files[@]} -gt 0 ]; then
          cp -r "\${files[@]}" figures/ || true
      fi
  fi

  Rscript ${projectDir}/Scripts/summarise_results.R \
      "${META_JSON}" \
      "${outdir_fallback}" \
      summary_core.html \
      --fig_dir figures \
      --rmd_out summary_core.Rmd
  """
}



/*
Step11: Build summary 
*/

process BUILD_SUMMARY_WRAPPER {
  tag "summary_wrapper"
  label 'small'
  publishDir "${params.resultsdirectory ?: 'results'}/summary", mode: 'copy'

  input:
    path CORE_HTML               // from BUILD_SUMMARY_CORE
    val  TS                      // params.ts
    val  TS2                     // params.ts2                         // params.ts2
  output:
    path "summary_report.html"

  script:
  """
  # Write a lightweight HTML that references existing artefacts by name
  cat > summary_report.html <<EOF
<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>Pipeline Summary (Final)</title>
<style>
.tabs{display:flex;gap:12px;margin:12px 0}
.tab{padding:6px 12px;border:1px solid #ccc;border-bottom:none;cursor:pointer;background:#f7f7f7}
.tab.active{background:#fff;font-weight:bold}
.panel{display:none;border:1px solid #ccc;padding:10px}
.panel.active{display:block}
iframe,embed{border:1px solid #ddd;width:100%;height:900px}
</style>
</head>
<body>
<h1>Pipeline Summary (Final)</h1>
<div class="tabs">
  <div class="tab active" data-target="core">Overview (core)</div>
  <div class="tab" data-target="timeline">Timeline</div>
  <div class="tab" data-target="report">Report</div>
  <div class="tab" data-target="dag">DAG</div>
</div>

<div id="core" class="panel active">
  <iframe src="./summary_core.html" width="100%" height="1200px"></iframe>
</div>

<div id="timeline" class="panel">
  <iframe src="../Reports/${TS}_timeline.html"></iframe>
  <p>Open directly: <a href="../Reports/${TS}_timeline.html">${TS}_timeline.html</a></p>
</div>

<div id="report" class="panel">
  <iframe src="../Reports/${TS}_report.html"></iframe>
  <p>Open directly: <a href="../Reports/${TS}_report.html">${TS}_report.html</a></p>
</div>

<div id="dag" class="panel">
  <iframe src="../Reports/${TS2}_dag.html"></iframe>
  <p>Open directly: <a href="../Reports/${TS2}_dag.html">${TS2}_dag.html</a></p>
</div>

<script>
document.querySelectorAll('.tab').forEach(t=>{
  t.addEventListener('click',()=>{
    document.querySelectorAll('.tab').forEach(x=>x.classList.remove('active'));
    document.querySelectorAll('.panel').forEach(x=>x.classList.remove('active'));
    t.classList.add('active');
    document.getElementById(t.dataset.target).classList.add('active');
  });
});
</script>
</body>
</html>
EOF
  """
}
