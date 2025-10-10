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
  tag "${file(fasta).baseName.replaceAll('-split','')}"
  label 'large'
  input:
    path fasta
    val ref_origin
    val genomename
  output:
    val "${blastdatabase}"

  script:
  fastadir    = new File("${ref_origin}").getParent()
  blastsubdir = new File(fastadir + "/BLAST/")

  if( !new File(fastadir).canWrite() || (blastsubdir.exists() && !blastsubdir.canWrite()) ) {
    System.err.println("Error: No write access to ${fastadir}/BLAST; falling back to ${projectDir}/BLAST/${genomename}")
    blastdb = new File("${projectDir}/BLAST/${genomename}").toString()
  } else {
    blastdb = blastsubdir.toString()
  }

  blastdatabase = "${file(fasta).getName()}".replaceAll("-split.*","")
  refgenomename = "${new File(ref_origin).getName()}"

  if( blastdatabase != refgenomename ) {
    blastdb = new File(blastdb + "/brioche-blastdb/").toString()
  }
  if( !new File(blastdb).exists() ) {
    new File(blastdb).mkdirs()
  }
  blastdatabase = blastdb.toString() + "/" + blastdatabase

  """
  set -euo pipefail
  if [[ ! -e "${blastdatabase}.ndb" ]]; then
    echo "creating database for ${fasta} -> ${blastdatabase}"
    makeblastdb -in "${fasta}" -out "${blastdatabase}" -dbtype nucl -title "${params.genomename}"
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
  tag "${file(blastdb).baseName}"
  label 'large'
  input:
    path probefasta
    val blastdb
    val blastoutformat
  output:
        path "${splitname}*"
  script:
  splitname=file("${blastdb}").baseName
  """
  blastn -db ${blastdb} -query ${probefasta} \\
         -out blastoutput.txt \\
         -num_threads $task.cpus \\
         -outfmt ${blastoutformat} \\
         -evalue ${params.evalue} ${params.otherblastoptions} \\
         -dust ${params.dust}
  split -l 50000 blastoutput.txt ${splitname}
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
    val targetdesign
    val minmismatches
    val extendablebps
    path genomefasta
    val blastoutformat
    val istarget3primeend
    val markerchar
  output:
    path "${rds}"
    path "${snpmapped}"
  script:
  snpmapped="${blastout.baseName}_samtools_coordinates_mapped.tsv"
  rds="SNPcoordinates/${blastout.baseName}_blast_results.RDS"
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
  mv  \${snpcoords}/processed_blast_results.RDS ${rds}
  samtools faidx ${genomefasta} -r "\${snpcoords}/samtools_marker_positions.txt"|awk '
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
  ' >"${snpmapped}"
  """
}


/*
Step4: Perform filtering and generate annotations
*/
process ADD_MARKERINFO{
    tag "${blast_results.baseName}"
    input:
        path blast_results
        path coordinates
        val probename
        val genomename
        val keepduplicates
        val coveragefilter
        val pidentfilter
        val istarget3primeend
    output:
        tuple path("${allmappings}"), path("${mappings}"), path("${rawsallmappings}")
    script:
    def random_factor = new Random().nextInt()
    allmappings="${blast_results.baseName}_"+random_factor+"_all_mappings.csv"
    mappings="${blast_results.baseName}_"+random_factor+"_mapping.tsv"
    rawsallmappings="${blast_results.baseName}_"+random_factor+"_complete_unfiltered_blast_results.csv"
    """
    R --version
    Rscript -e "briocheR::addSNPdetails(
               blast.path='${blast_results}',
               reference.bases='${coordinates}',
               probe.name='${probename}',
               genome.name='${genomename}',
               output.path = '\$(pwd)',
               pident.filter=${pidentfilter},
               coverage.filter= ${coveragefilter},
               is3prime = ${istarget3primeend})"
    mv "${probename}_with_${genomename}_all_mappings.csv" ${allmappings}
    mv "${probename}_with_${genomename}_mapping.tsv" ${mappings}
    mv "${probename}_with_${genomename}_complete_unfiltered_blast_results.csv" ${rawsallmappings}
    """
}
/*
Step5: Process to merge individual mapping results into a single file
*/
process MERGE_MAPPINGS {
    publishDir "${params.resultsdirectory}/",mode:'copy'
    label 'merge'
    tag 'Merge_all_CSVs'
    input:
        path allmappings
        path rawsallmappings
        path resultsdir
        val probename
        val genomename
    output:
        tuple path("${filteredmappings}"), path("${allmappings}"), path("${allmappingsunfiltered}")
    script:
    allmaps=allmappings[0]
    filteredmappings="${probename}_with_${genomename}_filtered_mappings.csv"
    allmappingsunfiltered="${probename}_with_${genomename}_blastmappings_unfiltered.csv"
    allmappings="${probename}_with_${genomename}_all_mappings.csv"
    """
    #-k\${SNP_ID},\${SNP_ID} -k\${coverage},\${coverage}r -k\${pident},\${pident}r -k\${chrom},\${chrom} -k\${pos},\${pos}
    header=\$(head -n1 $allmaps)
    chrom=\$(echo \$header|awk -F, '{for(i=1;i<=NF;i++){if(\$i=="saccver")print i}}')
    pos=\$(echo \$header|awk -F, '{for(i=1;i<=NF;i++){if(\$i=="sstart")print i}}')
    SNP_ID=\$(echo \$header|awk -F, '{for(i=1;i<=NF;i++){if(\$i=="qaccver")print i}}')
    pident=\$(echo \$header|awk -F, '{for(i=1;i<=NF;i++){if(\$i=="pident")print i}}')
    coverage=\$(echo \$header|awk -F, '{for(i=1;i<=NF;i++){if(\$i=="Coverage")print i}}')
    find -name "*all_mappings.csv"|split -l 2000 - subsetcsvs
    for ss in subsetcsvs*;
    do
        cat \$(awk '{printf("%s ",\$0)}' \$ss)|grep -v "qaccver" >> Mergedmappings.csv
    done

    sort -t, -k\${SNP_ID},\${SNP_ID} -k\${coverage},\${coverage}r -k\${pident},\${pident}r -k\${chrom},\${chrom} -k\${pos},\${pos} Mergedmappings.csv \\
    |awk -F, -v h="\$header" 'BEGIN{print "Alternate_SNP_ID,"h} {count[\$1]++;print \$1"."count[\$1]","\$0}' \\
    >allmaps.csv

    find -name "*unfiltered*.csv"|split -l 2000 - subsetcsvs2
    for ss in subsetcsvs2*;
    do
        cat \$(awk '{printf("%s ",\$0)}' \$ss)|grep -v "qaccver" >> Mergedmappingsunfilts.csv
    done

    sort -t, -k\${SNP_ID},\${SNP_ID} -k\${coverage},\${coverage}r -k\${pident},\${pident}r -k\${chrom},\${chrom} -k\${pos},\${pos} Mergedmappingsunfilts.csv \\
    |awk -F, -v h="\$header" 'BEGIN{print "Alternate_SNP_ID,"h} {count[\$1]++;print \$1"."count[\$1]","\$0}' \\
    >allmapsunfilt.csv 
    awk -F, -v OFS=, -v maxhit="${params.maximumhits}" 'NR==1{for(i=1;i<=NF;i++)if(\$i=="Hybridized")hyb=i; if(!hyb){print "ERROR: Hybridized column not found" > "/dev/stderr"; exit 1} header=\$0; next}{split(\$1,a,".");pref=a[1];cnt[pref]++;row[NR]=\$0;p[NR]=pref;h[NR]=\$hyb} END{print header; for(i=2;i<=NR;i++) if(cnt[p[i]]<=maxhit && h[i]!="No ") print row[i]}' allmaps.csv > map.csv
    #awk -F, -v maxhit="${params.maximumhits}" 'BEGIN{split(\$1,a,".");data[\$1]=\$0;IDs[a[1]]=0} {split(\$1,a,".");IDs[a[1]]++;data[\$1]=\$0} END{for(n in data){split(n,b,".");if(IDs[b[1]]<=maxhit)print data[n]}}' allmaps.csv > map.csv   
    LC_ALL=C tr -d ' \t\r' < allmapsunfilt.csv | sponge allmapsunfilt.csv
    LC_ALL=C tr -d ' \t\r' < allmaps.csv | sponge allmaps.csv
    LC_ALL=C tr -d ' \t\r' < map.csv | sponge map.csv
    echo "##brioche parameters" > header.txt
    echo "##evalue            = ${params.evalue} #The number of expected hits of similar quality (score) that could be found just by chance." >> header.txt
    echo "##minimumlength     = ${params.minlength} #minimum length of HSP to be considered fully hybridized" >> header.txt
    echo "##extendablebps     =${params.extendablebps}  #number of matching base pairs from the 3 prime end for a probe to be considered as extendable" >> header.txt
    echo "##maximumgaps       = ${params.maxgaps} #Maximum number of gaps allowed per HSP" >> header.txt
    echo "##coverage          = ${params.coverage} #Option to filter any hits with coverage less than the provided percentage" >> header.txt
    echo "##pident            = ${params.pident}  #Option to filter any hits with pident less than the provided percentage" >> header.txt
    echo "##maximumhits       = ${params.maximumhits} #Filter probes with hits more than maximumhits" >> header.txt
    cat header.txt map.csv >"${filteredmappings}"
    cat header.txt allmaps.csv>"${allmappings}"
    cat header.txt allmapsunfilt.csv>"${allmappingsunfiltered}"
    """
}


/*
Step5: Process to merge individual mapping results into a single file
*/
process MERGE_MIN_MAPPINGS {
    publishDir "${params.resultsdirectory}/",mode:'copy'
    label 'merge'
    tag 'Merge_all_TSVs'
    input:
        path mappings
        val resultsdir
        val probename
        val genomename
    output:
        path("${outfile}")
    script:
    allmaps=mappings[0]
    outfile="${probename}_with_${genomename}_mappings.tsv"
    """
    header=\$(head -n1 $allmaps)
    chrom=\$(echo \$header|awk '{for(i=1;i<=NF;i++){if(\$i=="Chr")print i}}')
    pos=\$(echo \$header|awk '{for(i=1;i<=NF;i++){if(\$i=="SNP_position")print i}}')
    SNP_ID=\$(echo \$header|awk '{for(i=1;i<=NF;i++){if(\$i=="SNP_ID")print i}}')
    pident=\$(echo \$header|awk '{for(i=1;i<=NF;i++){if(\$i=="pident")print i}}')
    header=\$(echo \$header|awk '{for(i=1;i<=NF;i++){if(i==1){printf \$i"\tAlternate_SNP_ID"}else{printf "\t"\$i}}}')
    coverage=\$(echo \$header|awk '{for(i=1;i<=NF;i++){if(\$i=="Coverage")print i}}')
    find -name "*.tsv"|split -l 2000 - subsettsvs
    for ss in subsettsvs*;
    do
        cat \$(awk '{printf("%s ",\$0)}' \$ss)|grep -v "SNP_ID" >> Mergedmappings.tsv
    done
    sort -k\${SNP_ID},\${SNP_ID} -k\${coverage},\${coverage}r -k\${pident},\${pident}r -k\${chrom},\${chrom} -k\${pos},\${pos} Mergedmappings.tsv \\
    |awk 'BEGIN{count[\$1]=0} {count[\$1]++;for(i=1;i<=NF;i++){if(i==1){printf\$i"\t"\$i"."count[\$1]}else{printf"\t"\$i}};print""}' \\
    |awk -v maxhit=${params.maximumhits} -v h="\$header" 'BEGIN{print h; split(\$2,a,".");data[\$2]=\$0;count[\$1]=0} {count[\$1]++;data[\$2]=\$0} END{for(n in data){split(n,a,".");if(count[a[1]]<=maxhit){print data[n]}}}' \\
    |awk '{if(\$0!=""){print \$0}}' \\
    > map.tsv
    echo "##brioche parameters" > header.txt
    echo "##evalue            = ${params.evalue} #The number of expected hits of similar quality (score) that could be found just by chance." >> header.txt
    echo "##minimumlength     = ${params.minlength} #minimum length of HSP to be considered fully hybridized" >> header.txt
    echo "##extendablebps     =${params.extendablebps}  #number of matching base pairs from the 3 prime end for a probe to be considered as extendable" >> header.txt
    echo "##maximumgaps       = ${params.maxgaps} #Maximum number of gaps allowed per HSP" >> header.txt
    echo "##coverage          = ${params.coverage} #Option to filter any hits with coverage less than the provided percentage" >> header.txt
    echo "##pident            = ${params.pident}  #Option to filter any hits with pident less than the provided percentage" >> header.txt
    echo "##maximumhits       = ${params.maximumhits} #Filter probes with hits more than maximumhits" >> header.txt
    cat header.txt map.tsv >"${outfile}"  
    """
}

/*
Step6: Perform filtering and generate annotations
*/
process INTERMEDIATE_FILTERING {
    tag "Apply_intermediate_filters"
    label 'large'
    publishDir "${params.resultsdirectory}/",mode:'copy'
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
    output:
        tuple path("${filteredmappretzelcsv}"), path("${intermediatefilteredmap}")
    script:
    filteredmappretzelcsv="${probename}_with_${genomename}_intermediate_filtering_hits.csv"
    intermediatefilteredmap="${probename}_with_${genomename}_intermediate_filtering_mappings.tsv"
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
            blast.hits='${filteredmappings}',
            mappings.file='${outfile}')"
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
    publishDir "${params.resultsdirectory}/",mode:'copy'
    input:
        path filteredmappretzelcsv
        path intermediatefilteredmap
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
    output:
        tuple path("${strictmappedcsv}"), path("${strictmappedctsv}"), path("${pretzelfile}")
    script:
    strictmappedcsv="${probename}_with_${genomename}_strict_filtering_hits.csv"
    strictmappedctsv="${probename}_with_${genomename}_strict_filtering_mappings.tsv"
    pretzelfile="${probename}_with_${genomename}_mappings-pretzel-alignment.xlsx"
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
            mappings.file='${intermediatefilteredmap}')"
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

/*
Step8: Gather summary stats for filtering
*/
process COLLECT_SUMSTATS {
    tag "Collect_results"
    label 'large'
    publishDir "${params.resultsdirectory}/",mode:'copy'
    input:
        path strictmappedcsv
        path strictmappedtsv
        val resultsdir
        val probename
        val genomename
        val targetdesign 
    output:
        path("${summaryfiltering}")
    script:
    summaryfiltering="${probename}_with_${genomename}_summary_filtering.csv"
    """
    R --version
    Rscript -e "briocheR::Dosumstats(
            probe.name='${probename}',
            genome.name='${genomename}',
            output.path='${resultsdir}',
            blast.hits='${strictmappedcsv}',
            mappings.file='${strictmappedtsv}',
            target.file='${targetdesign}')"
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
