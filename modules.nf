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

process BUILD_BLASTDb {
  tag "${file(fasta).baseName.replaceAll('-split','')}"
  label 'large'
  input:
    val fasta
    val refgenome
  output:
        val "${blastdatabase}"
  script:
  fastadir=new File("${refgenome}").getParent()
  blastsubdir=new File(fastadir+"/BLAST/")

  // Check: if FASTA dir not writable or BLAST subdir exists but not writable ? fallback
  if( !new File(fastadir).canWrite() || (blastsubdir.exists() && !blastsubdir.canWrite()) )
  {
      System.err.println("Error: No write access to ${fastadir}/BLAST; falling back to ${projectDir}/BLAST")
      blastdb=new File("${projectDir}/BLAST").toString()
  }
  else
  {
      blastdb=blastsubdir.toString()
  }

  blastdatabase="${file(fasta).getName()}".replaceAll("-split.*","")
  refgenomename="${file(refgenome).getName()}"

  if(blastdatabase!=refgenomename)
  {
      blastdb=new File(blastdb+"/brioche-blastdb/").toString()
  }
  if(!new File(blastdb).exists())
  {
      new File(blastdb).mkdirs()
  }
  blastdatabase=blastdb.toString()+"/"+blastdatabase

  """
  #check if blastdb exists
  if [[ ! -e "${blastdatabase}.ndb" ]]
  then
    echo creating database
    makeblastdb -in $fasta -out ${blastdatabase} -dbtype "nucl" -title "${params.genomename}"
  fi
  """
}
 

process SPLIT_REFGENOME{
    tag "${file(refgenome).baseName}"
    input:
        val refgenome
        val resultsdir
        val chromstoexclude
        val buildblastdbonly
    output:
        path "*split.fasta"
    script:
    path2briocheR=params.pathtobriocher
    genomefasta="${file(refgenome).getName()}-split.fasta"
    """
    mkdir -p ${resultsdir}

    #Check if reference dictionary exists, if not, build one
    dict=\$(echo ${refgenome}|sed -e "s/.[^.]*\$/.dict/")
    if [[ ! -e "\${dict}" ]]
    then
        picard -Xms4g -Xmx10g CreateSequenceDictionary -R ${refgenome} -O \${dict}
    else
        echo "Dictionary file found"
    fi

    exclude_array="DC"
    if [[ ! -z "${chromstoexclude}" ]];
    then
        exclude_array=\$(echo ${chromstoexclude}|sed 's/,/|/g')
    fi

    chromosomes=\$(grep "@SQ" \${dict}|cut -f2|sed -e "s/SN://g"|grep -v -P "^\${exclude_array}")
    split=\$(echo \${chromosomes}|awk -v maxchrom=30 'NF > maxchrom {print "no"}; NF < maxchrom {print "yes"}')
    if [[ "\${split}" == "yees" && "$buildblastdbonly" == "false" ]];
    then
        if [[ ! -z "\${chromosomes}" ]];
        then
            for chr in \${chromosomes};
            do
                samtools faidx $refgenome "\$chr" > "\${chr}-split.fasta"
            done
        else
            echo "No chromosomes found in fasta file"
            exit 1
        fi
    elif [[ ! -z "${params.restrict2chrom}" ]]
    then
        samtools faidx $refgenome "${params.restrict2chrom}" > "${params.restrict2chrom}-split.fasta"
    else
            ln -s ${refgenome} ${genomefasta}
    fi
    echo "Try and install briocheR that has been updated"
    #install briocheR if not installed
    Rscript -e "if(!require('briocheR')){
                install.packages('${path2briocheR}',repos=NULL);
                print('BriocheR installed')}"
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
    output:
        tuple path("${allmappings}"), path("${mappings}")
    script:
    def random_factor = new Random().nextInt()
    allmappings="${blast_results.baseName}_"+random_factor+"_all_mappings.csv"
    mappings="${blast_results.baseName}_"+random_factor+"_mapping.tsv"
    """
    R --version
    Rscript -e "briocheR::addSNPdetails(
               blast.path='${blast_results}',
               reference.bases='${coordinates}',
               probe.name='${probename}',
               genome.name='${genomename}',
               output.path = '\$(pwd)',
               pident.filter=${pidentfilter},
               coverage.filter= ${coveragefilter})"
    mv "${probename}_with_${genomename}_all_mappings.csv" ${allmappings}
    mv "${probename}_with_${genomename}_mapping.tsv" ${mappings}
    """
}
/*
Step5: Process to merge individual mapping results into a single file
*/
process MERGE_MAPPINGS {
    publishDir "${params.resultsdirectory}/",mode:'move'
    label 'merge'
    tag 'Merge_all_CSVs'
    input:
        path allmappings
        path resultsdir
        val probename
        val genomename
    output:
        tuple path("${filteredmappings}"), path("${allmappings}")
    script:
    allmaps=allmappings[0]
    filteredmappings="${probename}_with_${genomename}_filtered_mappings.csv"
    allmappings="${probename}_with_${genomename}_all_mappings.csv"
    """
    #-k\${SNP_ID},\${SNP_ID} -k\${coverage},\${coverage}r -k\${pident},\${pident}r -k\${chrom},\${chrom} -k\${pos},\${pos}
    header=\$(head -n1 $allmaps)
    chrom=\$(echo \$header|awk -F, '{for(i=1;i<=NF;i++){if(\$i=="saccver")print i}}')
    pos=\$(echo \$header|awk -F, '{for(i=1;i<=NF;i++){if(\$i=="sstart")print i}}')
    SNP_ID=\$(echo \$header|awk -F, '{for(i=1;i<=NF;i++){if(\$i=="qaccver")print i}}')
    pident=\$(echo \$header|awk -F, '{for(i=1;i<=NF;i++){if(\$i=="pident")print i}}')
    coverage=\$(echo \$header|awk -F, '{for(i=1;i<=NF;i++){if(\$i=="Coverage")print i}}')
    find -name "*.csv"|split -l 100 - subsetcsvs
    for ss in subsetcsvs*;
    do
        cat \$(awk '{printf("%s ",\$0)}' \$ss)|grep -v "qaccver" >> Mergedmappings.csv
    done

    sort -t, -k\${SNP_ID},\${SNP_ID} -k\${coverage},\${coverage}r -k\${pident},\${pident}r -k\${chrom},\${chrom} -k\${pos},\${pos} Mergedmappings.csv \\
    |awk -F, -v h="\$header" 'BEGIN{print "Alternate_SNP_ID,"h} {count[\$1]++;print \$1"."count[\$1]","\$0}' \\
    >allmaps.csv

    awk -F, -v maxhit="${params.maximumhits}" 'BEGIN{split(\$1,a,".");data[\$1]=\$0;IDs[a[1]]=0} {split(\$1,a,".");IDs[a[1]]++;data[\$1]=\$0} END{for(n in data){split(n,b,".");if(IDs[b[1]]<=maxhit)print data[n]}}' allmaps.csv > map.csv
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
    """
}


/*
Step5: Process to merge individual mapping results into a single file
*/
process MERGE_MIN_MAPPINGS {
    publishDir "${params.resultsdirectory}/",mode:'move'
    label 'merge'
    tag 'Merge_all_TSVs'
    input:
        path mappings
        val resultsdir
        val probename
        val genomename
    output:
        tuple path("${outfile}"),
        path("${pretzelfile}")
    script:
    allmaps=mappings[0]
    outfile="${probename}_with_${genomename}_mappings.tsv"
    pretzelfile="${probename}_with_${genomename}_mappings-pretzel-alignment.xlsx"
    """
    header=\$(head -n1 $allmaps)
    chrom=\$(echo \$header|awk '{for(i=1;i<=NF;i++){if(\$i=="Chr")print i}}')
    pos=\$(echo \$header|awk '{for(i=1;i<=NF;i++){if(\$i=="SNP_position")print i}}')
    SNP_ID=\$(echo \$header|awk '{for(i=1;i<=NF;i++){if(\$i=="SNP_ID")print i}}')
    pident=\$(echo \$header|awk '{for(i=1;i<=NF;i++){if(\$i=="pident")print i}}')
    header=\$(echo \$header|awk '{for(i=1;i<=NF;i++){if(i==1){printf \$i"\tAlternate_SNP_ID"}else{printf "\t"\$i}}}')
    coverage=\$(echo \$header|awk '{for(i=1;i<=NF;i++){if(\$i=="Coverage")print i}}')
    find -name "*.tsv"|split -l 100 - subsettsvs
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
    data<-read.delim('${outfile}',header=T,sep='\\\\t',comment.char = '#',stringsAsFactors=F,strip.white=T);
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
