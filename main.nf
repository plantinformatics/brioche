/*
 * Copyright (c) 2023
 */

/*
 * '' - A Nextflow pipeline for mapping markers to any reference genome
 *
 * This pipeline includes steps for mapping markers to any reference genome based on the
 * flanking sequences around that marker
 *
 *
 */

/*
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

log.info """\
Parameters used in brioche for $params.genomename
=======================================================================================================================
$params.loginfo
=======================================================================================================================
"""

/*
 * Import modules
 */
include {
  BUILD_TARGET;
  SPLIT_TARGET_FASTA;
  PREP_REFGENOME;
  BUILD_BLASTDb;
  RUN_BLAST;
  PROCESS_BLAST_RESULTS;
  ADD_MARKERINFO;
  MERGE_MAPPINGS;
  MERGE_MIN_MAPPINGS;
  INTERMEDIATE_FILTERING;
  ADVANCED_FILTERING;
  COLLECT_SUMSTATS;
  MAKE_RUN_META;
  BUILD_SUMMARY_CORE;
  BUILD_SUMMARY_WRAPPER;
  SPLIT_REFGENOME;
  } from './modules.nf'

/*
 * main pipeline logic
 */

workflow {
    // Part 1: Build blast DB
    PREP_REFGENOME( params.genomefasta )
    SPLIT_REFGENOME(
        PREP_REFGENOME.out, 
        params.resultsdirectory,
        params.chromstoexclude,
        params.buildblastdbonly
    )
    BUILD_BLASTDb( 
        SPLIT_REFGENOME.out,
        params.genomefasta,
        params.genomename
    )
    if(!params.buildblastdbonly)
    {
        BUILD_TARGET( params.targetdesign )

        // Channels: keep target table for later; split the fasta now
        targettable_ch = BUILD_TARGET.out.map{ it[0] }   // path to target.table
        fasta_ch       = BUILD_TARGET.out.map{ it[1] }   // path to target.fasta

        // Split target.fasta into 50k-sequence chunks or what user sets
        def query_chunk_size = params.query_chunk_size ?: 50000
        SPLIT_TARGET_FASTA( fasta_ch, query_chunk_size )
        split_queries = SPLIT_TARGET_FASTA.out.flatten()   // paths split/q_*.fa

        // Run BLAST once per split query FASTA channel
        RUN_BLAST(
          split_queries,           // one task per query chunk FASTA
          BUILD_BLASTDb.out,       // blastdb (used as `val blastdb` in the process)
          params.blastoutformat
        )

          //Part 2b: Process the output from blasta
        PROCESS_BLAST_RESULTS(
          RUN_BLAST.out,
          targettable_ch,
          params.minlength,
          params.extendablebps,
          PREP_REFGENOME.out,
          params.blastoutformat,
          params.istarget3primeend,
          params.markercharacter
        )

          // PART 3: Process to filter markers
        ADD_MARKERINFO(
          PROCESS_BLAST_RESULTS.out,
          params.probename,
          params.genomename,
          params.keepduplicates,
          params.coverage,
          params.pident,
          params.istarget3primeend,
          params.doorientation,
          params.forcebiallelic,
          params.orientationfile
        )

        // Split the 5-tuple output from ADD_MARKERINFO into separate channels (need to do this as more separately filtered individual files are present now before the concatenation
        //  0: *_all_mappings.csv
        //  1: *_filtered_mappings.csv
        //  2: *_mapping.tsv
        //  3: *_filtered_mappings.tsv
        //  4: *_complete_unfiltered_blast_results.csv

        ADD_MARKERINFO.out.map{ it[0] }.collect().set{ allmappings_ch }
        ADD_MARKERINFO.out.map{ it[1] }.collect().set{ filtered_allmappings_ch }
        ADD_MARKERINFO.out.map{ it[2] }.collect().set{ mappings_tsv_ch }
        ADD_MARKERINFO.out.map{ it[3] }.collect().set{ minmappings_ch }          // filtered TSVs
        ADD_MARKERINFO.out.map{ it[4] }.collect().set{ rawsallmappings_ch }

        // Part 4. Merge CSV outputs
        MERGE_MAPPINGS(
          allmappings_ch,          // *_all_mappings.csv
          filtered_allmappings_ch, // *_filtered_mappings.csv
          rawsallmappings_ch,      // *_complete_unfiltered_blast_results.csv
          params.resultsdirectory,
          params.probename,
          params.genomename
        )
        // Part 4. Merge TSVs
        MERGE_MIN_MAPPINGS(
          minmappings_ch,          // *_filtered_mappings.tsv
          params.resultsdirectory,
          params.probename,
          params.genomename
        )
          MERGE_MAPPINGS.out.map{it->it[0]}.set{filteredmappings}
          MERGE_MAPPINGS.out.map{it->it[0]}.set{filteredmappings}

          // PART 5. Intermediate filtering run proces 
          INTERMEDIATE_FILTERING(
            filteredmappings,
            MERGE_MIN_MAPPINGS.out,
            params.resultsdirectory,
            params.probename,
            params.genomename,
            params.usetargetchrom,
            params.chromchrommatch,
            params.markertargetsites,
            params.usesharedmarkersmap,
            params.similarmarkersmap,
            params.useldedgemap,
            params.ldedgemap,
            params.usegeneticmap,
            params.geneticmap 
          )
          INTERMEDIATE_FILTERING.out.map{it->it[0]}.set{filteredmappretzelcsv}
          INTERMEDIATE_FILTERING.out.map{it->it[1]}.set{intermediatefilteredmap}
          ADVANCED_FILTERING(
            filteredmappretzelcsv,
            intermediatefilteredmap,
            params.resultsdirectory,
            params.probename,
            params.genomename,
            params.usetargetchrom,
            params.chromchrommatch,
            params.markertargetsites,
            params.usesharedmarkersmap,
            params.similarmarkersmap,
            params.useldedgemap,
            params.ldedgemap,
            params.usegeneticmap,
            params.geneticmap,
            params.doorientation,
            params.orientationfile,
            projectDir
          )
          ADVANCED_FILTERING.out.map{ it -> it[0] }.set { strictmappedcsv }
          ADVANCED_FILTERING.out.map{ it -> it[1] }.set { strictmappedctsv }

          COLLECT_SUMSTATS(
            strictmappedcsv,
            strictmappedctsv,
            params.resultsdirectory,
            params.probename,
            params.genomename,
            params.targetdesign
          )

          COLLECT_SUMSTATS.out.collect().set { SUMSTATS_DONE_CH }

          // Build run metadata (for run_meta.json)
          def run_info = [
            runName         : workflow.runName,
            sessionId       : workflow.sessionId,
            profile         : workflow.profile,
            start           : workflow.start,
            complete        : workflow.complete,
            duration        : workflow.duration,
            success         : workflow.success,
            exitStatus      : workflow.exitStatus,
            containerEngine : workflow.containerEngine,
            nextflowVersion : workflow.nextflow.version,
            projectDir      : workflow.projectDir.toString(),
            launchDir       : workflow.launchDir.toString(),
            workDir         : workflow.workDir.toString(),
            outdir          : (params.resultsdirectory ?: "${workflow.projectDir}/results").toString(),
            commandLine     : workflow.commandLine
          ]

          Channel.value(run_info).set { RUN_INFO_CH }
          MAKE_RUN_META(RUN_INFO_CH)
          MAKE_RUN_META.out.set { META_JSON_CH }

          // Pass the Figs directory path (always exists; created in config)
          Channel.value( file("${params.resultsdirectory}/Figs") ).set { FIGS_DIR_CH }

          BUILD_SUMMARY_CORE(
            META_JSON_CH,
            FIGS_DIR_CH,
            SUMSTATS_DONE_CH
          )
          def SUMMARY_CORE_HTML_CH = BUILD_SUMMARY_CORE.out[0]
          BUILD_SUMMARY_WRAPPER(
            SUMMARY_CORE_HTML_CH,
            params.ts,                  
            params.ts2                 
          )
    }
}
