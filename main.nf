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
        probefasta=BUILD_TARGET.out.map{it->it[1]}
        targettable=BUILD_TARGET.out.map{it->it[0]}
          //Part 2: Run blast on probe fasta
          RUN_BLAST(
            probefasta,
            BUILD_BLASTDb.out,
            params.blastoutformat
          )

          //Part 2: Process the output from blasta
          PROCESS_BLAST_RESULTS(
                    RUN_BLAST.out.flatten(),
                    targettable,
                    params.minlength,
                    params.extendablebps,
                    PREP_REFGENOME.out,
                    params.blastoutformat,
                    params.istarget3primeend,
                    params.markercharacter)

          // PART 3: Process to filter markers
          ADD_MARKERINFO(
                    PROCESS_BLAST_RESULTS.out,
                    params.probename,
                    params.genomename,
                    params.keepduplicates,
                    params.coverage,
                    params.pident,
                    params.istarget3primeend)
          ADD_MARKERINFO.out.map{it->it[0]}.collect().set{allmappings}
          ADD_MARKERINFO.out.map{it->it[1]}.collect().set{minmappings}
          ADD_MARKERINFO.out.map{it->it[2]}.collect().set{rawsallmappings}
          MERGE_MAPPINGS(
             allmappings,
             rawsallmappings,
             params.resultsdirectory,
             params.probename,
             params.genomename
          )
          MERGE_MIN_MAPPINGS(
             minmappings,
             params.resultsdirectory,
             params.probename,
             params.genomename
          )
          // PART 3: Process to filter markers
          MERGE_MAPPINGS.out.map{it->it[0]}.set{filteredmappings}

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
            params.ldedgemap
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
            params.ldedgemap
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
