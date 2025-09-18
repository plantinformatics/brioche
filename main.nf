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
        params.genomefasta
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
                    params.pident)
          ADD_MARKERINFO.out.map{it->it[0]}.collect().set{allmappings}
          ADD_MARKERINFO.out.map{it->it[1]}.collect().set{minmappings}
          MERGE_MAPPINGS(
             allmappings,
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
    }
}