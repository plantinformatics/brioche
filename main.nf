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


// Adding git monitoring for version recording (probably should move this to the nextflow.config or run.config file sometime soon
// Git commit version, git url, added branch name, plus cast paramaters information for future incorporation into final 1to1 mappings file used for anchoring
// Utilises First, Nextflow parameters and secondary failsafe using .git config file
// Note: I've tested this part and it seems robust but it's mostly ChatGPT, I'm not familiar enough with Nextflow and nextflow standard parameters to break down every component of it 
def yn = { v ->
    if (v == null) return 'NO'
    def s = v.toString().trim().toLowerCase()
    (s in ['1','true','t','yes','y']) ? 'YES'
  : (s in ['0','false','f','no','n']) ? 'NO'
  : s.toUpperCase()
}

def shortShaFromGit = { File projDir ->
    try {
        def head = new File(projDir, '.git/HEAD')
        if (!head.exists()) return null
        def txt = head.text.trim()
        if (txt.startsWith('ref:')) {
            def refPath = txt.split(':',2)[1].trim()
            def refFile = new File(projDir, ".git/${refPath}")
            return refFile.exists() ? refFile.text.trim().take(7) : null
        } else {
            return txt.take(7)
        }
    } catch (ignored) { null }
}

def headBranchFromGit = { File projDir ->
    try {
        def head = new File(projDir, '.git/HEAD')
        if (!head.exists()) return null
        def txt = head.text.trim()
        if (!txt.startsWith('ref:')) return null
        def refPath = txt.split(':',2)[1].trim()
        return refPath.tokenize('/').last()
    } catch (ignored) { null }
}

def readGitConfig = { File projDir ->
    def url = null
    def branches = []
    try {
        def cfg = new File(projDir, '.git/config')
        if (!cfg.exists()) return [url: null, branches: []]
        def curSection = ''
        cfg.eachLine { line ->
            def ln = line.trim()
            if (!ln) return
            def mSec = (ln =~ /^\[(.+?)\]\s*$/)
            if (mSec.matches()) {
                curSection = mSec[0][1]
                return
            }
            if (curSection.toLowerCase().startsWith('remote "origin"')) {
                def mUrl = (ln =~ /^url\s*=\s*(.+)$/)
                if (mUrl.matches()) url = mUrl[0][1].trim()
            }
            if (curSection.toLowerCase().startsWith('branch "')) {
                def mBr = (curSection =~ /^branch\s+"(.+)"$/)
                if (mBr.matches()) branches << mBr[0][1]
            }
        }
    } catch (ignored) { /* noop */ }
    [url: url, branches: branches]
}

// SAFE normalizer (no regex backrefs)
def normalizeRepoUrl = { String s ->
    if (!s) return ''
    String t = s.toString()
    if (t.startsWith('git@') && t.contains(':')) {
        int colon = t.indexOf(':')
        String host = t.substring('git@'.length(), colon)
        String path = t.substring(colon + 1)
        t = "https://${host}/${path}"
    } else if (t.startsWith('ssh://git@')) {
        String rest = t.substring('ssh://git@'.length())
        int slash = rest.indexOf('/')
        if (slash > 0) {
            String host = rest.substring(0, slash)
            String path = rest.substring(slash + 1)
            t = "https://${host}/${path}"
        }
    }
    if (t.endsWith('.git')) t = t.substring(0, t.length() - 4)
    return t
}

final File   _projDir   = new File( (workflow.projectDir ?: '.').toString() )
final String _repoRaw   = (workflow.repository ?: '').toString()
final String _commitSha = (workflow.commitId ?: shortShaFromGit(_projDir) ?: '').toString()
final String _revision  = (workflow.revision ?: '').toString()

final def    _cfg       = readGitConfig(_projDir)
final String _cfgUrl    = (_cfg.url ?: '')
final String _cfgBranch = headBranchFromGit(_projDir) ?: (_cfg.branches ? _cfg.branches[0] : '')

final String brioche_repo_url = normalizeRepoUrl( _repoRaw ?: _cfgUrl )
final String brioche_branch   = (_revision ?: _cfgBranch ?: '')

final String brioche_version = _commitSha
    ? (brioche_branch ? "${brioche_branch}@${_commitSha.take(7)}" : "local@${_commitSha.take(7)}")
    : 'unknown@unknown'

log.info "[Brioche] version: ${brioche_version}"
log.info "[Brioche] repo:    ${brioche_repo_url ?: '(none)'}"
log.info "[Brioche] branch:  ${brioche_branch ?: '(none)'}"
log.info "[Brioche] flags (raw):  targetchrom=${params.usetargetchrom}, sharedmap=${params.usesharedmarkersmap}, ldedgemap=${params.useldedgemap}, geneticmap=${params.usegeneticmap}"
log.info "[Brioche] flags (cast): targetchrom=${yn(params.usetargetchrom)}, sharedmap=${yn(params.usesharedmarkersmap)}, ldedgemap=${yn(params.useldedgemap)}, geneticmap=${yn(params.usegeneticmap)}"

// end of adding git monitoring 


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
            params.geneticmap,
            params.localdupdist,
            params.keeplocalduppos
          )
          INTERMEDIATE_FILTERING.out.map{it->it[0]}.set{filteredmappretzelcsv}
          INTERMEDIATE_FILTERING.out.map{it->it[1]}.set{intermediatefilteredmap}
          INTERMEDIATE_FILTERING.out.map{it->it[2]}.set{dupmapinter}
          ADVANCED_FILTERING(
            filteredmappretzelcsv,
            intermediatefilteredmap,
            dupmapinter,
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
          )
          ADVANCED_FILTERING.out.map{ it -> it[0] }.set { strictmappedcsv }
          ADVANCED_FILTERING.out.map{ it -> it[1] }.set { strictmappedtsv }
          ADVANCED_FILTERING.out.map{ it -> it[3] }.set { priors_informatives_strict }

          COLLECT_SUMSTATS(
            strictmappedcsv,
            strictmappedtsv,
            dupmapinter,
            priors_informatives_strict,
            params.resultsdirectory,
            params.probename,
            params.genomename,
            params.targetdesign,
            brioche_version,
            params.coverage,
            params.pident,
            (params.otherblastoptions ?: ''),
            yn(params.usetargetchrom),
            yn(params.usesharedmarkersmap),
            yn(params.useldedgemap),
            yn(params.usegeneticmap),
            (brioche_repo_url   ?: ''),   
            (brioche_branch ?: '')    
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
