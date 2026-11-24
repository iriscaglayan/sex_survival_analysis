#!/usr/bin/env nextflow
/*
 * Sex differences in cancer survival across TCGA cohorts
 * -------------------------------------------------------
 * For each cancer (in parallel): DOWNLOAD -> REPRESENTATION_AUDIT -> CLEAN_CLINICAL -> COX_MODEL
 * Then all cancers are combined by MERGE_RESULTS and drawn by FIGURES.
 */

include { DOWNLOAD             } from './modules/download'
include { REPRESENTATION_AUDIT } from './modules/audit'
include { CLEAN_CLINICAL       } from './modules/clean'
include { COX_MODEL            } from './modules/cox'
include { MERGE_RESULTS        } from './modules/merge'
include { FIGURES              } from './modules/figures'

log.info """
    S E X - S U R V I V A L   P I P E L I N E
    =========================================
    cancers    : ${params.cancers}
    min_events : ${params.min_events}
    raw_dir    : ${params.raw_dir}
    outdir     : ${params.outdir}
    """.stripIndent()

workflow {
    // 'STAD,ESCA,COAD' -> channel emitting 'STAD', 'ESCA', 'COAD' (one item each)
    cancers_ch = channel.fromList( params.cancers.tokenize(',')*.trim() )

    DOWNLOAD(cancers_ch)                                   // one task per cancer
    REPRESENTATION_AUDIT(DOWNLOAD.out.tables)
    CLEAN_CLINICAL(REPRESENTATION_AUDIT.out.patients)
    COX_MODEL(CLEAN_CLINICAL.out.clean)

    // fan-in: gather every cancer's results into one task
    MERGE_RESULTS(
        REPRESENTATION_AUDIT.out.summary.collect(),
        CLEAN_CLINICAL.out.dropped.collect(),
        COX_MODEL.out.cox.collect()
    )

    FIGURES(
        MERGE_RESULTS.out.representation,
        MERGE_RESULTS.out.dropped,
        MERGE_RESULTS.out.cox,
        MERGE_RESULTS.out.sex_effect,
        CLEAN_CLINICAL.out.clean.map { cancer, table -> table }.collect()
    )
}

workflow.onComplete {
    log.info workflow.success
        ? "\nDone. Results in: ${params.outdir}\n"
        : "\nFailed. Check .nextflow.log and the work/ folder of the failed task.\n"
}
