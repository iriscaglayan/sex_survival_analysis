/*
 * Four summary figures across all cancers.
 */
process FIGURES {
    publishDir "${params.outdir}/figures", mode: 'copy'

    input:
    path representation
    path dropped
    path cox
    path sex_effect
    path clean_tables

    output:
    path "*.png"

    script:
    """
    make_figures.py --representation ${representation} --dropped ${dropped} \\
        --cox ${cox} --sex-effect ${sex_effect} --clean ${clean_tables} \\
        --min-events ${params.min_events} --outdir .
    """
}
