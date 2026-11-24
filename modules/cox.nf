/*
 * Cox models: all patients (sex effect), women only, men only.
 */
process COX_MODEL {
    tag "$cancer"
    publishDir "${params.outdir}/cox", mode: 'copy'

    input:
    tuple val(cancer), path(clean)

    output:
    path "${cancer}_cox.csv", emit: cox

    script:
    """
    cox_model.py --cancer ${cancer} --clean ${clean} \\
        --min-events ${params.min_events} --outdir .
    """
}
