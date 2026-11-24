/*
 * Count women vs men, deaths per sex and missing data per sex.
 */
process REPRESENTATION_AUDIT {
    tag "$cancer"
    publishDir "${params.outdir}/audit", mode: 'copy'

    input:
    tuple val(cancer), path(clinical), path(survival)

    output:
    tuple val(cancer), path("${cancer}_patients.csv"), emit: patients
    path "${cancer}_representation.csv",               emit: summary
    path "${cancer}_representation.png",               emit: plot

    script:
    """
    representation_audit.py --cancer ${cancer} \\
        --clinical ${clinical} --survival ${survival} --outdir .
    """
}
