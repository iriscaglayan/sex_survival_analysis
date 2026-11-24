/*
 * Fan-in: combine every cancer's tables into one table per step.
 */
process MERGE_RESULTS {
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path summaries
    path dropped
    path cox

    output:
    path "all_representation.csv", emit: representation
    path "all_dropped.csv",        emit: dropped
    path "all_cox.csv",            emit: cox
    path "sex_effect.csv",         emit: sex_effect

    script:
    """
    merge_results.py --representation ${summaries} \\
        --dropped ${dropped} --cox ${cox} --outdir .
    """
}
