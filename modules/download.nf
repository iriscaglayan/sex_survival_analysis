/*
 * Download clinical + survival tables for one TCGA cohort from the UCSC Xena GDC hub.
 * The files are made read-only right after download (raw data is never edited).
 * storeDir keeps the files in data/raw/ permanently: if they are already there,
 * the task is skipped, so the data is downloaded only once, ever.
 * A failed download is retried twice, then stops the run with a clear error.
 */
process DOWNLOAD {
    tag "$cancer"
    storeDir params.raw_dir
    errorStrategy 'retry'   // network hiccup? try again...
    maxRetries 2            // ...twice, then stop with an error

    input:
    val cancer

    output:
    tuple val(cancer),
          path("TCGA-${cancer}.clinical.tsv.gz"),
          path("TCGA-${cancer}.survival.tsv.gz"), emit: tables
    path "TCGA-${cancer}.md5",                    emit: md5

    script:
    """
    for f in clinical survival; do
        curl -fsSL --retry 3 -o TCGA-${cancer}.\$f.tsv.gz \\
            "${params.xena_hub}/TCGA-${cancer}.\$f.tsv.gz"
    done
    md5sum TCGA-${cancer}.*.tsv.gz > TCGA-${cancer}.md5 2>/dev/null \\
        || md5 -r TCGA-${cancer}.*.tsv.gz > TCGA-${cancer}.md5

    # raw data is read-only: nothing (script or person) can edit it by accident
    chmod a-w TCGA-${cancer}.*.tsv.gz TCGA-${cancer}.md5
    """
}
