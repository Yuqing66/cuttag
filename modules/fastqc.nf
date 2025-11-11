#!/usr/bin/env nextflow

process fastqc {
    tag "${sample_id}"

    publishDir "results/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc ${r1} ${r2} --outdir ./
    """
}
