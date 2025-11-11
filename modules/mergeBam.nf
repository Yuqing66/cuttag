#!/usr/bin/env nextflow

process mergeBam {
    tag "${sample_id}"
    label 'samtools'

    cpus params.merge_cpus ?: 4
    memory params.merge_memory ?: '32 GB'
    time params.merge_time ?: '4h'
    queue params.merge_queue ?: 'short'

    input:
    tuple val(sample_id), val(meta), path(bams)

    output:
    tuple val(sample_id), val(meta), path("merge/${sample_id}/${sample_id}.bam"), path("merge/${sample_id}/${sample_id}_sorted.bam"), path("merge/${sample_id}/${sample_id}_sorted.bam.bai")

    script:
    """
    set -euo pipefail

    outdir="merge/${sample_id}"
    mkdir -p \${outdir}

    merged_bam=\${outdir}/${sample_id}.bam
    sorted_bam=\${outdir}/${sample_id}_sorted.bam

    samtools merge -f \${merged_bam} ${bams.collect { "'${it}'" }.join(' ')}
    samtools sort -o \${sorted_bam} \${merged_bam}
    samtools index \${sorted_bam}
    """
}
