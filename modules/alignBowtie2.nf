#!/usr/bin/env nextflow

process alignBowtie2 {
    tag "${sample_id}:${chunk_id}"
    label 'bowtie2'

    cpus params.align_cpus ?: 4
    memory params.align_memory ?: '16 GB'
    time params.align_time ?: '6h'
    queue params.align_queue ?: 'short'
    container params.align_container ?: ''

    def bt2_index = params.bowtie2_index
    if( !bt2_index ) {
        throw new IllegalArgumentException('params.bowtie2_index must point to a Bowtie2 index prefix')
    }

    input:
    tuple val(sample_id), val(meta), val(chunk_id), path(r1), path(r2)

    output:
    tuple val(sample_id), val(meta), val(chunk_id), path("align/${sample_id}/${sample_id}_${chunk_id}.bam"), path("align/${sample_id}/${sample_id}_${chunk_id}.bow")

    script:
    """
    set -euo pipefail

    outdir="align/${sample_id}"
    mkdir -p ${outdir}

    sam=${outdir}/${sample_id}_${chunk_id}.sam
    bow=${outdir}/${sample_id}_${chunk_id}.bow
    bam=${outdir}/${sample_id}_${chunk_id}.bam

    bowtie2 \
        -x ${bt2_index} \
        -p ${task.cpus} \
        --no-unal \
        -1 ${r1} \
        -2 ${r2} \
        ${params.bowtie2_extra_args ?: ''} \
        -S ${sam} \
        > ${bow} 2>&1

    grep -v Warning ${bow} > ${bow}.tmp || true
    mv ${bow}.tmp ${bow}

    samtools view -bS ${sam} > ${bam}
    rm -f ${sam}
    """
}
