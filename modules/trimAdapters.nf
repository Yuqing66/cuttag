#!/usr/bin/env nextflow

process trimAdapters {
    tag "${sample_id}:${chunk_id}"
    label 'cutadapt'

    cpus params.trimAdapters_cpus ?: 2
    memory params.trimAdapters_memory ?: '8 GB'
    time params.trimAdapters_time ?: '2h'
    queue params.trimAdapters_queue ?: 'short'

    def adapter_r1 = params.adapter_r1 ?: 'CTGTCTCTTATACACATCTG'
    def adapter_r2 = params.adapter_r2 ?: 'CTGTCTCTTATACACATCTG'
    def min_length = params.adapter_min_length ?: 20

    input:
    tuple val(sample_id), val(meta), val(chunk_id), path(r1), path(r2)

    output:
    tuple val(sample_id), val(meta), val(chunk_id), path("cut/${sample_id}_${chunk_id}_R1_cut.fastq.gz"), path("cut/${sample_id}_${chunk_id}_R2_cut.fastq.gz")

    script:
    """
    set -euo pipefail

    mkdir -p cut

    cutadapt \
        -a ${adapter_r1} \
        -A ${adapter_r2} \
        --minimum-length=${min_length} \
        -o cut/${sample_id}_${chunk_id}_R1_cut.fastq.gz \
        -p cut/${sample_id}_${chunk_id}_R2_cut.fastq.gz \
        ${r1} ${r2}
    """
}
