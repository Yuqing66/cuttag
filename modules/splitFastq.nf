#!/usr/bin/env nextflow

/*
 * split fastq files into smaller chunks for parallel processing
 */

process splitFastq {
    tag "${sample_id}"

    cpus params.splitFastq_cpus ?: 2
    memory params.splitFastq_memory ?: '16 GB'
    time params.splitFastq_time ?: '4h'
    queue params.splitFastq_queue ?: 'short'

    input:
    tuple val(sample_id), val(meta), path(r1), path(i2), path(r3)

    output:
    tuple val(sample_id), val(meta), path("split/${sample_id}/${sample_id}_R1-*.fastq.gz"), path("split/${sample_id}/${sample_id}_R2-*.fastq.gz"), path("split/${sample_id}/${sample_id}_R3-*.fastq.gz")

    script:
    """
    set -euo pipefail

    chunk_dir="split/${sample_id}"
    mkdir -p "\${chunk_dir}"

    reads_per_chunk=${params.reads_per_chunk ?: 1000000}
    lines_per_chunk=\$(( reads_per_chunk * 4 ))

    pigz -dc ${r1} | split -d -l \${lines_per_chunk} -a 4 - \${chunk_dir}/${sample_id}_R1-
    pigz -dc ${i2} | split -d -l \${lines_per_chunk} -a 4 - \${chunk_dir}/${sample_id}_R2-
    pigz -dc ${r3} | split -d -l \${lines_per_chunk} -a 4 - \${chunk_dir}/${sample_id}_R3-

    for fq in \${chunk_dir}/${sample_id}_R1-* \${chunk_dir}/${sample_id}_R2-* \${chunk_dir}/${sample_id}_R3-*; do
        mv "\${fq}" "\${fq}.fastq"
        gzip "\${fq}.fastq"
    done
    """
}
