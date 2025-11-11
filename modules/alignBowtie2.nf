#!/usr/bin/env nextflow

process alignBowtie2 {
    tag "${sample_id}:${chunk_id}"
    label 'bowtie2'

    cpus params.align_cpus ?: 4
    memory params.align_memory ?: '16 GB'
    time params.align_time ?: '6h'
    queue params.align_queue ?: 'short'

    input:
    tuple val(sample_id), val(meta), val(chunk_id), path(r1), path(r2), path(bt2_index_dir), val(bt2_index_basename)

    output:
    tuple val(sample_id), val(meta), val(chunk_id), path("align/${sample_id}/${sample_id}_${chunk_id}.bam"), path("align/${sample_id}/${sample_id}_${chunk_id}.bow")

    script:
    def align_dir = "align/${sample_id}"
    def base_name = "${sample_id}_${chunk_id}"
    def sam_path = "${align_dir}/${base_name}.sam"
    def bow_path = "${align_dir}/${base_name}.bow"
    def bam_path = "${align_dir}/${base_name}.bam"
    def staged_bt2_prefix = "${bt2_index_dir}/${bt2_index_basename}"

    """
    set -euo pipefail

    mkdir -p ${align_dir}

    bowtie2 \
        -x ${staged_bt2_prefix} \
        -p ${task.cpus} \
        --no-unal \
        -1 ${r1} \
        -2 ${r2} \
        ${params.bowtie2_extra_args ?: ''} \
        -S ${sam_path} \
        > ${bow_path} 2>&1

    grep -v Warning ${bow_path} > ${bow_path}.tmp || true
    mv ${bow_path}.tmp ${bow_path}

    samtools view -bS ${sam_path} > ${bam_path}
    rm -f ${sam_path}
    """
}
