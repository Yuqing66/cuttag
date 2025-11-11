#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process extractValidReads {
    tag "${sample_id}:${chunk_id}"

    cpus params.extractValidReads_cpus ?: 1
    memory params.extractValidReads_memory ?: '8 GB'
    time params.extractValidReads_time ?: '4h'
    queue params.extractValidReads_queue ?: 'short'

    input:
    tuple val(sample_id), val(meta), val(chunk_id), path(r1), path(i), path(r2), path(extract_script), path(cell_barcode_reference)

    output:
    tuple val(sample_id), val(meta), val(chunk_id), path("valid/${sample_id}/${sample_id}_${chunk_id}_R1_valid.fastq.gz"), path("valid/${sample_id}/${sample_id}_${chunk_id}_R2_valid.fastq.gz")

    script:
    """
    set -euo pipefail

    outdir="valid/${sample_id}"
    mkdir -p \${outdir}

    input_basename_r1=\$(basename ${r1})
    script_output_r1="\${outdir}/\${input_basename_r1%.fastq.gz}_valid.fastq.gz"
    input_basename_r2=\$(basename ${r2})
    script_output_r2="\${outdir}/\${input_basename_r2%.fastq.gz}_valid.fastq.gz"
    final_r1="\${outdir}/${sample_id}_${chunk_id}_R1_valid.fastq.gz"
    final_r2="\${outdir}/${sample_id}_${chunk_id}_R2_valid.fastq.gz"

    b1_list=( ${meta.antibody_b1_list.collect { '"' + it + '"' }.join(' ')} )
    b2_list=( ${meta.antibody_b2_list.collect { '"' + it + '"' }.join(' ')} )

    python ${extract_script} \
        -r1 ${r1} \
        -index ${i} \
        -r2 ${r2} \
        -odir \${outdir} \
        -b1 \${b1_list[@]} \
        -b2 \${b2_list[@]} \
        -b ${cell_barcode_reference} \
        ${params.extract_extra_args ?: ''}

    mv \${script_output_r1} \${final_r1}
    mv \${script_output_r2} \${final_r2}
    """
}
