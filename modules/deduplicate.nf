#!/usr/bin/env nextflow

process deduplicateFragments {
    tag "${sample_id}"
    label 'dedup'

    cpus params.dedup_cpus ?: 2
    memory params.dedup_memory ?: '16 GB'
    time params.dedup_time ?: '6h'
    queue params.dedup_queue ?: 'short'
    container params.dedup_container ?: ''

    input:
    tuple val(sample_id), val(meta), path(merged_bam), path(sorted_bam), path(sorted_bai)

    output:
    tuple val(sample_id), val(meta), path("fragments/${sample_id}/${sample_id}_fragment_sorted.bam"), path("fragments/${sample_id}/${sample_id}_fragment_sorted.bam.bai"), path("dedup/${sample_id}/${sample_id}_fragment_freq.txt"), path("dedup/${sample_id}/${sample_id}_dedup.txt")

    script:
    """
    set -euo pipefail

    : \${params.dedup_script:?"Set --dedup_script to the path of cuttag_dedup_multiAb.py"}

    fragment_dir="fragments/${sample_id}"
    dedup_dir="dedup/${sample_id}"
    mkdir -p ${fragment_dir} ${dedup_dir}

    fragment_raw=${fragment_dir}/${sample_id}_fragment.bam
    fragment_sorted=${fragment_dir}/${sample_id}_fragment_sorted.bam
    fragment_index=${fragment_sorted}.bai
    freq_file=${dedup_dir}/${sample_id}_fragment_freq.txt
    dedup_report=${dedup_dir}/${sample_id}_dedup.txt

    samtools view -b -f 67 -q ${params.fragment_min_mapq ?: 30} ${sorted_bam} > ${fragment_raw}
    samtools sort -o ${fragment_sorted} ${fragment_raw}
    samtools index ${fragment_sorted}
    rm -f ${fragment_raw}

    samtools view -f 67 -q ${params.fragment_min_mapq ?: 30} ${sorted_bam} | \
        awk -v OFS=":" '{print $1,$3,$4,$8}' | \
        awk -F":" -v OFS="\t" '{print $8,$9,$10,$11,$12,$13}' | \
        awk -v OFS="\t" '{print $4":"$5":"$6":"$1":"$2":"$3}' | \
        sort | uniq -c > ${freq_file}

    if [[ -z "${meta.antibodies}" ]]; then
        echo "No antibody metadata found for sample ${sample_id}" >&2
        exit 1
    fi

    ab_args=( ${meta.antibodies.collect { anti -> "\"${anti.antibody}:${anti.b1_raw},${anti.b2_raw}\"" }.join(' ')} )
    if [[ ${#ab_args[@]} -eq 0 ]]; then
        echo "Antibody argument list is empty for sample ${sample_id}" >&2
        exit 1
    fi

    python ${params.dedup_script} \
        -odir ${dedup_dir} \
        -f ${freq_file} \
        -endbp ${params.dedup_endbp ?: 3} \
        -perc_cutoff ${params.dedup_perc_cutoff ?: 0.7} \
        -ab ${ab_args[@]} \
        ${params.dedup_extra_args ?: ''} \
        > ${dedup_report}
    """
}
