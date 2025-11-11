#!/usr/bin/env nextflow

process deduplicateFragments {
    tag "${sample_id}"
    label 'dedup'

    cpus params.dedup_cpus ?: 2
    memory params.dedup_memory ?: '16 GB'
    time params.dedup_time ?: '6h'
    queue params.dedup_queue ?: 'short'

    input:
    tuple val(sample_id), val(meta), path(merged_bam), path(sorted_bam), path(sorted_bai), path(dedup_script)

    output:
    tuple val(sample_id), val(meta), path("fragments/${sample_id}_fragment_sorted.bam"), path("fragments/${sample_id}_fragment_sorted.bam.bai"), path("dedup/${sample_id}_fragment_freq.txt"), path("dedup/${sample_id}_dedup.txt"), path("dedup/${sample_id}_fragment_freq_*.txt"), path("dedup/${sample_id}_fragment_freqDedup.txt")

    script:
    """
    set -euo pipefail

    fragment_dir="fragments"
    dedup_dir="dedup"
    mkdir -p \${fragment_dir} \${dedup_dir}

    fragment_raw=\${fragment_dir}/${sample_id}_fragment.bam
    fragment_sorted=\${fragment_dir}/${sample_id}_fragment_sorted.bam
    fragment_index=\${fragment_sorted}.bai
    freq_file=\${dedup_dir}/${sample_id}_fragment_freq.txt
    dedup_report=\${dedup_dir}/${sample_id}_dedup.txt

    samtools view -b -f 67 -q ${params.fragment_min_mapq ?: 30} ${sorted_bam} > \${fragment_raw}
    samtools sort -o \${fragment_sorted} \${fragment_raw}
    samtools index \${fragment_sorted}
    rm -f \${fragment_raw}

    samtools view -f 67 -q ${params.fragment_min_mapq ?: 30} ${sorted_bam} | \
        awk -v OFS=":" '{print \$1,\$3,\$4,\$8}' | \
        awk -F":" -v OFS="\t" '{print \$8,\$9,\$10,\$11,\$12,\$13}' | \
        awk -v OFS="\t" '{print \$4":"\$5":"\$6":"\$1":"\$2":"\$3}' | \
        sort | uniq -c > \${freq_file}

    ab_args=( ${meta.antibody_args.collect { "\"${it}\"" }.join(' ')} )

    python ${dedup_script} \
        -odir \${dedup_dir} \
        -f \${freq_file} \
        -endbp ${params.dedup_endbp ?: 3} \
        -perc_cutoff ${params.dedup_perc_cutoff ?: 0.7} \
        -ab \${ab_args[@]} \
        ${params.dedup_extra_args ?: ''} \
        > \${dedup_report}
    """
}
