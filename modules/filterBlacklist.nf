#!/usr/bin/env nextflow

process filterBlacklistPeaks {
    tag "${sample_id}:${antibody}"
    label 'bedtools'

    publishDir "results/peaks", mode: 'copy'

    cpus params.blacklist_peaks_cpus ?: 1
    memory params.blacklist_peaks_memory ?: '4 GB'
    time params.blacklist_peaks_time ?: '1h'
    queue params.blacklist_peaks_queue ?: 'short'

    input:
    tuple val(sample_id), val(meta), val(antibody), path(merged_bed), path(blacklist_bed)

    output:
    tuple val(sample_id), val(meta), val(antibody), path("peaks/${sample_id}_${antibody}_filtered.bed")

    script:
    def peak_dir = "peaks"
    def filtered_bed = "${peak_dir}/${sample_id}_${antibody}_filtered.bed"
    """
    set -euo pipefail

    mkdir -p ${peak_dir}

    bedtools intersect -v -a ${merged_bed} -b ${blacklist_bed} | \
        bedtools sort -i - > ${filtered_bed}
    """
}
