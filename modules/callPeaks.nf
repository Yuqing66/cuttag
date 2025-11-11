#!/usr/bin/env nextflow

process callAntibodyPeaks {
    tag "${sample_id}:${antibody}"
    label 'macs2'

    cpus params.call_peaks_cpus ?: 4
    memory params.call_peaks_memory ?: '32 GB'
    time params.call_peaks_time ?: '6h'
    queue params.call_peaks_queue ?: 'short'

    input:
    tuple val(sample_id), val(meta), val(antibody), path(cutsite_txt), path(chrom_sizes), val(control_bam)

    output:
    tuple val(sample_id), val(meta), val(antibody), path("peaks/${sample_id}_${antibody}_peaks.broadPeak"), path("peaks/${sample_id}_${antibody}_cutsites.bed")

    script:
    def peak_dir = "peaks"
    def workspace_dir = "${peak_dir}/workspace_${sample_id}_${antibody}"
    def cutsites_bed = "${peak_dir}/${sample_id}_${antibody}_cutsites.bed"
    def windows_bed = "${workspace_dir}/${sample_id}_${antibody}_windows.bed"
    def macs_prefix = "${sample_id}_${antibody}"
    """
    set -euo pipefail

    mkdir -p ${peak_dir} ${workspace_dir}

    cat ${cutsite_txt} | \
        awk -v OFS="\\t" '{print \$1,\$2,\$2,\$3}' | \
        bedtools sort -i - > ${cutsites_bed}

    bedtools slop -b ${params.call_peaks_slop_bp ?: 100} -i ${cutsites_bed} -g ${chrom_sizes} | \
        awk -v OFS="\\t" '{print \$1,\$2,\$3,\$4}' > ${windows_bed}

    macs2 callpeak \
        -t ${windows_bed} \
        -c ${control_bam} \
        -n ${macs_prefix} \
        --outdir ${workspace_dir} \
        -g mm \
        -q ${params.call_peaks_qvalue ?: 0.0001} \
        --nomodel \
        --broad \
        -s ${params.call_peaks_bandwidth ?: 200}

    mv ${workspace_dir}/${macs_prefix}_peaks.broadPeak ${peak_dir}/${sample_id}_${antibody}_peaks.broadPeak
    rm -rf ${workspace_dir}
    """
}
