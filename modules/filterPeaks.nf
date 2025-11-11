#!/usr/bin/env nextflow

process filterAntibodyPeaks {
    tag "${sample_id}:${antibody}"
    label 'bedtools'

    cpus params.filter_peaks_cpus ?: 2
    memory params.filter_peaks_memory ?: '8 GB'
    time params.filter_peaks_time ?: '2h'
    queue params.filter_peaks_queue ?: 'short'

    input:
    tuple val(sample_id), val(meta), val(antibody), path(broad_peak), path(cutsites_bed), val(cpbpmc), val(cutnum), val(merge_distance)

    output:
    tuple val(sample_id), val(meta), val(antibody), path("peaks/${sample_id}_${antibody}_m${merge_distance}.bed")

    script:
    def peak_dir = "peaks"
    def coverage_tmp = "${peak_dir}/${sample_id}_${antibody}_coverage_tmp.bed"
    def filtered_tmp = "${peak_dir}/${sample_id}_${antibody}_coverage_tmp_filtered.bed"
    def merged_bed = "${peak_dir}/${sample_id}_${antibody}_m${merge_distance}.bed"
    """
    set -euo pipefail

    mkdir -p ${peak_dir}

    totalcut=\$(wc -l < ${cutsites_bed} || echo 0)
    if [[ \${totalcut} -eq 0 ]]; then
        : > ${merged_bed}
        exit 0
    fi

    bedtools sort -i ${cutsites_bed} | \
        bedtools coverage -a ${broad_peak} -b stdin | \
        awk -v OFS="\\t" '{print \$1,\$2,\$3,\$4,\$6}' > ${coverage_tmp}

    awk -v cs=\${totalcut} -v OFS="\\t" -v cp=${cpbpmc} -v cut=${cutnum} \
        '{if(\$5>0 && (\$4/\$5*1000000000/cs)>cp && \$4>cut) print \$1,\$2,\$3}' \
        ${coverage_tmp} > ${filtered_tmp}

    if [[ -s ${filtered_tmp} ]]; then
        bedtools sort -i ${filtered_tmp} | \
            bedtools merge -d ${merge_distance} -i - > ${merged_bed}
    else
        : > ${merged_bed}
    fi

    rm -f ${coverage_tmp} ${filtered_tmp}
    """
}
