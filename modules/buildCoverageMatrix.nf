#!/usr/bin/env nextflow

process buildCoverageMatrix {
    tag "${sample_id}:${antibody}"
    label 'bedtools'

    publishDir "results", mode: 'copy', pattern: 'peak_bc_matrix_*'

    cpus params.coverage_matrix_cpus ?: 2
    memory params.coverage_matrix_memory ?: '8 GB'
    time params.coverage_matrix_time ?: '4h'
    queue params.coverage_matrix_queue ?: 'short'

    input:
    tuple val(sample_id), val(meta), val(antibody), path(filtered_peaks_bed), path(cutsite_txt)

    output:
    tuple val(sample_id), val(meta), val(antibody),
          path("peak_bc_matrix_${sample_id}_${antibody}.txt"),
          path("peak_bc_matrix_${sample_id}_${antibody}_colnames.txt"),
          path("peak_bc_matrix_${sample_id}_${antibody}_rownames.txt")

    script:
    def min_cut_sites = params.coverage_min_cut_sites ?: 50
    def split_size = params.coverage_matrix_split_size ?: 1000
    def tmp_dir = "${sample_id}_${antibody}_coverage_tmp"
    def matrix_file = "peak_bc_matrix_${sample_id}_${antibody}.txt"
    def colnames_file = "peak_bc_matrix_${sample_id}_${antibody}_colnames.txt"
    def rownames_file = "peak_bc_matrix_${sample_id}_${antibody}_rownames.txt"
    def barcode_freq = "${sample_id}_${antibody}_barcode_freq.txt"
    def selected_barcodes = "${sample_id}_${antibody}_selected_barcodes.txt"
    def cutsites_bed = "${sample_id}_${antibody}_cutsites_sorted.bed"
    """
    set -euo pipefail

    tmp_dir="${tmp_dir}"
    matrix_file="${matrix_file}"
    colnames_file="${colnames_file}"
    rownames_file="${rownames_file}"
    barcode_freq="${barcode_freq}"
    selected_barcodes="${selected_barcodes}"
    cutsites_bed="${cutsites_bed}"
    min_cut_sites=${min_cut_sites}
    split_size=${split_size}

    awk '{print \$3}' ${cutsite_txt} | sort | uniq -c > \${barcode_freq}
    awk -v min=\${min_cut_sites} '{if(\$1>min) print \$2}' \${barcode_freq} > \${selected_barcodes}
    cp \${selected_barcodes} \${colnames_file}

    if [[ ! -s ${filtered_peaks_bed} ]]; then
        : > \${matrix_file}
        : > \${rownames_file}
        rm -f \${barcode_freq}
        rm -rf \${tmp_dir}
        exit 0
    fi

    awk '{print \$1":"\$2"-"\$3}' ${filtered_peaks_bed} > \${rownames_file}

    if [[ ! -s \${selected_barcodes} ]]; then
        : > \${matrix_file}
        rm -f \${barcode_freq}
        rm -rf \${tmp_dir}
        exit 0
    fi

    mkdir -p \${tmp_dir}
    awk -v OFS="\\t" '{print \$1,\$2,\$2,\$3}' ${cutsite_txt} | bedtools sort -i - > \${cutsites_bed}

    coverage_list="\${tmp_dir}/coverage_files.list"
    : > \${coverage_list}

    while IFS= read -r barcode; do
        [[ -z "\${barcode}" ]] && continue
        coverage_file="\${tmp_dir}/\${barcode}_coverage.txt"
        awk -v bc="\${barcode}" -v OFS="\\t" '\$4==bc {print \$1,\$2,\$3}' \${cutsites_bed} | \
            bedtools sort -i - | \
            bedtools coverage -counts -a ${filtered_peaks_bed} -b stdin | \
            awk '{print \$4}' > "\${coverage_file}"
        echo "\${coverage_file}" >> \${coverage_list}
    done < \${selected_barcodes}

    if [[ ! -s \${coverage_list} ]]; then
        : > \${matrix_file}
        rm -f \${barcode_freq}
        rm -rf \${tmp_dir}
        exit 0
    fi

    split -l \${split_size} -d \${coverage_list} "\${tmp_dir}/lists"
    for list in \${tmp_dir}/lists*; do
        files=\$(tr '\\n' ' ' < \${list})
        paste \${files} > "\${tmp_dir}/merge\${list##*lists}"
    done
    paste \${tmp_dir}/merge* > \${matrix_file}

    rm -f \${barcode_freq} \${cutsites_bed}
    rm -rf \${tmp_dir}
    """
}
