#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.samplesheet              = params.samplesheet              ?: 'input_example.csv'
params.reads_per_chunk          = params.reads_per_chunk          ?: 1000000
params.adapter_r1               = params.adapter_r1               ?: 'CTGTCTCTTATACACATCTG'
params.adapter_r2               = params.adapter_r2               ?: 'CTGTCTCTTATACACATCTG'
params.adapter_min_length       = params.adapter_min_length       ?: 20
params.fragment_min_mapq        = params.fragment_min_mapq        ?: 30
params.antibody_sheet           = params.antibody_sheet           ?: 'input_antibody_barcode.csv'
params.coverage_min_cut_sites   = params.coverage_min_cut_sites   ?: 50

def extract_script_path = file(params.extract_script)
def cell_barcode_reference_path = file(params.cell_barcode_reference)
def dedup_script_path = file(params.dedup_script)
def chrom_size_path = file(params.chrom_size)
def blacklist_bed_path = file(params.blacklist_bed)
def control_bam_path = params.control_bam

def antibody_meta = file(params.antibody_sheet)
  .splitCsv()                                        // returns a List of Lists
  .drop(1)                                           // skip header row
  .with { rows -> [
      antibody_b1_list: rows*.getAt(1),
      antibody_b2_list: rows*.getAt(2),
      antibody_args   : rows.collect { r -> "${r[0]}:${r[1]},${r[2]}" }
  ]}

print(antibody_meta)
def bowtie2_index_file = new File(params.bowtie2_index as String)
def bowtie2_index_dir = file(bowtie2_index_file.parent, type: 'dir')
def bowtie2_index_basename = bowtie2_index_file.name

include { fastqc }              from './modules/fastqc.nf'
include { splitFastq }          from './modules/splitFastq.nf'
include { extractValidReads }   from './modules/extractValidReads.nf'
include { trimAdapters }        from './modules/trimAdapters.nf'
include { alignBowtie2 }        from './modules/alignBowtie2.nf'
include { mergeBam }             from './modules/mergeBam.nf'
include { deduplicateFragments } from './modules/deduplicate.nf'
include { callAntibodyPeaks }    from './modules/callPeaks.nf'
include { filterAntibodyPeaks }  from './modules/filterPeaks.nf'
include { filterBlacklistPeaks } from './modules/filterBlacklist.nf'
include { buildCoverageMatrix }  from './modules/buildCoverageMatrix.nf'

workflow scCUTTag_10x {

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map { row ->
            def sample_id = (row.sample_id ?: row.sample) as String
            tuple(
                sample_id,
                [:] + [sample: sample_id] + antibody_meta,
                file(row.fastq_r1),
                file(row.fastq_i),
                file(row.fastq_r2)
            )
        }
        .tap { split_input }
        .set { qc_input }

    fastqc(
        qc_input.map { sample_id, meta, r1, i, r2 ->
            tuple(sample_id, r1, r2)
        }
    )

    splitFastq(split_input)

    splitFastq.out
        .flatMap { sample_id, meta, r1_matches, i_matches, r2_matches ->
            def r1_files = ((r1_matches instanceof List) ? r1_matches : [r1_matches]).sort { it.fileName.toString() }
            def i_files  = ((i_matches instanceof List) ? i_matches : [i_matches]).sort { it.fileName.toString() }
            def r2_files = ((r2_matches instanceof List) ? r2_matches : [r2_matches]).sort { it.fileName.toString() }

            (0..<r1_files.size()).collect { idx ->
                def suffix = r1_files[idx].fileName.toString()
                    .replaceAll(/^.*_R1-/, '')
                    .replace('.fastq.gz', '')
                tuple(sample_id, meta, suffix, r1_files[idx], i_files[idx], r2_files[idx])
            }
        }
        .map { sample_id, meta, chunk_id, r1_chunk, i_chunk, r2_chunk ->
            tuple(sample_id, meta, chunk_id, r1_chunk, i_chunk, r2_chunk, extract_script_path, cell_barcode_reference_path)
        }
        .set { chunk_pairs }

    extractValidReads(chunk_pairs)

    trimAdapters(extractValidReads.out)

    trimAdapters.out
        .map { sample_id, meta, chunk_id, r1_cut, r2_cut ->
            tuple(sample_id, meta, chunk_id, r1_cut, r2_cut, bowtie2_index_dir, bowtie2_index_basename)
        }
        .set { align_input }

    alignBowtie2(align_input)

    alignBowtie2.out
        .map { sample_id, meta, chunk_id, bam, bow -> tuple(sample_id, meta, bam) }
        .groupTuple()
        .map { sample_id, metas, bams ->
            tuple(sample_id, metas[0], bams)
        }
        .set { merge_input }

    mergeBam(merge_input)

    mergeBam.out
        .map { sample_id, meta, merged_bam, sorted_bam, sorted_bai ->
            tuple(sample_id, meta, merged_bam, sorted_bam, sorted_bai, dedup_script_path)
        }
        .set { dedup_input }

    deduplicateFragments(dedup_input)

    deduplicateFragments.out
        .flatMap { sample_id, meta, fragment_sorted_bam, fragment_index, freq_file, dedup_report, antibody_tracks, dedup_track ->
            def tracks = antibody_tracks instanceof List ? antibody_tracks : [antibody_tracks]
            def freq_prefix = "${sample_id}_fragment_freq_"
            def freq_diag = "${sample_id}_fragment_freq_freq.txt"
            tracks
                // Downstream steps only use antibody-specific _freq_{antibody}.txt files (drop the diagnostic _freq_dedup.txt).
                .findAll { track ->
                    def track_name = track?.name
                    track_name &&
                    track_name.startsWith(freq_prefix) &&
                    track_name != freq_diag &&
                    !track_name.endsWith('_dedup.txt')
                }
                .collect { track ->
                    def antibody = track.name.substring(freq_prefix.length())
                    antibody = antibody?.endsWith('.txt') ? antibody[0..-5] : antibody
                    tuple(sample_id, meta, antibody, track)
                }
        }
        .tap { antibody_cutsite_for_peaks}
        .set { antibody_cutsite_for_matrix }

    antibody_cutsite_for_peaks
        .map { sample_id, meta, antibody, cutsite_file ->
            tuple(sample_id, meta, antibody, cutsite_file, chrom_size_path, control_bam_path)
        }
        .set { call_peaks_input }

    antibody_cutsite_for_matrix
        .map { sample_id, meta, antibody, cutsite_file ->
            tuple("${sample_id}:${antibody}", tuple(sample_id, meta, antibody, cutsite_file))
        }
        .set { coverage_cutsite_keyed }

    callAntibodyPeaks(call_peaks_input)

    callAntibodyPeaks.out
        .map { sample_id, meta, antibody, peaks_broad, cutsites_bed ->
            tuple(sample_id, meta, antibody, peaks_broad, cutsites_bed, params.peak_cpbpmc ?: 2, params.peak_cutnum ?: 2, params.peak_merge_distance ?: 3000)
        }
        .set { peak_filter_input }

    filterAntibodyPeaks(peak_filter_input)

    filterAntibodyPeaks.out
        .map { sample_id, meta, antibody, merge_distance, merged_bed ->
            tuple(sample_id, meta, antibody, merge_distance, merged_bed, blacklist_bed_path)
        }
        .set { blacklist_filter_input }

    filterBlacklistPeaks(blacklist_filter_input)

    filterBlacklistPeaks.out
        .map { sample_id, meta, antibody, merge_distance, filtered_bed ->
            tuple("${sample_id}:${antibody}", tuple(sample_id, meta, antibody, filtered_bed))
        }
        .join(coverage_cutsite_keyed)
        .map { key, peak_data, cutsite_data ->
            def (sample_id, meta, antibody, filtered_bed) = peak_data
            def cutsite_file = cutsite_data[3]
            tuple(sample_id, meta, antibody, filtered_bed, cutsite_file)
        }
        .set { coverage_matrix_input }

    buildCoverageMatrix(coverage_matrix_input)
}

workflow {
    scCUTTag_10x()
}
