#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.samplesheet              = params.samplesheet              ?: 'input_example.csv'
params.reads_per_chunk          = params.reads_per_chunk          ?: 1000000
params.adapter_r1               = params.adapter_r1               ?: 'CTGTCTCTTATACACATCTG'
params.adapter_r2               = params.adapter_r2               ?: 'CTGTCTCTTATACACATCTG'
params.adapter_min_length       = params.adapter_min_length       ?: 20
params.fragment_min_mapq        = params.fragment_min_mapq        ?: 30
params.antibody_sheet           = params.antibody_sheet           ?: 'input_antibody_barcode.csv'

def extract_script_path = file(params.extract_script)
def cell_barcode_reference_path = file(params.cell_barcode_reference)
def dedup_script_path = file(params.dedup_script)

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
include { mergeBam }            from './modules/mergeBam.nf'
include { deduplicateFragments } from './modules/deduplicate.nf'

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
}

workflow {
    scCUTTag_10x()
}
