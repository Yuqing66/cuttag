#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * scCUTTag_10x.nf
 *
 * End-to-end Nextflow implementation of the manual CUT&Tag workflow:
 *  1. FASTQ QC (FastQC)
 *  2. Split raw FASTQs into chunks for parallelism
 *  3. Extract valid reads per chunk and read (R1/R2)
 *  4. Trim residual adapters (cutadapt)
 *  5. Align each chunk with Bowtie2 + samtools
 *  6. Merge BAMs per sample, generate fragment-level BAMs, and deduplicate
 */

params.samplesheet              = params.samplesheet              ?: 'input_example.csv'
params.reads_per_chunk          = params.reads_per_chunk          ?: 1000000
params.adapter_r1               = params.adapter_r1               ?: 'CTGTCTCTTATACACATCTG'
params.adapter_r2               = params.adapter_r2               ?: 'CTGTCTCTTATACACATCTG'
params.adapter_min_length       = params.adapter_min_length       ?: 20
params.fragment_min_mapq        = params.fragment_min_mapq        ?: 30
params.antibody_sheet           = params.antibody_sheet           ?: 'input_antibody_barcode.csv'

def parseCsvLine(String line) {
    def values = []
    def current = new StringBuilder()
    boolean inQuotes = false

    for (int idx = 0; idx < line.size(); idx++) {
        char c = line.charAt(idx)
        if (c == '"') {
            if (inQuotes && idx + 1 < line.size() && line.charAt(idx + 1) == '"') {
                current.append('"')
                idx++
            } else {
                inQuotes = !inQuotes
            }
        } else if (c == ',' && !inQuotes) {
            values << current.toString()
            current.setLength(0)
        } else {
            current.append(c)
        }
    }

    values << current.toString()
    return values.collect { it.trim() }
}

def parseBarcodeField(String value) {
    (value ?: '')
        .split(/\s+/)
        .findAll { it }
}

def antibody_sheet_path = file(params.antibody_sheet)
if( !antibody_sheet_path.exists() ) {
    throw new IllegalArgumentException("Antibody sheet '${params.antibody_sheet}' was not found")
}

def antibody_rows = []
def antibody_reader = java.nio.file.Files.newBufferedReader(antibody_sheet_path)
try {
    def headerLine = antibody_reader.readLine()
    assert headerLine : "Antibody sheet '${params.antibody_sheet}' is empty"

    def headers = parseCsvLine(headerLine)
    def idxAntibody = headers.indexOf('antibody')
    def idxB1 = headers.indexOf('b1')
    def idxB2 = headers.indexOf('b2')
    assert idxAntibody >= 0 && idxB1 >= 0 && idxB2 >= 0 : "Antibody sheet must contain 'antibody', 'b1', and 'b2' columns"

    String line
    while( (line = antibody_reader.readLine()) != null ) {
        if( !line?.trim() ) {
            continue
        }

        def values = parseCsvLine(line)
        def antibody = values[idxAntibody]
        def b1_field = values[idxB1]
        def b2_field = values[idxB2]
        def b1_values = parseBarcodeField(b1_field)
        def b2_values = parseBarcodeField(b2_field)

        assert antibody : "Missing antibody value in '${params.antibody_sheet}'"
        assert b1_values : "Missing b1 barcode(s) for antibody ${antibody}"
        assert b2_values : "Missing b2 barcode(s) for antibody ${antibody}"

        antibody_rows << [
            antibody: antibody as String,
            b1      : b1_values,
            b2      : b2_values,
            b1_raw  : b1_field,
            b2_raw  : b2_field
        ]
    }
} finally {
    antibody_reader?.close()
}

assert antibody_rows : "No antibody definitions found in '${params.antibody_sheet}'"

def all_b1_barcodes = antibody_rows.collectMany { it.b1 }
def all_b2_barcodes = antibody_rows.collectMany { it.b2 }

include { fastqc }              from './modules/fastqc.nf'
include { splitFastq }          from './modules/splitFastq.nf'
include { extractValidReads }   from './modules/extractValidReads.nf'
include { trimAdapters }        from './modules/trimAdapters.nf'
include { alignBowtie2 }        from './modules/alignBowtie2.nf'
include { mergeBam }            from './modules/mergeBam.nf'
include { deduplicateFragments } from './modules/deduplicate.nf'

workflow scCUTTag_10x {

    def rows_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map { row ->
            def sample_id = row.sample_id ?: row.sample ?: null
            assert sample_id : "samplesheet must provide 'sample_id' column"

            def fastq_r1 = row.fastq_r1
            def fastq_i  = row.fastq_i
            def fastq_r2 = row.fastq_r2

            assert fastq_r1 && fastq_i && fastq_r2 : "samplesheet must provide fastq_r1, fastq_i, and fastq_r2 (read 2) for sample ${sample_id}"

            tuple(
                sample_id as String,
                [sample: sample_id],
                file(fastq_r1),
                file(fastq_i),
                file(fastq_r2)
            )
        }

    rows_ch.into { split_input; qc_input }

    fastqc(
        qc_input.map { sample_id, meta, r1, i, r2 ->
            tuple(sample_id, r1, r2)
        }
    )

    splitFastq(split_input)

    splitFastq.out
        .flatMap { sample_id, meta, r1_matches, i_matches, r2_matches ->
            def r1_files = (r1_matches instanceof List ? r1_matches : [r1_matches]).sort { it.getFileName().toString() }
            def i_files  = (i_matches instanceof List ? i_matches : [i_matches]).sort { it.getFileName().toString() }
            def r2_files = (r2_matches instanceof List ? r2_matches : [r2_matches]).sort { it.getFileName().toString() }

            assert r1_files.size() == i_files.size() && r1_files.size() == r2_files.size() : "Mismatched chunk counts for sample ${sample_id}"

            (0..<r1_files.size()).collect { idx ->
                def r1_file = r1_files[idx]
                def i_file  = i_files[idx]
                def r2_file = r2_files[idx]
                def suffix = r1_file.getFileName().toString()
                    .replaceAll(/^.*_R1-/, '')
                    .replace('.fastq.gz', '')
                def chunk_id = suffix
                tuple(sample_id, meta, chunk_id, r1_file, i_file, r2_file)
            }
        }
        .set { base_chunk_pairs }

    base_chunk_pairs
        .map { sample_id, meta, chunk_id, r1_chunk, i_chunk, r2_chunk ->
            def meta_with_antibody_info = meta + [
                antibodies      : antibody_rows,
                antibody_b1_list: all_b1_barcodes,
                antibody_b2_list: all_b2_barcodes
            ]

            tuple(sample_id, meta_with_antibody_info, chunk_id, r1_chunk, i_chunk, r2_chunk)
        }
        .set { chunk_pairs }

    extractValidReads(chunk_pairs)

    trimAdapters(extractValidReads.out)

    alignBowtie2(trimAdapters.out)

    alignBowtie2.out
        .map { sample_id, meta, chunk_id, bam, bow -> tuple(sample_id, meta, bam) }
        .groupTuple()
        .map { sample_id, metas, bams ->
            def meta = metas[0]
            tuple(sample_id, meta, bams)
        }
        .set { merge_input }

    mergeBam(merge_input)

    deduplicateFragments(mergeBam.out)
}

workflow {
    scCUTTag_10x()
}
