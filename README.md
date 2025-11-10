# scCUTTag 10x Nextflow Pipeline

This repository captures the manual multi-step single-cell CUT&Tag workflow in a reproducible Nextflow DSL2 pipeline. Each manual command from the notebook (split, extract valid reads, cutadapt, Bowtie2, samtools merge, fragment filtering, deduplication) is mapped to an isolated process so the entire workflow can be launched with a single `nextflow run` command.

## Input samplesheet

Provide a CSV file (default `input_example.csv`) where **each row represents the fully merged FASTQ trio for one sample**. Required columns:

| column      | description |
|-------------|-------------|
| `sample_id` | Logical sample name |
| `fastq_r1`  | Path to the merged raw R1 FASTQ (read 1) |
| `fastq_i`   | Path to the merged raw index FASTQ (10x cell barcode read, typically `R2`) |
| `fastq_r2`  | Path to the merged raw R2 FASTQ (read 2) |

All sequencing lanes for the same sample must be concatenated (or otherwise merged) before launching this pipeline. The samplesheet therefore contains exactly one row per sample. Antibody barcode definitions are tracked separately (see below).

## Antibody barcode sheet

Define all antibodies that were multiplexed together in `input_antibody_barcode.csv` (configurable via `--antibody_sheet`). The file must have a header `antibody,b1,b2` where:

- `antibody` – display name (e.g. `H3K27ac`).
- `b1` – **space-separated** list of the antibody barcodes to pass to `extractValidReads` via `-b1`.
- `b2` – space-separated list for the `-b2` argument.

Each row describes one antibody, and the workflow runs every sample once per row (so two antibodies double the per-sample analyses). Outputs are namespaced as `<sample>/<antibody>/...` to keep results distinct. If a barcode list contains whitespace wrap the field in quotes (see `input_antibody_barcode.csv`).

## Required parameters

Set the following parameters via `nextflow.config`, `-params-file`, or `-params` on the command line:

- `samplesheet`: path to the FASTQ samplesheet (defaults to `input_example.csv`).
- `extract_script`: path to `extractValidReads_multiCUTTag_nobccorrect_hamdist8.py`.
- `cell_barcode_reference`: path to the 10x whitelist (e.g. `737K-cratac-v1.txt`).
- `antibody_sheet`: path to the antibody barcode sheet (`input_antibody_barcode.csv`).
- `bowtie2_index`: Bowtie2 genome index prefix (e.g. `/.../mm10/Bowtie2Index/genome`).
- `dedup_script`: path to `cuttag_dedup_multiAb.py`.

Optional knobs: chunk size (`reads_per_chunk`), adapter sequences/length thresholds, Bowtie2 extra flags, dedup thresholds, queue/memory/cpu hints per process (see `modules/*.nf`).

## Workflow stages

1. **FastQC** – sanity check the biological read pair (R1/R2) for every merged sample.
2. **Split FASTQ** – chunk each sample’s R1/index/R2 reads into ~`reads_per_chunk` blocks to maximize parallelism.
3. **Extract Valid Reads** – run the provided Python script on each chunk (consuming R1/index/R2 and emitting valid R1/R2 pairs).
4. **Cutadapt** – trim residual Nextera adapters on chunk-level valid pairs.
5. **Bowtie2 + samtools** – align every chunk pair, collect `.bow` alignment summaries, emit unsorted BAMs.
6. **Merge BAMs** – samtools merge/sort/index per sample.
7. **Fragment filtering + deduplication** – reproduce the manual `samtools view -f 67 -q 30` + AWK counting + `cuttag_dedup_multiAb.py`, emitting fragment-level BAMs, frequency tables, and dedup reports.

All intermediate artefacts follow the original folder layout (`split/`, `valid/`, `cut/`, `align/`, `merge/`, `fragments/`, `dedup/`), making it easy to cross-validate with the manually generated files.

## Running the pipeline

```
nextflow run scCUTTag_10x.nf \
    -profile local \
    --samplesheet samples.csv \
    --antibody_sheet input_antibody_barcode.csv \
    --extract_script /path/to/extractValidReads_multiCUTTag_nobccorrect_hamdist8.py \
    --cell_barcode_reference /path/to/737K-cratac-v1.txt \
    --bowtie2_index /share/.../mm10/Bowtie2Index/genome \
    --dedup_script /path/to/cuttag_dedup_multiAb.py
```

Adjust profiles/resources in `nextflow.config` to match your environment (local, LSF, etc.).
