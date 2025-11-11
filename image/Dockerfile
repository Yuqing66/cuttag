# Base OS (arch is chosen at build time; on Apple Silicon, build with --platform linux/amd64)
FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
ENV MAMBA_ROOT_PREFIX=/opt/conda

# Minimal system deps + pigz for fast gz handling
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl ca-certificates bzip2 bash coreutils pigz \
    && rm -rf /var/lib/apt/lists/*

# Micromamba (lighter/faster than conda)
RUN curl -L https://micro.mamba.pm/api/micromamba/linux-64/latest \
  | tar -xvj -C /usr/local/bin/ bin/micromamba --strip-components=1

SHELL ["/bin/bash", "-lc"]

RUN micromamba create -y -n env -c conda-forge -c bioconda \
    "python=3.8" \
    "bowtie2=2.5.0" \
    "samtools=1.16.1" \
    "bedtools=2.30.0" \
    "macs2=2.2.7.1" \
    "fastqc" \
    "cutadapt=4.1" \
    "pysam" \
    "biopython" \
    && micromamba clean -a -y

# Put the env on PATH
ENV PATH=/opt/conda/envs/env/bin:$PATH

# Default workdir
WORKDIR /data

# Default shell
CMD ["bash"]