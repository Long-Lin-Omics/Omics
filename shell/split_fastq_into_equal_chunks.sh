#!/bin/bash

# Usage:
# ./split_fastq_sequential.sh input.fastq.gz 5 output_prefix

INPUT=$1
NSAMPLES=$2
PREFIX=$3

TOTAL_READS=$(zcat "$INPUT" | echo $(( $(wc -l) / 4 )))
READS_PER_SAMPLE=$(( TOTAL_READS / NSAMPLES ))

zcat "$INPUT" \
| split -l $(( READS_PER_SAMPLE * 4 )) -d -a 2 --additional-suffix=.fastq -

i=0
for f in x*.fastq; do
    gzip "$f"
    mv "$f.gz" "${PREFIX}_pseudo${i}.fastq.gz"
    i=$((i + 1))
done
