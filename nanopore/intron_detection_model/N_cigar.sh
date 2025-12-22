#!/bin/bash

# Usage: ./check_splicing_with_unmapped.sh input.bam output.txt

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.bam> <output.txt>"
    exit 1
fi

INPUT_BAM="$1"
OUTPUT_FILE="$2"

if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input BAM file '$INPUT_BAM' not found!"
    exit 1
fi

# 输出每条 read 的状态：UNMAPPED / SPLICED / NO_N
samtools view "$INPUT_BAM" | \
awk 'BEGIN{OFS="\t"}{
    flag = $2;
    cigar = $6;

    if (and(flag, 4)) {
        print $1, "UNMAPPED";
    } else {
        if (cigar ~ /N/) {
            print $1, "SPLICED";
        } else {
            print $1, "NO_N";
        }
    }
}' > "$OUTPUT_FILE"

echo "Output written to $OUTPUT_FILE"

cut  -f 2 $OUTPUT_FILE |sort |uniq -c