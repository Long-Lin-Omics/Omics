#!/bin/bash

# Set variables
THREADS=8
GENOME_INDEX=/path/to/mm10/genome_index  # <-- replace with your mm10 Bowtie2 index prefix
OUTDIR=./results
FASTQ_DIR=./fastq
TRIMGALORE_PATH=trim_galore
BOWTIE2_PATH=bowtie2
SAMTOOLS_PATH=samtools
MACS2_PATH=macs2

mkdir -p ${OUTDIR}/trimmed ${OUTDIR}/aligned ${OUTDIR}/macs2

# Process each sample pair
for SAMPLE in SRR11294129 SRR11294130
do
    echo "Processing sample: $SAMPLE"

    # Step 1: Trim adapters and low-quality bases
    ${TRIMGALORE_PATH} --paired --phred33 --quality 20 --stringency 1 -e 0.1 --length 20 \
        -o ${OUTDIR}/trimmed \
        ${FASTQ_DIR}/${SAMPLE}_1.fastq.gz ${FASTQ_DIR}/${SAMPLE}_2.fastq.gz

    # Get trimmed file names
    R1_TRIM="${OUTDIR}/trimmed/${SAMPLE}_1_val_1.fq.gz"
    R2_TRIM="${OUTDIR}/trimmed/${SAMPLE}_2_val_2.fq.gz"

    # Step 2: Align to mm10 genome with Bowtie2
    ${BOWTIE2_PATH} -p ${THREADS} -x ${GENOME_INDEX} -1 ${R1_TRIM} -2 ${R2_TRIM} \
        -I 50 -X 800 --fr -N 0 -L 22 -i 'S,1,1.15' --n-ceil 'L,0,0.15' --dpad 15 --gbar 4 \
        --end-to-end --score-min 'L,-0.6,-0.6' \
        -S ${OUTDIR}/aligned/${SAMPLE}.sam

    # Step 3: Convert SAM to BAM, filter MAPQ ≥20, sort, and index
    ${SAMTOOLS_PATH} view -@ ${THREADS} -b -q 20 ${OUTDIR}/aligned/${SAMPLE}.sam | \
        ${SAMTOOLS_PATH} sort -@ ${THREADS} -o ${OUTDIR}/aligned/${SAMPLE}.sorted.bam

    ${SAMTOOLS_PATH} index ${OUTDIR}/aligned/${SAMPLE}.sorted.bam

    # Cleanup intermediate SAM file
    rm ${OUTDIR}/aligned/${SAMPLE}.sam

done

# Step 4: Call peaks with MACS2 for DPPA2 (example sample SRR11294130, with IgG as control)
${MACS2_PATH} callpeak -t ${OUTDIR}/aligned/SRR11294130.sorted.bam \
    -c ${OUTDIR}/aligned/SRR11294129.sorted.bam \
    -f BAMPE -g mm --outdir ${OUTDIR}/macs2 -n DPPA2_peaks --qvalue 1e-5 --nomodel --extsize 200

echo "Pipeline complete!"

