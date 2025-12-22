
rule trimmed_reads_count:
    input:
        trimmed_fastq1="{output_dir}/trimmed/{sample}_R1_trimmed.fastq.gz",
        trimmed_fastq2="{output_dir}/trimmed/{sample}_R2_trimmed.fastq.gz"
    output:
        "{output_dir}/trimmed/{sample}.trimmed.fastq.count.txt"
    shell:
        "source /ddn/gs1/home/linl7/bin/scripts/shell/fastq.sh && pair_fq_count_check {input.trimmed_fastq1} {input.trimmed_fastq2} {output}"