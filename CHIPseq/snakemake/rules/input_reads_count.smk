
rule input_reads_count:
    input:
        fastq1="{output_dir}/input/{sample}_R1.fastq.gz",
        fastq2="{output_dir}/input/{sample}_R2.fastq.gz"
    output:
        "{output_dir}/input/{sample}.fastq.count.txt"
    shell:
        "source /ddn/gs1/home/linl7/bin/scripts/shell/fastq.sh && pair_fq_count_check {input.fastq1} {input.fastq2} {output}"
