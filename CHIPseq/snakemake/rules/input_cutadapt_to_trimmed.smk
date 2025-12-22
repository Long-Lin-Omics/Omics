
rule input_cutadapt_to_trimmed:
    input:
        fastq1="{output_dir}/input/{sample}_R1.fastq.gz",
        fastq2="{output_dir}/input/{sample}_R2.fastq.gz"
    output:
        trimmed_fastq1="{output_dir}/trimmed/{sample}_R1_trimmed.fastq.gz",
        trimmed_fastq2="{output_dir}/trimmed/{sample}_R2_trimmed.fastq.gz"
    threads: rule_threads_n
    shell:
        "cutadapt -j {threads} -q 20 -m 20 -a AGATCGGAAGAGC -A AGATCGGAAGAGC "
        "-o {output.trimmed_fastq1} -p {output.trimmed_fastq2} {input.fastq1} {input.fastq2}"
