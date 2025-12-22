
rule trimmed_fastqc:
    input:
        trimmed_fastq1="{output_dir}/trimmed/{sample}_R1_trimmed.fastq.gz",
        trimmed_fastq2="{output_dir}/trimmed/{sample}_R2_trimmed.fastq.gz"
    output:
        "{output_dir}/fastqc_trimmed/{sample}_R1_trimmed_fastqc.html"
    threads: rule_threads_n
    shell:
        "fastqc -t {threads} {input.trimmed_fastq1} {input.trimmed_fastq2} -o {output_dir}/fastqc_trimmed/"
