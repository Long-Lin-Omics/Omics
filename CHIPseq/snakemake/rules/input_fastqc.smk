
rule input_fastqc:
    input:
        fastq1="{output_dir}/input/{sample}_R1.fastq.gz",
        fastq2="{output_dir}/input/{sample}_R2.fastq.gz"        
    output:
        report="{output_dir}/fastqc_raw/{sample}_R1_fastqc.html"
    threads: rule_threads_n
    shell:
        "fastqc -t {threads} {input.fastq1} {input.fastq2} -o {output_dir}/fastqc_raw/"
