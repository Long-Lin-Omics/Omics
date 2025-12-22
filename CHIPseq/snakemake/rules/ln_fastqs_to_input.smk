
rule ln_fastqs_to_input:
    input:
        fastq1=lambda wildcards: all_samples[wildcards.sample]["fastq1"],
        fastq2=lambda wildcards: all_samples[wildcards.sample]["fastq2"]
    output:
        ####### this makes it only accept fastq.gz in config file. ##########
        fastq1="{output_dir}/input/{sample}_R1.fastq.gz",
        fastq2="{output_dir}/input/{sample}_R2.fastq.gz"
    shell:
        "ln -s {input.fastq1} {output.fastq1} && ln -s {input.fastq2} {output.fastq2}"
