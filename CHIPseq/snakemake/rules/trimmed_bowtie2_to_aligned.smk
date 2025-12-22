
rule trimmed_bowtie2_to_aligned:
    input:
        trimmed_fastq1="{output_dir}/trimmed/{sample}_R1_trimmed.fastq.gz",
        trimmed_fastq2="{output_dir}/trimmed/{sample}_R2_trimmed.fastq.gz"
    output:
        bam = "{output_dir}/aligned/{sample}.bam"
    threads: rule_threads_n
    params:
        rg_id = "{sample}",
        rg_lb = "Unknown",
        rg_pl = "ILLUMINA",
        rg_pu = "Unknown",
        rg_sm = "{sample}"
    shell:
        "bowtie2 -p {threads} --local -x {config[bowtie2_index]} -1 {input.trimmed_fastq1} -2 {input.trimmed_fastq2} --rg-id {params.rg_id} --rg LB:{params.rg_lb} --rg PL:{params.rg_pl} --rg PU:{params.rg_pu} --rg SM:{params.rg_sm}"
        "| samtools view --threads {threads} -Sb - | samtools sort --threads {threads} -o {output.bam}"

