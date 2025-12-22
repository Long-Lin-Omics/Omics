
rule trimmed_multiqc:
    input:
        expand("{output_dir}/fastqc_trimmed/{sample}_R1_trimmed_fastqc.html", sample=sample_names, output_dir=output_dir)
    output:
        "{output_dir}/multiqc/multiqc_trimmed_fastq_report.html"
    shell:
        "multiqc {output_dir}/fastqc_trimmed/ -o {output_dir}/multiqc -n multiqc_trimmed_fastq_report"
