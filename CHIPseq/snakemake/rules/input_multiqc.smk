
rule input_multiqc:
    input:
        expand("{output_dir}/fastqc_raw/{sample}_R1_fastqc.html", sample=sample_names, output_dir=output_dir)
    output:
        "{output_dir}/multiqc/multiqc_raw_fastq_report.html"
    shell:
        "multiqc {output_dir}/fastqc_raw/ -o {output_dir}/multiqc/ -n multiqc_raw_fastq_report"