
rule stats_from_fqs_trim_dup_peaks:
    input:  
        fqs=expand("{output_dir}/input/{sample}.fastq.count.txt", output_dir=output_dir,sample=sample_names),
        trim_fqs=expand("{output_dir}/trimmed/{sample}.trimmed.fastq.count.txt",output_dir=output_dir,sample=sample_names),
        dup_mets=expand("{output_dir}/aligned/{sample}_dup_metrics.txt",output_dir=output_dir,sample=sample_names),
        peaks=expand("{output_dir}/peaks/{sample}_peaks.narrowPeak.clean",output_dir=output_dir,sample=cases)
    output:
        "{output_dir}/report/{identifier}.stat.tsv"
    params:
        samples=expand(sorted(sample_names)) # {params} will show [] in shell
    shell:
        "python /ddn/gs1/home/linl7/bin/scripts/pipeline_report_stat.py --fastq {input.fqs} --trimmed_fastq {input.trim_fqs} --fastq_report --metrics {input.dup_mets} --peaks {input.peaks} --samples {params.samples} --output {output}"
