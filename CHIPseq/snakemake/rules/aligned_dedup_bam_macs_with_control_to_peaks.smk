
rule aligned_dedup_bam_macs_with_control_to_peaks:
    input:
        case = lambda wildcards: "{output_dir}/aligned/{case}_dedup.bam".format(output_dir=output_dir,case=wildcards.case),
        control = lambda wildcards: "{output_dir}/aligned/{control}_dedup.bam".format(output_dir=output_dir, control=config["cases"][wildcards.case]["control"])
    output:
        "{output_dir}/peaks/{case}_peaks.narrowPeak"
    shell:
        "macs2 callpeak -t {input.case} -c {input.control} -f BAMPE -g {config[macs2_genome_size]} "
        "--outdir {output_dir}/peaks/ -n {wildcards.case}"
