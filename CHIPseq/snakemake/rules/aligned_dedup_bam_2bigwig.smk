
rule aligned_dedup_bam_2bigwig:
    input:
        "{output_dir}/aligned/{sample}_dedup.bam"
    output:
        normalized="{output_dir}/ucsc/{sample}_normalized.bw",
        unnormalized="{output_dir}/ucsc/{sample}_unnormalized.bw"
    threads: rule_threads_n
    shell:
        "samtools index --threads {threads} {input} && bamCoverage --numberOfProcessors {threads} -b {input} -o {output.normalized} --normalizeUsing BPM"
        " &&  bamCoverage --numberOfProcessors {threads} -b {input} -o {output.unnormalized}"
