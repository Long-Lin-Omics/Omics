
rule aligned_markduplicate:
    input:
        bam="{output_dir}/aligned/{sample}.bam"
    output:
        dedup_bam = "{output_dir}/aligned/{sample}_dedup.bam",
        dup_mets = "{output_dir}/aligned/{sample}_dup_metrics.txt"
    shell:
        "java -Xmx12g -jar /ddn/gs1/biotools/picard/picard.jar MarkDuplicates I={output.bam} O={output.dedup_bam} M={output.dup_mets} REMOVE_DUPLICATES=true"
