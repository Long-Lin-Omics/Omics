
rule aligned_dedup_bam_plot_fragment:
    input:
        "{output_dir}/aligned/{sample}_dedup.bam"
    output:
       fragment_file = "{output_dir}/aligned/{sample}_fragment_lengths.txt",
       fragment_plot =  "{output_dir}/aligned/{sample}_fragment.hist.png"
    threads: rule_threads_n
    shell:
        "samtools view --threads {threads} {input} | awk '$9 > 0 {{print $9}}' > {output.fragment_file} && Rscript /ddn/gs1/home/linl7/bin/scripts/chip-seq/plot_fragment_length.R  {output.fragment_file}  {output.fragment_plot} {wildcards.sample}"
