

rule comparison_metaplot:
    input: 
        peaks=lambda wildcards: expand("{output_dir}/peaks/{case}_peaks.narrowPeak",output_dir=output_dir,case=config["comparisons"][wildcards.comparisons]),
        bws=lambda wildcards: expand("{output_dir}/ucsc/{case}_normalized.bw",output_dir=output_dir,case=config["comparisons"][wildcards.comparisons]),
        un_bws=lambda wildcards: expand("{output_dir}/ucsc/{case}_unnormalized.bw",output_dir=output_dir,case=config["comparisons"][wildcards.comparisons]) 
    output:
        xlsx="{output_dir}/metaplot/normalized/{identifier}.{comparisons}.xlsx",
        merged_bed="{output_dir}/metaplot/normalized/{identifier}.{comparisons}.merged.bed",
        matrix="{output_dir}/metaplot/normalized/{identifier}.{comparisons}.matrix.gz",
        tsv="{output_dir}/metaplot/normalized/{identifier}.{comparisons}.matrix.tsv",
        un_xlsx="{output_dir}/metaplot/unnormalized/{identifier}.{comparisons}.xlsx",
        un_matrix="{output_dir}/metaplot/unnormalized/{identifier}.{comparisons}.matrix.gz",
        un_tsv="{output_dir}/metaplot/unnormalized/{identifier}.{comparisons}.matrix.tsv"
    params:
        n=lambda wildcards: len(config["comparisons"][wildcards.comparisons])
    threads: rule_threads_n
    shell:
        "cat {input.peaks} | sort -k1,1 -k2,2n | bedtools merge -d 100 -i - > {output.merged_bed} && "
        "computeMatrix reference-point -p {threads} --referencePoint center -S {input.bws} -R {output.merged_bed} "
        "--beforeRegionStartLength 2000 --afterRegionStartLength 2000 --binSize 40 -o {output.matrix} --outFileNameMatrix {output.tsv} && "
        "python /ddn/gs1/home/linl7/bin/scripts/chip-seq/metaplot.py {output.tsv} {output.xlsx} {params.n} && "
        "computeMatrix reference-point -p {threads} --referencePoint center -S {input.un_bws} -R {output.merged_bed} "
        "--beforeRegionStartLength 2000 --afterRegionStartLength 2000 --binSize 40 -o {output.un_matrix} --outFileNameMatrix {output.un_tsv} && "
        "python /ddn/gs1/home/linl7/bin/scripts/chip-seq/metaplot.py {output.un_tsv} {output.un_xlsx} {params.n}"