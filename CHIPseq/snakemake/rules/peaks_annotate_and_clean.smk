
rule peaks_annotate_and_clean:
    input:
        "{output_dir}/peaks/{case}_peaks.narrowPeak"
    output:
        filt = "{output_dir}/peaks/{case}_peaks.narrowPeak.filtered",
        clean = "{output_dir}/peaks/{case}_peaks.narrowPeak.clean",
        annotated = "{output_dir}/peaks/{case}_peaks.narrowPeak.clean.annotated"
    shell:
        "awk '$9 >= 4' {input}  > {output.filt} && bedtools subtract -a {output.filt} -b {blacklist} > {output.clean} && "
        "perl /ddn/gs1/home/linl7/bin/homer/bin/annotatePeaks.pl {output.clean} {ref_to_annotate} > {output.annotated}"
