
rule peaks_clean_peaks_2bed:
    input:
        "{output_dir}/peaks/{case}_peaks.narrowPeak.clean"
    output:
        "{output_dir}/ucsc/{case}.peaks.bb"
    shell:
        "Rscript ~/bin/scripts/chip-seq/peak_to_bigbed.R {input} {output} {fai}"
