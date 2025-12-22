import os

#configfile: "config.yaml"

cases = list(config["cases"].keys())
output_dir = config["output_dir"]
identifier = config['identifier']
fai=config['fai']
ref_to_annotate=config['ref_to_annotate']


rule all:
    input:
        expand("{output_dir}/peaks/{sample}_peaks.narrowPeak",sample=cases,output_dir=output_dir),
        expand("{output_dir}/ucsc/{sample}_normalized.bw", sample=cases, output_dir=output_dir),
        "/data/wade/linl7/{identifier}/hub.txt".format(identifier=identifier)
        
rule link_input:
    input:
        fastq1=lambda wildcards: config["cases"][wildcards.sample]["fastq1"],
        fastq2=lambda wildcards: config["cases"][wildcards.sample]["fastq2"]
    output:
        ####### this makes it only accept fastq.gz in config file. ##########
        fastq1="{output_dir}/input/{sample}_R1.fastq.gz",
        fastq2="{output_dir}/input/{sample}_R2.fastq.gz"
    shell:
        "ln -s {input.fastq1} {output.fastq1} && ln -s {input.fastq2} {output.fastq2} "

rule trim_galore:
    input:
        fastq1="{output_dir}/input/{sample}_R1.fastq.gz",
        fastq2="{output_dir}/input/{sample}_R2.fastq.gz"
    output:
        trimmed_fastq1="{output_dir}/trimmed/{sample}_R1_val_1.fq.gz",
        trimmed_fastq2="{output_dir}/trimmed/{sample}_R2_val_2.fq.gz"
    shell:
        "trim_galore  --paired --phred33 --quality 20 --stringency 1 -e 0.1 --length 20 -o {output_dir}/trimmed/ {input.fastq1} {input.fastq2}"

rule bowtie2:
    input:
        trimmed_fastq1="{output_dir}/trimmed/{sample}_R1_val_1.fq.gz",
        trimmed_fastq2="{output_dir}/trimmed/{sample}_R2_val_2.fq.gz"
    output:
        bam = "{output_dir}/aligned/{sample}.sorted.bam"
    params:
        rg_id = "{sample}",
        rg_lb = "Unknown",
        rg_pl = "ILLUMINA",
        rg_pu = "Unknown",
        rg_sm = "{sample}"
    shell:
        "bowtie2 -p 8 --local -x {config[bowtie2_index]} -1 {input.trimmed_fastq1} -2 {input.trimmed_fastq2} --rg-id {params.rg_id} --rg LB:{params.rg_lb} --rg PL:{params.rg_pl} --rg PU:{params.rg_pu} --rg SM:{params.rg_sm} "
        "-I 50 -X 800 --fr -N 0 -L 22 -i 'S,1,1.15' --n-ceil 'L,0,0.15' --dpad 15 --gbar 4 --end-to-end --score-min 'L,-0.6,-0.6'"
        "| samtools view -@ 8 -q 20 -Sb - | samtools sort -@ 8 -o {output.bam} && samtools index  -@ 8 {output.bam}"

rule ucsc_bam2bigwig:
    input:
        case="{output_dir}/aligned/{sample}.sorted.bam"
    output:
        normalized="{output_dir}/ucsc/{sample}_normalized.bw",
        unnormalized="{output_dir}/ucsc/{sample}_unnormalized.bw"
    shell:
        "bamCoverage -b {input} -o {output.normalized} --normalizeUsing BPM &&  bamCoverage -b {input} -o {output.unnormalized}"

rule macs2:
    input:
        case="{output_dir}/aligned/{sample}.sorted.bam"
    output:
        "{output_dir}/peaks/{sample}_peaks.narrowPeak"
    shell:
        "macs2 callpeak -t {input.case} -f BAMPE -g {config[macs2_genome_size]} "
        "--outdir {output_dir}/peaks/ -n {wildcards.sample} --pvalue 1e-5 --nomodel --extsize 200"

rule uscs_peak2bed:
    input:
        "{output_dir}/peaks/{case}_peaks.narrowPeak"
    output:
        "{output_dir}/ucsc/{case}.peaks.bb"
    shell:
        "Rscript ~/bin/scripts/chip-seq/peak_to_bigbed.R {input} {output} {fai}"

rule ucsc_hub:
    input:
        bw=expand("{output_dir}/ucsc/{sample}_normalized.bw", sample=cases,output_dir=output_dir),
        peaks=expand("{output_dir}/ucsc/{case}.peaks.bb", case=cases,output_dir=output_dir)
    output:
        hub_dir=directory("/data/wade/linl7/{identifier}/mm10/"),
        hub_file="/data/wade/linl7/{identifier}/hub.txt",
        genomes_file="/data/wade/linl7/{identifier}/genomes.txt",
        trackdb_file="/data/wade/linl7/{identifier}/mm10/trackDb.txt"
    params:
        sorted_sample_names = sorted(cases),
        sorted_cases_names = sorted(cases)
    shell:
        """
        #mkdir -p {output.hub_dir}
        cp {input.bw} {input.peaks} {output.hub_dir}
    
        # Create the main hub file
        echo "hub ChIPseq_Hub" > {output.hub_file}
        echo "shortLabel ChIP-seq Hub" >> {output.hub_file}
        echo "longLabel ChIP-seq Data Hub" >> {output.hub_file}
        echo "genomesFile genomes.txt" >> {output.hub_file}
        echo "email long.lin@nih.gov" >> {output.hub_file}
        
        # Specify the genome version
        echo "genome mm10" > {output.genomes_file}
        echo "trackDb mm10/trackDb.txt" >> {output.genomes_file}
        
        # Create trackDb.txt with all samples
        echo -e "track cont\ncontainer multiWig\nshortLabel {identifier}\nlongLabel {identifier}\ntype bigWig\nvisibility full\nalwaysZero on\nautoScale on\nshowSubtrackColorOnUi on\naggregate none\ncolor 1,1,1\npriority 1\n " > {output.trackdb_file}
        
        index=1


        for sample in {params.sorted_sample_names};
        do 
            echo "track $sample " >> {output.trackdb_file}
            echo "bigDataUrl ${{sample}}_normalized.bw" >> {output.trackdb_file}
            echo "shortLabel $sample signal" >> {output.trackdb_file}
            echo "longLabel ChIP-seq $sample signal track" >> {output.trackdb_file}
            echo "type bigWig" >> {output.trackdb_file}
            echo "visibility full" >> {output.trackdb_file}
            echo "alwaysZero on" >> {output.trackdb_file}
            echo -e "autoScale on\nmaxHeightPixels 1:60:9999\ncolor 20,20,255\nparent cont" >> {output.trackdb_file}
            ((index++))
            echo "priority $index" >> {output.trackdb_file} 
            echo "" >> {output.trackdb_file}
        done

        for case in {params.sorted_cases_names}
        do
            echo "track ${{case}}_peaks" >> {output.trackdb_file}
            echo "bigDataUrl ${{case}}.peaks.bb" >> {output.trackdb_file}
            echo "shortLabel $case peaks" >> {output.trackdb_file}
            echo "longLabel ChIP-seq $case peak calls" >> {output.trackdb_file}
            echo "type bigBed" >> {output.trackdb_file}
            echo "visibility dense" >> {output.trackdb_file}
            echo -e "color 200,0,0" >> {output.trackdb_file}
            ((index++))
            echo "priority $index" >> {output.trackdb_file}
            echo "" >> {output.trackdb_file}
        done
        """

