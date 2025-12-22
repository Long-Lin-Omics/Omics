import os

#configfile: "config.yaml"

cases = list(config["cases"].keys())
controls = list(config["controls"].keys())
all_samples = config["cases"] | config["controls"]
sample_names = list(all_samples.keys())
comparisons = list(config["comparisons"].keys())
output_dir = config["output_dir"]
identifier = config['identifier']
fai=config['fai']
ref_to_annotate=config['ref_to_annotate']
blacklist=config['blacklist']
scripts_folder = config['scripts_folder']

# rule_threads_n = config['rule_threads_n']

rule all:
    input:
        #expand("{output_dir}/fastqc_raw/{sample}_R1_fastqc.html", sample=sample_names, output_dir=output_dir),
        # expand("{output_dir}/fastqc_trimmed/{sample}_trimmed_fastqc.html", sample=sample_names, output_dir=output_dir),
        # expand("{output_dir}/trimmed/{sample}_R1_trimmed.fastq.gz", sample=sample_names, output_dir=output_dir),
        # expand("{output_dir}/aligned/{sample}.bam", sample=sample_names, output_dir=output_dir),
        # expand("{output_dir}/motif/{case}_homer_results.txt", case=cases, output_dir=output_dir),
        expand("/data/wade/linl7/{identifier}/mm10/trackDb.txt",identifier=identifier),
        # "{output_dir}/multiqc/multiqc_raw_fastq_report.html".format(output_dir=output_dir),
        # "{output_dir}/multiqc/multiqc_trimmed_fastq_report.html".format(output_dir=output_dir),
        expand("{output_dir}/peaks/{case}_peaks.narrowPeak",case=cases,output_dir=output_dir),
        expand("{output_dir}/peaks/{case}_peaks.narrowPeak.clean.annotated",case=cases,output_dir=output_dir),
        expand("{output_dir}/aligned/{sample}_fragment_lengths.txt",sample=sample_names,output_dir=output_dir),
        expand("{output_dir}/ucsc/{sample}_normalized.bw", sample=sample_names, output_dir=output_dir),
        expand("{output_dir}/metaplot/normalized/{comparisons}.xlsx",comparisons=comparisons,output_dir=output_dir),
        expand("{output_dir}/metaplot/unnormalized/{comparisons}.xlsx",comparisons=comparisons,output_dir=output_dir),
        "{output_dir}/report/stat.tsv".format(output_dir=output_dir)

#        expand("{output_dir}/ucsc/{sample}_unnormalized.bw", sample=sample_names, output_dir=output_dir),
        # "{output_dir}/diffbind/differential_binding.csv".format(output_dir=output_dir),
        # "{output_dir}/go_kegg/go_analysis.csv".format(output_dir=output_dir),
        # expand("{output_dir}/comparisons/{comparison}_heatmap.png", comparison=[c["name"] for c in comparisons], output_dir=output_dir),
        # expand("{output_dir}/comparisons/{comparison}_profile.png", comparison=[c["name"] for c in comparisons], output_dir=output_dir)


### --- FastQC Before Trimming --- ###
rule fastqc_raw:
    input:
        fastq1=lambda wildcards: all_samples[wildcards.sample]["fastq1"],
        fastq2=lambda wildcards: all_samples[wildcards.sample]["fastq2"]
    output:
        ####### this makes it only accept fastq.gz in config file. ##########
        fastq1="{output_dir}/input/{sample}_R1.fastq.gz",
        fastq2="{output_dir}/input/{sample}_R2.fastq.gz",
        report="{output_dir}/fastqc_raw/{sample}_R1_fastqc.html"
    shell:
        "ln -s {input.fastq1} {output.fastq1} && ln -s {input.fastq2} {output.fastq2} && {scripts_folder}/CHIPseq/softwares/fastqc {output.fastq1} {output.fastq2} -o {output_dir}/fastqc_raw/"

rule report_fastq:
    input:
        fastq1="{output_dir}/input/{sample}_R1.fastq.gz",
        fastq2="{output_dir}/input/{sample}_R2.fastq.gz"
    output:
        "{output_dir}/input/{sample}.fastq.count.txt"
    shell:
        "source {scripts_folder}/shell/fastq.sh && pair_fq_count_check {input.fastq1} {input.fastq2} {output}"

rule multiqc_raw:
    input:
        expand("{output_dir}/fastqc_raw/{sample}_R1_fastqc.html", sample=sample_names, output_dir=output_dir),
    output:
        "{output_dir}/multiqc/multiqc_raw_fastq_report.html"
    shell:
        "{scripts_folder}/CHIPseq/softwares/multiqc {output_dir}/fastqc_raw/ -o {output_dir}/multiqc/ -n multiqc_raw_fastq_report"

### --- Cutadapt for Adapter & Quality Trimming --- ###
rule cutadapt:
    input:
        fastq1=lambda wildcards: all_samples[wildcards.sample]["fastq1"],
        fastq2=lambda wildcards: all_samples[wildcards.sample]["fastq2"]
    output:
        trimmed_fastq1="{output_dir}/trimmed/{sample}_R1_trimmed.fastq.gz",
        trimmed_fastq2="{output_dir}/trimmed/{sample}_R2_trimmed.fastq.gz"
    shell:
        "{scripts_folder}/CHIPseq/softwares/cutadapt -q 20 -m 20 -a AGATCGGAAGAGC -A AGATCGGAAGAGC "
        "-o {output.trimmed_fastq1} -p {output.trimmed_fastq2} {input.fastq1} {input.fastq2}"

### --- FastQC After Trimming --- ###
rule fastqc_trimmed:
    input:
        trimmed_fastq1="{output_dir}/trimmed/{sample}_R1_trimmed.fastq.gz",
        trimmed_fastq2="{output_dir}/trimmed/{sample}_R2_trimmed.fastq.gz"
    output:
        "{output_dir}/fastqc_trimmed/{sample}_R1_trimmed_fastqc.html"
    shell:
        "{scripts_folder}/CHIPseq/softwares/fastqc {input.trimmed_fastq1} {input.trimmed_fastq2} -o {output_dir}/fastqc_trimmed/"

rule report_trim_fastq:
    input:
        trimmed_fastq1="{output_dir}/trimmed/{sample}_R1_trimmed.fastq.gz",
        trimmed_fastq2="{output_dir}/trimmed/{sample}_R2_trimmed.fastq.gz"
    output:
        "{output_dir}/trimmed/{sample}.trimmed.fastq.count.txt"
    shell:
        "source {scripts_folder}/shell/fastq.sh && pair_fq_count_check {input.trimmed_fastq1} {input.trimmed_fastq2} {output}"

rule multiqc_trimmed:
    input:
        expand("{output_dir}/fastqc_trimmed/{sample}_R1_trimmed_fastqc.html", sample=sample_names, output_dir=output_dir)
    output:
        "{output_dir}/multiqc/multiqc_trimmed_fastq_report.html"
    shell:
        "{scripts_folder}/CHIPseq/softwares/multiqc {output_dir}/fastqc_trimmed/ -o {output_dir}/multiqc -n multiqc_trimmed_fastq_report"


### --- Continue with Alignment and Analysis --- ###
rule bowtie2:
    input:
        trimmed_fastq1="{output_dir}/trimmed/{sample}_R1_trimmed.fastq.gz",
        trimmed_fastq2="{output_dir}/trimmed/{sample}_R2_trimmed.fastq.gz"
    output:
        bam = "{output_dir}/aligned/{sample}.bam",
        sortbam =  "{output_dir}/aligned/{sample}_dedup.bam",
        dup_mets = "{output_dir}/aligned/{sample}_dup_metrics.txt"
    params:
        rg_id = "{sample}",
        rg_lb = "Unknown",
        rg_pl = "ILLUMINA",
        rg_pu = "Unknown",
        rg_sm = "{sample}"
    shell:
        "{scripts_folder}/CHIPseq/softwares/bowtie2 -p 8 --local -x {config[bowtie2_index]} -1 {input.trimmed_fastq1} -2 {input.trimmed_fastq2} --rg-id {params.rg_id} --rg LB:{params.rg_lb} --rg PL:{params.rg_pl} --rg PU:{params.rg_pu} --rg SM:{params.rg_sm}"
        "| {scripts_folder}/CHIPseq/softwares/samtools view -Sb - | {scripts_folder}/CHIPseq/softwares/samtools sort -o {output.bam}"
        " && {scripts_folder}/CHIPseq/softwares/java -Xmx90g -jar {scripts_folder}/CHIPseq/softwares/picard.jar MarkDuplicates I={output.bam} O={output.sortbam} M={output.dup_mets} REMOVE_DUPLICATES=true"


# rule mark_duplicates:
#     input:
#         "{output_dir}/aligned/{sample}.bam"
#     output:
#         "{output_dir}/aligned/{sample}_dedup.bam"
#     shell:
#         "java -Xmx2g -jar /ddn/gs1/biotools/picard/picard.jar MarkDuplicates I={input} O={output} M={output_dir}/aligned/{wildcards.sample}_dup_metrics.txt REMOVE_DUPLICATES=true"

rule plot_fragment:
    input:
        "{output_dir}/aligned/{sample}_dedup.bam"
    output:
       fragment_file = "{output_dir}/aligned/{sample}_fragment_lengths.txt",
       fragment_plot =  "{output_dir}/aligned/{sample}_fragment.hist.png"
    shell:
        "{scripts_folder}/CHIPseq/softwares/samtools view {input} | awk '$9 > 0 {{print $9}}' > {output.fragment_file} && Rscript {scripts_folder}/R/plot_fragment_length.R  {output.fragment_file}  {output.fragment_plot} {wildcards.sample}"

rule ucsc_bam2bigwig:
    input:
        "{output_dir}/aligned/{sample}_dedup.bam"
    output:
        normalized="{output_dir}/ucsc/{sample}_normalized.bw",
        unnormalized="{output_dir}/ucsc/{sample}_unnormalized.bw"
    shell:
        "{scripts_folder}/CHIPseq/softwares/samtools index {input} && {scripts_folder}/CHIPseq/softwares/bamCoverage -b {input} -o {output.normalized} --normalizeUsing BPM &&  {scripts_folder}/CHIPseq/softwares/bamCoverage -b {input} -o {output.unnormalized}"

rule macs2:
    input:
        case = lambda wildcards: "{output_dir}/aligned/{case}_dedup.bam".format(output_dir=output_dir,case=wildcards.case),
        control = lambda wildcards: "{output_dir}/aligned/{control}_dedup.bam".format(output_dir=output_dir, control=config["cases"][wildcards.case]["control"])
    output:
        "{output_dir}/peaks/{case}_peaks.narrowPeak"
    shell:
        "{scripts_folder}/CHIPseq/softwares/macs2 callpeak -t {input.case} -c {input.control} -f BAMPE -g {config[macs2_genome_size]} "
        "--outdir {output_dir}/peaks/ -n {wildcards.case}"

rule clean_and_annotate_peak:
    input:
        "{output_dir}/peaks/{case}_peaks.narrowPeak"
    output:
        filt = "{output_dir}/peaks/{case}_peaks.narrowPeak.filtered",
        clean = "{output_dir}/peaks/{case}_peaks.narrowPeak.clean",
        annotated = "{output_dir}/peaks/{case}_peaks.narrowPeak.clean.annotated"
    shell:
        "awk '$9 >= 4' {input}  > {output.filt} && {scripts_folder}/CHIPseq/softwares/bedtools subtract -a {output.filt} -b {blacklist} > {output.clean} && "
        "perl {scripts_folder}/CHIPseq/softwares/annotatePeaks.pl {output.clean} {ref_to_annotate} > {output.annotated}"



rule uscs_peak2bed:
    input:
        "{output_dir}/peaks/{case}_peaks.narrowPeak.clean"
    output:
        "{output_dir}/ucsc/{case}.peaks.bb"
    shell:
        "Rscript {scripts_folder}/CHIPseq/peak_to_bigbed.R {input} {output} {fai}"


# rule diffbind:
#     input:
#         expand("{output_dir}/aligned/{comp_case}_dedup.bam", comp_case=config["comparisons"][wildcards.comparisons], output_dir=output_dir),
#     output:
#         "{output_dir}/comparisons/{comp_case}_peaks.narrowPeak"
#     output:
#         expand("{output_dir}/comparisons/{comparison}_heatmap.png", comparison=[c["name"] for c in comparisons], output_dir=output_dir),


rule report_stat:
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
        "python {scripts_folder}/python/pipeline_report_stat.py --fastq {input.fqs} --trimmed_fastq {input.trim_fqs} --fastq_report --metrics {input.dup_mets} --peaks {input.peaks} --samples {params.samples} --output {output}"

rule ucsc_hub:
    input:
        bw=expand("{output_dir}/ucsc/{sample}_normalized.bw", sample=sample_names,output_dir=output_dir),
        peaks=expand("{output_dir}/ucsc/{case}.peaks.bb", case=cases,output_dir=output_dir)
    output:
        hub_dir=directory("/data/wade/linl7/{identifier}/mm10/"),
        hub_file="/data/wade/linl7/{identifier}/hub.txt",
        genomes_file="/data/wade/linl7/{identifier}/genomes.txt",
        trackdb_file="/data/wade/linl7/{identifier}/mm10/trackDb.txt"
    params:
        sorted_sample_names = sorted(sample_names),
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


rule mataplot:
    input: 
        peaks=lambda wildcards: expand("{output_dir}/peaks/{case}_peaks.narrowPeak.clean",output_dir=output_dir,case=config["comparisons"][wildcards.comparisons]),
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
    shell:
        "cat {input.peaks} | sort -k1,1 -k2,2n | bedtools merge -d 100 -i - > {output.merged_bed} && "
        "{scripts_folder}/CHIPseq/softwares/computeMatrix reference-point -p 8 --referencePoint center -S {input.bws} -R {output.merged_bed} "
        "--beforeRegionStartLength 2000 --afterRegionStartLength 2000 --binSize 40 -o {output.matrix} --outFileNameMatrix {output.tsv} && "
        "python {scripts_folder}/python/metaplot.py {output.tsv} {output.xlsx} {params.n} && "
        "{scripts_folder}/CHIPseq/softwares/computeMatrix reference-point -p 8 --referencePoint center -S {input.un_bws} -R {output.merged_bed} "
        "--beforeRegionStartLength 2000 --afterRegionStartLength 2000 --binSize 40 -o {output.un_matrix} --outFileNameMatrix {output.un_tsv} && "
        "python {scripts_folder}/python/metaplot.py {output.un_tsv} {output.un_xlsx} {params.n}"
        



# rule homer_motif:
#     input:
#         "{output_dir}/peaks/{case}_peaks.narrowPeak"
#     output:
#         "{output_dir}/motif/{case}_homer_results.txt"
#     shell:
#         "findMotifsGenome.pl {input} {config[genome]} {output_dir}/motif/{wildcards.case}/"

# rule diffbind:
#     input:
#         expand("{output_dir}/peaks/{case}_peaks.narrowPeak", case=cases, output_dir=output_dir)
#     output:
#         "{output_dir}/diffbind/differential_binding.csv"
#     shell:
#         "Rscript scripts/run_diffbind.R {output_dir}/peaks/ {output}"

# rule go_kegg:
#     input:
#         "{output_dir}/diffbind/differential_binding.csv"
#     output:
#         "{output_dir}/go_kegg/go_analysis.csv"
#     shell:
#         "Rscript scripts/run_go_kegg.R {input} {output}"


# rule compute_matrix:
#     input:
#         case_bw="{output_dir}/ucsc/{case}.bigwig",
#         control_bw=lambda wildcards: "{output_dir}/ucsc/{control}.bigwig".format(output_dir=output_dir,control=cases[wildcards.case]["control"]),
#         regions=config["regions_bed"]
#     output:
#         "{output_dir}/comparisons/{comparison}_matrix.gz"
#     params:
#         before=3000, after=3000, bin_size=50
#     shell:
#         "computeMatrix reference-point --referencePoint center -S {input.case_bw} {input.control_bw} "
#         "-R {input.regions} -o {output} --beforeRegionStartLength {params.before} --afterRegionStartLength {params.after} "
#         "--binSize {params.bin_size}"

# rule plot_heatmap:
#     input:
#         "{output_dir}/comparisons/{comparison}_matrix.gz"
#     output:
#         "{output_dir}/comparisons/{comparison}_heatmap.png"
#     shell:
#         "plotHeatmap -m {input} -out {output}"

# rule plot_profile:
#     input:
#         "{output_dir}/comparisons/{comparison}_matrix.gz"
#     output:
#         "{output_dir}/comparisons/{comparison}_profile.png"
#     shell:
#         "plotProfile -m {input} -out {output}"


