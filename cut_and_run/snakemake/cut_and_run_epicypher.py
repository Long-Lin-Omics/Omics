import os
import re
from collections import defaultdict
import subprocess

#configfile: "config.yaml"

cases = list(config["cases"].keys())
comparisons = list(config["comparisons"].keys())
output_dir = config["output_dir"]
identifier = config['identifier']
fai=config['fai']
ref_to_annotate=config['ref_to_annotate']
blacklist=config['blacklist']
scripts_folder = config['scripts_folder']

if "controls" in config:
    controls = list(config["controls"].keys())
    all_samples = config["cases"] | config["controls"]
else:
    all_samples = config["cases"]


# merge_reps = config.get("merge_fq", False)

# if merge_reps:
#     to_merge_fqs = {}
#     merged_samples = {}
#     group_dict = defaultdict(list)
#     for sample_name in all_samples:
#         m = re.match(r"(.*)-rep-\d+$", sample_name)
#         if m:
#             group_name = m.group(1)
#             group_dict[group_name].append(sample_name)
#     for group_name, reps in group_dict.items():
#         to_merged_fastq1 = [all_samples[r]["fastq1"] for r in reps]
#         to_merged_fastq2 = [all_samples[r]["fastq2"] for r in reps]
#         to_merge_fqs[group_name] = {
#             "fastq1": to_merged_fastq1,
#             "fastq2": to_merged_fastq2
#         }
#         merged_samples[group_name]={
#             "fastq1": output_dir + '/input/merged/' + group_name + '_R1.fastq.gz',
#             "fastq2": output_dir + '/input/merged/' + group_name + '_R2.fastq.gz'
#         }
#     all_samples.update(merged_samples)

sample_names = list(all_samples.keys())

rule all:
    input:
        # expand("{output_dir}/fastqc_raw/{sample}_R1_fastqc.html", sample=sample_names, output_dir=output_dir),
        # # expand("{output_dir}/fastqc_trimmed/{sample}_trimmed_fastqc.html", sample=sample_names, output_dir=output_dir),
        # # expand("{output_dir}/trimmed/{sample}_R1_trimmed.fastq.gz", sample=sample_names, output_dir=output_dir),
        expand("{output_dir}/aligned/{sample}_dedup.bam", sample=sample_names, output_dir=output_dir),
        # expand("{output_dir}/motif/{case}_homer_results.txt", case=cases, output_dir=output_dir),
        expand("/data/wade/linl7/{identifier}/mm10/trackDb.txt",identifier=identifier),
        "{output_dir}/multiqc/multiqc_raw_fastq_report.html".format(output_dir=output_dir),
        "{output_dir}/multiqc/multiqc_trimmed_fastq_report.html".format(output_dir=output_dir),
        expand("{output_dir}/peaks/{case}_peaks.narrowPeak",case=cases,output_dir=output_dir),
        expand("{output_dir}/peaks/{case}_peaks.narrowPeak.clean.annotated",case=cases,output_dir=output_dir),
        expand("{output_dir}/aligned/{sample}_fragment_lengths.txt",sample=sample_names,output_dir=output_dir),
        expand("{output_dir}/ucsc/{sample}_normalized.bw", sample=sample_names, output_dir=output_dir),
        "{output_dir}/trimmed/target.barcodes.heatmap.png".format(output_dir=output_dir),
        expand("{output_dir}/metaplot/normalized/{comparisons}.xlsx",comparisons=comparisons,output_dir=output_dir),
        # expand("{output_dir}/metaplot/unnormalized/{comparisons}.xlsx",comparisons=comparisons,output_dir=output_dir),
        "{output_dir}/report/stat.tsv".format(output_dir=output_dir),
        "{output_dir}/aligned/scaleFactor.overview.txt".format(output_dir=output_dir)

        # expand("{output_dir}/ucsc/{sample}_unnormalized.bw", sample=sample_names, output_dir=output_dir),
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
        "ln -s {input.fastq1} {output.fastq1} && ln -s {input.fastq2} {output.fastq2} && {scripts_folder}/cut_and_run/softwares/fastqc {output.fastq1} {output.fastq2} -o {output_dir}/fastqc_raw/"
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
        "{scripts_folder}/cut_and_run/softwares/multiqc {output_dir}/fastqc_raw/ -o {output_dir}/multiqc/ -n multiqc_raw_fastq_report"

### --- Cutadapt for Adapter & Quality Trimming --- ###
rule cutadapt:
    input:
        fastq1=lambda wildcards: all_samples[wildcards.sample]["fastq1"],
        fastq2=lambda wildcards: all_samples[wildcards.sample]["fastq2"]
    output:
        trimmed_fastq1="{output_dir}/trimmed/{sample}_R1_trimmed.fastq.gz",
        trimmed_fastq2="{output_dir}/trimmed/{sample}_R2_trimmed.fastq.gz"
    shell:
        "{scripts_folder}/cut_and_run/softwares/cutadapt -q 20 -m 20 -a AGATCGGAAGAGC -A AGATCGGAAGAGC "
        "-o {output.trimmed_fastq1} -p {output.trimmed_fastq2} {input.fastq1} {input.fastq2}"

### --- FastQC After Trimming --- ###
rule fastqc_trimmed:
    input:
        trimmed_fastq1="{output_dir}/trimmed/{sample}_R1_trimmed.fastq.gz",
        trimmed_fastq2="{output_dir}/trimmed/{sample}_R2_trimmed.fastq.gz"
    output:
        "{output_dir}/fastqc_trimmed/{sample}_R1_trimmed_fastqc.html"
    shell:
        "{scripts_folder}/cut_and_run/softwares/fastqc {input.trimmed_fastq1} {input.trimmed_fastq2} -o {output_dir}/fastqc_trimmed/"

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
        "{scripts_folder}/cut_and_run/softwares/multiqc {output_dir}/fastqc_trimmed/ -o {output_dir}/multiqc -n multiqc_trimmed_fastq_report"


### --- Continue with Alignment and Analysis --- ###
rule bowtie2:
    input:
        trimmed_fastq1="{output_dir}/trimmed/{sample}_R1_trimmed.fastq.gz",
        trimmed_fastq2="{output_dir}/trimmed/{sample}_R2_trimmed.fastq.gz"
    output:
        bam = "{output_dir}/aligned/{sample}.bam",
        mappedBam = "{output_dir}/aligned/{sample}_mapped.bam",
        sortbam =  "{output_dir}/aligned/{sample}_dedup.bam",
        noBlacklistBam = "{output_dir}/aligned/{sample}_noblacklist.bam",
        dup_mets = "{output_dir}/aligned/{sample}_dup_metrics.txt",
        read_count = "{output_dir}/aligned/{sample}_uniq_readCount.txt"
    params:
        rg_id = "{sample}",
        rg_lb = "Unknown",
        rg_pl = "ILLUMINA",
        rg_pu = "Unknown",
        rg_sm = "{sample}"
    shell:
        # "bowtie2 -p 8 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -x {config[bowtie2_index]} -1 {input.trimmed_fastq1} -2 {input.trimmed_fastq2} --rg-id {params.rg_id} --rg LB:{params.rg_lb} --rg PL:{params.rg_pl} --rg PU:{params.rg_pu} --rg SM:{params.rg_sm}"
        "{scripts_folder}/cut_and_run/softwares/bowtie2 -p 8 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -x {config[bowtie2_index]} -1 {input.trimmed_fastq1} -2 {input.trimmed_fastq2} --rg-id {params.rg_id} --rg LB:{params.rg_lb} --rg PL:{params.rg_pl} --rg PU:{params.rg_pu} --rg SM:{params.rg_sm}"
        "| {scripts_folder}/cut_and_run/softwares/samtools view -@ 8 -Sb - | {scripts_folder}/cut_and_run/softwares/samtools sort -@ 8 -o {output.bam}"
        " && {scripts_folder}/cut_and_run/softwares/samtools view -@ 8 -b -F 4 {output.bam} > {output.mappedBam}"
        " && {scripts_folder}/cut_and_run/softwares/bedtools intersect -v -abam {output.mappedBam} -b {blacklist} > {output.noBlacklistBam}"
        " && {scripts_folder}/cut_and_run/softwares/java -Xmx90g -jar {scripts_folder}/cut_and_run/softwares/picard.jar MarkDuplicates I={output.noBlacklistBam} O={output.sortbam} M={output.dup_mets} REMOVE_DUPLICATES=false"
        " && {scripts_folder}/cut_and_run/softwares/samtools view -@ 8 -c {output.sortbam} > {output.read_count}"

rule ecoli_bowtie2:
    input:
        trimmed_fastq1="{output_dir}/trimmed/{sample}_R1_trimmed.fastq.gz",
        trimmed_fastq2="{output_dir}/trimmed/{sample}_R2_trimmed.fastq.gz"
    output:
        bam = "{output_dir}/aligned/{sample}.ecoli.bam",
        mappedBam = "{output_dir}/aligned/{sample}_mapped_ecoli.bam",
        sortbam =  "{output_dir}/aligned/{sample}_dedup_ecoli.bam",
        dup_mets = "{output_dir}/aligned/{sample}_dup_metrics_ecoli.txt",
        read_count = "{output_dir}/aligned/{sample}_uniq_readCount_ecoli.txt"
    params:
        rg_id = "{sample}",
        rg_lb = "Unknown",
        rg_pl = "ILLUMINA",
        rg_pu = "Unknown",
        rg_sm = "{sample}"
    shell:
        "{scripts_folder}/cut_and_run/softwares/bowtie2 -p 8 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -x {config[ecoli_bowtie2_index]} -1 {input.trimmed_fastq1} -2 {input.trimmed_fastq2} --rg-id {params.rg_id} --rg LB:{params.rg_lb} --rg PL:{params.rg_pl} --rg PU:{params.rg_pu} --rg SM:{params.rg_sm}"
        "| {scripts_folder}/cut_and_run/softwares/samtools view -@ 8 -Sb - | {scripts_folder}/cut_and_run/softwares/samtools sort -@ 8 -o {output.bam}"
        " && {scripts_folder}/cut_and_run/softwares/samtools view -@ 8 -b -F 4 {output.bam} > {output.mappedBam}"
        " && {scripts_folder}/cut_and_run/softwares/java -Xmx90g -jar {scripts_folder}/cut_and_run/softwares/picard.jar MarkDuplicates I={output.mappedBam} O={output.sortbam} M={output.dup_mets} REMOVE_DUPLICATES=true"
        " && {scripts_folder}/cut_and_run/softwares/samtools view -@ 8 -c {output.sortbam} > {output.read_count}"

rule plot_fragment:
    input:
        "{output_dir}/aligned/{sample}_dedup.bam"
    output:
       fragment_file = "{output_dir}/aligned/{sample}_fragment_lengths.txt",
       fragment_plot =  "{output_dir}/aligned/{sample}_fragment.hist.png"
    shell:
        "{scripts_folder}/cut_and_run/softwares/samtools view {input} | awk '$9 > 0 {{print $9}}' > {output.fragment_file} && Rscript {scripts_folder}/R/plot_fragment_length.R  {output.fragment_file}  {output.fragment_plot} {wildcards.sample}"

rule get_barcodes_count:
    input:
        fastq1 = "{output_dir}/trimmed/{sample}_R1_trimmed.fastq.gz",
        fastq2 = "{output_dir}/trimmed/{sample}_R2_trimmed.fastq.gz",
        read_count = "{output_dir}/aligned/{sample}_uniq_readCount.txt"
    output:
        "{output_dir}/trimmed/{sample}.scalor.txt"
    shell:
        "python {scripts_folder}/cut_and_run/epicypher.scalor.py {input.fastq1} {input.fastq2} {input.read_count} {output} {wildcards.sample}"

rule gather_barcodes_count:
    input:
        expand("{output_dir}/trimmed/{sample}.scalor.txt", sample=sample_names, output_dir=output_dir)
    output:
        "{output_dir}/trimmed/target.barcodes.table.txt"
    script:
        "{scripts_folder}/cut_and_run/snakemake/gather_barcodes.py"

rule plot_barcodes_heatmap:
    input:
        "{output_dir}/trimmed/target.barcodes.table.txt"
    output:
        "{output_dir}/trimmed/target.barcodes.heatmap.png"
    script:
        "{scripts_folder}/cut_and_run/snakemake/plot_barcodes_heatmap.py"


rule get_scalefactor:
    input:
        mm = "{output_dir}/aligned/{sample}_uniq_readCount.txt",
        ecoli = "{output_dir}/aligned/{sample}_uniq_readCount_ecoli.txt"
    output:
        "{output_dir}/aligned/{sample}_scaleFactor_raw.txt"
    run:
        mm_count = int(open(input.mm).read().strip())
        ecoli_count = int(open(input.ecoli).read().strip())
        
        spike_fraction = ecoli_count / mm_count
        raw_factor = 1 / spike_fraction

        with open(output[0], "w") as out:
            out.write(str(raw_factor) + "\n")

rule normalize_scalefactors:
    input:
        expand("{output_dir}/aligned/{sample}_scaleFactor_raw.txt",
               sample=sample_names, output_dir=output_dir)
    output:
        expand("{output_dir}/aligned/{sample}_scaleFactor.txt",
               sample=sample_names, output_dir=output_dir)
    run:
        # 读取所有 raw_factor
        raw = {wildcard: float(open(f).read())
               for wildcard, f in zip(sample_names, input)}

        max_raw = max(raw.values())

        for s in sample_names:
            final = raw[s] / max_raw
            with open(f"{output_dir}/aligned/{s}_scaleFactor.txt", "w") as out:
                out.write(str(final) + "\n")

rule scale_factor_overview:
    input:
        scale_factor = expand("{output_dir}/aligned/{sample}_scaleFactor.txt", sample=sample_names, output_dir=output_dir)
    output:
        "{output_dir}/aligned/scaleFactor.overview.txt"
    run:
        with open(output[0], "w") as out:
            out.write('sample\tuniq_read_count\tuniq_spikeIn_count\traw_factor\tscale_factor\n')
            for s in sample_names:
                read_count = int(open(f"{output_dir}/aligned/{s}_uniq_readCount.txt").read())
                ecoli_count = int(open(f"{output_dir}/aligned/{s}_uniq_readCount_ecoli.txt").read())
                raw_factor = float(open(f"{output_dir}/aligned/{s}_scaleFactor_raw.txt").read())
                scale_factor = float(open(f"{output_dir}/aligned/{s}_scaleFactor.txt").read())
                out.write(s + '\t' + str(read_count) + '\t' + str(ecoli_count) + '\t' + str(raw_factor) + '\t' + str(scale_factor) + '\n')


rule ucsc_bam2bigwig:
    input:
        bam = "{output_dir}/aligned/{sample}_dedup.bam",
        scale_factor="{output_dir}/aligned/{sample}_scaleFactor.txt"
    output:
        normalized="{output_dir}/ucsc/{sample}_normalized.bw",
        unnormalized="{output_dir}/ucsc/{sample}_unnormalized.bw"
    params:
        scale_factor=lambda wildcards, input: open(input.scale_factor).readlines()[0].strip().split('\t')[0]
    shell:
        "{scripts_folder}/cut_and_run/softwares/samtools index {input.bam} && {scripts_folder}/cut_and_run/softwares/bamCoverage -b {input.bam} -o {output.normalized} --scaleFactor {params.scale_factor} &&  {scripts_folder}/cut_and_run/softwares/bamCoverage -b {input.bam} -o {output.unnormalized}"

rule macs2_with_control:
    input:
        case = lambda wildcards: "{output_dir}/aligned/{case}_dedup.bam".format(output_dir=output_dir,case=wildcards.case),
        control = lambda wildcards: "{output_dir}/aligned/{control}_dedup.bam".format(output_dir=output_dir, control=config["cases"][wildcards.case]["control"])
    output:
        "{output_dir}/peaks/{case}_peaks.narrowPeak"
    shell:
        "{scripts_folder}/cut_and_run/softwares/macs2 callpeak -t {input.case} -c {input.control} -f BAMPE -g {config[macs2_genome_size]} --keep-dup all"
        " --outdir {output_dir}/peaks/ -n {wildcards.case} "


rule macs2_no_control:
    input:
        bam = lambda wc: f"{output_dir}/aligned/{wc.case}_dedup.bam"
    output:
        "{output_dir}/peaks/{case}_peaks.narrowPeak"
    shell:
        "EXTSIZE=$({scripts_folder}/cut_and_run/softwares/samtools stats -@ 8 {input.bam} | grep '^IS' | awk '{{sum+=$2*$3; count+=$3}} END {{print int(sum/count)}}')"
        " && {scripts_folder}/cut_and_run/softwares/macs2 callpeak -t {input.bam} -f BAMPE -g {config[macs2_genome_size]} --keep-dup all "
        " --outdir {output_dir}/peaks/ -n {wildcards.case} --pvalue 1e-5 --nomodel --extsize $EXTSIZE"



rule clean_and_annotate_peak:
    input:
        "{output_dir}/peaks/{case}_peaks.narrowPeak"
    output:
        filt = "{output_dir}/peaks/{case}_peaks.narrowPeak.filtered",
        clean = "{output_dir}/peaks/{case}_peaks.narrowPeak.clean",
        annotated = "{output_dir}/peaks/{case}_peaks.narrowPeak.clean.annotated"
    shell:
        "awk '$9 >= 4' {input}  > {output.filt} && {scripts_folder}/cut_and_run/softwares/bedtools subtract -a {output.filt} -b {blacklist} > {output.clean} && "
        "perl {scripts_folder}/cut_and_run/softwares/annotatePeaks.pl {output.clean} {ref_to_annotate} > {output.annotated}"



rule uscs_peak2bed:
    input:
        "{output_dir}/peaks/{case}_peaks.narrowPeak"
    output:
        "{output_dir}/ucsc/{case}.peaks.bb"
    shell:
        "Rscript {scripts_folder}/../CHIPseq/peak_to_bigbed.R {input} {output} {fai}"


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
        peaks=expand("{output_dir}/peaks/{sample}_peaks.narrowPeak",output_dir=output_dir,sample=cases)
    output:
        "{output_dir}/report/stat.tsv"
    params:
        samples=expand(sorted(sample_names)) # {params} will show [] in shell
    shell:
        "python {scripts_folder}/python/pipeline_report_stat.py --fastq {input.fqs} --trimmed_fastq {input.trim_fqs} --fastq_report --metrics {input.dup_mets} --peaks {input.peaks} --samples {params.samples} --output {output}"

rule ucsc_hub:
    input:
        bw=expand("{output_dir}/ucsc/{sample}_normalized.bw", sample=sample_names,output_dir=output_dir),
        peaks=expand("{output_dir}/ucsc/{case}.peaks.bb", case=cases, output_dir=output_dir, allow_missing=True)
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
        cp {input.bw} {output.hub_dir}
        for peak in {input.peaks}; do
            if [ -f $peak ]; then
                cp $peak {output.hub_dir}/
            fi
        done

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
            if [ -f "{output_dir}/ucsc/$case.peaks.bb" ]; then
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
            fi
        done

        """


rule mataplot:
    input: 
        peaks=lambda wildcards: expand("{output_dir}/peaks/{case}_peaks.narrowPeak",output_dir=output_dir,case=config["comparisons"][wildcards.comparisons]),
        bws=lambda wildcards: expand("{output_dir}/ucsc/{case}_normalized.bw",output_dir=output_dir,case=config["comparisons"][wildcards.comparisons]),
        un_bws=lambda wildcards: expand("{output_dir}/ucsc/{case}_unnormalized.bw",output_dir=output_dir,case=config["comparisons"][wildcards.comparisons]) 
    output:
        xlsx="{output_dir}/metaplot/normalized/{comparisons}.xlsx",
        merged_bed="{output_dir}/metaplot/normalized/{comparisons}.merged.bed",
        matrix="{output_dir}/metaplot/normalized/{comparisons}.matrix.gz",
        tsv="{output_dir}/metaplot/normalized/{comparisons}.matrix.tsv",
        un_xlsx="{output_dir}/metaplot/unnormalized/{comparisons}.xlsx",
        un_matrix="{output_dir}/metaplot/unnormalized/{comparisons}.matrix.gz",
        un_tsv="{output_dir}/metaplot/unnormalized/{comparisons}.matrix.tsv"
    params:
        n=lambda wildcards: len(config["comparisons"][wildcards.comparisons])
    shell:
        "cat {input.peaks} | sort -k1,1 -k2,2n | bedtools merge -d 100 -i - > {output.merged_bed} && "
        "{scripts_folder}/cut_and_run/softwares/computeMatrix reference-point -p 8 --referencePoint center -S {input.bws} -R {output.merged_bed} "
        "--beforeRegionStartLength 2000 --afterRegionStartLength 2000 --binSize 40 -o {output.matrix} --outFileNameMatrix {output.tsv} && "
        "python {scripts_folder}/python/metaplot.py {output.tsv} {output.xlsx} {params.n} && "
        "{scripts_folder}/cut_and_run/softwares/computeMatrix reference-point -p 8 --referencePoint center -S {input.un_bws} -R {output.merged_bed} "
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
# #         "plotProfile -m {input} -out {output}"


