
import re


#configfile: "config.yaml"

cases = list(config["cases"].keys())
comparisons = list(config["comparisons"].keys())
output_dir = config["output_dir"]
identifier = config['identifier']
SPIKEIN = config.get("spike_in", False)
genome = config['genome']
gtf = config['gtf']
bit = config['bit']
tx2gene =  config['tx2gene']
transcriptome = config['transcriptome']
fragment_length = config['fragment_length']
effectGenomeSize = config['effectGenomeSize']
scripts_folder = config['scripts_folder']


if "controls" in config:
    controls = list(config["controls"].keys())
    all_samples = config["cases"] | config["controls"]
else:
    all_samples = config["cases"]

sample_names = list(all_samples.keys())


rule all:
    input:
        # expand("{output_dir}/fastqc_raw/{sample}_R1_fastqc.html", sample=sample_names, output_dir=output_dir),
        # # expand("{output_dir}/fastqc_trimmed/{sample}_trimmed_fastqc.html", sample=sample_names, output_dir=output_dir),
        # # expand("{output_dir}/trimmed/{sample}_R1_trimmed.fastq.gz", sample=sample_names, output_dir=output_dir),
        "{output_dir}/multiqc/multiqc_raw_fastq_report.html".format(output_dir=output_dir),
        "{output_dir}/multiqc/multiqc_trimmed_fastq_report.html".format(output_dir=output_dir),
        # expand("{output_dir}/align/{sample}.Aligned.out.bam", sample=sample_names, output_dir=output_dir),
        # expand("{output_dir}/align/{sample}_fragment_lengths.txt",sample=sample_names,output_dir=output_dir,
        "/data/wade/linl7/{identifier}_BPM_normalised/hub.txt".format(identifier=identifier),
        expand("{output_dir}/DEG_by_DESeq2_no_spike/{comparisons}/sample.list",comparisons=comparisons,output_dir=output_dir),
        "{output_dir}/report/stat.tsv".format(output_dir=output_dir),
        ("{output_dir}/align/scaleFactor.overview.txt".format(output_dir=output_dir) if SPIKEIN else []),
        (expand("{output_dir}/venn/{comparisons}.venn.png",comparisons=comparisons,output_dir=output_dir) if SPIKEIN else [])


### --- FastQC Before Trimming --- ###
rule fastqc_raw:
    input:
        fastq1=lambda wildcards: all_samples[wildcards.sample]["fastq1"]
    output:
        ####### this makes it only accept fastq.gz in config file. ##########
        fastq1="{output_dir}/input/{sample}_R1.fastq.gz",
        report="{output_dir}/fastqc_raw/{sample}_R1_fastqc.html"
    shell:
        "ln -s {input.fastq1} {output.fastq1} && {scripts_folder}/RNAseq/softwares/fastqc {output.fastq1} -o {output_dir}/fastqc_raw/ -t 8"
rule report_fastq:
    input:
        fastq1="{output_dir}/input/{sample}_R1.fastq.gz"
    output:
        "{output_dir}/input/{sample}.fastq.count.txt"
    shell:
        "source {scripts_folder}/shell/fastq.sh && count_reads {input.fastq1} > {output}"

rule multiqc_raw:
    input:
        expand("{output_dir}/fastqc_raw/{sample}_R1_fastqc.html", sample=sample_names, output_dir=output_dir),
    output:
        "{output_dir}/multiqc/multiqc_raw_fastq_report.html"
    shell:
        "{scripts_folder}/RNAseq/softwares/multiqc {output_dir}/fastqc_raw/ -o {output_dir}/multiqc/ -n multiqc_raw_fastq_report"

### --- Cutadapt for Adapter & Quality Trimming --- ###
rule cutadapt:
    input:
        fastq1=lambda wildcards: all_samples[wildcards.sample]["fastq1"]
    output:
        trimmed_fastq1="{output_dir}/trimmed/{sample}_R1_trimmed.fastq.gz"
    shell:
        "{scripts_folder}/RNAseq/softwares/cutadapt -q 20 -m 20 -a AGATCGGAAGAGC "
        "-o {output.trimmed_fastq1} {input.fastq1} -j 8"

### --- FastQC After Trimming --- ###
rule fastqc_trimmed:
    input:
        trimmed_fastq1="{output_dir}/trimmed/{sample}_R1_trimmed.fastq.gz"
    output:
        "{output_dir}/fastqc_trimmed/{sample}_R1_trimmed_fastqc.html"
    shell:
        "{scripts_folder}/RNAseq/softwares/fastqc {input.trimmed_fastq1} -o {output_dir}/fastqc_trimmed/ -t 8"

rule report_trim_fastq:
    input:
        trimmed_fastq1="{output_dir}/trimmed/{sample}_R1_trimmed.fastq.gz"
    output:
        "{output_dir}/trimmed/{sample}.trimmed.fastq.count.txt"
    shell:
        "source {scripts_folder}/shell/fastq.sh && count_reads {input.trimmed_fastq1} > {output}"

rule STAR:
    input:
        trimmed_fastq1="{output_dir}/trimmed/{sample}_R1_trimmed.fastq.gz"
    output:
        bam = "{output_dir}/align/{sample}.Aligned.out.bam",
        transcriptome_bam = "{output_dir}/align/{sample}.Aligned.toTranscriptome.out.bam",
        sortbam =  "{output_dir}/align/{sample}.Aligned.out.sorted.bam",
        flagstat = "{output_dir}/align/{sample}.Aligned.out.sorted.flagstat.txt",
        dup_mets = "{output_dir}/align/{sample}_dup_metrics.txt",
        read_count = "{output_dir}/align/{sample}_uniq_readCount.txt"
    params:
        out_prefix = "{output_dir}/align/{sample}."
    shell:
        "{scripts_folder}/RNAseq/softwares/STAR --runThreadN 8 --genomeDir {genome} --readFilesIn {input.trimmed_fastq1}"
        " --outFileNamePrefix {params.out_prefix} --readFilesCommand zcat --outSAMtype BAM Unsorted --quantTranscriptomeBan Singleend --outFilterType BySJout"
        " --alignSJoverhangMin 8 --outFilterMultimapNmax 20"
        " --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999"
        " --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20"
        " --alignIntronMax 1000000 --alignMatesGapMax 1000000"
        " --quantMode TranscriptomeSAM --outSAMattributes NH HI AS NM MD"
        " && {scripts_folder}/RNAseq/softwares/samtools sort -@ 8 -o {output.sortbam} {output.bam} && {scripts_folder}/RNAseq/softwares/samtools index -@ 8 {output.sortbam}"
        " && {scripts_folder}/RNAseq/softwares/samtools flagstat -@ 8 {output.sortbam} > {output.flagstat}"
        " && {scripts_folder}/RNAseq/softwares/samtools view -@ 8 -F 0x100 -c {output.sortbam} > {output.read_count}"
        " && {scripts_folder}/RNAseq/softwares/java -Xmx90g -jar {scripts_folder}/RNAseq/softwares/picard.jar MarkDuplicates I={output.sortbam} O=/dev/null M={output.dup_mets} REMOVE_DUPLICATES=false CREATE_INDEX=false"

rule qualimap:
    input: 
        sortbam =  "{output_dir}/align/{sample}.Aligned.out.sorted.bam"
    output:
        qc = "{output_dir}/align/{sample}.rnaseq.qualimap.report/rnaseq_qc_results.txt"
    params:
        bamqc_output_prefix = "{output_dir}/align/{sample}.bamqc.qualimap.report/",
        rnaseq_output_prefix = "{output_dir}/align/{sample}.rnaseq.qualimap.report/"
    shell:
        "{scripts_folder}/RNAseq/softwares/qualimap bamqc -nt 8 -bam {input.sortbam} -gff {gtf} -outdir {params.bamqc_output_prefix} --java-mem-size=90G"
        " && {scripts_folder}/RNAseq/softwares/qualimap rnaseq -bam {input.sortbam} -gtf {gtf} -outdir {params.rnaseq_output_prefix} --java-mem-size=90G"

rule multiqc_trimmed:
    input:
        trim = expand("{output_dir}/fastqc_trimmed/{sample}_R1_trimmed_fastqc.html", sample=sample_names, output_dir=output_dir),
        qc =  expand("{output_dir}/align/{sample}.rnaseq.qualimap.report/rnaseq_qc_results.txt", sample=sample_names, output_dir=output_dir)
    output:
        "{output_dir}/multiqc/multiqc_trimmed_fastq_report.html"
    shell:
        "{scripts_folder}/RNAseq/softwares/multiqc {output_dir}/fastqc_trimmed/ {output_dir}/align/ -o {output_dir}/multiqc -n multiqc_trimmed_fastq_report"

rule GCbias:
    input:
        sortbam =  "{output_dir}/align/{sample}.Aligned.out.sorted.bam"
    output:
        GC_freq = "{output_dir}/align/{sample}.GC.freq.txt",
        biasplot = "{output_dir}/align/{sample}.GC.biasPlot.pdf"
    shell:
        "{scripts_folder}/RNAseq/softwares/computeGCBias -p 8 -b {input.sortbam} -l {fragment_length} --effectiveGenomeSize {effectGenomeSize} -g {bit} --GCbiasFrequenciesFile {output.GC_freq}  --biasPlot  {output.biasplot}"
        # " && correctGCBias -p 8 -b {input.sortbam} --effectiveGenomeSize {effectGenomeSize}  -g {bit} --GCbiasFrequenciesFile {output.GC_freq} -o $output_dir/align/${sample_name}.gc_corrected.bam"

rule plot_fragment:
    input:
        "{output_dir}/align/{sample}.Aligned.out.sorted.bam"
    output:
       fragment_file = "{output_dir}/align/{sample}_fragment_lengths.txt",
       fragment_plot =  "{output_dir}/align/{sample}_fragment.hist.png"
    shell:
        "{scripts_folder}/RNAseq/softwares/samtools view -@ 8 {input} | awk '$9 > 0 {{print $9}}' > {output.fragment_file} && Rscript {scripts_folder}/R/plot_fragment_length.R  {output.fragment_file}  {output.fragment_plot} {wildcards.sample}"


rule ercc_count:
    input:
        "{output_dir}/align/{sample}.Aligned.out.sorted.bam"
    output:
        read_count = "{output_dir}/align/{sample}_uniq_readCount_ERCC.txt"
    shell:
       # "samtools idxstats {input} | awk '$1 ~ /ERCC/ {{sum += $3}} END {{print sum}}' > {output.read_count}"
       "{scripts_folder}/RNAseq/softwares/samtools view -@ 8 -F 0x100  {input} | awk '$3~/ERCC/' |wc -l > {output.read_count}"

rule get_scalefactor:
    input:
        mm = "{output_dir}/align/{sample}_uniq_readCount.txt",
        ecoli = "{output_dir}/align/{sample}_uniq_readCount_ERCC.txt"
    output:
        "{output_dir}/align/{sample}_scaleFactor_raw.txt"
    run:
        mm_count = int(open(input.mm).read().strip())
        ecoli_count = int(open(input.ecoli).read().strip())
        
        spike_fraction = ecoli_count / mm_count
        raw_factor = 1 / spike_fraction

        with open(output[0], "w") as out:
            out.write(str(raw_factor) + "\n")

rule normalize_scalefactors:
    input:
        expand("{output_dir}/align/{sample}_scaleFactor_raw.txt",
               sample=sample_names, output_dir=output_dir)
    output:
        expand("{output_dir}/align/{sample}_scaleFactor.txt",
               sample=sample_names, output_dir=output_dir)
    run:
        # 读取所有 raw_factor
        raw = {wildcard: float(open(f).read())
               for wildcard, f in zip(sample_names, input)}

        max_raw = max(raw.values())

        for s in sample_names:
            final = raw[s] / max_raw
            with open(f"{output_dir}/align/{s}_scaleFactor.txt", "w") as out:
                out.write(str(final) + "\n")

rule scale_factor_overview:
    input:
        scale_factor = expand("{output_dir}/align/{sample}_scaleFactor.txt", sample=sample_names, output_dir=output_dir)
    output:
        "{output_dir}/align/scaleFactor.overview.txt"
    run:
        with open(output[0], "w") as out:
            out.write('sample\tuniq_read_count\tuniq_spikeIn_count\traw_factor\tscale_factor\n')
            for s in sample_names:
                read_count = int(open(f"{output_dir}/align/{s}_uniq_readCount.txt").read())
                ecoli_count = int(open(f"{output_dir}/align/{s}_uniq_readCount_ERCC.txt").read())
                raw_factor = float(open(f"{output_dir}/align/{s}_scaleFactor_raw.txt").read())
                scale_factor = float(open(f"{output_dir}/align/{s}_scaleFactor.txt").read())
                out.write(s + '\t' + str(read_count) + '\t' + str(ecoli_count) + '\t' + str(raw_factor) + '\t' + str(scale_factor) + '\n')



rule ucsc_bam2bigwig:
    input:
        bam = "{output_dir}/align/{sample}.Aligned.out.sorted.bam",
        scale_factor = lambda wc: f"{wc.output_dir}/align/{wc.sample}_scaleFactor.txt" if SPIKEIN else []
    output:
        unnormalized = "{output_dir}/ucsc/{sample}_unnormalized.bw",
        BPMnormalized = "{output_dir}/ucsc/{sample}_BPM_normalized.bw",
        SpikeINnormalized = "{output_dir}/ucsc/{sample}_spikeIn_normalized.bw" if SPIKEIN else []
    params:
        scale_factor = lambda wc, input: open(input.scale_factor).readline().strip().split("\t")[0] if SPIKEIN else ''
    shell:
        r"""
        set -euo pipefail

        {scripts_folder}/RNAseq/softwares/bamCoverage -p 8 -b {input.bam} \
            -o {output.unnormalized} \
            --samFlagExclude 1804

        {scripts_folder}/RNAseq/softwares/bamCoverage -p 8 -b {input.bam} \
            -o {output.BPMnormalized} \
            --normalizeUsing BPM \
            --samFlagExclude 1804

        if [ "{SPIKEIN}" = "True" ]; then
            {scripts_folder}/RNAseq/softwares/bamCoverage -p 8 -b {input.bam} \
                -o {output.SpikeINnormalized} \
                --scaleFactor {params.scale_factor} \
                --normalizeUsing None \
                --samFlagExclude 1804
        fi
        """

rule ucsc_hub:
    input:
        unnormalized = expand("{output_dir}/ucsc/{sample}_unnormalized.bw",sample=sample_names,output_dir=output_dir),
        bpm_normalized = expand("{output_dir}/ucsc/{sample}_BPM_normalized.bw",sample=sample_names,output_dir=output_dir),
        spikeIn_normalized = expand("{output_dir}/ucsc/{sample}_spikeIn_normalized.bw",sample=sample_names,output_dir=output_dir) if SPIKEIN else []
    output:
        un_hub = "/data/wade/linl7/" + identifier + "_unnormalised/hub.txt",
        bpm_hub = "/data/wade/linl7/" + identifier + "_BPM_normalised/hub.txt",
        spike_hub = ("/data/wade/linl7/" + identifier + "_spikeIn_normalised/hub.txt" if SPIKEIN else [])
    params:
        unnormalized_hub = identifier + "_unnormalised",
        bpm_normalized_hub = identifier + "_BPM_normalised",
        spikeIn_normalized_hub = identifier + "_spikeIn_normalised" if SPIKEIN else ''
    shell:
        "mkdir -p {output_dir}/bws_unnormalised  &&"
        " for bw in {input.unnormalized}; do ln -s $bw {output_dir}/bws_unnormalised/; done &&"
        " bash {scripts_folder}/shell/make_ucsc_hub.sh {output_dir}/bws_unnormalised {params.unnormalized_hub} &&"
        " mkdir -p {output_dir}/bws_BPM_normalised  &&"
        " for bw in {input.bpm_normalized}; do ln -s $bw {output_dir}/bws_BPM_normalised/; done &&"
        " bash {scripts_folder}/shell/make_ucsc_hub.sh {output_dir}/bws_BPM_normalised {params.bpm_normalized_hub} &&"
        " if [ \"{SPIKEIN}\" = \"True\" ]; then "
        " mkdir -p {output_dir}/bws_spikeIn_normalised &&"
        " for bw in {input.spikeIn_normalized}; do ln -s $bw {output_dir}/bws_spikeIn_normalised/; done &&"
        " bash {scripts_folder}/shell/make_ucsc_hub.sh {output_dir}/bws_spikeIn_normalised {params.spikeIn_normalized_hub}; "
        " fi"

rule salmon:
    input:
        bam = "{output_dir}/align/{sample}.Aligned.toTranscriptome.out.bam"
    output:
        salmon_sf = "{output_dir}/align/{sample}.Aligned.toTranscriptome.salmon_quant/quant.sf"
    params:
        out_prefix = "{output_dir}/align/{sample}.Aligned.toTranscriptome.salmon_quant"
    shell:
        "{scripts_folder}/RNAseq/softwares/salmon quant -p 8 -t {transcriptome} --fldMean 300 --libType A -a {input.bam} -o {params.out_prefix}  --gcBias --seqBias"



rule report_stat:
    input:  
        fqs=expand("{output_dir}/input/{sample}.fastq.count.txt", output_dir=output_dir,sample=sample_names),
        trim_fqs=expand("{output_dir}/trimmed/{sample}.trimmed.fastq.count.txt",output_dir=output_dir,sample=sample_names),
        dup_mets=expand("{output_dir}/align/{sample}_dup_metrics.txt",output_dir=output_dir,sample=sample_names)
        # peaks=expand("{output_dir}/peaks/{sample}_peaks.narrowPeak",output_dir=output_dir,sample=cases)
    output:
        "{output_dir}/report/stat.tsv"
    params:
        samples=expand(sorted(sample_names)) # {params} will show [] in shell
    shell:
        "python {scripts_folder}/python/pipeline_report_stat.py --fastq {input.fqs} --trimmed_fastq {input.trim_fqs} --fastq_report --metrics {input.dup_mets} --samples {params.samples} --output {output}"

def needs_rev_ref(comparison):
    A, B = comparison.split("_vs_")
    return A < B   # A alphabetically comes before B → need reverse

def get_condition(sample):
    return re.sub(r"-rep-\d+$", "", sample)



rule deseq2:
    input: 
        salmon_sf=lambda wildcards: expand("{output_dir}/align/{case}.Aligned.toTranscriptome.salmon_quant/quant.sf",output_dir=output_dir,case=config["comparisons"][wildcards.comparisons])
    output:
        sample_list = "{output_dir}/DEG_by_DESeq2_no_spike/{comparisons}/sample.list",
        degs_no_spike = "{output_dir}/DEG_by_DESeq2_no_spike/{comparisons}/{comparisons}.count10.groupsize3.deg.txt",
        degs_spike =  ("{output_dir}/DEG_by_DESeq2_with_spikeIn/{comparisons}/{comparisons}.count10.groupsize3.deg.txt" if SPIKEIN else [])
    params:
        no_spike_out_prefix = "{output_dir}/DEG_by_DESeq2_no_spike/{comparisons}/{comparisons}.count10.groupsize3",
        spike_out_prefix = "{output_dir}/DEG_by_DESeq2_with_spikeIn/{comparisons}/{comparisons}.count10.groupsize3",
        revref = lambda wc: "--rev_ref" if needs_rev_ref(wc.comparisons) else ""
    shell:
        """
        echo -e "samples\tpath\tcondition" > {output.sample_list}
        for sf in {input.salmon_sf}; do
            # sf: output_dir/align/sample.Aligned.toTranscriptome.salmon_quant/quant.sf
            sm=$(basename $(dirname ${{sf}}))
            sm=${{sm%.Aligned.toTranscriptome.salmon_quant}}

            cond=$(echo $sm | sed -E 's/-rep-[0-9]+$//')

            echo -e "$sm\t$sf\t$cond" >> {output.sample_list}
        done

        # 运行 DESeq2
        Rscript {scripts_folder}/RNAseq/DESeq2/run.deseq.R \
            --samples {output.sample_list} \
            --tx2gene {tx2gene} \
            --output_prefix {params.no_spike_out_prefix} \
            --count 10 --group_size 3 \
            {params.revref}

        if [ "{SPIKEIN}" = "True" ]; then  
            Rscript {scripts_folder}/RNAseq/DESeq2/run.deseq.R \
                --samples {output.sample_list} \
                --tx2gene {tx2gene} \
                --output_prefix {params.spike_out_prefix} \
                --count 10 --group_size 3 --ercc_norm \
                {params.revref}
        fi
        """

rule venn_plot:
    input:
        degs_no_spike = "{output_dir}/DEG_by_DESeq2_no_spike/{comparisons}/{comparisons}.count10.groupsize3.deg.txt",
        degs_spike = "{output_dir}/DEG_by_DESeq2_with_spikeIn/{comparisons}/{comparisons}.count10.groupsize3.deg.txt"
    output:
        "{output_dir}/venn/{comparisons}.venn.png"
    params:
        title = "{comparisons}"
    shell:
        "Rscript {scripts_folder}/RNAseq/DESeq2/venn_plot_deseq2_genes.R {input.degs_no_spike} {input.degs_spike} {output} no_spike with_spike {params.title}"
