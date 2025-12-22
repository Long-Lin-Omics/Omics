# following https://github.com/CebolaLab/RNA-seq?tab=readme-ov-file
fastq1=$1
genome=$2
gtf=$3
bit=$4
transcriptome=$5
output_dir=$6
sample_name=$7
fragment_length=$8
effectGenomeSize=$9

# fastq1=/ddn/gs1/project/nextgen/post/hug4/LongLin/rnaseq/test/WT-D2-1.single.fastq.gz
# genome=/ddn/gs1/project/nextgen/post/hug4/LongLin/rnaseq/reference/STAR_mm10_75bp_ucsc/
# gtf=/ddn/gs1/project/nextgen/post/hug4/LongLin/rnaseq/reference/mm10.ensGene.gtf
# bit=/ddn/gs1/project/nextgen/post/hug4/LongLin/rnaseq/reference/mm10.2bit
# transcriptome=/ddn/gs1/project/nextgen/post/hug4/LongLin/rnaseq/reference/mm10.transcripts.fa
# output_dir=./WT-D2
# sample_name=WT-D2




# echo $fastq1
# echo $genome
# echo $output_dir
# echo $sample_name

make_folder() {
  local folder="$1"
  if [ ! -d "$folder" ]; then
    mkdir -p "$folder"
  fi
}

make_folder $output_dir/data
make_folder $output_dir/fastqc_raw
make_folder $output_dir/trimmed
make_folder $output_dir/fastqc_trimmed
make_folder $output_dir/align



ln -s $fastq1 $output_dir/data/$sample_name.fastq.gz
fastqc $output_dir/data/$sample_name.fastq.gz -o $output_dir/fastqc_raw -t 8 
cutadapt -q 20 -m 20 -a AGATCGGAAGAGC -o $output_dir/trimmed/${sample_name}_single_trimmed.fastq.gz $output_dir/data/$sample_name.fastq.gz -j 8 
fastqc $output_dir/trimmed/${sample_name}_single_trimmed.fastq.gz -o $output_dir/fastqc_trimmed -t 8

STAR --runThreadN 8 --genomeDir $genome --readFilesIn $output_dir/trimmed/${sample_name}_single_trimmed.fastq.gz \
--outFileNamePrefix $output_dir/align/${sample_name}. --readFilesCommand zcat --outSAMtype BAM Unsorted --quantTranscriptomeBan Singleend --outFilterType BySJout \
--alignSJoverhangMin 8 --outFilterMultimapNmax 20 \
--alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 \
--alignIntronMax 1000000 --alignMatesGapMax 1000000 \
--quantMode TranscriptomeSAM --outSAMattributes NH HI AS NM MD \

# assess lane effects
# deeptools plotCorrelation

samtools sort --threads 8 -o $output_dir/align/${sample_name}.Aligned.out.sorted.bam $output_dir/align/${sample_name}.Aligned.out.bam  
samtools index --threads 8 $output_dir/align/${sample_name}.Aligned.out.sorted.bam
samtools flagstat --threads 8 $output_dir/align/${sample_name}.Aligned.out.sorted.bam > $output_dir/align/${sample_name}.Aligned.out.sorted.flagstat

qualimap bamqc -nt 8 -bam $output_dir/align/${sample_name}.Aligned.out.sorted.bam -gff $gtf -outdir $output_dir/align/${sample_name}.bamqc.qualimap.report --java-mem-size=90G
qualimap rnaseq -bam $output_dir/align/${sample_name}.Aligned.out.sorted.bam -gtf $gtf -outdir $output_dir/align/${sample_name}.rnaseq.qualimap.report --java-mem-size=90G
echo "aligned to exons%" 
cat $output_dir/align/${sample_name}.rnaseq.qualimap.report/rnaseq_qc_results.txt | grep exonic | cut -d '(' -f 2 | cut -d ')' -f1
echo "The total number of reads mapped"
cat $output_dir/align/${sample_name}.Aligned.out.sorted.flagstat | grep mapped | head -n1 | cut -d ' ' -f1
echo "The total number of properly paired reads"
cat $output_dir/align/${sample_name}.Aligned.out.sorted.flagstat | grep 'properly paired' | head -n1 | cut -d ' ' -f1


# qualimap multi-bamqc sample.txt # can do pca
                        # sample_name full_path_bamqc group

# multiqc to collect fastqc and qualimap

# not remove duplicates unless using UMIs

computeGCBias -p 8 -b $output_dir/align/${sample_name}.Aligned.out.sorted.bam -l $fragment_length --effectiveGenomeSize $effectGenomeSize -g $bit --GCbiasFrequenciesFile  $output_dir/align/${sample_name}.GC.freq.txt  --biasPlot  $output_dir/align/${sample_name}.GC.biasPlot.pdf
correctGCBias -p 8 -b $output_dir/align/${sample_name}.Aligned.out.sorted.bam --effectiveGenomeSize $effectGenomeSize  -g $bit --GCbiasFrequenciesFile $output_dir/align/${sample_name}.GC.freq.txt -o $output_dir/align/${sample_name}.gc_corrected.bam

# bamCoverage -p 8 -b $output_dir/align/${sample_name}.gc_corrected.bam -o $output_dir/align/${sample_name}.gc_corrected.bw --normalizeUsing BPM --samFlagExclude 512
bamCoverage -p 8 -b $output_dir/align/${sample_name}.Aligned.out.sorted.bam -o $output_dir/align/${sample_name}.BPM.normalised.bw --normalizeUsing BPM --samFlagExclude 1804
spikeIn=$(samtools idxstats $output_dir/align/${sample_name}.Aligned.out.sorted.bam | awk '$1 ~ /ERCC/ {sum += $3} END {print sum}')
scaleFactor=$(echo "scale=10; 1/$spikeIn" | bc)
bamCoverage -p 8 -b $output_dir/align/${sample_name}.Aligned.out.sorted.bam -o $output_dir/align/${sample_name}.spikeIn.scaled.bw --normalizeUsing None --samFlagExclude 1804 --scaleFactor $scaleFactor

#salmon can accept multiple bams as a list
#salmon has its own in-built method to correct for GC-bias
# --fldMean ?? --fldSD ?? 
# No --gcBias if using cqn to correct for sample-specific biases in DEseq2
salmon quant -p 8 -t $transcriptome --fldMean 300 --libType A -a $output_dir/align/${sample_name}.Aligned.toTranscriptome.out.bam -o $output_dir/align/${sample_name}.Aligned.toTranscriptome.salmon_quant  --gcBias --seqBias

