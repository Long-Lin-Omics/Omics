### Unzip reads ###
zcat $MATE1_FASTQ > out/$S/raw_reads.1.fastq
zcat $MATE2_FASTQ > out/$S/raw_reads.2.fastq
 
### Filter reads ###
/ddn/gs1/home/bennettb/bin/trim_and_filter_fastq --paired -i1 out/$S/raw_reads.1.fastq -i2 out/$S/raw_reads.2.fastq -q 20 -o1 out/$S/filtered_reads.1.fastq -o2 out/$S/filtered_reads.2.fastq -r out/$S/filtered_reads.FilterStats.txt
 
### Remove adapter ###
cutadapt --match-read-wildcards -m $MINIMUM_FRAGMENT_SIZE -a $ADAPTER_SEQUENCE_MATE1 -A $ADAPTER_SEQUENCE_MATE2 -o out/$S/trimmed_reads.1.fastq -p out/$S/trimmed_reads.2.fastq out/$S/filtered_reads.1.fastq out/$S/filtered_reads.2.fastq > out/$S/adapter_trimming_report.txt
 
### Alignment ###
bowtie2 -p $NUM_THREADS -I $MINIMUM_FRAGMENT_SIZE -X $MAXIMUM_FRAGMENT_SIZE --local -x $BOWTIE2_INDEX -1 out/$S/trimmed_reads.1.fastq -2 out/$S/trimmed_reads.2.fastq 2> out/$S/alignment_report.txt | samtools view -Sb -q5 -f2 - > out/$S/aligned_reads.bam
 
### Deduplicate ###
java -Xmx2g -jar /ddn/gs1/home/bennettb/tools/picard-tools-3.1.1/picard.jar MarkDuplicates --REMOVE_DUPLICATES TRUE --ASSUME_SORT_ORDER queryname --VALIDATION_STRINGENCY SILENT --VERBOSITY WARNING -I out/$S/aligned_reads.bam -O out/$S/aligned_reads_deduplicated.bam -M out/$S/deduplication_report.txt --TMP_DIR out/$S/temp; rm -rf out/$S/temp
 
### Extract fragments ###
bedtools bamtobed -bedpe -i out/$S/aligned_reads_deduplicated.bam > out/$S/hits_extended.bedpe
cut -f1,2,6 out/$S/hits_extended.bedpe > out/$S/hits_extended.bed
 
### Make fragment length plot ###
perl calculate_fragment_lengths.pl out/$S/hits_extended.bed out/$S/fragment_length.txt
R432script plot_histogram.R out/$S/fragment_length.txt $S out/$S/fragment_length.png
 
### Generate browser tracks ###
sort -k1,1 -k2,2n out/$S/hits_extended.bed > out/$S/hits_extended_sorted.bed
bedtools genomecov -bg -g $CHROM_SIZES -i out/$S/hits_extended_sorted.bed > out/$S/coverage_extended.bedGraph
/ddn/gs1/home/bennettb/tools/ucsc_18.07.24/bedGraphToBigWig out/$S/coverage_extended.bedGraph $CHROM_SIZES out/$S/coverage_extended.bigWig
 
### Normalize browser tracks ###
NUM_READS=`wc -l out/$S/hits_extended.bed | awk '{print $1}'`
perl normalize_bedgraph.pl out/$S/coverage_extended.bedGraph $NUM_READS $DEPTH_NORMALIZATION_FACTOR out/$S/coverage_extended_normalized.bedGraph
/ddn/gs1/home/bennettb/tools/ucsc_18.07.24/bedGraphToBigWig out/$S/coverage_extended_normalized.bedGraph $CHROM_SIZES out/$S/coverage_extended_normalized.bigWig