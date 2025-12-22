
merge_multiple_peaks() {
    output_merged="$1"
    shift  # 把第一个参数移掉，后面的都是 input files
    for file in "$@"; do
        if [ ! -e "$file" ]; then
            echo "File not exists: $file"
            return 1
        else
            echo "peak number of $file is $(wc -l $file |cut -d ' ' -f 1 )"
        fi
    done
    cat "$@" | sort -k1,1 -k2,2n | bedtools merge -d 100 -i - > "$output_merged"
    if [ $? -eq 0 ]; then
        echo "Merge done: $output_merged with $(wc -l $output_merged |cut -d ' ' -f 1 ) peaks"
    else
        echo "Merging failed!"
        return 1
    fi
}

plot_peak_heatmap() {
    bw_files="$1"        # bigwig 文件，空格分隔
    bed_file="$2"        # bed 文件
    selected_tag="$3"    # 标签名
    computeMatrix reference-point -S $bw_files -R $bed_file --skipZeros --binSize 40 -p 8 -a 3000 -b 3000 --referencePoint center -o "${selected_tag}.mat.gz"
    plotHeatmap -m "${selected_tag}.mat.gz" -o ${selected_tag}.png --colorMap YlGnBu --missingDataColor white --legendLocation none
    rm "${selected_tag}.mat.gz"
}
#plot_peak_heatmap "data/sample1.bw data/sample2.bw" regions.bed H3K9me3

plot_peak_near_gene_heatmap() {
    bw_files="$1"        # bigwig 文件，空格分隔
    peak_file="$2"        # bed 文件
    gene_bed_file="$3"
    selected_tag="$4"    # 标签名
    intersect_bed='nearGene.bed'
    bedtools window -a $peak_file -b $gene_bed_file -w 5000 | cut -f1-3 |sort |uniq > $intersect_bed 
    computeMatrix reference-point -S $bw_files -R $intersect_bed --skipZeros --binSize 40 -p 8 -a 3000 -b 3000 --referencePoint center -o "${selected_tag}.mat.gz"
    plotHeatmap -m "${selected_tag}.mat.gz" -o ${selected_tag}.nearGene.png --colorMap YlGnBu --missingDataColor white --legendLocation none
    rm "${selected_tag}.mat.gz"
    echo "Number of peaks reamining: $(wc -l $intersect_bed | cut -d ' '  -f 1)"
    rm $intersect_bed     
}

plot_nearest_peaks_of_genes_heatmap() {
    bw_files="$1"        # bigwig 文件，空格分隔
    peak_file="$2"        # bed 文件
    gene_bed_file="$3"
    selected_tag="$4"    # 标签名
    nearest_bed='nearest.bed'
    cut -f 1-3 $gene_bed_file | sort -k1,1 -k2,2n  > a.bed
    cut -f 1-3 $peak_file | sort -k1,1 -k2,2n  > b.bed  
    bedtools closest -a a.bed -b b.bed -d > $selected_tag.closest.peaks
    cut -f 4-6 $selected_tag.closest.peaks > $nearest_bed
    rm a.bed b.bed
    echo "Quantiles of $selected_tag:"
    awk '{print $7}' $selected_tag.closest.peaks | Rscript -e 'x=scan("stdin"); print(quantile(x))'
    computeMatrix reference-point -S $bw_files -R $nearest_bed --skipZeros --binSize 40 -p 8 -a 3000 -b 3000 --referencePoint center -o "${selected_tag}.mat.gz"
    plotHeatmap -m "${selected_tag}.mat.gz" -o ${selected_tag}.gene.nearest_peaks.png --colorMap YlGnBu --missingDataColor white --legendLocation none
    rm "${selected_tag}.mat.gz" 
    echo "Number of peaks reamining: $(wc -l $nearest_bed | cut -d ' '  -f 1)"
    rm $nearest_bed     
}


function location_for_gene(){
    gene_list=$1
    gtf=$2
    output=$3
    >$output
    cat $gene_list |while read line
    do
        # no p
        regex=$line'(?!-ps)'
        content=$(grep -Pw $regex $gtf | awk '$3 == "gene"')
        if [ -z "$content" ]; then
            echo $line "not in the annotation gtf"
        else
            count=$(echo "$content" | wc -l)
            if [ "$count" -ne 1 ]; then
                gene_version=$(echo "$content" |awk -F'gene_version "' '{print $2}' | awk -F'"' '{print $1}' |sort -nr | head -n 1)
                selected_line=$(echo "$content" | grep "gene_version \"$gene_version\"")
                if [ $(echo "$selected_line" | wc -l ) -ne 1 ]; then
                    echo $line ": Count is not 1" 
                else
                    echo "$selected_line" |awk '{print $1"\t"$4"\t"$5}' >> $output
                fi  
            else
                echo "$content" |awk '{print $1"\t"$4"\t"$5}' >> $output 
            fi
        fi
    done
}

#echo "PGC genes"
#location_for_gene /ddn/gs1/project/nextgen/post/hug4/LongLin/Custom/PGC_metaplot/data/PGC.list /ddn/gs1/project/nextgen/post/hug4/LongLin/rnaseq/reference/Mus_musculus.GRCm38.102.gtf /ddn/gs1/project/nextgen/post/hug4/LongLin/Custom/PGC_metaplot/data/PGC.bed

call_peaks_with_seacr_from_bws() {
  local target_bw="$1"
  local control_bw="$2"
  local threshold="${3:-0.01}"  # default to 0.01 if not provided
  local out_prefix="$4"

  if [[ -z "$target_bw" || -z "$control_bw" || -z "$out_prefix" ]]; then
    echo "Usage: call_peaks_with_seacr <target.bigWig> <control.bigWig|'none'> [threshold] <output_prefix>"
    return 1
  fi

  # Convert bigWig to bedGraph
  bigWigToBedGraph "$target_bw" "${out_prefix}_target.bedGraph"
  sort -k1,1 -k2,2n "${out_prefix}_target.bedGraph" > "${out_prefix}_target.sorted.bedGraph"
  rm "${out_prefix}_target.bedGraph"

  if [[ "$control_bw" != "none" ]]; then
    bigWigToBedGraph "$control_bw" "${out_prefix}_control.bedGraph"
    sort -k1,1 -k2,2n "${out_prefix}_control.bedGraph" > "${out_prefix}_control.sorted.bedGraph"
    rm "${out_prefix}_control.bedGraph"
    echo "Calling peaks with control..."
    bash /ddn/gs1/home/linl7/bin/SEACR_1.3/SEACR_1.3.sh "${out_prefix}_target.sorted.bedGraph" "${out_prefix}_control.sorted.bedGraph" norm stringent "${out_prefix}_peaks_with_control.bed"
  else
    echo "Calling peaks without control..."
    bash /ddn/gs1/home/linl7/bin/SEACR_1.3/SEACR_1.3.sh "${out_prefix}_target.sorted.bedGraph" "$threshold" non stringent "${out_prefix}_peaks_no_control.bed"
  fi

  echo "Done. Peaks saved with prefix: ${out_prefix}_peaks*.bed"
  
}


call_peaks_with_seacr_from_bws_non() {
  local target_bw="$1"
  local control_bw="$2"
  local threshold="${3:-0.01}"  # default to 0.01 if not provided
  local out_prefix="$4"

  if [[ -z "$target_bw" || -z "$control_bw" || -z "$out_prefix" ]]; then
    echo "Usage: call_peaks_with_seacr <target.bigWig> <control.bigWig|'none'> [threshold] <output_prefix>"
    return 1
  fi

  # Convert bigWig to bedGraph
  if [ ! -e "${out_prefix}_target.sorted.bedGraph" ]; then
    bigWigToBedGraph "$target_bw" "${out_prefix}_target.bedGraph"
    sort -k1,1 -k2,2n "${out_prefix}_target.bedGraph" > "${out_prefix}_target.sorted.bedGraph"
    rm "${out_prefix}_target.bedGraph"
  fi

  if [[ "$control_bw" != "none" ]]; then
    if [ ! -e "${out_prefix}_control.sorted.bedGraph" ]; then
      bigWigToBedGraph "$control_bw" "${out_prefix}_control.bedGraph"
      sort -k1,1 -k2,2n "${out_prefix}_control.bedGraph" > "${out_prefix}_control.sorted.bedGraph"
      rm "${out_prefix}_control.bedGraph"
    fi
    echo "Calling peaks with control..."
    bash /ddn/gs1/home/linl7/bin/SEACR_1.3/SEACR_1.3.sh "${out_prefix}_target.sorted.bedGraph" "${out_prefix}_control.sorted.bedGraph" non stringent "${out_prefix}_peaks_with_control.bed"
  else
    echo "Calling peaks without control..."
    bash /ddn/gs1/home/linl7/bin/SEACR_1.3/SEACR_1.3.sh "${out_prefix}_target.sorted.bedGraph" "$threshold" non stringent "${out_prefix}_peaks_no_control.bed"
  fi

  echo "Done. Peaks saved with prefix: ${out_prefix}_peaks*.bed"
}
