

dorado_basecall_mapping_polya_kit_and_custom_primer(){
    pod_folder="$1"
    reference="$2"
    kit_name="$3"
    polya_config="$4"
    output_bam="$5"
    ~/bin/dorado-1.0.0-linux-x64/bin/dorado basecaller hac $pod_folder --reference $reference --kit-name $kit_name --estimate-poly-a --poly-a-config $polya_config > $output_bam 2>$output_bam.dorodo.log
}

dorado_basecall_mm10_polya_kit_and_custom_primer(){
    pod_folder="$1"
    reference="/ddn/gs1/project/nextgen/post/hug4/LongLin/rnaseq/reference/mm10-ont.mmi"
    kit_name="$2"
    polya_config="$3"
    output_bam="$4"
    ~/bin/dorado-1.0.0-linux-x64/bin/dorado basecaller hac $pod_folder --reference $reference --kit-name $kit_name --estimate-poly-a --poly-a-config $polya_config > $output_bam 2>$output_bam.dorodo.log
    awk '/PolyA tails called/ {
        called=$8; not_called=$11;
        total=called+not_called;
        printf "Called: %.2f%%\nNot called: %.2f%%\n", (called/total)*100, (not_called/total)*100
    }' $output_bam.dorodo.log
    plot_polyA_distribution_from_bam $output_bam $output_bam
}


dorado_basecall_mm10_polya_kit_and_custom_primer_no_trim(){
    pod_folder="$1"
    reference="/ddn/gs1/project/nextgen/post/hug4/LongLin/rnaseq/reference/mm10-ont.mmi"
    kit_name="$2"
    polya_config="$3"
    output_bam="$4"
    ~/bin/dorado-1.0.0-linux-x64/bin/dorado basecaller hac $pod_folder --reference $reference --kit-name $kit_name --estimate-poly-a --poly-a-config $polya_config --no-trim > $output_bam 2>$output_bam.dorodo.log
    awk '/PolyA tails called/ {
        called=$8; not_called=$11;
        total=called+not_called;
        printf "Called: %.2f%%\nNot called: %.2f%%\n", (called/total)*100, (not_called/total)*100
    }' $output_bam.dorodo.log
    plot_polyA_distribution_by_barcode_from_bam $output_bam $output_bam
}




plot_polyA_distribution_from_bam() {
  local bam_file=$1
  local output_prefix=$2

  if [[ ! -f "$bam_file" ]]; then
    echo "Error: BAM file '$bam_file' not found."
    return 1
  fi

  echo "Extracting pt:i tag from $bam_file ..."

  # Extract pt:i tag values
  samtools view -F 0x100 -F 0x800 "$bam_file" | awk -F'\t' '
  {
    if (!seen[$1]++) {
        for (i=12; i<=NF; i++) {
          if ($i ~ /^pt:i:/) {
            split($i, a, ":")
            print a[3]
            break
          }
        }
    }
  }' > "${output_prefix}_polyA_lengths.txt"

  echo "Extracted tail lengths saved to ${output_prefix}_polyA_lengths.txt"

  # Generate Python plotting script
  cat << EOF > "${output_prefix}_plot_polyA.py"
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Load data
polyA_lengths = pd.read_csv("${output_prefix}_polyA_lengths.txt", header=None, names=['length'])

# Filter out unestimated reads (-1)
polyA_lengths_filtered = polyA_lengths[polyA_lengths['length'] > 0]
print(polyA_lengths_filtered.describe())

# Plot
plt.figure(figsize=(10,6))
sns.histplot(polyA_lengths_filtered['length'], bins=50, kde=True)
plt.xlabel('PolyA Tail Length (nt)')
plt.ylabel('Read Count')
plt.title('Distribution of PolyA Tail Lengths')
plt.grid(True)

# Save to PNG
plt.savefig("${output_prefix}_polyA_length_distribution.png", dpi=300, bbox_inches='tight')
EOF

  # Run Python script
  python "${output_prefix}_plot_polyA.py"

  echo "Plot saved as ${output_prefix}_polyA_length_distribution.png"
  rm ${output_prefix}_polyA_lengths.txt "${output_prefix}_plot_polyA.py"

}

plot_polyA_distribution_by_barcode_from_bam() {
  local bam_file=$1
  local output_prefix=$2

  if [[ ! -f "$bam_file" ]]; then
    echo "Error: BAM file '$bam_file' not found."
    return 1
  fi

  echo "Extracting pt:i and BC:Z tags from $bam_file ..."

  # Extract pt:i and BC:Z tag values per read
  samtools view -F 0x100 -F 0x800 "$bam_file" | awk -F'\t' '
  {
    if (!seen[$1]++) {
      for (i=12; i<=NF; i++) {
        if ($i ~ /^BC:Z:/) {
          split($i, b, ":")
          printf b[3] "\t"
        }
        if ($i ~ /^pt:i:/) {
          split($i, a, ":")
          print a[3]
        }
      }
    }
  }' > "${output_prefix}_barcode_polyA_lengths.tsv"

  echo "Extracted tail lengths saved to ${output_prefix}_barcode_polyA_lengths.tsv"

  # Generate Python plotting script
   cat << EOF > "${output_prefix}_plot_polyA_by_barcode.py"
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math

# Load data
df = pd.read_csv("${output_prefix}_barcode_polyA_lengths.tsv", sep='\\t', header=None, names=['barcode', 'length'])

# Filter out unestimated (-1)
df = df[df['length'] > 0]

# Get unique barcodes and number
barcodes = sorted(df['barcode'].unique())
num_barcodes = len(barcodes)
cols = 2
rows = math.ceil(num_barcodes / cols)

# Set plot style
sns.set(style="whitegrid")
fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 4*rows))
axes = axes.flatten()

# Plot histogram for each barcode
for i, bc in enumerate(barcodes):
    ax = axes[i]
    subset = df[df['barcode'] == bc]
    count = subset.shape[0]
    print(bc)
    print(subset['length'].quantile([0, 0.25, 0.5, 0.75, 1]))
    sns.histplot(subset['length'], bins=50, kde=True, ax=ax)
    ax.set_title(f"{bc} (n={count})")
    ax.set_xlabel('PolyA Tail Length (nt)')
    ax.set_ylabel('Read Count')
    ax.grid(True)

# Hide unused subplots
for j in range(i+1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.savefig("${output_prefix}_polyA_length_distribution_by_barcode.png", dpi=300)
EOF

  # Run Python script
  python "${output_prefix}_plot_polyA_by_barcode.py"

  echo "Plot saved as ${output_prefix}_polyA_length_distribution_by_barcode.png"
  #rm "${output_prefix}_barcode_polyA_lengths.tsv" "${output_prefix}_plot_polyA_by_barcode.py"
}

dorado_downstream_pipeline() {
    bamfile="$1"
    tts_bed="/ddn/gs1/project/nextgen/post/hug4/LongLin/rnaseq/reference/mm10.tts.bed"
    gene_bed="/ddn/gs1/project/nextgen/post/hug4/LongLin/rnaseq/reference/mm10.gene.bed"
    output_prefix="${bamfile%.bam}"

    # check BAM exists or not
    if [ ! -e "$bamfile" ]; then
        echo "BAM not exist: $bamfile"
        return 1
    fi

    # check .bai 
    if [ ! -e "${bamfile}.bai" ]; then
        echo "no .bai file. sorting nad indexing..."
        sorted_bam="${output_prefix}.sorted.bam"
        samtools sort -@ 80 -o "$sorted_bam" "$bamfile"
        samtools index -@ 80 "$sorted_bam"
        bamfile="$sorted_bam"
    fi

    if [ ! -e "${output_prefix}.tsv" ]; then
      echo "extreating polyA length from bam"
      python /ddn/gs1/project/nextgen/post/hug4/LongLin/nanopore/051625-PAS_WT_sue/polyA.alignment.py "$bamfile" "${output_prefix}.tsv"
    fi

    if [ ! -e "${output_prefix}.intersect.tts.txt" ]; then
      echo "calculating distace to tts"
      bedtools intersect -a "${output_prefix}.tsv" -b "$gene_bed" -wa -wb | \
      awk '{
          if($5=="+" && $13 == "+") 
              print $0, $3-$10
          else if($5=="-" && $13=="-") 
              print $0, ($2-$9)*-1
            else print $0, "NA"
      }' > "${output_prefix}.intersect.tts.txt"
    fi

    # 画图
    echo "plotting distributions "
    Rscript /ddn/gs1/project/nextgen/post/hug4/LongLin/nanopore/051625-PAS_WT_sue/plota.R "${output_prefix}.intersect.tts.txt"

    echo "getting stats"
    Rscript /ddn/gs1/project/nextgen/post/hug4/LongLin/nanopore/051625-PAS_WT_sue/stata.R "${output_prefix}.tsv" "${output_prefix}.intersect.tts.txt" "${output_prefix}"


    echo "done"
}




