#!/bin/bash
# Usage: ./run_IR_polyA.sh input.bam dorado.tsv output_dir

set -euo pipefail

BAM="$1"
DORADO_TSV="$2"
OUTDIR="$3"

# 创建输出目录
mkdir -p "$OUTDIR"

# 1️⃣ 建立 BAM 链接到输出目录
BAM_LINK="$OUTDIR/input.bam"
ln -sf "$(realpath "$BAM")" "$BAM_LINK"

# 2️⃣ 建立 BAM 索引
if [ ! -f "$BAM_LINK.bai" ]; then
    echo "Indexing BAM..."
    samtools index "$BAM_LINK"
fi

# 3️⃣ 运行 IR_from_miminap.py
IR_OUT="$OUTDIR/IR.txt"
echo "Running IR detection..."
gtf=/ddn/gs1/project/nextgen/post/hug4/LongLin/rnaseq/reference/Mus_musculus.GRCm38.102.gtf
python /ddn/gs1/home/linl7/bin/scripts/nanopore/intron_detection_model/IR_from_miminap.py "$BAM_LINK" $gtf "$IR_OUT" 80 5000000

# 4️⃣ 画 polyA length histogram by IR
OUT_PNG="$OUTDIR/IR.hist.png"
echo "Plotting polyA length histogram..."
Rscript /ddn/gs1/home/linl7/bin/scripts/nanopore/intron_detection_model/polya_length_hist_by_IR.R "$DORADO_TSV" "$IR_OUT" "$OUT_PNG"

echo "Done! All results are in $OUTDIR"
