#!/bin/bash
set -euo pipefail

#########################
# 默认参数
#########################
THREADS=8
REFERENCE_POINT="center"   # center | TSS | TES
BED=""
PEAKS=""
BW_FILES=""
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_METAPLOT="$DIR/../python/metaplot.py"

# 修法1：-m 为开关；-M 为脚本路径
METAPLOT=0
METAPLOT_SCRIPT=""

#########################
# 用法函数
#########################
usage() {
    echo "Usage:"
    echo "  $0 -p <peaks> [-b <bed>] -s <bigWig(s)> -o <out_prefix> [-t <threads>] [-r <center|TSS|TES>] [-m] [-M <metaplot_script>]"
    echo ""
    echo "Required:"
    echo "  -p  peaks file(s): one or multiple files (space-separated, quoted) OR"
    echo "  -b  bed file (optional; if provided, uses it directly)"
    echo "  -s  bigWig files (space-separated, quoted)"
    echo "  -o  output prefix"
    echo ""
    echo "Optional:"
    echo "  -t  threads (default 8)"
    echo "  -m  use Python metaplot (if -M not provided, uses default: $DEFAULT_METAPLOT)"
    echo "  -M  path to Python metaplot script (works only when -m is set)"
    echo "  -r  referencePoint for computeMatrix reference-point: center|TSS|TES (default center)"
    echo ""
    echo "Behavior:"
    echo "  * Without -m: plotHeatmap output PDF (with shortened sample labels)."
    echo "  * With -m: run Python metaplot and output XLSX."
    exit 1
}

#########################
# 名称缩短函数：basename + 去扩展名 + 去常见 suffix
#########################
shorten_label() {
    local f="$1"
    local name
    name="$(basename "$f")"

    # 去掉扩展名
    name="$(echo "$name" | sed -E 's/\.(bw|bigWig)(\.gz)?$//')"

    # 去掉“多种可能的尾部后缀”（按需自行追加）
    # 规则：只去末尾一次/多次都行（用循环更激进；这里用一次性多后缀匹配更简单）
    name="$(echo "$name" | sed -E '
        s/(_BPM_normalized|_RPM_normalized|_CPM_normalized|_RPKM_normalized)$//;
        s/(_BPM|_RPM|_CPM|_RPKM)$//;
        s/(_normalized|_norm|_normalised)$//;
        s/(_signal|_coverage|_track|_wiggle|_bw|_bigwig)$//;
        s/(_clean|_filtered)$//;
    ')"

    # 再做一次：有些文件名可能叠了多个后缀（例如 *_BPM_normalized_clean）
    name="$(echo "$name" | sed -E '
        s/(_BPM_normalized|_RPM_normalized|_CPM_normalized|_RPKM_normalized)$//;
        s/(_BPM|_RPM|_CPM|_RPKM)$//;
        s/(_normalized|_norm|_normalised)$//;
        s/(_signal|_coverage|_track|_wiggle|_bw|_bigwig)$//;
        s/(_clean|_filtered)$//;
    ')"

    echo "$name"
}

#########################
# 解析命令行参数（修法1）
#########################
while getopts "p:b:s:o:t:r:mM:" opt; do
    case $opt in
        p) PEAKS="$OPTARG" ;;
        b) BED="$OPTARG" ;;
        s) BW_FILES="$OPTARG" ;;
        o) OUT_PREFIX="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        r) REFERENCE_POINT="$OPTARG" ;;
        m) METAPLOT=1 ;;
        M) METAPLOT_SCRIPT="$OPTARG" ;;
        *) usage ;;
    esac
done

#########################
# 参数检查
#########################
if [ -z "${PEAKS:-}" ] && [ -z "${BED:-}" ]; then
    echo "Error: You must provide either peaks (-p) or BED (-b)"
    usage
fi
[ -z "${BW_FILES:-}" ] && { echo "Error: -s is required"; usage; }
[ -z "${OUT_PREFIX:-}" ] && { echo "Error: -o is required"; usage; }

# 如果用户开了 -m 但没给 -M，就用默认脚本
if [ "$METAPLOT" -eq 1 ] && [ -z "$METAPLOT_SCRIPT" ]; then
    METAPLOT_SCRIPT="$DEFAULT_METAPLOT"
fi

case "$REFERENCE_POINT" in
  center|TSS|TES) ;;
  *)
    echo "Error: -r must be one of: center, TSS, TES (got: $REFERENCE_POINT)"
    usage
    ;;
esac


#########################
# 解析 peaks / bw 为数组（更稳）
#########################
# 用户一般会把它们用引号括起来：-s "a.bw b.bw"
read -r -a BW_ARR <<< "$BW_FILES"

# PEAKS 可为空（如果走 -b），不为空就解析
PEAKS_ARR=()
if [ -n "${PEAKS:-}" ]; then
    read -r -a PEAKS_ARR <<< "$PEAKS"
fi

#########################
# 文件定义
#########################
MERGED_BED="${OUT_PREFIX}.merged.bed"
MATRIX="${OUT_PREFIX}.matrix.gz"
TSV="${OUT_PREFIX}.matrix.tsv"
XLSX="${OUT_PREFIX}.xlsx"
HEATMAP="${OUT_PREFIX}.heatmap.pdf"

#########################
# Step 1: 准备 BED
#########################
if [ -n "$BED" ] && [ -f "$BED" ]; then
    echo "Using provided BED: $BED"
    REGION_FILE="$BED"
elif [ "${#PEAKS_ARR[@]}" -gt 0 ]; then
    if [ ! -f "$MERGED_BED" ]; then
        echo "Generating merged BED from peaks..."
        cat "${PEAKS_ARR[@]}" | sort -k1,1 -k2,2n | bedtools merge -d 100 -i - > "$MERGED_BED"
    else
        echo "Merged BED already exists: $MERGED_BED, skipping generation"
    fi
    REGION_FILE="$MERGED_BED"
else
    echo "Error: No BED or peaks file found!"
    exit 1
fi

#########################
# Step 2: 生成 samplesLabel（用于 plotHeatmap 显示更短样本名）
#########################
SAMPLE_LABELS=()
for f in "${BW_ARR[@]}"; do
    SAMPLE_LABELS+=("$(shorten_label "$f")")
done

echo "Samples:"
for i in "${!BW_ARR[@]}"; do
    echo "  ${BW_ARR[$i]}  ->  ${SAMPLE_LABELS[$i]}"
done

#########################
# Step 3: computeMatrix
#########################
if [ ! -f "$MATRIX" ]; then
    echo "Generating matrix with computeMatrix..."
    computeMatrix reference-point \
        -p "$THREADS" \
        --referencePoint "$REFERENCE_POINT" \
        -S "${BW_ARR[@]}" \
        -R "$REGION_FILE" \
        --beforeRegionStartLength 2000 \
        --afterRegionStartLength 2000 \
        --binSize 40 \
        --samplesLabel "${SAMPLE_LABELS[@]}" \
        -o "$MATRIX" \
        --outFileNameMatrix "$TSV"
else
    echo "Matrix already exists: $MATRIX, skipping generation"
fi

#########################
# Step 4: 统计 sample 数
#########################
N=$(computeMatrixOperations info -m "$MATRIX" | \
    awk '
    /^Samples:/ {flag=1; next}
    /^Groups:/ {flag=0}
    flag && NF {count++}
    END {print count}
    ')

echo "Detected $N samples."

#########################
# Step 5: 绘图
#########################
if [ "$METAPLOT" -eq 1 ]; then
    # 使用 Python 脚本绘图（输出 XLSX）
    if [ ! -f "$METAPLOT_SCRIPT" ]; then
        echo "Error: metaplot script not found: $METAPLOT_SCRIPT"
        exit 1
    fi
    if [ ! -f "$XLSX" ]; then
        echo "Running Python metaplot script: $METAPLOT_SCRIPT"
        python "$METAPLOT_SCRIPT" "$TSV" "$XLSX" "$N"
    else
        echo "Output XLSX already exists: $XLSX, skipping plotting"
    fi
else
    # 默认使用 plotHeatmap（输出 PDF）
    if [ ! -f "$HEATMAP" ]; then
        echo "No -m provided, using plotHeatmap..."
        plotHeatmap \
            -m "$MATRIX" \
            -out "$HEATMAP" \
            --colorMap RdBu \
            --whatToShow 'plot, heatmap and colorbar'
            # 如果你想带曲线：把上一行改成
            # --whatToShow 'heatmap and profile and colorbar'
    else
        echo "Heatmap already exists: $HEATMAP, skipping plotting"
    fi
fi

echo "All done."