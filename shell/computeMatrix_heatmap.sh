#!/bin/bash
set -euo pipefail

#########################
# 默认参数
#########################
THREADS=8
BED=""
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
    echo "Usage: $0 -p <peaks> [-b <bed>] -s <bigWig(s)> -o <out_prefix> [-t <threads>] [-m] [-M <metaplot_script>]"
    echo ""
    echo "Required:"
    echo "  -p  peaks file(s). Can be one file or multiple files (space-separated, quoted)."
    echo "  -s  bigWig files (space-separated, quoted)."
    echo "  -o  output prefix"
    echo ""
    echo "Optional:"
    echo "  -b  bed file (if provided, use it directly)"
    echo "  -t  threads (default 8)"
    echo "  -m  use Python metaplot (if -M not provided, uses default: $DEFAULT_METAPLOT)"
    echo "  -M  path to Python metaplot script (requires -m to take effect)"
    echo ""
    echo "Behavior:"
    echo "  * If -m is NOT set: use plotHeatmap to draw heatmap PDF."
    echo "  * If -m is set: run Python metaplot script to produce XLSX."
    exit 1
}

#########################
# 解析命令行参数
#########################
while getopts "p:b:s:o:t:mM:" opt; do
    case $opt in
        p) PEAKS="$OPTARG" ;;
        b) BED="$OPTARG" ;;
        s) BW_FILES="$OPTARG" ;;
        o) OUT_PREFIX="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        m) METAPLOT=1 ;;
        M) METAPLOT_SCRIPT="$OPTARG" ;;
        *) usage ;;
    esac
done

#########################
# 参数检查
#########################
# peaks 或 bed 至少提供一个（但你 usage 里把 -p 设成 required；这里仍兼容你原逻辑）
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
elif [ -n "${PEAKS:-}" ]; then
    if [ ! -f "$MERGED_BED" ]; then
        echo "Generating merged BED from peaks..."
        # PEAKS 可能是多个文件（空格分隔），这里故意不加引号以让 shell 展开成多个参数
        # 一般 narrowPeak/bed 文件名不会带空格；若会带空格，建议改用数组传参。
        cat $PEAKS | sort -k1,1 -k2,2n | bedtools merge -d 100 -i - > "$MERGED_BED"
    else
        echo "Merged BED already exists: $MERGED_BED, skipping generation"
    fi
    REGION_FILE="$MERGED_BED"
else
    echo "Error: No BED or peaks file found!"
    exit 1
fi

#########################
# Step 2: computeMatrix
#########################
if [ ! -f "$MATRIX" ]; then
    echo "Generating matrix with computeMatrix..."
    computeMatrix reference-point \
        -p "$THREADS" \
        --referencePoint center \
        -S $BW_FILES \
        -R "$REGION_FILE" \
        --beforeRegionStartLength 2000 \
        --afterRegionStartLength 2000 \
        --binSize 40 \
        -o "$MATRIX" \
        --outFileNameMatrix "$TSV"
else
    echo "Matrix already exists: $MATRIX, skipping generation"
fi

#########################
# Step 3: 统计 sample 数
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
# Step 4: 绘图
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
    else
        echo "Heatmap already exists: $HEATMAP, skipping plotting"
    fi
fi

echo "All done."