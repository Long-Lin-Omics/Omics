#!/usr/bin/env bash
set -euo pipefail

FA=$1
GTF=$2

echo "===================================="
echo " FA vs GTF Compatibility Check (gz OK)"
echo "===================================="

# ---------- 工具函数 ----------
read_file() {
    local f=$1
    if [[ $f == *.gz ]]; then
        gunzip -c "$f"
    else
        cat "$f"
    fi
}

# ---------- 检查输入 ----------
if [[ ! -f "$FA" ]]; then
    echo "ERROR: FASTA not found: $FA"
    exit 1
fi

if [[ ! -f "$GTF" ]]; then
    echo "ERROR: GTF not found: $GTF"
    exit 1
fi

if ! command -v samtools &> /dev/null; then
    echo "ERROR: samtools not found."
    exit 1
fi

# ---------- FASTA index ----------
echo "[1] Preparing FASTA index..."

if [[ $FA == *.gz ]]; then
    echo "Decompressing FASTA temporarily..."
    TMP_FA=$(mktemp --suffix=.fa)
    gunzip -c "$FA" > "$TMP_FA"
    samtools faidx "$TMP_FA"
    FAI="${TMP_FA}.fai"
else
    samtools faidx "$FA"
    FAI="${FA}.fai"
    TMP_FA=""
fi

# ---------- 提取 contig ----------
echo "[2] Extracting contig names..."

cut -f1 "$FAI" | sort > fa.chr
read_file "$GTF" | cut -f1 | sort | uniq > gtf.chr

echo "FASTA contigs: $(wc -l < fa.chr)"
echo "GTF contigs:   $(wc -l < gtf.chr)"

# ---------- GTF-only ----------
echo "[3] Checking GTF-only contigs..."
comm -23 gtf.chr fa.chr > gtf_only.chr

if [[ -s gtf_only.chr ]]; then
    echo " Contigs in GTF but NOT in FASTA:"
    cat gtf_only.chr
else
    echo "No GTF-only contigs"
fi

# ---------- FASTA-only ----------
echo "[4] Checking FASTA-only contigs..."
comm -13 gtf.chr fa.chr > fa_only.chr

if [[ -s fa_only.chr ]]; then
    echo " Contigs in FASTA but NOT in GTF:"
    cat fa_only.chr
else
    echo "No FASTA-only contigs"
fi

# ---------- 长度 ----------
echo "[5] Reading FASTA lengths..."
awk '{print $1"\t"$2}' "$FAI" > fa.len

# ---------- GTF 最大坐标 ----------
echo "[6] Computing GTF max positions..."
read_file "$GTF" | awk '{
    if ($5 > max[$1]) max[$1] = $5
}
END {
    for (chr in max)
        print chr"\t"max[chr]
}' | sort > gtf.max

# ---------- 越界检查 ----------
echo "[7] Checking coordinate bounds..."
join -1 1 -2 1 fa.len gtf.max > joined.txt || true

OUTBOUND=0
awk '{
    chr=$1
    fa_len=$2
    gtf_max=$3

    if (gtf_max > fa_len) {
        print "OUT-OF-BOUND:", chr, "GTF:", gtf_max, "FASTA:", fa_len
        bad=1
    }
}
END {
    if (bad != 1) print "✅ All coordinates within bounds"
}' joined.txt | tee bound_check.txt

grep -q "OUT-OF-BOUND" bound_check.txt && OUTBOUND=1 || true

# ---------- 总结 ----------
echo "===================================="
echo " SUMMARY"
echo "===================================="

if [[ -s gtf_only.chr ]]; then
    echo "Incompatible: GTF has contigs not in FASTA"
elif [[ $OUTBOUND -eq 1 ]]; then
    echo " Incompatible: Coordinate overflow detected"
else
    echo " Basic compatibility OK"
fi

echo "===================================="
echo "Done."

# ---------- 清理 ----------
if [[ -n "${TMP_FA:-}" ]]; then
    rm -f "$TMP_FA" "$FAI"
fi


