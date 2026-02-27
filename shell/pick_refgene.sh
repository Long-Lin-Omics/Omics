#!/bin/bash
set -euo pipefail

#####################################
# Defaults
#####################################
LIST_COL=1            # gene.list ID column（1-based）
KEY="transcript_id"   # refGene key to match in ID column 
ATTR_COL=4            # refGene ID column 
SEP=";"               # seperator
OUT_PREFIX="picked"   # output suffix

usage() {
  cat <<EOF
Usage:
  $0 -l gene.list -r mm10.refGene.gene [-c LIST_COL] [-k KEY] [-a ATTR_COL] [-o OUT_PREFIX]

Required:
  -l   gene.list file
  -r   refGene.gene file

Optional:
  -c   which column in gene.list to use as ID (1-based). Default: ${LIST_COL}
  -k   which attribute KEY in refGene to match (e.g. transcript_id, gene_id, gene_name). Default: ${KEY}
  -a   which column in refGene contains attributes string. Default: ${ATTR_COL}
  -o   output prefix. Default: ${OUT_PREFIX}

Outputs:
  \${OUT_PREFIX}.matches.refGene.gene   matched lines from refGene
  \${OUT_PREFIX}.not_found.ids          IDs from gene.list not found (keeps gene.list order)

Examples:
  $0 -l gene.list -r mm10.refGene.gene -c 1 -k transcript_id -o out
  $0 -l gene.list -r mm10.refGene.gene -c 2 -k gene_id -o out_gene
EOF
  exit 1
}

LIST=""
REF=""

while getopts "l:r:c:k:a:o:h" opt; do
  case "$opt" in
    l) LIST="$OPTARG" ;;
    r) REF="$OPTARG" ;;
    c) LIST_COL="$OPTARG" ;;
    k) KEY="$OPTARG" ;;
    a) ATTR_COL="$OPTARG" ;;
    o) OUT_PREFIX="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

[[ -z "${LIST}" || -z "${REF}" ]] && usage
[[ -f "${LIST}" ]] || { echo "Error: gene.list not found: ${LIST}" >&2; exit 1; }
[[ -f "${REF}"  ]] || { echo "Error: refGene not found: ${REF}" >&2; exit 1; }

MATCH_OUT="${OUT_PREFIX}.matches.refGene.gene"
MISS_OUT="${OUT_PREFIX}.not_found.ids"
: > "$MATCH_OUT"
: > "$MISS_OUT"

awk -F'\t' -v list_col="$LIST_COL" -v key="$KEY" -v attr_col="$ATTR_COL" -v sep="$SEP" \
    -v match_out="$MATCH_OUT" -v miss_out="$MISS_OUT" '
function trim(s) { gsub(/^[ \t]+|[ \t]+$/, "", s); return s }

NR==FNR {
    id = $list_col
    id = trim(id)
    if (id == "" || id ~ /^#/) next

    # 去掉版本号：NM_000000.1 -> NM_000000
    sub(/\.[0-9]+$/, "", id)

    want[id]=1
    order[++n]=id
    next
}
{
    attrs = $attr_col
    if (attrs == "") next

    m = split(attrs, a, sep)
    found = ""

    for (i=1; i<=m; i++) {
        s = trim(a[i])

        if (s ~ ("^" key "[ \t]*")) {
            gsub(("^" key "[ \t]*"), "", s)
            gsub(/"/, "", s)
            s = trim(s)
            sub(/\.[0-9]+$/, "", s)
            found = s
            break
        }
    }

    if (found != "" && want[found]) {
        if (! hit[found]) { 
            print > match_out
            hit[found]=1
         }
    }
}
END {
    for (i=1; i<=n; i++) {
        id = order[i]
        if (want[id] && !hit[id]) print id > miss_out
    }
}
' "$LIST" "$REF"

echo "Wrote: $MATCH_OUT"
echo "Wrote: $MISS_OUT"
echo "Done."