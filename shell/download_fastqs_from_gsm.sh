#!/bin/bash

# Usage:
# ./fetch_fastqs_from_gsms.sh gsm_ids.txt

if [ $# -ne 1 ]; then
    echo "Usage: $0 gsm_ids.txt"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

gsm_file="$1"
srx_file="srx_ids.txt"

# Step 1: Convert GSMs to SRXs
python3 $SCRIPT_DIR/../python/gsm_to_srx.py "$gsm_file" "$srx_file"

# Step 2: Download FASTQs
bash $SCRIPT_DIR/download_fastqs_from_ebi.sh "$srx_file"
