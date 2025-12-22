#!/bin/bash

# Usage:
# ./download_fastqs_from_ebi.sh sra_ids.txt

if [ $# -ne 1 ]; then
    echo "Usage: $0 sra_ids.txt"
    exit 1
fi

input_file="$1"

> downloaded.fastqs.list

while IFS= read -r sra_id; do
    echo "Fetching FASTQ links for $sra_id..."

    # Get fastq FTP URLs via ENA API
    urls=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${sra_id}&result=read_run&fields=fastq_ftp&format=tsv" \
        | tail -n +2)

    if [[ -z "$urls" ]]; then
        echo "  No FASTQ found for $sra_id."
        continue
    fi

    url_string=$(echo "$urls" | cut -f2- )

    echo -n $sra_id " " >> downloaded.fastqs.list

    # Split semicolon-separated URLs and download
    IFS=';' read -ra url_array <<< "$url_string"
    for url in "${url_array[@]}"; do
        echo "  Downloading  ftp://$url..."
        echo -n $(basename $url) " " >> downloaded.fastqs.list
        wget -nc ftp://"$url"
    done

    echo >> downloaded.fastqs.list

done < "$input_file"

echo " All downloads complete."
