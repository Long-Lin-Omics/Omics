#!/bin/bash

# Usage:
# ./make_ucsc_hub.sh /path/to/bw_dir  "hub_identifier"

# Check inputs
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <bw_directory> <hub_identifier>"
    exit 1
fi

BW_DIR="$1"
IDENTIFIER="$2"
HUB_DIR=/data/wade/linl7/$IDENTIFIER


# Create hub directory
mkdir -p "$HUB_DIR/mm10"

# Copy bigWig files to hub directory
cp "$BW_DIR"/*.bw "$HUB_DIR/mm10/"

# Create hub.txt
cat <<EOF > "$HUB_DIR/hub.txt"
hub $IDENTIFIER
shortLabel $IDENTIFIER
longLabel ChIP-seq Data Hub
genomesFile genomes.txt
email your.email@domain.com
EOF

# Create genomes.txt
cat <<EOF > "$HUB_DIR/genomes.txt"
genome mm10
trackDb mm10/trackDb.txt
EOF

# Create trackDb.txt
TRACKDB="$HUB_DIR/mm10/trackDb.txt"
echo -e "track cont\ncontainer multiWig\nshortLabel $IDENTIFIER\nlongLabel $IDENTIFIER\ntype bigWig\nvisibility full\nalwaysZero on\nautoScale on\nshowSubtrackColorOnUi on\naggregate none\ncolor 1,1,1\npriority 1\n" > "$TRACKDB"

index=1

# Loop over all bigWig files
for bw_file in "$HUB_DIR/mm10/"*.bw; do
    sample=$(basename "$bw_file" .bw)
    cat <<EOF >> "$TRACKDB"
track $sample
bigDataUrl $(basename "$bw_file")
shortLabel $sample signal
longLabel ChIP-seq $sample signal track
type bigWig
visibility full
alwaysZero on
autoScale on
maxHeightPixels 1:60:9999
color 20,20,255
parent cont
priority $index

EOF
    ((index++))
done

echo "UCSC Track Hub created at: ${HUB_DIR/\/data\/wade/https:\/\/orio.niehs.nih.gov\/ucscview}/hub.txt" | tee -a track.path.txt
