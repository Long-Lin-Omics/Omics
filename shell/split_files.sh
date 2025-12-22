#!/bin/bash

# Define variables
SRC_DIR=$1  # Change this to your source folder
FILES_PER_DIR=$2         # Number of files per subdirectory
PREFIX="batch_"          # Prefix for subdirectories
COUNT=0                  # Counter for subdirectory naming
FILE_COUNT=0             # Counter for files added

# Create destination directories and move files
mkdir -p "$SRC_DIR"  # Ensure the source directory exists
for file in "$SRC_DIR"/*; do
    # Skip if not a regular file
    [ -f "$file" ] || continue  

    # Create a new subdirectory every FILES_PER_DIR files
    if (( FILE_COUNT % FILES_PER_DIR == 0 )); then
        COUNT=$((COUNT + 1))
        DIR_NAME="$SRC_DIR/$PREFIX$COUNT"
        mkdir -p "$DIR_NAME"
    fi

    # Move the file into the current subdirectory
    mv "$file" "$DIR_NAME/"
    FILE_COUNT=$((FILE_COUNT + 1))
done

echo "Files have been organized into subdirectories."
