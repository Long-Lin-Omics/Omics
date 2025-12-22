#!/usr/bin/env python3

import gzip
import random
import argparse

def downsample_paired_fastq(fastq1, output1, subsample_fraction):
    """
    Downsamples paired-end FASTQ files to a specified fraction.
    
    Parameters:
    - fastq1: Path to input R1 FASTQ file (can be gzipped).
    - fastq2: Path to input R2 FASTQ file (can be gzipped).
    - output1: Path to output downsampled R1 FASTQ file.
    - output2: Path to output downsampled R2 FASTQ file.
    - subsample_fraction: Fraction of reads to keep (0 < fraction ≤ 1).
    """
    assert 0 < subsample_fraction <= 1, "subsample_fraction must be between 0 and 1"

    def open_file(filename, mode="rt"):
        return gzip.open(filename, mode) if filename.endswith(".gz") else open(filename, mode)
    
    with open_file(fastq1) as f1, \
         open_file(output1, "wt") as out1:
        
        while True:
            # Read 4 lines per read from both FASTQ files
            r1 = [f1.readline() for _ in range(4)]
            
            # Stop when reaching the end of file
            if not r1[0]:
                break
            
            # Randomly decide whether to keep the read pair
            if random.random() < subsample_fraction:
                out1.writelines(r1)

def main():
    parser = argparse.ArgumentParser(description="Downsample paired-end FASTQ files.")
    parser.add_argument("fastq1", help="Input FASTQ file for Read 1 (R1)")
    parser.add_argument("output1", help="Output downsampled FASTQ file for R1")
    parser.add_argument("fraction", type=float, help="Fraction of reads to keep (0 < fraction ≤ 1)")
    
    args = parser.parse_args()
    
    downsample_paired_fastq(args.fastq1, args.output1, args.fraction)

if __name__ == "__main__":
    main()
