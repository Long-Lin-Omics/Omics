#!/usr/bin/env python3

import os
import re
import argparse
from collections import defaultdict

def find_paired_fastq(folder):
    """
    Identifies paired-end FASTQ files in a given folder.
    Returns:
        A dictionary where keys are sample names and values are tuples (R1, R2).
    """
    fastq_files = [f for f in os.listdir(folder) if f.endswith((".fastq", ".fastq.gz"))]
    pairs = defaultdict(dict)
    for f in fastq_files:
        match = re.match(r"(.+?)(?:_R?1|\.1)(?:\.sanger)?\.fastq(?:\.gz)?$", f)  # Match R1 files
        if match:
            sample_name = match.group(1)
            pairs[sample_name]["R1"] = os.path.join(folder, f)
        else:
            match = re.match(r"(.+?)(?:_R?2|\.2)(?:\.sanger)?\.fastq(?:\.gz)?$", f)  # Match R2 files
            if match:
                sample_name = match.group(1)
                pairs[sample_name]["R2"] = os.path.join(folder, f)
    paired_samples = {
        sample: (data["R1"], data["R2"]) 
        for sample, data in pairs.items() if "R1" in data and "R2" in data
    }
    return paired_samples

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find paired-end FASTQ files in a folder.")
    parser.add_argument("-d", "--directory", type=str,required=True, help="Folder containing FASTQ files")
    parser.add_argument("-s", "--sampleName", type=str,default='(.*)', help="regular pattern to catch sample name")
    parser.add_argument("-m", "--mergeFastqs", action="store_true",  help="print out shell script to merge paired fastqs")
    parser.add_argument("-md", "--megrdDirectory", type=str, default='./', help="Folder to put merged FASTQ files")
    args = parser.parse_args()
    
    paired_fastq_files = find_paired_fastq(args.directory)

    count = {}
    sample_fqs = {}
    total_fqs = 0
    matched_fqs = 0
    if paired_fastq_files:
        print("\nDetected Paired-End FASTQ Files:")
        for sample, (r1, r2) in paired_fastq_files.items():
            total_fqs += 2 
            sampleNameMatch = re.search(args.sampleName, sample)
            if sampleNameMatch:
                sample = sampleNameMatch.groups()[0]
                matched_fqs += 2
                if sample not in count:
                    count[sample] = 1
                    sample_fqs[sample] = {}
                    sample_fqs[sample]['r1'] = [r1]
                    sample_fqs[sample]['r2'] = [r2]    
                else:
                    count[sample] +=1
                    sample_fqs[sample]['r1'].append(r1)
                    sample_fqs[sample]['r2'].append(r2)
                sample = sample + '-rep-' + str(count[sample]) 
            else:
                continue
            print(f"{sample}:\n  fastq1: \"{r1}\"\n  fastq2: \"{r2}\"")
            print(f"  control: \"\"")
    else:
        print("No paired FASTQ files found.")

    print("\nSummary Table")
    for k,v in sorted(count.items(), key=lambda item: item[1], reverse=True):
        print(k + ': ' + str(v))

    if matched_fqs == total_fqs:
        print("All fastqs captured!")
    else:
        print("The number of matched fastqs is " + str(matched_fqs) + ', different from the number of all fastq files: ' + str(total_fqs)) 

    if args.mergeFastqs:
        merged_fqs_dir = args.megrdDirectory + '/merged_fqs' 
        os.mkdir(merged_fqs_dir)
        merge_sh = open(merged_fqs_dir + '/merge.sh', 'wt')
        for sample in sample_fqs.keys():
            merged_fq1 = merged_fqs_dir + '/' + sample + '.merged.1.fastq.gz'
            merged_fq2 = merged_fqs_dir + '/' + sample + '.merged.2.fastq.gz'
            merge_sh.write("cat " + ' '.join(sample_fqs[sample]['r1']) + ' > ' + merged_fq1 + '\n')
            merge_sh.write("cat " + ' '.join(sample_fqs[sample]['r2']) + ' > ' + merged_fq2 + '\n')
            print(f"{sample}:\n  fastq1: \"{merged_fq1}\"\n  fastq2: \"{merged_fq2}\"")
    else:
        print(sample_fqs)

