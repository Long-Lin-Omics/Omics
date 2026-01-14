#!/usr/bin/env python3

import os
import re
import argparse
from collections import defaultdict

def find_single_fastq(folder):
    fastq_files = [f for f in os.listdir(folder) if f.endswith((".fastq", ".fastq.gz"))]
    sinles = defaultdict(dict)
    for f in fastq_files:
        match = re.match(r"(.+?)(?:_R?1|\.1|single)(?:\.sanger)?\.fastq(?:\.gz)?$", f)
        if match:
            sample_name = match.group(1)
            sinles[sample_name] = os.path.join(folder, f)
    return sinles

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find paired-end FASTQ files in a folder.")
    parser.add_argument("-d", "--directory", type=str,required=True, help="Folder containing FASTQ files")
    parser.add_argument("-s", "--sampleName", type=str,default='(.*)', help="regular pattern to catch sample name")
    parser.add_argument("-m", "--mergeFastqs", action="store_true",  help="print out shell script to merge paired fastqs")
    parser.add_argument("-md", "--megrdDirectory", type=str, default='./', help="Folder to put merged FASTQ files")
    args = parser.parse_args()
    
    single_fastqs = find_single_fastq(args.directory)

    count = {}
    sample_fqs = {}
    total_fqs = 0
    matched_fqs = 0
    named_rep = {}
    if single_fastqs:
        print("\nDetected Single-End FASTQ Files:")
        for sample, r1 in single_fastqs.items():
            total_fqs += 1 
            sampleNameMatch = re.search(args.sampleName, sample)
            if sampleNameMatch:
                sample = sampleNameMatch.groups()[0]
                matched_fqs += 1
                if sample not in count:
                    count[sample] = 1
                    sample_fqs[sample] = {}
                    sample_fqs[sample] = [r1]    
                else:
                    count[sample] +=1
                    sample_fqs[sample].append(r1)
                rep = sample + '-rep-' + str(count[sample])
                if sample not in named_rep:
                    named_rep[sample] = {} 
                named_rep[sample]["\"" + rep + "\""] = 1 
                sample = rep
            else:
                continue
            print(f"{sample}:\n  fastq1: \"{r1}\"")
    else:
        print("No Single-End FASTQ files found.")

    print("\nSummary Table")
    for k,v in sorted(count.items(), key=lambda item: item[1], reverse=True):
        print(k + ': ' + str(v) + " " + ','.join(named_rep[k].keys()))

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

