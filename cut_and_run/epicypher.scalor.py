import sys
import gzip
import re
from collections import Counter


barcodes = {
    'Unmodified': ['TTCGCGCGTAACGACGTACCGT','CGCGATACGACCGCGTTACGCG'],
    'H3K4me1': ['CGACGTTAACGCGTTTCGTACG', 'CGCGACTATCGCGCGTAACGCG'],
    'H3K4me2': ['CCGTACGTCGTGTCGAACGACG', 'CGATACGCGTTGGTACGCGTAA'],
    'H3K4me3': ['TAGTTCGCGACACCGTTCGTCG', 'TCGACGCGTAAACGGTACGTCG'],
    'H3K9me1': ['TTATCGCGTCGCGACGGACGTA', 'CGATCGTACGATAGCGTACCGA'],
    'H3K9me2': ['CGCATATCGCGTCGTACGACCG', 'ACGTTCGACCGCGGTCGTACGA'],
    'H3K9me3': ['ACGATTCGACGATCGTCGACGA', 'CGATAGTCGCGTCGCACGATCG'],
    'H3K27me1': ['CGCCGATTACGTGTCGCGCGTA', 'ATCGTACCGCGCGTATCGGTCG'],
    'H3K27me2': ['CGTTCGAACGTTCGTCGACGAT', 'TCGCGATTACGATGTCGCGCGA'],
    'H3K27me3': ['ACGCGAATCGTCGACGCGTATA', 'CGCGATATCACTCGACGCGATA'],
    'H3K36me1': ['CGCGAAATTCGTATACGCGTCG', 'CGCGATCGGTATCGGTACGCGC'],
    'H3K36me2': ['GTGATATCGCGTTAACGTCGCG', 'TATCGCGCGAAACGACCGTTCG'],
    'H3K36me3': ['CCGCGCGTAATGCGCGACGTTA', 'CCGCGATACGACTCGTTCGTCG'],
    'H4K20me1': ['GTCGCGAACTATCGTCGATTCG', 'CCGCGCGTATAGTCCGAGCGTA'],
    'H4K20me2': ['CGATACGCCGATCGATCGTCGG', 'CCGCGCGATAAGACGCGTAACG'],
    'H4K20me3': ['CGATTCGACGGTCGCGACCGTA', 'TTTCGACGCGTCGATTCGGCGA'],
}

def find_marker_from_filename(fastq_files):
    for marker in barcodes:
        if any(marker in f for f in fastq_files):
            return marker
    return 'None'


def count_barcodes(fastq_files, markers):
    counts = Counter({m: 0 for m in markers})

    for fq_file in fastq_files:
        open_func = gzip.open if fq_file.endswith('.gz') else open
        with open_func(fq_file, 'rt') as f:
            for i, line in enumerate(f):
                if i % 4 == 1:  # sequence line in fastq
                    for m in markers:
                        for bc in barcodes[m]:
                            counts[m] += line.count(bc)
    return counts



def main():
    if len(sys.argv) != 6:
        print("Usage: python count_barcodes.py <fastq_file1> <fastq_file2> <sample_read_count_file> <output_file> <sample_name>")
        sys.exit(1)

    fastq_files = [sys.argv[1], sys.argv[2]]
    sample_read_count_file = sys.argv[3]
    output_file = sys.argv[4]
    sample_name = sys.argv[5]

    sample_read_count = int(open(sample_read_count_file).read().strip())

    marker = find_marker_from_filename(fastq_files)

    counts = count_barcodes(fastq_files, barcodes.keys())

    with open(output_file, 'w') as out:
        out.write(sample_name + '\n')
        out.write( "Marker: " + marker +'\n' )
        for bc, count in sorted(counts.items()):
            result = count / counts[marker]
            out.write(f"{bc}\t{count}\t{result}\n")
        out.write("Uniq align reads: " + str(sample_read_count) + '\n')
        count_barcodes_total = sum(counts.values())
        out.write("Total barcode reads: " + str(count_barcodes_total) + '\n')
        bar_rate = count_barcodes_total/sample_read_count*100 
        out.write("% total barcode reads: " + f"{bar_rate:.3g}" + '%\n')

if __name__ == '__main__':
    main()


# a=count_barcodes(['/ddn/gs1/project/nextgen/post/hug4/LongLin/cut_and_tag/NS20170/input/test/2i-wt-H3K9me3-rep-1_R1.fastq','/ddn/gs1/project/nextgen/post/hug4/LongLin/cut_and_tag/NS20170/input/test/2i-wt-H3K9me3-rep-1_R2.fastq'],'H3K9me3')
# print(a)
