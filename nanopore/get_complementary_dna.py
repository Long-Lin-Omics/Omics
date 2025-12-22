def get_complementary_dna(dna_strand):
    # DNA base pair mappings
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # Generate complementary strand
    complementary_strand = ''.join(complement[base] for base in dna_strand.upper())
    
    return complementary_strand[::-1]


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        for i in sys.argv[1:]:
            print(i + '\t' + get_complementary_dna(i) + '\n')
    else:
        A1 = 'GGCTTCTTCTTGCTCTTAGGTAGTAGGTTC'
        A2 = 'GTGATTCTCGTCTTTCTGCGTAGTAGGTTC'
        A3 = 'GTACTTTTCTCTTTGCGCGGTAGTAGGTTC'
        A4 = 'GGTCTTCGCTCGGTCTTATTTAGTAGGTTC'
        print('A1\t'+A1 + '\t' + get_complementary_dna(A1) +'\n')
        print('A2\t'+A2 + '\t' + get_complementary_dna(A2) +'\n')
        print('A3\t'+A3 + '\t' + get_complementary_dna(A3) +'\n')
        print('A4\t'+A4 + '\t' + get_complementary_dna(A4) +'\n')

