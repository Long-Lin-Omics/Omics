#!/usr/bin/bash
Rscript -e "df <- tailfindr::find_tails(fast5_dir='$IN', save_dir='$OUT/A1', csv_filenames='n3p_tails.csv', num_cores=30, dna_datatype='custom-cdna', end_primer='GAACCTACTACCTAAGAGCAAGAAGAAGCC', front_primer='')"
Rscript -e "df <- tailfindr::find_tails(fast5_dir='$IN', save_dir='$OUT/A2', csv_filenames='n3p_tails.csv', num_cores=30, dna_datatype='custom-cdna', end_primer='GAACCTACTACGCAGAAAGACGAGAATCAC', front_primer='')"
Rscript -e "df <- tailfindr::find_tails(fast5_dir='$IN', save_dir='$OUT/A3', csv_filenames='n3p_tails.csv', num_cores=30, dna_datatype='custom-cdna', end_primer='GAACCTACTACCGCGCAAAGAGAAAAGTAC', front_primer='')"
Rscript -e "df <- tailfindr::find_tails(fast5_dir='$IN', save_dir='$OUT/A4', csv_filenames='n3p_tails.csv', num_cores=30, dna_datatype='custom-cdna', end_primer='GAACCTACTAAATAAGACCGAGCGAAGACC', front_primer='')"
