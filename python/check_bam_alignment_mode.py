import pysam
import sys

bam = sys.argv[1]
bamfile = pysam.AlignmentFile(bam, "rb")
has_secondary = False
has_supplementary = False

for i, read in enumerate(bamfile.fetch(until_eof=True)):
    if read.is_secondary:
        has_secondary = True
    if read.is_supplementary:
        has_supplementary = True
    if has_secondary and has_supplementary:
        break
    
print("Has secondary:", has_secondary)
print("Has supplementary:", has_supplementary)