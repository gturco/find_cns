# write a genomic fasta file with all sequences covered by features
# in the specified Bed file masked to N.
from flatfeature import Bed
import sys
# b = Bed(sys.argv[1], sys.argv[2])
b = Bed("/Users/gturco/data/rice_v6.bed", "/Users/gturco/data/rice_v6.fasta")

for seqid, seq in b.mask_cds():
    seqids =  []
    seq.tostring()
