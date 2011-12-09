import sys
sys.path.insert(0,".")
import blast_misc
import time
import operator

blast_file = 'data/t.blast'
b = blast_misc.blast_array(blast_file, best_hit=0, maxkeep=999999, dopickle=0)
print b
print b.shape
