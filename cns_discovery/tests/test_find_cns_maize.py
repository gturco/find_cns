import unittest
import collections

import sys
import itertools

sys.path.append("../scripts")
from flatfeature import Bed
from find_cns_maize import parse_blast
from find_cns_maize import main

# 
# handle = open('/Users/gturco/code/freeling_lab/find_cns_gturco/pipeline/tests/test_blast_2.txt')
# fh =handle.readlines()
# blast_str= ' , '.join(fh)
# 
# 
# rice_bed = Bed('/Users/gturco/rice_v6_cns_res/04_08_10/test_org/rice_v6.bed'); rice_bed.fill_dict()
# maize_bed = Bed('/Users/gturco/maize/maize_v2.bed', '/Users/gturco/maize/maize_v2.fasta'); maize_bed.fill_dict()
# # sfeat = maize_bed.accn('GRMZM2G039586') #5 
# # qfeat = rice_bed.accn('Os06g37450')
# sfeat = maize_bed.accn('GRMZM2G159142') #5 
# qfeat = rice_bed.accn('Os02g01720')
# 
# cns = parse_blast(blast_str, -1 , qfeat, sfeat, rice_bed, maize_bed, 12000, 26000)
# print cns
#maybe print if used opp or same for both when running it...
from pyfasta import Fasta

class TestMaize(unittest.TestCase):
    def setUp(self):
        handle = open('/Users/gturco/code/freeling_lab/find_cns_gturco/pipeline/tests/blast_3.txt')
        fh =handle.readlines()
        self.blast_str = ' , '.join(fh)
        self.unmasked_fasta = Fasta('/Users/gturco/find_cns/maize_v2_UM.fasta')
        
        self.qbed = Bed('/Users/gturco/rice_maize/rice_v6.bed'); self.qbed.fill_dict()
        self.sbed = Bed('/Users/gturco/maize/maize_v2.bed', '/Users/gturco/maize/maize_v2.fasta'); self.sbed.fill_dict()
        self.sfeat = self.sbed.accn('GRMZM2G086714') 
        self.qfeat = self.qbed.accn('Os09g27050')
        
    def test_get_cmd(self):
        sfasta = 'data/rice_v6_maize_v2/maize_v2_split/2.fasta'
        qfasta = 'data/rice_v6_maize_v2/rice_v6_split/4.fasta'
    
    def test_parse_balse(self):
        orientaion = -1
        cns = parse_blast(self.blast_str, orientaion, self.qfeat, self.sfeat, self.qbed, self.sbed, 12000, 26000,self.unmasked_fasta )
        print cns
        # assertEqual(cns, cns_same_strand)
        
        


    # def test_map(self):
    #     num_list = [(1,2,3,4), (3,5,4,6)]
    #     result = map(num_list, multiply)
    #     print result
    # def multiply(x):
    #     return (x[0], x[1], x[2] * -1, x[3] * -1)


if __name__ == '__main__':
    unittest.main()