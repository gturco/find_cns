import unittest
import collections

import sys
import itertools

sys.path.append("../scripts")
from flatfeature import Bed
from find_cns import main

class TestFind_CNS(unittest.TestCase):
    def setUp(self):
        self.qfasta = "/Users/gturco/code/freeling_lab/find_cns_gturco/pipeline/tests/data/rice_v6_sorghum_v1/rice_v6.fasta"
        self.sfasta = "/Users/gturco/code/freeling_lab/find_cns_gturco/pipeline/tests/data/rice_v6_sorghum_v1/sorghum_v1.fasta"
        self.qbed = "/Users/gturco/code/freeling_lab/find_cns_gturco/pipeline/tests/data/rice_v6_sorghum_v1/rice_v6.bed"
        self.sbed = "/Users/gturco/code/freeling_lab/find_cns_gturco/pipeline/tests/data/rice_v6_sorghum_v1/sorghum_v1.bed"
        self.pairs = "/Users/gturco/code/freeling_lab/find_cns_gturco/pipeline/tests/data/rice_v6_sorghum_v1/rice_v6_sorghum_v1.pairs.txt"
        self.blast_path = "/Users/gturco/blast-2.2.25/bin/bl2seq"
    def test_main(self):
        """test for test_get_cns_dict"""
        qbed = Bed(self.qbed, self.qfasta); qbed.fill_dict()
        sbed = Bed(self.sbed, self.sfasta); sbed.fill_dict()
        main(qbed, sbed, self.pairs, 12000,12000, "pair", self.blast_path, "T",2)

if __name__ == '__main__':
    unittest.main()
