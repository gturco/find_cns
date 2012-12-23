import unittest
import collections

import sys
import itertools

sys.path.append("../scripts")
from flatfeature import Bed
from assign import get_cns_dict
from assign import assign
from assign import make_pair_maps
from assign import cns_fmt_dict
from assign import main


class TestAssign(unittest.TestCase):
    def setUp(self):
        self.cns_filename = "data/rice_v6_sorghum_v1/rice_v6_sorghum_v1.cns.txt"
        self.pairsfile = "data/rice_v6_sorghum_v1/rice_v6_sorghum_v1.pairs.txt"
        self.qbed = Bed("data/rice_v6_sorghum_v1/rice_v6.bed") ;self.qbed.fill_dict()
        self.sbed = Bed("data/rice_v6_sorghum_v1/sorghum_v1.bed") ;self.sbed.fill_dict()
        self.cns_dict, self.evalue_dict = get_cns_dict(self.cns_filename)
        self.qpair_map, self.spair_map = make_pair_maps(self.pairsfile, "pair", self.qbed, self.sbed)

    def test_get_cns_dict(self):
        """test for test_get_cns_dict"""
        #print self.cns_dict.keys()
        print "keys!",  self.evalue_dict.keys()


    def test_assign(self):
      assign(self.cns_dict, self.qbed, self.sbed, self.qpair_map, self.spair_map)
    
    def test_cns_fmt_dict(self):
      for cns, qfeat, sfeat in assign(self.cns_dict, self.qbed, self.sbed, self.qpair_map, self.spair_map):
        d = cns_fmt_dict(cns, qfeat, sfeat, self.evalue_dict)
        print "dddddddd", d
        
    def test_main(self):
      pass
      #main(self.cns_filename, self.qbed, self.sbed,  self.pairsfile, "pair", 123, 123, 12000,12000)

if __name__ == '__main__':
    unittest.main()
