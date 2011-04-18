import unittest
import collections

import sys
import itertools

sys.path.append("../scripts")
from flatfeature import Bed
from assign import get_cns_dict
from assign import assign
from assign import nearest_feat

class TestAssign(unittest.TestCase):
    def setUp(self):
        self.cns_filename = "resources/cnsfile.csv"
        self.cns_dict = get_cns_dict(self.cns_filename)

    def test_get_cns_dict(self):
        """test for test_get_cns_dict"""
        
        self.assertEqual(type(collections.defaultdict(list)), type(self.cns_dict))
        self.assertEqual(2, len(self.cns_dict))
        
        test_key = ('1', '5', (3093668, 3093683, 3714671, 3714686))
        expected_value = [('Os01g06550', '25', '[Os05g07040', 'Os05g07060]'), ('Os01g06560', '26', '[Os05g07040', 'Os05g07060]')]
        self.assertEqual(expected_value, self.cns_dict[test_key])
        
        test_key = ('2', '4', (25504628,25504678,26128206,26128256))
        expected_value = [('Os02g42406', '3613', '[Os04g44440','Os04g44500]'), ('Os02g42400', '3612', '[Os04g44440','Os04g44500]')]
        self.assertEqual(expected_value, self.cns_dict[test_key])

    def test_nearest_feat(self):
        
        postion_of_nearest_gene = nearest_feat(((10,12),(20,21)), 18, 19)
        self.assertEqual(1,postion_of_nearest_gene)
        
    def test_assign(self):
        """test for test_assign"""
        qbed_file = "/Users/gturco/rice_v6_cns_res/04_08_10/test_org/rice_v6.bed"
        qbed = Bed(qbed_file)
        qbed.fill_dict()
        
        filtered_cns_values = list(assign(self.cns_dict, qbed, {}))
        
        self.assertEqual(2, len(filtered_cns_values))
        
        first_filtered_cns_value = filtered_cns_values[0]
        self.assertTrue(5, first_filtered_cns_value)
        #cns, saccn, saccn_l, saccn_r, qfeat
        #print first_filtered_cns_value
        
if __name__ == '__main__':
    unittest.main()
