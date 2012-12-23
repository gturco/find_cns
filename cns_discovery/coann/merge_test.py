import unittest
import sys
sys.path.append("../scripts")
from flatfeature import Bed
from merge2 import parse_missed_genes,group_genes_in_bed,update_locs,merge_hits,no_intervening_genes

class TestMask(unittest.TestCase):
    def setUp(self):
        self.old_bed = Bed("data/rice_t_sorghum_v1/sorghum_v1.bed")
        self.missed_bed = Bed("data/rice_t_sorghum_v1/missed_sorghum_v1_from_rice_b.bed")
        self.matches = "data/rice_t_sorghum_v1/missed_sorghum_v1_from_rice_b.matches.txt"
        self.missed_genes = parse_missed_genes(self.matches)
        self.missed_genes_grouped, self.missed_genes_dict = group_genes_in_bed(self.missed_genes,self.old_bed,self.missed_bed)
    
    def test_group_genes_in_bed(self):
        missed_genes_grouped, missed_genes_dict = group_genes_in_bed(self.missed_genes,self.old_bed,self.missed_bed)
        ### adding to old bed example
        self.assertEqual(missed_genes_dict['Sb01g039400']['locs'], [(62821196, 62822809), (62822899, 62823011)])
        #### example with more then one hit for os03
        self.assertEqual(missed_genes_grouped["Os03g06330"],[('1', 1035243, 1035376, 'sorghum_v1_1_1035243_1035376'), ('1', 43157679, 43159029, 'sorghum_v1_1_43157679_43159029')])
    def test_update_locs(self):
        pass
        #make sure gets samllest and largerst correct
    def test_merge_hits(self):
        hits = tuple(self.missed_genes_grouped["Os03g49400"])
        hits = list(hits)
        merged_hit = merge_hits(hits,self.old_bed, self.missed_genes_dict)
        self.assertEqual(len(merged_hit),1)
        self.assertEqual(len(merged_hit['sorghum_v1_1_9901089_9901320']['locs']),3)
        self.assertEqual(merged_hit['sorghum_v1_1_9901089_9901320']['start'],9896924)
        self.assertEqual(merged_hit['sorghum_v1_1_9901089_9901320']['end'],9901320)
        ### find a large pair of hits that should be merged
    def test_merge_hits_none(self):
        ## none of the hits should merge
        hits = (self.missed_genes_grouped["Os01g58037"])
        merged_hit = merge_hits(hits,self.old_bed,self.missed_genes_dict)
        self.assertEqual(len(merged_hit),5)

    def test_no_intervening_genes(self):
        hit =  self.missed_genes_grouped["Os01g58037"][0]
        b_hit = ('1',10267337,10267500, 'example')
        booln = no_intervening_genes(hit,b_hit,self.old_bed)
        self.assertEqual(booln,True)
        ### confrim that if within an intervinging gene does not find it

if __name__ == '__main__':
    unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSequenceFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)
