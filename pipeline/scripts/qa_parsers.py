from bed_utils import Bed

class ParsePairs(object):
    """takes line from pairs returns the qaccn and saccn"""
    ##TODO add dag and cluster

    def __init__(self,line):
        self.line = line.strip().split('\t')

    def dag(self):
        assert len(self.line) > 5, self.line
        pairs = self.line[1], self.line[5]
        return pairs

    def cluster(self):
        pass
        #assert len(self.line) == 5, self.line
        #pair = self.line[1], self.line[3]


    def qa(self):
        pairs = self.raw()
        return pairs
        #### add in parsing part

    def raw(self,qnolocaldups,snolocaldups):
        pair_pos = self.cluster()
        pair_pos = map(int,pair_pos)
        qoder = get_order(qnolocaldups)
        sorder = get_order(snolocaldups)
        pairs = qoder[pair_pos[0]], sorder[pair_pos[1]]
        return pairs

    def pair(self):
        if len(self.line) == 1:
            line = self.line[0].split(",")
        assert len(self.line) >= 2, "dont know how to handle %s" % self.line
        pair = self.line[0], self.line[1]
        return pair

