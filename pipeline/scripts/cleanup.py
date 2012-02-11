
class LocalDups(object):
    def __init__(self):
        args = line.strip().split("\t")
        self.parent = args[0]
        self.children = args[1:]
        self.dups = parent + children

    def __str__(self):
        return "\t".join(self.dups)

    def get_interving_genes(self,bed):
        sorted_dups = self.get_order(bed)
        dup_start = all[0]
        dup_end = all[-1]
        intervening = flat.get_features_in_region(dup_start['seqid'],
                dup_start['start'], dup_end['end'])
        intervening_dups = [i['accn'] for i in intervening if i not in self.dups]
        return intervening_dups


    def get_order(self,bed):
        all = [bed.accn(dup) for dup in self.dups]
        all.sort(key=operator.itemgetter('start'))
        return all

