from flatfeature import Bed
import operator

class LocalDups(object):
    def __init__(self,filename,bed):
        self.filename = filename
        self.bed = Bed(bed)
        self.bed.fill_dict()

    def get_order_dups(self):
        d = {}
        for line in open(self.filename):
            dupline = DupLine(line)
            dups = dupline.get_order(self.bed)
            d[dups[0]['accn']] = "P"
            for dup in dups[1:]:
                d[dup['accn']] = dups[0]['accn']
            intervening = dupline.get_interving_genes(self.bed)
            for i in intervening:
                if i in d.keys():continue
                d[i] = "I"
        self.filename.close()
        return d

    def write_ordered(self,out_fh):
        """write localdups to outfile"""
        localdup_fh = open(out_fh, "w")
        d = {}
        for line in open(self.filename):
            dupline = DupLine(line)
            dups = dupline.get_order(self.bed)
            line = "{0}\n".format("\t".join(dups))
            localdup_fh.write(line)
        localdup_fh.close()


    def get_dups(self):
        d = {}
        for  line in open(self.filename):
            dupline = DupLine(line)
            d[dupline.parent] = 'P'
            for dup in dupline.children:
                d[dup] = dupline.parent
            intervening = dupline.get_interving_genes(self.bed)
            for i in intervening:
                if i in d.keys(): continue
                d[i] = "I"
        self.filename.close()
        return d

class DupLine(object):
    def __init__(self,line):
        args = line.strip().split("\t")
        self.parent = args[0]
        self.children = list(set(args[1:]))
        self.dups = list(set([self.parent] + self.children))

    def __str__(self):
        return "\t".join(self.dups)

    def get_interving_genes(self,bed):
        sorted_dups = self.get_order(bed)
        dup_start = sorted_dups[0]
        dup_end = sorted_dups[-1]
        intervening = bed.get_features_in_region(str(dup_start['seqid']),
                int(dup_start['start']), int(dup_end['end']))
        intervening_dups = [i['accn'] for i in intervening if i not in self.dups]
        return intervening_dups


    def get_order(self,bed):
        dups = [bed.row_to_dict(bed.d[dup]) for dup in self.dups]
        dups.sort(key=operator.itemgetter('start'))
        return dups

