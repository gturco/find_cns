import datetime
import os.path as op
PATH=op.dirname(__file__)
import collections


class CNS(object):
    def __init__(self, cns_key,cns_info):
        self.qseqid, self.qaccn, self.sseqid, self.saccn = cns_key
        self.qstart, self.qstop, self.sstart, self.sstop = map(int,cns_info[:4])
        self.evalue = cns_info[-1]
        self.cns_id = "{0}|{1}|{2}|{3}|{4}|{5}|{6}".format(self.qseqid,self.qstart,self.qstop,self.sseqid,self.sstart,self.sstop,self.evalue)

    @classmethod
    def parse_raw_line(self,lines):
        for line in open(lines):
            if line[0] == "#": continue
            cnss = CNS.get_cns_info(line)
            for cns in cnss:
                yield cns

    @classmethod
    def get_cns_info(self,line):
        args = line.rstrip().split(",")
        cns_key = args[:4]
        cnss = args[4:]
        cns = [CNS(cns_key,cnss[cns:cns + 5]) for cns in range(0,len(cnss),5)]
        return cns
    
    def to_dict(self):
         return {'qseqid':self.qseqid,'qaccn':self.qaccn,'sseqid':self.sseqid,'saccn':self.saccn,
             'qstart':self.qstart,'qend':self.qstop,'sstart':self.sstart,'send':self.sstop,'eval':self.evalue}

class BlastLine(object):
    __slots__ = ('query', 'subject', 'pctid', 'hitlen', 'nmismatch', 'ngaps', \
                 'qstart', 'qstop', 'sstart', 'sstop', 'eval', 'score')

    def __init__(self, sline):
        args = sline.split("\t")
        self.query  =args[0]
        self.subject  = args[1]
        self.pctid =float(args[2])
        self.hitlen =int(args[3])
        self.nmismatch =int(args[4])
        self.ngaps =int(args[5])
        self.qstart =int(args[6])
        self.qstop =int(args[7])
        self.sstart =int(args[8])
        self.sstop =int(args[9])
        self.eval =float(args[10])
        self.score =float(args[11])

    def __repr__(self):
        return "BLine('%s' to '%s', eval=%.3f, score=%.1f)" % (self.query, self.subject, self.eval, self.score)


    def to_blast_line(self):
        #def g(attr):
        #    return getattr(self, attr)
        return "\t".join(map(str, (getattr(self, attr) for attr in BlastLine.__slots__)))

def parse_at_description(desc_file):
    desc = collections.defaultdict(str)
    for line in open(desc_file):
        line = line.split("\t")
        name = line[0][:line[0].find(".")]
        if name in desc:
            desc[name] += ";;" + line[-1].rstrip()
        else:
            desc[name] = line[-1].rstrip()
    return desc

def read_cns_to_rna(outdir):
    d = {}
    for line in open(outdir + "/cns_to_rna.csv"):
        hash, evalue, desc = line.rstrip("\n").split("\t")
        d[hash] = evalue + "|" +  desc
    return d

def read_cns_to_protein_exons(outdir):
    d = {}
    for sline in open(outdir + "/cns_to_protein_exons.csv"):
        hash, accns = sline.rstrip("\n").split("\t")
        d[hash] = "#".join(accns.split(","))
    return d
