class CNS(object):
    def __init__(self, cns_key,cns_info):
        self.qseqid, self.qaccn, self.sseqid, self.saccn = cns_key
        self.qstart, self.qstop, self.sstart, self.sstop, self.evalue = cns_info
        self.cns_id = "{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}".format(self.qseqid,self.qaccn,
                self.sseqid,self.saccn,self.qstart,self.qstop,self.sstart,self.sstop)

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




