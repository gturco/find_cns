import datetime
import os.path as op
PATH=op.dirname(__file__)
import collections

def get_sheet_from_date(flanker_tmpl):
    d = datetime.date.today()
    one_day = datetime.timedelta(days=1)
    while not op.exists(flanker_tmpl.replace('DATE', str(d))):
        d -= one_day
    flanker_file = flanker_tmpl.replace('DATE', str(d))
    return flanker_file

def parse_cns_datasheet(cns_datasheet):
    headers = None
    for line in open(cns_datasheet):
        if headers is None:
            headers = line.lstrip("#").rstrip().split(",")
            continue
        line = line.rstrip().split(",")
        line = dict(zip(headers, line))
        for h in ('qstart', 'qstop', 'sstart', 'sstop'): line[h] = int(line[h])
        key = "%(qchr)s|%(qstart)i|%(qstop)i|%(schr)s|%(sstart)i|%(sstop)i"

        yield key % line, line

#Bd1,Bradi1g59970,1,Sb01g032360,59025115,59025142,55266333
#Bd1,Bradi1g54400,Bd1,Bradi1g05050,52780740,52780900,3350383,3350222
#Bd1,Bradi1g54690,Bd1,Bradi1g05040,53020891,53020906,3342856,3342871

def parse_raw_cns(raw_cns):
    for line in open(raw_cns):
        if line[0] == "#": continue
        line = line.rstrip().split(",")
        qseqid, qaccn, sseqid, saccn = line[:4]
        locs = map(float, line[4:])

        key = "%(qseqid)s|%(qstart)i|%(qend)i|%(sseqid)s|%(sstart)i|%(send)i|%(eval)s"
        for i in range(0, len(locs), 5):
            qstart, qend, sstart, send, evalue = locs[i:i + 5]
            d = {'qseqid': qseqid, 'sseqid':sseqid,
                   'qaccn': qaccn, 'saccn': saccn, 'qstart': int(qstart),
                   'qend': int(qend), 'sstart': int(sstart), 'send': int(send),"eval": evalue}
            yield key % d, d


def parse_orthology(ortho_file, qbed, sbed):
    orthos = {}
    for line in open(ortho_file):
        if line[0] == "#": continue
        line = line.split("\t")
        qseqid, qi, sseqid, si = line[:4]
        qi, si = int(qi), int(si)

        q, s = qbed[qi], sbed[si]
        qaccn, saccn = q['accn'], s['accn']
        orthos[(qaccn, saccn)] = {'qseqid': qseqid, 'qpos': qi,
                                  'sseqid': sseqid, 'spos': si}

        if qbed.filename == sbed.filename:
            orthos[(saccn, qaccn)] = {'qseqid': sseqid, 'qpos': si,
                                  'sseqid': qseqid, 'spos': qi}
    return orthos


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

def parse_os_description(desc_file=op.join(PATH, "../data/sativa_v6.1.description")):
    d = {}
    for line in open(desc_file):
        line = line.replace('LOC_', '').split()
        d[line[0]] = ";".join(line[1:])
    return d

def parse_at_description(desc_file=op.join(PATH, "../data/thaliana_v9.description")):
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
