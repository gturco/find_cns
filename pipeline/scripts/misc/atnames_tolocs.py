import sys

import gtpym
# python data/athaliana_athaliana.genes.orthology data/athaliana_athaliana/athaliana.gff
names = sys.argv[1]
gff = sys.argv[2]

def get_loc(g, key):
    firstkey = "".join(list(key))
    while True:
        try:
            return g[key]
        except KeyError:
            i = int(key[-1])
            if i == 0:
                key = key[:-2] + str(int(key[-2]) - 1) + '9'
            else:
                key = (key[:-1] + str(i - 1))
            assert not '-' in key, (firstkey, key)

PAD = 2000

g = gtpym.FeatureIndexMemory(gff)
for line in open(names):
    assert not '-' in line, line
    line = [x.strip('-').upper() for x in line.rstrip('-\n\r').split(",")]
    alpha, qstart, qstop, sstart, sstop = line
    qgstart = get_loc(g, qstart)
    qgstop = get_loc(g, qstop)
    sgstart = get_loc(g, sstart)
    sgstop = get_loc(g, sstop)
    assert qgstart.seqid == qgstop.seqid
    assert sgstart.seqid == sgstop.seqid, (sgstart.seqid, sgstop.seqid, line)
    nline = [qgstart.seqid, sgstart.seqid, qgstart.start - PAD, qgstop.end + PAD,
                                sgstart.start - PAD, sgstop.end + PAD]
    print ",".join(map(str, nline))
