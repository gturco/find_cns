fblast = 'data/athaliana_athaliana/athaliana_athaliana.blast'
fgff = 'data/athaliana_athaliana/athaliana.gff'
forthology = 'data/athaliana_athaliana.orthology'

import sys
import gtpym
sys.path.insert(0, "scripts2/")
from make_genelist import is_ortho, parse_orthos
from biostuff import BlastLine

ortho = parse_orthos(forthology, is_same=True)
gff = gtpym.FeatureIndexMemory(fgff)

xs = []
ys = []
QSEQ = '2'
SSEQ = '4'

print >>sys.stderr, ortho[(QSEQ, SSEQ)]

for line in open(fblast):
    b = BlastLine(line)
    q = gff[b.query]
    s = gff[b.subject]
    if not (q.seqid == QSEQ and s.seqid == SSEQ): continue

    b.qstart += q.start
    b.qstop += q.start
    b.sstart += s.start
    b.sstop += s.start
    print b.to_blast_line()

    xs.append(b.qstart)
    ys.append(b.sstart)

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


f = plt.figure()
ax = f.add_subplot(1, 1, 1)

for xmin, xmax, ymin, ymax in ortho[('2', '4')]:
    r = Rectangle((xmin, ymin), width=xmax-xmin, height=ymax-ymin, ec='r', fc='none')
    ax.add_patch(r)

ax.plot(xs, ys, 'b,')
ax.set_xlabel(QSEQ)
ax.set_ylabel(SSEQ)
print >>sys.stderr, ax.get_xlim()
print >>sys.stderr, ax.get_ylim()
plt.show()

