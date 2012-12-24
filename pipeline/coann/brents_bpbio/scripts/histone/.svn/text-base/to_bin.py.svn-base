import sys
import math
from pyfasta import Fasta
import numpy as np

"""
NOTE: gaps in coverage are filled with the average
of the value above and the value below.
"""

hists = (
    "GSM310840_Cy5_Cy3.txt",
    "GSM310841_Cy5_Cy3.txt",
    "GSM310842_Cy5_Cy3.txt",
    "GSM310843_Cy5_Cy3.txt",
)

posns = (
    "GPL7143_2007-06-19_ATH1_chr_all_meth01.pos",
    "GPL7143_2007-06-19_ATH1_chr_all_meth02.pos",
    "GPL7143_2007-06-19_ATH1_chr_all_meth03.pos",
)


"""
SEQ_ID  CHROMOSOME  PROBE_ID    POSITION    LENGTH  COUNT
CHR1v01212004   chr1v01212004   CHR1V01212004FS008761566    8761566 46  1
CHR1v01212004   chr1v01212004   CHR1V01212004FS013192017    13192017    52  1
"""
def parse_posns(pfiles):

    lookups = {}
    for f in pfiles:
        header = None
        for line in open(f):
            line = line.rstrip().split("\t")
            if header is None:
                header = line 
                continue
            lookups[line[2]] = int(line[3]), int(line[4])
    return lookups


"""
ID_REF  Cy5 Cy3
CHR1V01212004FS000000001    39524.11    31755.55
CHR1V01212004FS000000061    3463.56 3682.89
CHR1V01212004FS000000121    1434.22 1427.44
"""
import collections
def parse_hists(fhists):

    vals_by_id = collections.defaultdict(list)
    for f in fhists:
        header = None
        for line in open(f):
            line = line.rstrip().split("\t")
            if header is None:
                header = line 
                continue
            id, cy5, cy3 = line[0], float(line[1]), float(line[2])
            vals_by_id[id].append((cy5, cy3))
    # TODO: see zscore here: for confidence in measurement...
    # http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM310842
    def avg(li): return sum(li)/len(li)
    for id, c53li in vals_by_id.iteritems():
        cy5 = avg([x[0] for x in c53li])
        cy3 = avg([x[1] for x in c53li])
        yield id, math.log(1 + (cy5 / cy3))

def fill(a):
    """
    fill zeros with the preceding value. (dont want moving avg).

    >>> a = np.zeros((12,))
    >>> a[:3] = 12
    >>> a[6:9] = 15
    >>> a[11:] = 13

    >>> fill(a)
    >>> a
    array([ 12. ,  12. ,  12. ,  13.5,  13.5,  13.5,  15. ,  15. ,  15. ,
            14. ,  14. ,  13. ])

    >>> a = np.arange(5)
    >>> a[0] = 0
    >>> fill(a)
    >>> a
    array([1, 1, 2, 3, 4])
    >>> a[-1] = 0
    >>> fill(a)
    >>> a
    array([1, 1, 2, 3, 3])
    """
    s0, = np.where((a[1:] == 0) & (a[:-1] != 0))
    s1, = np.where((a[:-1] == 0) & (a[1:] != 0))
    s0 += 1
    s1 += 1
    if len(s0) < len(s1):
        s0 = np.concatenate(([0], s0))
    if len(s1) < len(s0):
        s1 = np.concatenate((s1, [0]))
        
    idxs = np.transpose((s0, s1))

    for start, end in idxs:
        if start == 0:
            a[start:end] = a[end]
        elif end == 0:
            a[start] = a[start - 1]
        else:
            a[start:end] = (a[start - 1] + a[end]) / 2


if __name__ == "__main__":
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    plookups = parse_posns(posns)
    fasta = Fasta("/opt/src/flatfeature/data/thaliana_v9.fasta")

    arrs = {}
    counts = {}
    for seqid, seq in fasta.iteritems():
        arrs[seqid] = np.zeros((len(seq),), dtype=np.float32)
        counts[seqid] = np.zeros((len(seq),), dtype=np.float32)
    print arrs.keys()

    success = 0
    for id, val in parse_hists(hists):
        start, length = plookups[id]
        assert val != 0
        #print id, val, start, length
        seqid = id[3]
        if not seqid in "12345": continue
        a = arrs[seqid][start - 1: start + length]
        counts[seqid][start - 1: start + length] += 1
        c = counts[seqid][start -1: start + length]
    
        # keep the average. since there are overlaps, this 
        # weights by the number of existing measurements taht have
        # already contributed to the value.
        arrs[seqid][start - 1: start + length] = ((c - 1) * a + val) / c
        success += 1
    print "ADDED:", success
    for seqid, arr in arrs.iteritems():
        arr = fill(arr) 
        arr.tofile("histone.%s.bin" % seqid)
