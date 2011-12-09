"""
This module provides a class BlastTree() which is a wrapper
over Rtree. It allows a full organismal blast to be saved into
a searchable structure. It assumes either that the blast is 
chromosomal sequence or that it is gene sequence--in which
case the qgff, and sgff files are used to update the local 
blast postions to the global positions.
It can also handle mixed, where the query is chromosomal positions
and the subject is gene positions or vice-versa. DO NOT specify
the gff file if the positions are already chromsomal as this
will give false positons.
"""

from cblastline import BlastLine
from rtree import Rtree
from gff_reader import GFFDict


class BlastTree(dict):
    def __init__(self, blast_file, qgff=None, sgff=None):
        self.blast_file = blast_file
        if qgff:
            qgff = GFFDict(qgff)
        if sgff:
            sgff = GFFDict(sgff)
        self.qgff = qgff
        self.sgff = sgff

        self.tree = blast_to_tree(blast_file, self.qgff, self.sgff)

    def __getitem__(self, k):
        return self.tree[k]
    def keys():
        return self.tree.keys()
    def iterkeys():
        return self.tree.iterkeys()
    def values():
        return self.tree.values()
    def itervalues():
        return self.tree.itervalues()
    def contains(self, k):
        return k in self.tree
    __contains__ = contains


    def find(self, qseqid, sseqid, bounds):
        """specifiy the query and subject chromsomes and bounds
        in the form: (xmin, ymin, xmax, ymax)
        """
        t = self.tree[qseqid][sseqid]
        idxs = t['tree'].intersection(bounds)
        if len(idxs) == 0: return idxs
        return [t['blasts'][i] for i in idxs]


def blast_to_tree(blast_file, qgff=None, sgff=None):
    """
    create a series of rtree's from a blast file
    the gff files are used to adjust the blast positions
    from local to global.
    if the blast position are chromosome based, the gff files
    should not be specified.
    """
    if qgff and isinstance(qgff, str):
        qgff = GFFDict(qgff)
    if sgff and isinstance(sgff, str):
        sgff = GFFDict(sgff)

    r = {}
    counts = {}

    for line in open(blast_file):
        b = BlastLine(line)

        if not b.query in r:
            r[b.query] = {}
            counts[b.query] = {}
        if not b.subject in r[b.query]:
            r[b.query][b.subject] = {'blasts':[], 'tree': Rtree()}
            counts[b.query][b.subject] = 0
        d = r[b.query][b.subject]
        d['blasts'].append(b)

        smin, smax = b.sstart, b.sstop
        if smin > smax: smax, smin = smin, smax

        qmin, qmax = b.qstart, b.qstop
        query = b.query
        subject = b.subject

        if sgff is not None:
            s = sgff[b.subject]
            smin += s.start
            smax += s.start
            subject = s.seqid
        if qgff is not None:
            q = qgff[b.query]
            qmin += q.start
            qmax += q.start
            query = q.seqid

        #r[b.query][b.subject]['tree'].add(counts[b.query][b.subject], (b.qstart, smin, b.qstop, smax))
        d['tree'].add(counts[query][subject], (qmin, smin, qmax, smax))
        counts[query][subject] += 1
    return r


