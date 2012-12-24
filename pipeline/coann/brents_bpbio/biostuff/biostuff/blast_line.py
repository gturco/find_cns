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
        #return "\t".join(map(str, map(g, BlastLine.__slots__)))
