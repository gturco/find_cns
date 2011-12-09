cdef extern from *:
    ctypedef char* const_char_star "const char*"

cimport stdlib

cdef extern from "Python.h":
    char *PyString_AsString(object)


cdef const_char_star blast_format = "%s\t%s\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\t%f"
cdef const_char_star blast_format_line = "%s\t%s\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\t%f\n"

cdef extern from "stdio.h":
    ctypedef struct FILE:
        pass
    FILE * fopen(char *, char *)
    int fscanf(FILE * i, char * fmt, ...)
    void rewind(FILE * f)
    cdef int fclose(FILE *)
    cdef Py_ssize_t strlen(char *)
    int EOF

    int sscanf(char* astr, const_char_star format, ...)
    char *fgets(char *line, int maxline, FILE *fp)



cdef class BlastFile:
    cdef FILE* fh
    cdef object filename
    def __cinit__(self, char* filename):
        self.fh = fopen(filename, 'r')
        self.filename = filename

    def __iter__(self):
        rewind(self.fh)
        return self
        
    def __next__(self):
        cdef float pct = 0.0, evalue = 0.0, bit = 0.0
        cdef char qname[512], sname[512]
        cdef int hlen, nmiss, ngap, qstart, qstop, sstart, sstop
        cdef char *tmp
        cdef int success
        success = fscanf(self.fh, blast_format_line, qname, sname, \
                         &pct, &hlen, &nmiss, &ngap, &qstart, &qstop,\
                         &sstart, &sstop, &evalue, &bit )
        if success == EOF:
            raise StopIteration
        return create_blast_line(qname, sname, pct, hlen, nmiss, ngap,
                        qstart, qstop, sstart, sstop, evalue, bit)

    def __dealloc__(self):
        fclose(self.fh)

    def __repr__(self):
        return "BlastFile('%s')" % (self.filename, )

import sys

cdef class BlastLine:
    r"""
    given a string of tab-delimited (-m 8) blast output, parse it and create
    an object with the usual attrs:

    >>> from biostuff import BlastLine 
    >>> b = BlastLine("Os09g11510	Os08g13650	92.31	39	3	0	2273	2311	3237	3199	0.001	54.0")
    >>> b.query
    'Os09g11510'

    >>> attrs = ('query', 'subject', 'pctid', 'hitlen', 'nmismatch', 'ngaps', \
    ...  'qstart', 'qstop', 'sstart', 'sstop', 'evalue', 'score')
    >>> [getattr(b, attr) for attr in attrs]# doctest: +ELLIPSIS
    ['Os09g11510', 'Os08g13650', 92.3..., 39, 3, 0, 2273, 2311, 3237, 3199, 0.001..., 54.0]
        
    """
    cdef public int hitlen, nmismatch, ngaps, qstart, qstop, sstart, sstop
    cdef public float pctid, score
    cdef public double evalue
    cdef char _cquery[128], _csubject[128]
    cdef object _pysubject, _pyquery
    attrs = ('query', 'subject', 'pctid', 'hitlen', 'nmismatch', 'ngaps', \
                 'qstart', 'qstop', 'sstart', 'sstop', 'evalue', 'score')

    property query:
        def __get__(self):
            if self._pyquery is None:
                return self._cquery
            return self._pyquery
        def __set__(self, val):
            self._pyquery = val

    property subject:
        def __get__(self):
            if self._pysubject is None:
                return self._csubject
            return self._pysubject
        def __set__(self, val):
            self._pysubject = val

    def __init__(self, char *sline):
        if sline != NULL:
            sscanf(sline, blast_format, self._cquery, self._csubject,
                &self.pctid, &self.hitlen, &self.nmismatch, &self.ngaps,
                &self.qstart, &self.qstop,
                &self.sstart, &self.sstop,
                &self.evalue, &self.score)

    def __richcmp__(BlastLine self, BlastLine other, size_t op):
        if op == 2: # ==
            if self.query != other.query and self.qstart != other.qstart: 
                return False
            return self.subject == other.subject and \
                    self.qstop == other.qstop and \
                    self.sstop == other.sstop and \
                    self.evalue == other.evalue and \
                    self.hitlen == other.hitlen

        elif op == 3: # !=
            return not self.__richcmp__(other, 2)
        else:
            raise Exception("that comparison not implemented")


    def __repr__(self):
        return ("BlastLine(%s[%i:%i]-%s[%i:%i], " +\
                   "ptcid=%.3f, eval=%.3f)") \
            % (self.query, self.qstart, self.qstop, self.subject, self.sstart,\
               self.sstop, self.pctid, 
               self.evalue)

    """
    @classmethod
    def yield_lines(cls, filename):
        cdef float pct = 0.0, evalue = 0.0, bit = 0.0
        cdef char qname[512], sname[512]
        cdef char *tmp
        cdef FILE *fh
        fh = fopen(filename, "r")
        while fscanf(fh, blast_format_line, qname, sname, &pct, &hlen, &nmiss, &ngap, &qstart, &qstop, &sstart, &sstop, &evalue, &bit ) != EOF:
            yield create_blast_line(qname, sname, ptc, hlen, nmiss, ngap,
                                    qstart, qstop, sstart, sstop, evalue, bit)
    """
    def to_blast_line(self, as_str=True):
        if as_str:
            s = [getattr(self, attr) for attr in BlastLine.attrs]
            s[2] = "%.2f" % (s[2], ) # pctid
            s[-1] = "%.1f" % (s[-1], ) # bitscore
            return "\t".join(map(str, s))
        else:
            return [getattr(self, attr) for attr in BlastLine.attrs]


    def __reduce__(self):
        return create_blast_line, (
            self.query, self.subject, self.pctid, self.hitlen, self.nmismatch,
            self.ngaps, self.qstart, self.qstop, self.sstart, self.sstop, 
            self.evalue, self.score)

    #def __getstate__(self):
        #    return tuple([getattr(self, k) for k in BlastLine.attrs])


cdef extern from "pnew.h":
        cdef BlastLine NEW_BLASTLINE "PY_NEW" (object t)

cpdef BlastLine create_blast_line(char *query, char*subject, float pctid, int hitlen, 
                       int nmismatch, int ngaps, int qstart, int qstop, 
                       int sstart, int sstop, float evalue, float score):
    cdef BlastLine b = NEW_BLASTLINE(BlastLine)
    b.query = query
    b.subject = subject
    b.pctid = pctid
    b.hitlen = hitlen
    b.nmismatch = nmismatch
    b.ngaps = ngaps
    b.qstart = qstart
    b.qstop = qstop
    b.sstart = sstart
    b.sstop = sstop
    b.evalue = evalue
    b.score = score
    return b
