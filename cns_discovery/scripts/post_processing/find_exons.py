import sys
import os
sys.path.insert(0, os.path.dirname(__file__))
from cns_utils import BlastLine, CNS
import collections
"""
Each CNS was blastx's to all arabidopsis protein, and the hit's with
an e-val < 0.01 were listed.  The score , CNS length, % identity of
hit, hit length, CNS length/ hit length (protein x 3)= coverge,
definition (organism, gene name), sequence, etc.

CNS= exon is any hit has a score of of > or =50.0  (this e-val is
always < 0.0001 and often more significant.)

OR

score > or =45.0 AND
a coverge of > or = 0.90.  Coverge= length hit in aa x 3/ CNS length.
(this is usually e-val < 0.001).

Almost all CNSs in readout are bonifide exons.  Even the ones
excluded have a good chance of not being exons.  These are clearly
the big ones.
"""


def main(blast_files, out_dir,raw_cns):
    """empty docstring"""
    cns_by_id = {}  
    for cns in CNS.parse_raw_line(raw_cns):
	cns_by_id[cns.cns_id] = cns

    exons = collections.defaultdict(dict)
    for blast_file in blast_files:
        for line in open(blast_file):
            b = BlastLine(line)
            # chop the q__ and s__
            key = b.query[3:]
            #assert key == b.subject[3:], (key, b.subject[3:])
            # convert piped rice names to short canonical names.
            subject = b.subject.split("|")[0] if "|" in b.subject else b.subject
            # chop At2g26540.1 to At2g26540
            subject = subject[:-2] if subject[-2] == "." else subject
            subject = subject.replace('LOC_', '')
           
            if b.score > 50:
                if not subject in exons[key]:
                    exons[key][subject] = [b.eval]
                else:
                    exons[key][subject].append(b.eval)
                continue

            if b.score < 45: continue
            cns = cns_by_id[key]

            # qstart?
            qlen = cns.qstop - cns.qstart
            coverage = (b.hitlen * 3.) / qlen
            #print >>sys.stderr, coverage
            if coverage < 0.90: continue
            if not subject in exons[key]:
                exons[key][subject] = [b.eval]
            else:
                exons[key][subject].append(b.eval)


    exons = dict(exons)
    write_exons(exons, out_dir)
    #for cns_hash, at_exons in exons.iteritems():
    print >>sys.stderr, "%i total unique cnss are exons" % (len(exons), )
    return exons

def write_exons(exons, out_dir):
    of = out_dir + "/cns_to_protein_exons.csv"
    print >>sys.stderr, "writing exons to %s" %of
    fh = open(of, "wb")
    for exon, subjects in exons.iteritems():
        line_str = []
        for s, evalues in subjects.iteritems():
            evalues.sort()
            line_str.append("%s|%.3g" % (s, evalues[0]))
        print >>fh, "%s\t%s" % (exon, ",".join(line_str))


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser("%prog [options] [blast_files]")
    parser.add_option("-o", "--out_dir", dest="out_dir", help="")
    parser.add_option("--cns", dest="cns", help=" raw cns file cns.txt or cns.txt.local")
    (opts, blast_files) = parser.parse_args()

    main(blast_files, opts.out_dir, opts.cns)
