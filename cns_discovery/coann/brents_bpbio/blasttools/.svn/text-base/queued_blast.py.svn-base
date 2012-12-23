"""
this will parallelize the blast on a directory of fasta files. keeps
the number of CPUs on the given machine full until all jobs are done.

run as file.py --help to see useage options.
"""

import commands
import os
import sys
import re
import pp
import glob

format_db = "/usr/bin/formatdb -p F -i %s"
# dont need to keep anything with > 80 hits as we mask at 50 anyway.
blast_command = "/usr/bin/blastall -p blastn -K 80 -i %s -d %s -e 0.001 -m 8 -o %s/%schr%s_vs_%schr%s.blast "

s = pp.Server()

def gen_command(query, subject, directory, outdir):
    seen = {}
    for q_fasta in glob.glob("%s/*%s*.fasta" % (directory, query)):
        (qchr,) = re.search("chr(\d+)", q_fasta).groups(0)
        for s_fasta in glob.glob("%s/*%s*.fasta" % (directory, subject)):
            (schr,) = re.search("chr(\d+)", s_fasta).groups(0)

            if not os.path.exists("%s.nin" % s_fasta):
                print >>sys.stderr, "formatting %s" % s_fasta
                commands.getoutput(format_db % s_fasta)

            com = blast_command % (q_fasta, s_fasta, outdir, query, qchr, subject, schr)
            if com in seen:continue
            seen[com] = 1
            yield com



def consume(command):
    return command, commands.getoutput(command)

def main(query, subject, directory="fasta", outdir="blast"):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    bp = open(os.path.join(outdir, "00blast.params"), "w")
    print >> bp, blast_command
    bp.close()

    jobs = []
    for c in gen_command(query, subject, directory, outdir):
        jobs.append(s.submit(consume, (c,), (), ("commands",)))
    for j in jobs:
        command, r = j()
        print command
        if r: print r


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    parser.add_option("-q", "--query",   dest="query",   help="the query organism (used in creating the output file names)")
    parser.add_option("-s", "--subject", dest="subject", help="the subject organism (used in creating the output file names)")
    parser.add_option("-d", "--directory", dest="dir", help="directory containing the fasta sequences", default="fasta")
    parser.add_option("-o", "--outdir", dest="outdir", help="directory to send the blasts", default="blast")

    (options, _) = parser.parse_args()
    if options.query is not None and options.subject is not None:
        print options.query, options.subject, options.dir
        main(options.query, options.subject, options.dir, options.outdir)

    else:
        parser.print_help()

