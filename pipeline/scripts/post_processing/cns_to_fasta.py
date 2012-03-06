import sys
from pyfasta import Fasta
from cns_utils import CNS

def main(cnsfile, qfasta_file, sfasta_file, qorg, sorg, min_len):
    """empty docstring"""
    lens = []
    qfasta = Fasta(qfasta_file)
    sfasta = Fasta(sfasta_file)

    seen = {}

    lens_append = lens.append
    qseq, sseq = None, None
    # so we only read a new fasta file as needed.
    last_qchr, last_schr = None, None

    seen = {}
    for cns in CNS.parse_raw_line(cnsfile):
        qseq = qfasta[str(cns.qseqid)]
        sseq = sfasta[str(cns.sseqid)]


        sstart, send = sorted((cns.sstart, cns.sstop))
        qkey = (cns.qseqid, cns.qstart, cns.qstop)
        skey = (cns.sseqid, cns.sstart, cns.sstop)

        assert sstart < send

        if cns.qstop - cns.qstart < min_len: continue
        if send - sstart < min_len: continue


        if not (qkey in seen and skey in seen):
            print ">q__" + cns.cns_id
            seqstr = str(qseq[cns.qstart - 1: cns.qstop]).replace('R', 'N').replace('W', 'N').replace('M', 'N')
            assert set(seqstr.lower()).issubset("actgnx"), ('q', 'q__' + cns.cns_id, seqstr)
            print seqstr.upper()

            print ">s__" + cns.cns_id
            seqstr = str(sseq[sstart - 1: send]).replace('R', 'N').replace('W', 'N').replace('M', 'N')
            assert set(seqstr.lower()).issubset("actgnx"), ('s', 's__' + cns.cns_id, seqstr)
            print seqstr.upper()

        seen[qkey] = 1
        seen[skey] = 1


if __name__ == "__main__":
    import optparse
    import os
    parser = optparse.OptionParser()
    parser.add_option("-c", "--cnsfile", dest="cnsfile", help="template path to the cns.txt from find_cns.py")
    parser.add_option("--qfasta", dest="qfasta", help="query fasta with CDS masked")
    parser.add_option("--sfasta", dest="sfasta", help="subject fastas with CDS masked")
    parser.add_option("--sorg", dest="sorg", help="subject organism name")
    parser.add_option("--qorg", dest="qorg", help="query organism name")
    parser.add_option("--min_len", dest="min_len", help="skip cnss with len < than this", default=0)

    (options, _) = parser.parse_args()
    if not (options.cnsfile and options.qfasta and options.sfasta):
        sys.exit(parser.print_help())

    sys.path.insert(0, os.path.dirname(__file__))
    main(options.cnsfile, options.qfasta, options.sfasta, options.qorg, options.sorg, int(options.min_len))

