#import pylab
import sys
from pyfasta import Fasta
from common import parse_raw_cns

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
    for cns_id, cns_dict in parse_raw_cns(cnsfile):
        cns = cns_dict

        qseq = qfasta[cns['qseqid']]
        sseq = sfasta[cns['sseqid']]


        sstart, send = sorted((cns['sstart'], cns['send']))
        qkey = (cns['qseqid'], cns['qstart'], cns['qend'])
        skey = (cns['sseqid'], cns['sstart'], cns['send'])

        assert sstart < send

        if cns['qend'] - cns['qstart'] < min_len: continue
        if send - sstart < min_len: continue


        if not (qkey in seen and skey in seen):
            print ">q__" + cns_id
            seqstr = str(qseq[cns['qstart'] - 1: cns['qend']]).replace('R', 'N').replace('W', 'N').replace('M', 'N')
            assert set(seqstr.lower()).issubset("actgnx"), ('q', 'q__' + cns_id, seqstr)
            print seqstr

            print ">s__" + cns_id
            seqstr = str(sseq[sstart - 1: send]).replace('R', 'N').replace('W', 'N').replace('M', 'N')
            assert set(seqstr.lower()).issubset("actgnx"), ('s', 's__' + cns_id, seqstr)
            print seqstr

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
    from common import get_sheet_from_date
    #cnsfile = get_sheet_from_date(options.cnsfile)
    main(options.cnsfile, options.qfasta, options.sfasta, options.qorg, options.sorg, int(options.min_len))

