from pyfasta import Fasta
import os
import sys
from flatfeature import Flat, Bed
import numpy as np

import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(os.path.basename(__file__))


def is_up_to_date(m, g, f):
    """
    >>> is_up_to_date('a.1.t', 'a.2.t', 'a.3.t')
    False
    >>> is_up_to_date(__file__, __file__, __file__)
    False
    """
    ospe = os.path.exists
    if not ospe(m): return False
    if not (ospe(g) and ospe(f)): return False
    mtime = os.stat(m).st_mtime
    return mtime > os.stat(g).st_mtime \
            and mtime > os.stat(f).st_mtime

def get_fastas(fasta_file, genes=False):
    """
    >>> get_fastas('a.fasta')
    ('a.genomic.masked.fasta', 'a.features.fasta')

    >>> get_fastas('a.fasta', genes=True)
    ('a.genomic.masked.genes.fasta', False)
    """

    if not genes:
        genomic_masked = fasta_file.replace('.fa', '.genomic.masked.fa')
        features_fasta = fasta_file.replace('.fa', '.features.fa')
    else:
        genomic_masked = fasta_file.replace('.fa', '.genomic.masked.genes.fa')
        features_fasta = False #fasta_file.replace('.fa', '.genes.fa')

    return genomic_masked, features_fasta

def write_feature(fh, feat, fasta, seen_ids):
    r"""
    >>> from cStringIO import StringIO
    >>> fh = StringIO()
    >>> feat =  gtpym.FeatureNode('chr2', 'CDS', 2, 4, '+')
    >>> feat.add_attribute('ID', 'some_id')
    >>> fasta = Fasta('tests/data/t.fasta')
    >>> write_feature(fh, feat, fasta, {}, 'CDS')
    >>> fh.getvalue()
    '>some_id\nCTC\n'

    """
    name = feat['accn']
    if name in seen_ids:
        print >>sys.stderr, "WARNING:", name, "used more than 1X:", feat
        # probably it appeared on + and - strand.
        return
    seen_ids[name] = None
    fh.write('>' + name + '\n')
    fdict = {'chr': feat['seqid'], 'start': feat['start'], 'stop': feat['end'], 'strand': feat['strand'] }
    fh.write(fasta.sequence(fdict) + "\n")


def check_exists(f, raw_input=raw_input):
    """
    >>> raw_input = lambda a: 'yes'
    >>> check_exists('a.ttttttttt')
    >>> check_exists(__file__, raw_input=raw_input)
    """
    if f and os.path.exists(f):
        r = raw_input("%s exists, do you wish to overwrite? (y/(n)) " % (f,))
        if r.lower()[0] != 'y':
            raise Exception("wont overwrite existing file without permission")

class NotMaskedException(Exception): pass

def mask_seq_with_locs(sequence, locs, N):
    """
    >>> s = np.array('AAAAAAAAAAAAAAAAAAAAAA', dtype='c')
    >>> N = np.array('N', dtype='|S1')

    >>> s.tostring()
    'AAANNNNNAAAAAAAAAAAAAA'

    """
    n_masked = 0
    for start, end in locs:
        sequence[start - 1: end] = N
        n_masked += end - start + 1
    return n_masked

#TODO: use more flatfeatures stuff in here.
def main(flat_file, fasta_file, inverse=False):
    genomic_masked, features_fasta = get_fastas(fasta_file)

    if is_up_to_date(genomic_masked, flat_file, fasta_file) \
            and is_up_to_date(features_fasta, flat_file, fasta_file):
        log.debug("%s is up-to-date." % (genomic_masked, ))
        return False

    Klass = Flat if flat_file.endswith(".flat") else Bed
    flat = Klass(flat_file, fasta_file)
    fasta = Fasta(fasta_file)


    N = np.array('@' if inverse else 'N', dtype='S1').astype('c')
    if inverse:
        genomic_masked = genomic_masked.replace('genomic.masked', 'genomic.nonfeat.masked')

    seen_ids = {}
    for f in (genomic_masked, features_fasta):
        check_exists(f)

    try:
        genomic_masked  = open(genomic_masked, "w")
        features_fasta = f and open(features_fasta, "w")

        fkeys = sorted(fasta.keys())
        if len(fkeys) >= 100:
            log.debug("beginning masking for %i sequences" % (len(fkeys,)))

        for achr in fkeys:
            features = flat[flat['seqid'] == achr]
            sequence = np.array(fasta[achr])

            tot_masked = 0
            for feat in features:
                try:
                    tot_masked += mask_seq_with_locs(sequence, feat['locs'], N)
                except NotMaskedException:
                    print >>sys.stderr, feat
                    cleanup(genomic_masked, features_fasta)
                    raise NotMaskedException("the sequence %s was not masked" \
                                             % (feat.get_attribute("ID"), ))


                if features_fasta:
                    write_feature(features_fasta, feat, fasta, seen_ids)

            if inverse:
                nseq = np.array(fasta[achr])
                nseq[sequence != '@'] = np.array('N', dtype='c')
                sequence = nseq.copy()
                assert sequence[sequence != 'N'].shape[0] > 0, (achr, sequence[:100])

            genomic_masked.write('>' + achr + "\n")
            genomic_masked.write(sequence.tostring() + "\n")
            if len(fkeys) < 100:
                log.debug("%i basepairs of %i (%.1f%%) masked for %s" % \
                        (tot_masked, sequence.shape[0],
                            100. * float(tot_masked)/sequence.shape[0], achr))
    except:
        cleanup(genomic_masked, features_fasta)
        raise


def cleanup(*files):
    """
    >>> fh = open('ttt.test', 'w')
    >>> fh.close()
    >>> cleanup(fh)
    """
    for f in files:
        if hasattr(f, 'name'): f = f.name
        if os.path.exists(f): os.unlink(f)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print >>sys.stderr, """\
   usage: %s path/to.flat path/to.fasta --inverse

        feature types can be a list of names like: CDS mRNA gene
        in which case, the first of those that is found will be used.
        to mask the sequence in the fasta file. 
        if the list is empty, the entire range of each feature is masked.
        if the last argument is '--inverse', everythign _but_ the feature types
        is masked.
                """ % sys.argv[0]
        sys.exit()

    flat = sys.argv[1].rstrip()
    fasta = sys.argv[2].rstrip()

    inverse = sys.argv[-1] == '--inverse'
    if len(sys.argv) > 3:
        mask_types = [x.strip() for x in sys.argv[3:]]
        main(flat, fasta, inverse)
    else:
        main(flat, fasta, inverse)
