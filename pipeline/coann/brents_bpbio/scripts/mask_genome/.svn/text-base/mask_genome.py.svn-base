"""
take a genomic `fasta` and the `blast` output which resulted from the self-self
blast of that fasta. and generate a new fasta in which all basepairs occuring
more then `cutoff` times are masked. e.g.:
    $ python mask_genome.py -b rice_rice.blast -o rice -f rice.fasta -c 50
will create a new fasta file rice.masked.50.fasta
"""
import numpy as np
from pyfasta import Fasta
import tables
import numexpr
import sys, os

def cache_clear(cache, node, qchr, schr):
    """ keep some of the (hdf5) arrays in memory
    if there are too many, put one back into
    the h5 file
    """
    if len(cache) < 50: return
    #print "updating cache: %s, %s" % (qchr, schr)
    rmkey = [k for k in cache if not k in (qchr, schr)][0]
    getattr(node, 'c' + rmkey)[:] = cache[rmkey]
    del cache[rmkey]

def update_cache(achr, node, clen, h5, cache):
    """
    put a new (hdf5) array in memory, and make sure it's
    in the h5 as well
    """
    if not 'c' + achr in node:
        h5.createArray(node, 'c' + achr, np.zeros((clen,), dtype=np.uint32))

    cache[achr] = getattr(node, 'c' + achr)[:]

H5 = 'copy_count.h5'
def get_node(org, mode):
    # get the parent group node in the h5 file.
    if mode == 'w':
        h5 = tables.openFile(H5, mode='a')
        if org in h5.root:
            action = raw_input(\
               """%s copy counts exist in %s. what to do [d/a/u]?
                   'd': delete them and create new copy-counts
                   'a': abort
                   'u': use the existing copy-counts
                you can use the existing counts if the blast is unchanged."""
                               % (org, H5))[0].lower()
            if action == 'd':
                getattr(h5.root, org)._f_remove(recursive=True)
                h5.flush()
            elif action == 'u':
                h5.close()
                return None, None
            else:
                print('ABORT: %s already exists in %s' % (org, H5))
                h5.close(); sys.exit()
        return h5, h5.createGroup(h5.root, org, org)
    else:
        h5 = tables.openFile(H5, mode='r')
        return h5, getattr(h5.root, org)

def count_freq(blast_file, fasta, org, count_subject=True):
    """one large blast file """
    h5, node = get_node(org, 'w')

    # use existing counts.
    if (h5, node) == (None, None): return
    f = Fasta(fasta)

    print "counting..."
    cache = {}
    for sline in open(blast_file):
        line = sline.split("\t")
        qchr, schr = line[:2]

        qstart, qstop, sstart, sstop = map(int, line[6:10])

        if not qchr in cache:
            update_cache(qchr, node, len(f[qchr]), h5, cache)
            cache_clear(cache, node, qchr, schr)
        # convert to 0-based indexes:
        # 1 8 => 0 8, but range doesnt include upper boud.
        cache[qchr][qstart - 1: qstop] += 1

        if count_subject:
            if sstart > sstop: sstart, sstop = sstop, sstart
            if not schr in cache:
                update_cache(schr, node, len(f[schr]), h5, cache)
                cache_clear(cache, node, qchr, schr)
                cache[schr][sstart - 1: sstop] += 1


    for achr in cache:
        getattr(node, 'c' + achr)[:] = cache[achr]

    h5.close()

def mask(fasta_file, org, cutoff, mask_value='X'):
    h5, node = get_node(org, 'r')

    outfile = fasta_file[:fasta_file.rfind(".")] + (".masked.%i" % cutoff) \
                         + fasta_file[fasta_file.rfind("."):]

    print "> masking sequence to file:", outfile
    out = open(outfile ,'w')

    fasta = Fasta(fasta_file)

    soft_mask = mask_value.lower() == 'soft'
    for seqid in sorted(fasta.iterkeys()): 
        masked = 0
        if soft_mask:
            seq = str(fasta[seqid])
            # mask is the lowercase sequence.
            mask_value = np.array(seq.lower(), dtype='c')
            seq = np.array(seq.upper(), dtype='c')
        else:
            fasta[seqid].tostring = False
            seq = fasta[seqid][:] # a


        if not 'c' + seqid in node:
            print >>sys.stderr, seqid,\
                '! not found in masked, writing unchanged\n' \
                '  this means that no section of this sequence appeared\n' \
                '  more than %i times' % cutoff
            out.write('>' + seqid + '\n')
            out.write(seq.tostring() + '\n')
            continue
        
        hit_counts = getattr(node, 'c' + seqid)[:]
        masked_seq = np.where(numexpr.evaluate("hit_counts > %i" % cutoff)
                              , mask_value, seq).tostring() 

        l = len(masked_seq)
        print >>sys.stderr, "! seq:%s len:%i %%masked:%.3f" % (seqid, l, 
                                   100.0 * masked_seq.count(mask_value) / l)
        assert len(seq) == l
        out.write('>' + seqid + '\n')
        out.write(masked_seq + '\n')

    out.close()
    # write out a file .fasta.version containing
    # the svnversion (if available of this script
    # that was used to create the file.
    path = os.path.dirname(__file__)
    os.system('svnversion %s > %s.version' % (path, outfile))
    h5.close()

if __name__ == "__main__":
    import optparse
    p = optparse.OptionParser(__doc__)
    p.add_option("-b", dest="blast", help="path to self-self genomic blast")
    p.add_option("-f", dest="fasta", help="path to the fasta file (which was\n"
                 "used as query and subject in the blast")
    p.add_option("-o", dest="org", help="name of the organism. e.g. 'rice'")
    p.add_option("--h5", dest="h5", help="path to the hdf5 file to use.", default=H5)
    p.add_option("-c", dest="cutoff", 
                 help="cutoff value, bp locations appearing in the blast more\n"
                 "than this many times are masked", type='int', default=50)
    p.add_option("-S", dest="no_subject", default=False, action="store_true",
             help="do NOT count subject hits (only query) when tabulating")
    p.add_option("-m", dest="mask", help=\
         "mask sequence with this letter (usually 'X' or 'N'). if == 'SOFT',"
         "then the sequence subjected to soft-masking where all repetitive"
         "values are lower-cased. and all other sequence is upper-cased"
         "regardless of its case in the original fasta", default='X')

    options, _ = p.parse_args()
    if not (options.blast and options.fasta and options.org):
        sys.exit(p.print_help())
    ospe = os.path.exists
    if not (ospe(options.blast) and ospe(options.fasta)):
            print "make sure blast:%s and fasta:%s exist:" \
                                % (options.blast, options.fasta)
            sys.exit()

    if options.h5: H5=options.h5

    assert len(options.mask) == 1 or options.mask.lower() == 'soft'
    count_freq(options.blast, options.fasta, options.org, not options.no_subject)
    print "> done counting..."
    mask(options.fasta, options.org, options.cutoff, options.mask)
