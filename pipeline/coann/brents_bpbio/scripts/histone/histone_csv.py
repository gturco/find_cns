"""
given the binary files created with to_bin.py and a .gff
create a csv by gene. (goes to stdout when this script
is run. adjust the paths at the bottom before running. requires 
flatfeature/fatfeature.py in bpbio repository
NOTE: values created with to_bin.py (and reported here) are log(1 + cy5/cy3)
"""

from fatfeature import Fat
import numpy as np
import sys
import os


def pairs_to_slice(pairs):
    """
    given a list of tuples (like a list of CDS start, stops), return
    the numpy array that will work as a slice for those tuples
    """
    return np.concatenate([np.arange(s0-1, s1) for s0, s1 in pairs])


def calc_stats(mfiles_pat, gff):
    fat = Fat(gff)

    hists = {}
    for i in range(1, 6):
        i = str(i)
        mf = mfiles_pat % i
        assert os.path.exists(mf)
        hists[i] = np.fromfile(mf, dtype=np.float32)
    header = "accn,gene,cds,intron,up10,up100,up1000,down10,down100,down1000"

    print header

    for accn, f in sorted(fat.iteritems()):
        
        data = [accn]
        if not f.seqid in hists: continue # C, G chrs
        hist = hists[f.seqid]

        for locs in ([[f.start, f.end]], 
                     getattr(f, 'CDS', None), 
                     fat.introns(f),
                     fat.upstream(f, 10, noncoding=True), 
                     fat.upstream(f, 100, noncoding=True), 
                     fat.upstream(f, 1000, noncoding=True),
                     fat.downstream(f, 10, noncoding=True), 
                     fat.downstream(f, 100, noncoding=True), 
                     fat.downstream(f, 1000, noncoding=True)
        ):
            if locs is None:
                # occurs when there's no CDS.
                data.append("na")
                continue

            slicer = pairs_to_slice(locs)
            try:
                m = hist[slicer] # this context.
            except IndexError: 
                # difference between fasta and features due to version
                slicer = slicer[slicer < hist.shape[0]]
                m = hist[slicer]

            data.append("%.5f" % m.mean())
    
        print ",".join(data)



if __name__ == "__main__":

    # adjust these accordingly.
    histpat = "/tmp/GSE/histone.%s.bin"
    gff = "data/thaliana_v9_genes.gff"

    calc_stats(histpat, gff)
