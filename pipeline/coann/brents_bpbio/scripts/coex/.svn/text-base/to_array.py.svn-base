import numpy as np
import os
import sys

def get_names(directory):
    f = os.listdir(directory)[0]
    accns = [f.upper()] + [line.split("\t")[0].upper() for line in open(os.path.join(directory, f))]
    accns = sorted(a.replace("_", "") for a in accns)
    return dict((accn, i) for i, accn in enumerate(accns))

def fill_array(directory, names2idx, arr):
    L = len(names2idx)
    for f in sorted(os.listdir(directory)):
        i = names2idx[f.upper().replace("_", "")]
        a = arr[i]
        for line in open(os.path.join(directory, f)):
            accn, mr, pcc = line.split("\t")
            j = names2idx[accn.upper().replace("_", "")]
            #a[j] = mr
            a[j] = 1 - float(mr) / L



if __name__ == "__main__":

    patt = "data/Ath.coex.c4-1.g%i"

    names2idx = get_names(patt % 1)
    print >>sys.stderr, "creating"
    arr = np.zeros((len(names2idx), len(names2idx)), dtype='float32')
    print >>sys.stderr, "created"
    for seqid in range(1, 6):
        print >>sys.stderr, seqid
        fill_array(patt % seqid, names2idx, arr)
    arr.tofile('coexp.bin')

    order = open('names.order', 'w')
    for a in sorted(names2idx):
        print >>order, a
