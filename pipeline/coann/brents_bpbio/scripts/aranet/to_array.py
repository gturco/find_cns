import numpy as np

def get_names(f):
    names = []
    for line in open(f):
        line = line.split("\t")
        a, b = map(str.upper, line[:2])
        names.append(a)
        names.append(b)

    names = sorted(frozenset(names))
    return dict((accn, i) for i, accn in enumerate(names))


def fill_array(f, names2idx):
    L = len(names2idx)
    arr = np.zeros((len(names2idx), len(names2idx)), dtype='float32')
    for line in open(f):
        line = line.split("\t")
        a, b = map(str.upper, line[:2])
        i = names2idx[a]
        j = names2idx[b]
        arr[i, j] = float(line[-1])
    return arr


if __name__ == "__main__":

    f = "AraNet.v1.join.txt"
    names2idx = get_names(f)

    arr = fill_array(f, names2idx)
    arr.tofile('aranet.bin')
    order = open('names.order', 'w')
    for a in sorted(names2idx):
        print >>order, a
