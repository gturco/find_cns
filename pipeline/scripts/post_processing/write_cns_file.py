import os.path as op
import sys
sys.path.insert(0, op.dirname(__file__))
from common import get_sheet_from_date

def read_cns_to_outgroup_synteny_score(outdir):
    d = {}
    for sline in open(outdir + "/cns_to_outgroup_synteny_score.csv"):
        if sline[0] == '#': continue
        hash, oseqid, oloc, synteny_score = sline.rstrip("\n").split(" ")
        d[hash] = "%s:%i:%.4f" % (oseqid, int(oloc), float(synteny_score))
    return d


def main(datasheet_path, orthos):
    """empty docstring"""
    columns = ["outgroup", "orthologous"]
    of = op.splitext(datasheet_path)[0] + ".attrs.csv"
    print >>sys.stderr, "writing attrs to %s" % of
    fh = open(of, "w")
    outdir = op.dirname(datasheet_path)
    outgroup = read_cns_to_outgroup_synteny_score(outdir)

    in_fh = open(datasheet_path, "r")
    headers = in_fh.readline().strip().strip("#").split(",") + columns

    print >>fh, ",".join(headers)
    count = 0
    while True:
        line = in_fh.readline().strip()
        if not line: break
        d = dict(zip(headers, line.split(",")))
        cnsid = d['cns_id']
        og = outgroup.get(cnsid, "0")
        ortho = 'T' if (d['qaccn'], d['saccn']) in orthos else 'F'
        count += 1 if ortho == 'T' else 0
        print >>fh, line + ",%s,%s" % (og, ortho)
    print "orthologous:", count

def read_orthos(fortho):
    fh = open(fortho)
    header = fh.readline().rstrip().split(",")
    orthos = {}
    while True:
        line = fh.readline().strip()
        if not line: break
        if line[0] == "#": continue
        q,s = line.split(",")[:2]
        orthos[(q, s)] = True
    return orthos

if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("-o", "--orthology", dest="orthology", help="path to orthologies")
    parser.add_option("-d", "--datasheet", dest="datasheet", help="cns datasheet")
    (opts, _) = parser.parse_args()

    if not (opts.datasheet and opts.orthology):
        sys.exit(parser.print_help())

    orthos = read_orthos(opts.orthology)

    dpath = get_sheet_from_date(opts.datasheet)
    main(dpath, orthos)

