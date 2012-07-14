import sys
import os.path as op
sys.path.insert(0, op.dirname(__file__))
sys.path.append("../")
from cns_utils import BlastLine, parse_at_description
from pyfasta import Fasta
import operator
from collections import defaultdict
from BCBio.GFF import GFFParser

help = """\
All CNSs were blasted to Arabidopsis RNAs (TRNA, SNORNA, RIRNA, MIRNA). [Any arabidopsis accn from 
ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR9_genome_release/TAIR9_gff3/TAIR9_GFF3_genes.gff
that doesn't have a CDS.].
"""


def join_feat(key,seq_feat):
    feats = seq_feat[key]
    if len(feats) == 1: return feats[0]
    i,main_feat = [(i,feat) for i,feat in enumerate(feats) if feat.id == key][0]
    for fi,f in enumerate(feats):
        if i == fi: continue
        for subfeat in f.sub_features:
            main_feat.sub_features.append(subfeat)
    return main_feat


def conden_transcripts(seq_features):
    ids = set([])
    new_seq_feat = []
    seq_feat = defaultdict(list)
    for feat in seq_features:
        feat_id = feat.id.split('.')[0]
        seq_feat[feat_id].append(feat)
    for key in seq_feat.keys():
        new_feat = join_feat(key,seq_feat)
        new_seq_feat.append(new_feat)
    return new_seq_feat


def main(gff_file, outdir, th_fasta):
    """empty docstring"""
    parser = GFFParser()
    seqids = parser.parse(gff_file,None)

    fasta = Fasta(th_fasta, flatten_inplace=True)
    out_fasta = open(outdir + "/at_no_cds.fasta", "w")
    for seqid in seqids:
        seq_features = conden_transcripts(seqid.features)
        for feat in seq_features:
            has_cds = False
            ids = []
            ids.append(feat.id)
            for subf in feat.sub_features:
                if subf.type == 'CDS' or subf.type == 'chromosome':
                    has_cds = True
            if has_cds: continue
            #non_cds_feats.append(feat) 
            print >>out_fasta, ">%s" % ids[0]
            print >>out_fasta, fasta[seqid.id.lower()][int(feat.location.start) : int(feat.location.end)]


def read_descriptions(desc_file):
    desc = collections.defaultdict(str)
    for line in open(desc_file):
        line = line.split("\t")
        name = line[0][:line[0].find(".")]
        desc[name] += ";;" + line[-1].rstrip()
    return desc

def make_cns_to_at_map(des_file, blast_file, gff, query, subject, outdir):
    """
    take the cns vs at blast file, find the best at hit for the cns,
    and find teh tair desc for that hit, create a datastructure like:
        {'cns_id' : ('at_name', 'at_desc') ... } 
    e.g.
        {'q1|5766644|5766724|3|2158237|2158317' : ('AT2G38030', 'pre-tRNA')}

    this can then be used in   make_better_datasheet.
    """
    #os_desc = parse_os_description()
    at_desc = parse_at_description(des_file)

    blast_array = {}
    for line in open(blast_file):
        b = BlastLine(line)
        key = b.query, b.subject
        if not key in blast_array or b.score > blast_array[key].score:
            blast_array[key] = b
    blast_array = sorted(blast_array.values(), key=operator.attrgetter('score'), reverse=True)

    seen = {}
    updated = 0
    for row in blast_array:
        key = row.query[3:]
        #if key in seen: continue
        if key in seen:
            if len(seen[key][1]) < 10:
                seen[key][1] = at_desc[row.subject]
        else:
            seen[key] = [row.eval, at_desc[row.subject]]

        updated += 1
    of = outdir + "/cns_to_rna.csv"
    print >>sys.stderr, "writing to %s" % of
    fh = open(of, "w")
    for key, (eval, desc) in seen.iteritems():
        print >>fh, "\t".join(map(str, (key, eval, desc)))
    fh.close()

    print >>sys.stderr, "%i cns's tagged as hitting at_rna"  % updated

if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser(help)
    parser.add_option("-g", "--gff",     dest="gff", help="arabidopsis gff")
    parser.add_option("-b", "--blast",   dest="blast", help="the path of the file to send the blast output")
    parser.add_option("-d", "--desc",    dest="desc", help="TAIR descriptions file")
    parser.add_option("-o", "--outdir",  dest="outdir", help="output directory")
    parser.add_option("-q", "--query",   dest="query", help="name of query organism (rice)")
    parser.add_option("-s", "--subject", dest="subject", help="name of subject organism (sorghum)")
    parser.add_option("-f", "--fasta", dest="fasta", help="path to the thaliana fasta")
    parser.add_option("--blastpath", dest="blastpath", help="path to the blast dir")
    
    (options, _) = parser.parse_args()
    if not (options.gff and options.blast and options.desc and options.query and options.subject, options.outdir):
        sys.exit(parser.print_help())

    #print >>sys.stderr, "skipping main for testing"
    main(options.gff, options.outdir, options.fasta)

    import commands
    cmd = "bblast.py -b {0} -p blastn -e 0.001 -m 8 -W 7 -a 6 -i {1}/{2}_{3}.cns_test.fasta -d {1}/at_no_cds.fasta -o {4}".format(options.blastpath,options.outdir, options.query, options.subject, options.blast)
    print >>sys.stderr, "executing\n %s" % cmd
    print commands.getoutput(cmd)


    make_cns_to_at_map(options.desc, options.blast, options.gff, options.query,
                       options.subject, options.outdir)

