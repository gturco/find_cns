import sys
import os.path as op
sys.path.insert(0, op.dirname(__file__))
sys.path.append("../")
from common import BlastLine, parse_at_description
import collections
from pyfasta import Fasta
import operator
import gt
gt.warning_disable()


help = """\
All CNSs were blasted to Arabidopsis RNAs (TRNA, SNORNA, RIRNA, MIRNA). [Any arabidopsis accn from 
ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR9_genome_release/TAIR9_gff3/TAIR9_GFF3_genes.gff
that doesn't have a CDS.].
"""

def main(gff_file, outdir, th_fasta):
    """empty docstring"""
    fi = gt.FeatureIndexMemory()
    fi.add_gff3file(gff_file)

    #non_cds_feats = {}
    fasta = Fasta(th_fasta, flatten_inplace=True)
    out_fasta = open(outdir + "/at_no_cds.fasta", "w")
    for seqid in fi.get_seqids():
        for feat in fi.get_features_for_seqid(seqid):
            has_cds = False
            ids = []
            for subf in feat:
                if "ID" in feat.attribs: ids.append(feat.attribs["ID"])
                if subf.type == 'CDS': 
                    has_cds = True
                    break
            if has_cds: continue
            #non_cds_feats.append(feat) 
            print >>out_fasta, ">%s" % ids[0]
            print >>out_fasta, fasta[feat.seqid.lower()][feat.start - 1: feat.end]


def read_descriptions(desc_file):
    desc = collections.defaultdict(str)
    for line in open(desc_file):
        line = line.split("\t")
        name = line[0][:line[0].find(".")]
        desc[name] += ";;" + line[-1].rstrip()
    return desc

def make_cns_to_at_map(blast_file, gff, query, subject, outdir):
    """
    take the cns vs at blast file, find the best at hit for the cns,
    and find teh tair desc for that hit, create a datastructure like:
        {'cns_id' : ('at_name', 'at_desc') ... } 
    e.g.
        {'q1|5766644|5766724|3|2158237|2158317' : ('AT2G38030', 'pre-tRNA')}

    this can then be used in   make_better_datasheet.
    """
    #os_desc = parse_os_description()
    at_desc = parse_at_description()

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
    (options, _) = parser.parse_args()
    if not (options.gff and options.blast and options.desc and options.query and options.subject, options.outdir):
        sys.exit(parser.print_help())

    #print >>sys.stderr, "skipping main for testing"
    main(options.gff, options.outdir, options.fasta)

    import commands
    cmd = "bblast.py -p blastn -e 0.001 -m 8 -W 7 -a 6 -i %s/%s_%s.cns.fasta -d %s/at_no_cds.fasta -o %s" \
                          % (options.outdir, options.query, options.subject, options.outdir, options.blast)
    print >>sys.stderr, "executing\n %s" % cmd
    print commands.getoutput(cmd)


    make_cns_to_at_map(options.blast, options.gff, options.query,
                       options.subject, options.outdir)

