import sys
sys.path.append("scripts")
from find_exons import get_shelve
import os
import re
import collections
from genedex.misc.gff import read_gff
from pyfasta import Fasta
import operator
import cPickle

"""
Chr1    TAIR7   chromosome  1   30432563    .   .   .   ID=Chr1;Name=Chr1
Chr1    TAIR7   gene    3631    5899    .   +   .   ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
Chr1    TAIR7   five_prime_UTR  3631    3759    .   +   .   Parent=AT1G01010.1
Chr1    TAIR7   mRNA    3631    5899    .   +   .   ID=AT1G01010.1;Parent=AT1G01010;Name=AT1G01010.1
Chr1    TAIR7   exon    3631    3913    .   +   .   Parent=AT1G01010.1
Chr1    TAIR7   CDS 3760    3913    .   +   0   Parent=AT1G01010.1;Name=AT1G01010.1-CDS
Chr1    TAIR7   CDS 3996    4276    .   +   2   Parent=AT1G01010.1;Name=AT1G01010.1-CDS
Chr1    TAIR7   exon    3996    4276    .   +   .   Parent=AT1G01010.1
Chr1    TAIR7   CDS 4486    4605    .   +   0   Parent=AT1G01010.1;Name=AT1G01010.1-CDS
Chr1    TAIR7   exon    4486    4605    .   +   .   Parent=AT1G01010.1
Chr1    TAIR7   CDS 4706    5095    .   +   0   Parent=AT1G01010.1;Name=AT1G01010.1-CDS
"""

help = """\
All CNSs were blasted to Arabidopsis RNAs (TRNA, SNORNA, RIRNA, MIRNA). [Any arabidopsis accn from 
ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR7_genome_release/TAIR7_gff3/TAIR7_GFF3_genes.gff
that doesn't have a CDS.].
"""

def main(gff_file, outdir):
    """empty docstring"""
    name = re.compile("parent=([^.;]+)", re.I)

    feats = {}
    non_cds_feats = collections.defaultdict(list)
    for line in open(gff_file):
        line = line.split("\t")
        match = re.search(name, line[-1])
        if not match: continue
        fname = match.groups(0)[0]
        non_cds_feats[fname].append(line)
        if line[2].upper() == 'CDS':
            feats[fname] = True
            continue
        if fname in feats: continue
        feats[fname] = None
    i = 0 
    for k, v in sorted(feats.items()):
        if not v is None: del non_cds_feats[k]

    seen = {}
    RNA = open(outdir + '/at_non_cds.gff', 'w')
    for k, feat_list in sorted(non_cds_feats.items()):
        for feat in feat_list:
            if feat[0] in ('ChrC', 'ChrM'): continue
            if feat[2] == 'exon': continue
            key = (feat[0], feat[3], feat[4])
            if key in seen: continue
            feat[0] = feat[0].upper().replace('CHR', '')
            seen[key] = True
            feat[-1] = k
            print >> RNA, "\t".join(feat)
    RNA.close()

    
    gff = read_gff(outdir + '/at_non_cds.gff')
    fasta = Fasta('/home/gturco/src/find_cns_gturco/pipeline/data/arabidopsis.fasta')
    ftypes = {}
    FA = open(outdir + '/at_rnas.fasta','w')
    for chr, feature_list in gff.iteritems():
        for fname, feature in feature_list.iteritems():
            seq = fasta.sequence(feature)
            print >>FA, ">", feature['name']
            print >>FA, seq
    FA.close()

def read_descriptions(desc_file):
    desc = collections.defaultdict(str)
    for line in open(desc_file):
        line = line.split("\t")
        name = line[0][:line[0].find(".")]
        desc[name] += ";;" + line[-1].rstrip()
    return desc

def make_cns_to_at_map(blast_file, tair_desc, gff, query, subject, outdir):
    """
    take the cns vs at blast file, find the best at hit for the cns,
    and find teh tair desc for that hit, create a datastructure like:
        {'cns_id' : ('at_name', 'at_desc') ... } 
    e.g.
        {'q1|5766644|5766724|3|2158237|2158317' : ('AT2G38030', 'pre-tRNA')}

    this can then be used in   make_better_datasheet.
    """
    shelve = get_shelve(outdir, query, subject)
    from blast_misc import blast_array
    b = blast_array(blast_file, best_hit=True)
    seen = {}
    updated = 0
    for row in b:
        if str(row['query']) in seen: continue
        seen[str(row['query'])] = True
        besthit = sorted(b[b['query'] == row['query']], key=operator.itemgetter('eval'), reverse=True)[0]
        # 1: get rid of starting SB or OS
        key = str(row['query'])[3:]

        current = shelve[key]
        current["at_rna"] = str(row['subject']) + ';;eval(' + str(row['eval']) + ')' + tair_desc[str(row['subject'])]
        updated += 1
        shelve[key] = current
    print >>sys.stderr, "%i cns's tagged as hitting at_rna"  % updated
    shelve.close()

if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser(help)
    parser.add_option("-g", "--gff",     dest="gff", help="arabidopsis gff")
    parser.add_option("-b", "--blast",   dest="blast", help="the path of the file to send the blast output")
    parser.add_option("-d", "--desc",    dest="desc", help="TAIR descriptions file")
    parser.add_option("-o", "--outdir",  dest="outdir", help="output directory")
    parser.add_option("-q", "--query",   dest="query", help="name of query organism (rice)")
    parser.add_option("-s", "--subject", dest="subject", help="name of subject organism (sorghum)")
    (options, _) = parser.parse_args()
    if not (options.gff and options.blast and options.desc and options.query and options.subject, options.outdir):
        sys.exit(parser.print_help())


    main(options.gff, options.outdir)

    import commands
    cmd = "bblast.py -p blastn -e 0.001 -m 8 -W 7 -a 6 -i %s/%s_%s_cns.fasta -d %s/at_rnas.fasta -o %s" \
                          % (options.outdir, options.query, options.subject, options.outdir, options.blast)

    print >>sys.stderr, "executing\n %s" % cmd
    print commands.getoutput("formatdb -p F -i %s/at_rnas.fasta" % options.outdir)
    print commands.getoutput(cmd)


    descriptions = read_descriptions(options.desc)

    make_cns_to_at_map(options.blast, descriptions, options.gff, options.query,
                       options.subject, options.outdir)

