from flatfeature import Fasta
from BCBio.GFF import GFFParser

def main(gff_file,th_fasta):
    parser = GFFParser()
    seqids = parser.parse(gff_file, None)
    fasta = Fasta(th_fasta, flatten_inplace=True)
    out_fasta = open('this_is_a_test','w')
    for seqid in seqids:
        has_cds = False
        for feat in seqid.features:
            #print feat
            subf = feat.sub_features
            ids = []
            for f in subf:
                if 'ID' in f.qualifiers: ids.append(f.qualifiers['ID'])
                #print subf.attribs
                #print feat.attribs
                if f.type == 'CDS':
                    has_cds = True
                    break
            if has_cds: continue
            print >>out_fasta, '>%s' %ids[0]
            #print 'subF', subf
            print 'FEAT' ,feat.location.start
            print 'f',type(f.location.start)
            print >>out_fasta, fasta[seqid.id.lower()][int(feat.location.start)-1:int(feat.location.end)]


main('data/thaliana_v9_test.gff','data/thaliana_v9.fasta')
