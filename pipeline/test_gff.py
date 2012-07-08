from flatfeature import Fasta
from BCBio.GFF import GFFParser,GFFExaminer
from collections import defaultdict

def join_feat(key,seq_feat):
    feats = seq_feat[key]
    if len(feats) == 1: return feats[0]
    i,main_feat = [i,feat for i,feat in enumerate(feats) if feat == key][0]
    for fi,f in enumerate(feats):
        if i == fi: continue
        for subfeat in f.sub_features:
            main_feat.sub_features.append(subfeat)
    return main_feat

def condens_transcript(seq_features):
    ids = set([])
    new_seq_feat = []
    seq_feat = defaultdict(list)
    for feat in seq_features:
        print type(feat)
        feat_id = feat.id.split('.')[0]
        seq_feat[feat_id].append(feat)i
    for key in seq_feat.keys():
        new_feat = join_feat(key,seq_feat)
        new_seq_feat.append(new_feat)
    return new_seq_feat
    
    

def main(gff_file,th_fasta):
    parser = GFFParser()
    #parser = GFFExaminer()
    seqids = parser.parse(gff_file, None)
    #seqids = parser.parent_child_map(gff_file)
    fasta = Fasta(th_fasta, flatten_inplace=True)
    out_fasta = open('this_is_a_test','w')
    for i,seqid in enumerate(seqids):
        ss= condens_transcript(seqid.features)
        for i,feat in enumerate(ss):
            #print feat
            ids = []
            has_cds = False
            ids.append(feat.id)
            for subf in feat.sub_features:
                if str(feat.type) == 'CDS' or feat.type == 'gene'  or feat.type == 'protein':
                    has_cds = True
            if has_cds: continue
            print >>out_fasta, '>%s' %ids[0]
            print >>out_fasta, fasta[seqid.id.lower()][int(feat.location.start):int(feat.location.end)]


#main('data/thaliana_v9_test2.gff','data/thaliana_v9.fasta')
main('data/thaliana_v9_test.gff','data/thaliana_v9.fasta')
#main('data/thaliana_v9.gff','data/thaliana_v9.fasta')


#{'a':[,1,2,3,4],'a.1':[5,6,7],'b':[7,8]}


