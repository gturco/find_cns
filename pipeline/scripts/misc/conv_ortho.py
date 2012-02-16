from genedex.misc.gff import gff_line
genes = {}
for line in open("data/rice_rice/rice.gff"):
    if not 'gene' in line: continue
    feat = gff_line(line)
    if feat['type'] != 'gene': continue
    genes[feat['name'].upper()] = feat


for sline in open('data/rice_rice.orthology'):
    line = sline.rstrip().split(",")
    try:
        print ",".join(map(str, [line[0], line[1],
                   genes.get(line[2])['start'],
                   genes.get(line[3])['stop'],
                   genes.get(line[4])['start'],
                   genes.get(line[5])['stop']]))

    except:
        print sline.rstrip()
