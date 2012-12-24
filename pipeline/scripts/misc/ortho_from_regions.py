import re
import collections

regions = 'data/thaliana.alpha.regions.csv'

region_re = re.compile("([SA]\d{2,3})\w")

qs = collections.defaultdict(list)
ss = collections.defaultdict(list)
for line in open(regions):
    q, s, alpha = [x.strip() for x in line.upper().split("\t")]
    if alpha[0] == 'B': continue
    q, s = sorted([q, s])
    
    region = region_re.search(alpha).groups(0)[0]

    ss[region].append(s)
    qs[region].append(q)

for region in ss:
    slist = sorted(ss[region])
    qlist = sorted(qs[region])
    print ",".join((region, qlist[0], qlist[-1], slist[0], slist[-1]))
