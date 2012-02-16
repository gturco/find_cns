from cns_utils import CNS

x = "../data/rice_j_sorghum_n/rice_j_sorghum_n.cns.txt"
for line in open(x):
    if line[0] == "#": continue
    x = CNS(line)
    y= x.get_cns()
    print y.sstart, x.qaccn
