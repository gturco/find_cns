import csv
def parse_cns_pair(cns_pairs):
    handle = open("/Users/gturco/rice_v6_rice_v6.txt", "w")
    for line in cns_pairs:
        cnss = line.split(',')
        if len(cnss) > 11:
            info = cnss[:6]
            res = cnss[6:]
            urls = res[(len(res)/5)  * 4:]
            for number in range(len(res)/5):
                cns = res[(number * 4):((number+1) * 4)]
                url = urls[(number)]
                joined = ','.join(info) , ','.join(cns) , url
                row = ','.join(joined) + '\n'
                handle.write(row)
        else:
            handle.write(line)
    handle.close()
             
handle = open('/Users/gturco/find_regions_data/output/3_23_11/rice_v6_rice_v6_cns/rice_v6_rice_v6.cns.txt','r')
l = handle.readlines()        
x = l[:10]
m= parse_cns_pair(l)
