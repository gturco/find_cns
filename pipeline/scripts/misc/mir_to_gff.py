
import sys
import collections


seen = {}
chrs = collections.defaultdict(list)
for line in open(sys.argv[1]):
    if line[0] == '#': continue
    name, rfam, chr, start, stop, strand, notes = [x.strip() for x in line.split("\t")]
    chr, start, stop = map(int, (chr, start, stop))
    if start > stop:
        start, stop = stop, start
    chr, start, stop = map(str, (chr, start, stop))
    notes = notes.strip()

    strand = strand == '-1' and '-' or '+'
    if name in seen:
        1/0

    attrs = "ID=%s;rfam=%s;rname=%s" % (name, rfam, name)
    #if notes:
    #    attrs += ";notes=\"" + notes + '"'

    line = "\t".join((chr, "ucb", "MIR", start, stop, ".", strand, ".", attrs))
    chrs[chr].append(line)


print "##gff-version 3"
for chr, lines in sorted(chrs.items()):
    for line in lines:
        print line

