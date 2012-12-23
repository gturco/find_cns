import sys

import logging
logging.basicConfig(level=logging.DEBUG)

def _get_non_comment_line(fh):
    while True:
        line = fh.readline()
        if line and line[0] == '#': continue
        return line

def _get_lines_until_next(fh, parent_ids):
    """
    so this just continues to append lines to lines[]
    until Parent attribute of a line does not match
    any of the values in parent_ids
    any time a line is found with a new ID attribute
    (whose Parent attr matches the current parent_ids list),
    that line's own ID is added to the parent_ids list"""

    lines = [_get_non_comment_line(fh)]
    while True:
        new_parent = False
        for parent_id in parent_ids:
            if 'Parent=' + parent_id in lines[-1]:
                if not lines[-1]: return lines # end of file
                lines.append(_get_non_comment_line(fh))
                if 'ID=' in lines[-1]:
                    new_parent=True
                break
        else:
            break
        if new_parent:
            parent_ids.append(GFFLine(lines[-1]).attribs["ID"])
            parent_ids = list(set(parent_ids))

    return lines

class GFFNode(object):
    __slots__ = ('start', 'stop', 'end', 'parent', 'nodes', 'seqid')
    def __init__(self, node_list):
        self.start = min(n.start for n in node_list)
        self.stop = max(n.stop for n in node_list)
        self.end = self.stop
        assert "ID" in node_list[0].attribs, (node_list[0], node_list[0].attribs)
        self.parent = node_list[0]
        self.nodes = node_list[1:]
        self.seqid = self.parent.seqid

        if self.parent.start != self.start:
            logging.debug(("the start of the parent != the start of the item:" + \
                        ", ".join(map(str, (self, self.start, self.parent.start)))))
        if self.parent.end != self.end:
            logging.debug("the end of the parent != the end of the item:" + \
                    ", ".join(map(str, (self, self.end, self.parent.end))))
    
    def __repr__(self):
        return "GFFNode(%s: %i .. %i, %i sub-nodes)" % \
                   (self.parent.attrs.get('ID', ''), self.start, self.end,
                    len(self.nodes))


    @classmethod
    def yield_nodes(cls, fh):
        close = False
        if not hasattr(fh, 'read'):
            close = True
            fh = open(fh, 'r')

        next_gene_line = _get_non_comment_line(fh)
        while next_gene_line:
            if "rname=" in next_gene_line and not "ID=" in next_gene_line:
                next_gene_line = next_gene_line.replace("rname=", "ID=")
            assert "ID=" in next_gene_line,\
                    ("should have an id to be a parent feature", next_gene_line)

            parent = GFFLine(next_gene_line)
            block = _get_lines_until_next(fh, [parent.attribs["ID"]])
            block.insert(0, next_gene_line)
            next_gene_line = block.pop()

            yield GFFNode([GFFLine(l) for l in block])


        if close: fh.close()


class GFFLine(object):
    __slots__ = ('seqid', 'com', 'type', 'start', 'stop', 'end', 'strand', 'other',
                    'attrs', 'attribs', 'sattrs', 'orig')
    def __init__(self, sline):
        line = sline.rstrip().split("\t")
        self.seqid = line[0]
        self.com  = line[1]
        self.type = line[2]
        self.start = int(line[3])
        self.stop = self.end = int(line[4])
        self.orig = line[5]
        self.strand = line[6] in ('-', '-1') and -1 or 1
        self.other = line[7]
        self.sattrs = line[8]
        self.attrs = self.attribs = self._parse_attrs(line[8])
    
    def _parse_attrs(self, sattrs):
        attrs = {}
        if "=" in sattrs:
            for pair in sattrs.split(';'):
                if not pair: continue
                pair = pair.split('=')
                attrs[pair[0]] = pair[1]
        if attrs == {}: attrs["ID"] = sattrs.rstrip("\r\n; ")
        return attrs

    @classmethod
    def yield_lines(cls, fh):
        close = False
        if not hasattr(fh, 'read'):
            fh = open(fh, 'r')
            close = True
        line = _get_non_comment_line(fh)
        while line:
            yield GFFLine(line)
            line = _get_non_comment_line(fh)
    
        if close: fh.close()

    def __repr__(self):
        return "GFFLine(%s %s:%i .. %i)" % (self.seqid,
                          self.attrs.get("ID", self.attrs.get("Parent", "")), self.start, self.end)

    def to_line(self):
        s =  [getattr(self, s) for s in ('seqid', 'com', 'type', 'start', 'end', 'orig', 'strand',
                    'other', 'sattrs')]
        s[6] = s[6] in ('+', 1, '1') and '+' or '-'
        return "\t".join(map(str, s))

class GFFDict(dict):
    def __init__(self, gff_path):
        self.gff_path = gff_path
        self.load_gff(gff_path)
        self.chrs = {}


    def load_gff(self, gff_path):

        for gffnode in GFFNode.yield_nodes(gff_path):
            s = gffnode.seqid
            if not s in self.chrs: self.chrs[s] = {}
            self.chrs[s][gffnode.parent.attribs['ID']] = gffnode

    def __getitem__(self, k):
        for achr in self.chrs:
            try: return self.chrs[achr][k]
            except KeyError: continue
        raise KeyError
