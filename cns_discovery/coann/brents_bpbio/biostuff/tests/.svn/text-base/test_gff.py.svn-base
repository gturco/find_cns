from biostuff.gff_reader import GFFLine, GFFNode
import sys

GFF = "tests/data/th.gff" 

def test_gff_line_to_string():
    l = [ll for ll in open(GFF) if ll[0] != '#'][0]
    g = GFFLine(l)

    ga = str(g.to_line())[:30]
    gb = l[:30]
    assert ga == gb

def test_gff_node():

    node_gen = GFFNode.yield_nodes(GFF)
    n = node_gen.next()
    assert len(n.nodes) == 0

    n = node_gen.next()
    assert n.parent.start == 3631, n.parent
    assert len(list(n.nodes)) == 12, list(n.nodes)

    n = node_gen.next()
    assert n.parent.attrs["ID"] == "thaliana_1_6160_6269", n.parent

    n = node_gen.next()
    assert n.parent.end == 8737, n.parent



def test_gff_line_all():
    for f in GFFLine.yield_lines(GFF): pass

def test_gff_all():
    for f in GFFNode.yield_nodes(GFF): pass

def test_gff_line():
    line_gen = GFFNode.yield_nodes(GFF)
