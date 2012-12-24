from cblastline import BlastLine, BlastFile
from gff_reader import GFFNode, GFFLine
#from blasttree import blast_to_tree, BlastTree
import cblastline


def main():
    import sys
    progs = ['cblastline']
    if len(sys.argv) == 1 or not sys.argv[1] in progs:
        print """Usage:
    $ biostuff prog [opts]
or 
    $ biostuff prog -h

where current prog's are: %s
        """ % "\n".join(progs)
    else:
        prog = sys.argv.pop(1)
        print prog
        mod = __import__('biostuff.' + prog, fromlist=['main'])
        print mod
        mod.main()

if __name__ == "__main__":
    main()
