__author__ = 'james'

import sys


def usage():
    print "USAGE: " + sys.argv[0] + "<ref-genome file> <reads file>"


################################# START OF SCRIPT ###################################
if len(sys.argv) < 3:
    usage()
    sys.exit(1)

ref_genome = []
reads = []