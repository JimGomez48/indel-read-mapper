__author__ = 'james'

import sys


alleles = "ACGT"


def usage():
    print "USAGE: " + sys.argv[0] + " <ref-genome file> <reads file>"


def load_genome(filename):
    """
    Loads the reference genome into main memory from the specified file
    """

    genome_string = ""
    print "Loading reference genome..."
    try:
        with open(filename) as ref_file:
            for line in ref_file:
                if ">" in line:
                    continue
                line = line.rstrip('\n')
                genome_string += line
            genome = genome_string.split()
    except IOError as e:
        print "Couldn't open file \"" + filename + "\". Exiting"
        sys.exit(1)

    return genome


def create_lookup_table():
    """
    Creates a hash table to look up positions of allele sub-sequences of length 10
    within the genome
    """

    lookup_table = {}
    key = list("AAAAAAAAAA")
    for c1 in alleles:
        key[0] = c1
        for c2 in alleles:
            key[1] = c2
            for c3 in alleles:
                key[2] = c3
                for c4 in alleles:
                    key[3] = c4
                    for c5 in alleles:
                        key[4] = c5
                        for c6 in alleles:
                            key[5] = c6
                            for c7 in alleles:
                                key[6] = c7
                                for c8 in alleles:
                                    key[7] = c8
                                    for c9 in alleles:
                                        key[8] = c9
                                        for c10 in alleles:
                                            key[9] = c10
                                            lookup_table[''.join(key)] = ""

    print "Size of table: " + str(len(lookup_table))
    print "Creating sequence lookup hash table..."

    return lookup_table


################################# START OF SCRIPT ###################################
if len(sys.argv) < 3:
    usage()
    sys.exit(1)

ref_filename = sys.argv[1]
reads_file_name = sys.argv[2]

ref_genome = load_genome(ref_filename)
hash_table = create_lookup_table()

reads = []

print "DONE"