__author__ = 'james'

import sys


alleles = "ACGT"


def usage():
    print "USAGE: " + sys.argv[0] + " <ref-genome file> <reads file>"


def load_genome(filename):
    """
    Loads a genome into main memory from the specified file.
    Returns the genome as a list of alleles
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
            genome = list(genome_string)
    except IOError:
        print "Couldn't open file \"" + filename + "\". Exiting"
        sys.exit(1)

    return genome


def create_lookup_table(genome, seq_length):
    """
    Creates a hash table to look up positions of allele sub-sequences of length 10
    within the genome. Returns a dictionary as the lookup hash-table.
    """

    print "Creating sequence lookup hash table..."
    lookup_table = {}
    # store all positions of 10'mers occuring in the genome
    for i in range(0, len(genome) - seq_length):
        sequence = ''.join(genome[i: i + seq_length])
        lookup_table.setdefault(sequence, []).append(i)
    print "Created table of size: " + str(len(lookup_table))

    return lookup_table


################################# START OF SCRIPT ###################################
if len(sys.argv) < 3:
    usage()
    sys.exit(1)

ref_filename = sys.argv[1]
reads_file_name = sys.argv[2]

ref_genome = load_genome(ref_filename)
hash_table = create_lookup_table(ref_genome, 10)  # table with sequence length 10

reads = []

print "DONE"