__author__ = 'james'

import sys


alleles = "ACGT"


def usage():
    print "USAGE: " + sys.argv[0] + " <ref-genome file> <reads file>"


def error_die(msg):
    sys.stderr.write("ERROR: " + msg + "\n")
    sys.exit(1)


def load_genome(filename):
    """
    Loads a genome into main memory from the specified file.
    Returns genome as a list of alleles and the name of the genome
    """
    genome_string = ""
    try:
        with open(filename) as ref_file:
            line = ref_file.next()
            if not ">" in line:
                raise Exception("First line of genome does not begin with \'>\'")
            else:
                name = str(line.translate(None, ">\r\n"))
                print "Loading reference genome \'" + name + "\'..."

                for line in ref_file:
                    if ">" in line:
                        chr_name = line.translate(None, ">\r\n")
                        print "\tLoading chromosome \'" + chr_name + "\'..."
                        continue
                    line = line.rstrip("\r\n")
                    genome_string += line
                genome = list(genome_string)
    except IOError:
        sys.stderr.write("Couldn't open file \"" + filename + "\". Exiting\n")
        sys.exit(1)
    except Exception as e:
        error_die(e.message)

    return name, genome


def create_lookup_table(genome, seq_length):
    """
    Creates a hash table to look up positions of allele sub-sequences of specified
    length within the genome. Returns a dictionary as the lookup hash-table, with
    string sub-sequences as keys and a list of positions as the values
    """
    print "Creating sequence lookup hash table..."
    lookup_table = {}
    # store all positions of all 10'mers occuring in the genome
    for i in range(0, len(genome) - seq_length + 1):
        sequence = ''.join(genome[i: i + seq_length])
        lookup_table.setdefault(sequence, []).append(i)
        if len(sequence) != 10:
            print ("ERROR: sequence " + str(i) + "not of length 10")
            print "Exiting..."
            sys.exit(1)
    print "Created table of size: " + str(len(lookup_table))

    return lookup_table


def load_reads(filename):
    """
    Loads the reads file into main memory as a list. Returns the list of single reads
    and the name of the genome corresponding to the reads
    """
    reads = []
    try:
        with open(filename, "r") as reads_file:
            line = reads_file.next()
            if not ">" in line:
                raise Exception("First line of genome does not begin with \'>\'")

            name = str(line.translate(None, ">\r\n"))
            print "Loading reads from genome \'" + name + "\'..."

            for line in reads_file:
                if ">" in line:
                    chr_name = line.translate(None, ">\r\n")
                    print "\tChromosome \'" + chr_name + "\'..."
                    continue
                line = line.rstrip("\r\n")
                lines = line.split(",")
                reads.append(lines[0])
                reads.append(lines[1])
    except Exception as e:
        error_die(e.message)

    return name, reads


def main():
    if len(sys.argv) < 3:
        usage()
        sys.exit(1)

    ref_filename = sys.argv[1]
    reads_filename = sys.argv[2]

    ref_name, ref_genome = load_genome(ref_filename)
    reads_name, reads = load_reads(reads_filename)
    lookup_table = create_lookup_table(ref_genome, 10)  # table with sequence length 10

    print "Mapping reads to reference..."
    for read in reads:
        subsequences = []
        try:
            subsequences.append(read[0:10])
            subsequences.append(read[10:20])
            subsequences.append(read[20:30])
            subsequences.append(read[30:40])
            subsequences.append(read[40:50])
            for s in subsequences:
                try:
                    positions = lookup_table[s]
                    print s + ":" + str(positions)
                except KeyError:
                    print "KEY NOT FOUND"
                    continue
        except Exception as e:
            error_die(e.message)

    print "DONE"


# RUN MAIN
if __name__ == '__main__':
    main()