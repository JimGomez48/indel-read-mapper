__author__ = 'James Gomez'

import sys


alleles = "ACGT"


def usage():
    print "USAGE: " + sys.argv[0] + " <ref-genome file> <reads file> <match-thresh>"


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


def get_diff_list(sequence, ref_genome, position):
    """
    Returns a list containing the mismatches between the reads and ref genome.
    Specifically, original allele, snp, position in reference
    """
    read_list = list(sequence)
    ref_list = list(ref_genome[position:position+len(sequence)])
    mismatches = list()
    try:
        if len(read_list) != len(ref_list):
            return mismatches
        for i in range(0, len(sequence)):
            if read_list[i] != ref_list[i]:
                mismatches.append([ref_list[i], read_list[i], position + i])
    except Exception as e:
        error_die("num_mismatches(): " + e.message)

    return mismatches


def get_num_mismatches(sequence, ref_genome, position):
    characters = list(sequence)
    num_mismatches = 0
    for i in range(0, len(characters)):
        if characters[i] != ref_genome[position + i]:
            num_mismatches += 1

    return num_mismatches

def get_best_read_position(ref_genome, read, positions, thresh):
    least = 100
    best_pos = None
    for p in positions:
        num_mismatches = get_num_mismatches(read, ref_genome, p)
        if num_mismatches < thresh:
            if num_mismatches < least:
                least = num_mismatches
                best_pos = p

    return best_pos


def main():
    if len(sys.argv) < 4:
        usage()
        sys.exit(1)

    ref_filename = sys.argv[1]
    reads_filename = sys.argv[2]
    threshold = sys.argv[3]

    ref_name, ref_genome = load_genome(ref_filename)
    reads_name, reads = load_reads(reads_filename)
    if ref_name != reads_name:
        error_die("Reference genome id and reads genome id do not match")
    lookup_table = create_lookup_table(ref_genome, 10)  # table with seq length 10

    read_map = []
    snps = {}

    print "Mapping reads..."
    for read in reads:
        sub_sequences = [
            read[0:10],
            read[10:20],
            read[20:30],
            read[30:40],
            read[40:50]
        ]

        positions = []
        for s in sub_sequences:
            try:
                positions = lookup_table[s]  # throws key error if not in table
            except KeyError:
                continue
            if len(positions) > 0:
                best_pos = get_best_read_position(ref_genome, read, positions, threshold)
                if best_pos is None:
                    continue
                for i in range(best_pos, best_pos + len(read)):
                    read_map[best_pos] = read[i]

    ans_key_name = "myanswers.txt"
    print "writing answer key to \'" + ans_key_name + "\'..."
    with open(ans_key_name, "w") as my_answer_file:
        my_answer_file.write(">" + ref_name + "\n")
        my_answer_file.write(">SNP\n")
        keys = snps.keys()
        for k in keys:
            m = snps[k]
            my_answer_file.write("1," + str(m[0]) + "," + str(m[1]) + "," + str(k) + "\n")

    print "DONE"


# RUN MAIN
if __name__ == '__main__':
    main()