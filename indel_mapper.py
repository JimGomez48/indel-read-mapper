__author__ = 'James Gomez'

import sys
from math import *


allele_alphabet = "ACGT"


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
        lookup_table.setdefault(sequence, []).append(int(i))
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


def find_snps(sequence, ref_genome, position):
    """
    Returns a list containing the mismatches between the reads and ref genome, or
    None if no SNPs are found
    Specifically: [original allele, snp, position in reference]
    """
    seq_list = list(sequence)
    ref_list = list(ref_genome[position:position+len(sequence)])
    mismatches = list()
    try:
        if len(seq_list) != len(ref_list):
            return None
        for i in range(0, len(seq_list)):
            if seq_list[i] != ref_list[i]:
                mismatches.append([ref_list[i], seq_list[i], position + i])
    except Exception as e:
        error_die("num_mismatches(): " + e.message)

    return mismatches


def get_num_mismatches(sequence, ref_genome, position):
    """
    Returns the number of mismatches between the reference genome starting at the
    specified position and the given sequence
    """
    characters = list(sequence)
    num_mismatches = 0
    for i in range(0, len(characters)):
        if position + i >= len(ref_genome):
            break
        if characters[i] != ref_genome[position + i]:
            num_mismatches += 1

    return num_mismatches


def get_best_read_position(ref_genome, read, positions, thresh):
    """
    Maps the given read to the best position in the reference genome, if possible.
    Returns the best matching position in the reference genome, or None if no match
    is found
    """
    least = 100
    best_pos = None
    for p in positions:
        num_mismatches = get_num_mismatches(read, ref_genome, p)
        if num_mismatches < thresh and num_mismatches < least:
                least = num_mismatches
                best_pos = p

    return best_pos


def get_consensus_allele(alleles):
    """
    Returns the allele with the highest count
    """
    allele_counts = {}
    for a in alleles:
        assert a in allele_alphabet
        try:
            allele_counts[a] += 1
        except KeyError:
            allele_counts[a] = 1

    keys = allele_counts.keys()
    win_count = 0
    winner = None
    for k in keys:
        if allele_counts[k] > win_count:
            winner = k
            win_count = allele_counts[k]

    return winner


def main():
    # ALGORITHM VARIABLES
    hash_key_length = 10


    # ensure correct number of arguments
    if len(sys.argv) < 4:
        usage()
        sys.exit(1)

    ref_filename = sys.argv[1]
    reads_filename = sys.argv[2]
    thresh = int(sys.argv[3])

    ref_name, ref_genome = load_genome(ref_filename)
    reads_name, reads = load_reads(reads_filename)
    if ref_name != reads_name:
        error_die("Reference genome id and reads genome id do not match")
    lookup_table = create_lookup_table(ref_genome, hash_key_length)  # table with seq length 10

    print "Mapping reads..."
    read_map = {}
    for read in reads:
        # subdivide the read into smaller sub-sequences for perfect-match hashing
        # then collect the candidate positions of the read in a list
        positions = []
        for start_pos in range(0, len(read), 10):
            # don't use non-full-length sub seqs from end of read
            if start_pos + hash_key_length > len(read):
                continue

            sub_sequence = read[start_pos: start_pos + hash_key_length]

            try:
                temp_pos = lookup_table[sub_sequence]  # throws key error if no key
                for p in temp_pos:
                    positions.append(int(p) - int(start_pos))
            except KeyError:
                continue

        # find the best matching position for the read from the list of candidate
        # positions, i.e. find the position that minimizes the error between the read
        # and the corresponding sub-sequence in the reference genome. Then store it
        # in the read map, or continue if no position satisfies the threshold
        read_position = get_best_read_position(ref_genome, read, positions, thresh)
        if not read_position is None:
            for i in range(0, len(read)):
                # if read_map[i] is None:
                #     read_map[read_position + i] = []
                try:
                    read_map[read_position + i].append(str(read[i]))
                except KeyError:
                    read_map[read_position + i] = [str(read[i])]

    # Use the consensus algorithm to determine SNPs relative to the reference genome
    # Write the SNPs to the answer file
    print "Finding SNPs..."
    answer_file_name = "myanswers.txt"
    with open(answer_file_name, "w") as answer_file:
        answer_file.write(">" + ref_name + "\n")
        answer_file.write(">SNP" + "\n")
        for i in range(0, len(ref_genome)):
            ref_allele = ref_genome[i]
            try:
                # count how many of each allele appears in the current position of the
                # read map
                read_alleles = read_map[i]  # throws key error if bad key
                winner = get_consensus_allele(read_alleles)
                print winner
                print read_alleles
                if winner != ref_allele:  # if not the same, it's a SNP
                    answer_file.write("1," + str(ref_allele) + "," + str(winner) + "," + str(i) + "\n")
            except KeyError:
                continue

    print "DONE"


# RUN MAIN
if __name__ == '__main__':
    main()