__author__ = 'James Gomez'

import sys
from math import *

from paired_end_read import PairedEndRead
import seq_align
import eval


allele_alphabet = "ACGT"
answer_file_name = "myanswers.txt"
answer_key_name = "ans_"


def usage():
    print "USAGE: " + sys.argv[0] + " <genome ID>"


def error_die(msg):
    sys.stderr.write("ERROR: " + msg + "\n")
    sys.exit(1)


def load_genome(filename):
    """
    Loads a genome into main memory from the specified file.
    Returns genome as a list of alleles and the name of the genome
    """
    genome = ""
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
                genome += line

    return name, genome


def load_reads(filename, paired_end):
    """
    Loads the reads file into main memory as a either list of single reads or a list
    of paired-end reads, depending on the 'paired-end' flag. Returns the list
    and the name of the genome corresponding to the reads.

    :param filename: A string. The name of the file to load the reads from
    :param paired_end: A boolean. Set to True to load reads as paired-end reads, or
    False to load the reads as single reads
    """
    reads = []
    with open(str(filename), "r") as reads_file:
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
            if paired_end:
                reads.append(PairedEndRead(lines[0], lines[1]))
            else:
                reads.append(lines[0])
                reads.append(lines[1])

    return name, reads


def create_lookup_table(genome, seq_length):
    """
    Creates a hash table to look up positions of allele sub-sequences of specified
    length within the genome. Returns a dictionary as the lookup hash-table, with
    string sub-sequences as keys and a list of positions as the values.
    Complexity: O(l/m) where l is the length of the reference genome, m is seq length

    :param genome: the genome for which the lookup table will be created
    :param seq_length: the length to use for the sub-sequence hash keys
    """
    print "Creating sequence lookup hash table..."
    lookup_table = {}
    # store all positions of all 10'mers occuring in the genome
    for i in range(0, len(genome) - seq_length + 1):
        sequence = ''.join(genome[i: i + seq_length])
        lookup_table.setdefault(sequence, []).append(int(i))
        # if len(sequence) != 10:
        #     raise Exception("sub-sequence " + str(i) + " is not of length 10")

    print "\tCreated table of size: " + str(len(lookup_table))
    return lookup_table


def map_reads(ref_genome, reads, lookup_table, subseq_length, min_score):
    print "Mapping reads..."
    read_length = 50
    read_map = {}
    for read in reads:
        if read_length % subseq_length != 0:
            raise Exception("Can't equally subdivide 50 by subseq_length")

        # subdivide the read into smaller sub-sequences for perfect-match hashing
        # then collect the candidate positions of the read in a list
        # positions = []
        positions = set()
        for start_pos in range(0, read_length, subseq_length):
            # don't use non-full-length sub seqs from end of read
            if start_pos + subseq_length > read_length:
                continue

            sub_sequence = read[start_pos: start_pos + subseq_length]

            try:
                temp_pos = lookup_table[sub_sequence]  # throws key error if no key
                for p in temp_pos:
                    # positions.append(int(p) - int(start_pos))
                    positions.add(int(p) - int(start_pos))
            except KeyError:
                continue

        # Find the best matching position for the read (max align score)
        REF = 0
        TEST = 1
        SCORE = 2
        best = (None, 0, 0)
        best_pos = None
        for p in positions:
            current = seq_align.align(
                ref_seq=ref_genome[p:p + len(read)],
                test_seq=read,
                local=True
            )
            if current[SCORE] > best[SCORE] and current[SCORE] > min_score:
                best = current
                best_pos = p

        if not best_pos is None:
            for i in range(0, len(best[TEST])):
                if best[TEST][i] == '-':  # skip gaps
                    continue
                try:
                    read_map[best_pos + i].append(str(best[TEST][i]))
                except KeyError:
                    read_map[best_pos + i] = [str(best[TEST][i])]

    return read_map


# def get_num_mismatches(sequence, ref_genome, position):
#     """
#     Returns the number of mismatches between the reference genome starting at the
#     specified position and the given sequence.
#     Complexity: O(len(sequence))
#
#     :param sequence: the sub-sequence to test against the reference
#     :param ref_genome: the reference genome to be tested against  the sequence
#     :param position: the position in the reference to compare with the sequence
#     """
#     num_mismatches = 0
#     for i in range(0, len(sequence)):
#         if position + i >= len(ref_genome):
#             break
#         if sequence[i] != ref_genome[position + i]:
#             num_mismatches += 1
#
#     return num_mismatches


# def get_best_read_position(ref_genome, read, positions, max_mismatches):
#     """
#     Maps the given read to the best position in the reference genome, if possible.
#     Returns the best matching position in the reference genome, or None if no match
#     is found.
#     Complexity: O(m*n) m is number of positions, n is the length of the read
#
#     :param ref_genome: the reference genome
#     :param read: the read to be mapped to the reference
#     :param positions: the positions in the reference genome to test the read against
#     :param max_mismatches: the max allowed mismatches
#     """
#     least = 100
#     best_pos = None
#     for p in positions:
#         num_mismatches = get_num_mismatches(read, ref_genome, p)
#         if num_mismatches < max_mismatches and num_mismatches < least:
#             least = num_mismatches
#             best_pos = p
#
#     return best_pos


def get_consensus_allele(alleles):
    """
    Returns the allele with the highest count.
    Complexity: O(len(alleles))
    :rtype : str
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


def find_snps(answer_file, ref_genome, ref_name, read_map, thresh=0.6):
    """
    Finds SNPs with respect to the reference genome and writes the SNPs to the answer
    file.
    Complexity: O(n) where n is the length of the reference genome

    :param answer_file: The file to write the SNPs to
    :param ref_genome: the reference genome to compare against
    :param ref_name: the name of the reference genome
    :param read_map: the read_map to use hashing sub-sequences
    """
    # Use the consensus algorithm to determine SNPs relative to the reference genome
    # Write the SNPs to the answer file
    # TODO use the thresh param
    print "Finding SNPs..."
    answer_file.write(">" + ref_name + "\n")
    answer_file.write(">SNP" + "\n")
    for i in range(0, len(ref_genome)):
        ref_allele = ref_genome[i]
        try:
            # count how many of each allele appears in the current position of
            # the read map
            read_alleles = read_map[i]  # throws key error if bad key
            winner = get_consensus_allele(read_alleles)
            if winner != ref_allele:  # if not the same, it's a SNP
                answer_file.write(
                    "1," + str(ref_allele) + "," + str(winner) + "," + str(
                        i) + "\n")
        except KeyError:
            continue


def find_indels(seq1, seq2):
    inserts = []
    deletes = []

    return inserts, deletes


def main():
    # ALGORITHM VARIABLES
    global answer_key_name
    seq_length = 10
    min_align_score = 0.95
    snp_thresh = 0.6
    paired_end = False

    # ensure correct number of arguments
    if len(sys.argv) < 2:
        usage()
        sys.exit(1)

    print
    ref_filename = "ref_" + str(sys.argv[1]) + ".txt"
    reads_filename = "reads_" + str(sys.argv[1]) + ".txt"

    ref_name, ref_genome = load_genome(ref_filename)
    reads_name, reads = load_reads(reads_filename, paired_end=paired_end)
    if ref_name != reads_name:
        error_die("Reference genome id and reads genome id do not match")
    answer_key_name += (str(ref_name) + ".txt")

    lookup_table = create_lookup_table(ref_genome, seq_length)
    read_map = map_reads(
        ref_genome,
        reads,
        lookup_table,
        subseq_length=seq_length,
        min_score=min_align_score
    )

    with open(answer_file_name, "w") as answer_file:
        find_snps(
            answer_file,
            ref_genome,
            ref_name,
            read_map,
            thresh=snp_thresh
        )
        # TODO find indels

    print "DONE\n"


def run_eval():
    try:
        with open(answer_file_name, "r") as student_ans:
            with open(answer_key_name, "r") as answer_key:
                eval.Eval(answer_key, student_ans)
    except IOError as e:
        sys.stderr.write("Couldn't open answer key \'" + answer_key_name + "\'\n")
        sys.stderr.write(e.message + "\n")


# RUN MAIN
if __name__ == '__main__':
    # sequence_aligner.align(ref_seq="TACTGGATGA", test_seq="TTGGATGCTA")
    main()
    run_eval()