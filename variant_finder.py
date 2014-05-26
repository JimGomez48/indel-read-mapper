__author__ = 'James Gomez'

import sys
from math import *

from paired_end_read import PairedEndRead
from sequence_aligner import SequenceAligner
import Eval



allele_alphabet = "ACGT"
answer_file_name = "myanswers.txt"
seq_aligner = SequenceAligner()


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

    print "Created table of size: " + str(len(lookup_table))
    return lookup_table


def create_read_map(ref_genome, reads, lookup_table, subseq_length, thresh):
    print "Mapping reads..."
    read_length = 50
    read_map = {}
    for read in reads:
        if read_length % subseq_length != 0:
            raise Exception("Can't equally subdivide read by subseq_length")

        # subdivide the read into smaller sub-sequences for perfect-match hashing
        # then collect the candidate positions of the read in a list
        positions = []
        for start_pos in range(0, read_length, subseq_length):
            # don't use non-full-length sub seqs from end of read
            if start_pos + subseq_length > read_length:
                continue

            sub_sequence = read[start_pos: start_pos + subseq_length]

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
            for i in range(0, read_length):
                try:
                    read_map[read_position + i].append(str(read[i]))
                except KeyError:
                    read_map[read_position + i] = [str(read[i])]

    return read_map


def get_num_mismatches(sequence, ref_genome, position):
    """
    Returns the number of mismatches between the reference genome starting at the
    specified position and the given sequence.
    Complexity: O(len(sequence))

    :param sequence: the sub-sequence to test against the reference
    :param ref_genome: the reference genome to be tested against  the sequence
    :param position: the position in the reference to compare with the sequence
    """
    num_mismatches = 0
    for i in range(0, len(sequence)):
        if position + i >= len(ref_genome):
            break
        if sequence[i] != ref_genome[position + i]:
            num_mismatches += 1

    return num_mismatches


def get_best_read_position(ref_genome, read, positions, max_mismatches):
    """
    Maps the given read to the best position in the reference genome, if possible.
    Returns the best matching position in the reference genome, or None if no match
    is found.
    Complexity: O(m*n) m is number of positions, n is the length of the read

    :param ref_genome: the reference genome
    :param read: the read to be mapped to the reference
    :param positions: the positions in the reference genome to test the read against
    :param max_mismatches: the max allowed mismatches
    """
    least = 100
    best_pos = None
    for p in positions:
        num_mismatches = get_num_mismatches(read, ref_genome, p)
        if num_mismatches < max_mismatches and num_mismatches < least:
            least = num_mismatches
            best_pos = p

    return best_pos


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


def find_snps(answer_file, ref_genome, ref_name, read_map):
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


def find_insertions():
    # TODO
    pass


def find_deletions():
    # TODO
    pass


def needleman_wunsch_align(seq1, seq2):
    """
    Complexity: O(len(seq1*len(seq2))

    :param seq1:
    :param seq2:
    :return:
    """
    gap_score = -2
    match_score = 1
    mismatch_score = -1

    rows = len(seq1) + 1
    cols = len(seq2) + 1
    matrix = []

    for row in range(rows):
        matrix.append(list())
        for col in range(cols):
            matrix[row].append(0)

    #initialize first column
    for row in range(rows):
        matrix[row][0] = row * gap_score

    #initalize first row
    for col in range(cols):
        matrix[0][col] = col * gap_score

    for row in range(1, rows):
        for col in range(1, cols):
            diag = matrix[row - 1][col - 1]     # MATCH
            top = matrix[row - 1][col]          # DELETE
            left = matrix[row][col - 1]         # INSERT
            if seq1[row - 1] == seq2[col - 1]:
                diag_score = match_score
            else:
                diag_score = mismatch_score
            max_val = max(top + gap_score, diag + diag_score, left + gap_score)
            matrix[row][col] = max_val

    for row in range(1, rows):
        for col in range(1, cols):
            pass

    seq_aligner.align(seq1, seq2)

    # for i in range (0, len(seq1)):
    #   for j in range (0, len(seq2)):
    #     Match = matrix[i-1][j-1] + S(Ai, Bj)
    #     Delete = matrix[i-1][j] + d
    #     Insert = matrix[i][j-1] + d
    #     matrix[i][j] = max(Match, Insert, Delete)

    return matrix[rows - 1][cols - 1]


def main():
    # ALGORITHM VARIABLES
    seq_length = 10
    max_mismatches = 2
    paired_end = False

    # ensure correct number of arguments
    if len(sys.argv) < 3:
        usage()
        sys.exit(1)

    ref_filename = sys.argv[1]
    reads_filename = sys.argv[2]

    ref_name, ref_genome = load_genome(ref_filename)
    reads_name, reads = load_reads(reads_filename, paired_end=paired_end)
    if ref_name != reads_name:
        error_die("Reference genome id and reads genome id do not match")

    lookup_table = create_lookup_table(ref_genome, seq_length)
    read_map = create_read_map(
        ref_genome,
        reads,
        lookup_table,
        subseq_length=seq_length,
        thresh=max_mismatches
    )

    with open(answer_file_name, "w") as answer_file:
        find_snps(answer_file, ref_genome, ref_name, read_map)
        # TODO find indels


def eval():
    studentAns = open(answer_file_name, "r")
    answerKey = open("ans_genomeW1.txt", "r")
    test = Eval.Eval(answerKey, studentAns)
    for key in test:
        print key + ' grade: ' + str(test[key])
    studentAns.close()
    answerKey.close()


# RUN MAIN
if __name__ == '__main__':
    needleman_wunsch_align(
        "ACTGG",
        "GCCTGG"
    )
    main()
    eval()