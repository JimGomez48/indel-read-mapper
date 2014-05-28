"""
The main module used to re-sequence and find variants in a donor genome
"""

__author__ = 'James Gomez'

import sys
import math

from paired_end_read import PairedEndRead
import indel
import snp
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


def main():
    global answer_key_name
    if len(sys.argv) < 2:  # ensure correct args
        usage()
        sys.exit(1)

    ##################################### SETUP #####################################
    print
    ref_filename = "ref_" + str(sys.argv[1]) + ".txt"
    reads_filename = "reads_" + str(sys.argv[1]) + ".txt"

    ref_name, ref_genome = load_genome(ref_filename)
    reads_name, reads = load_reads(reads_filename, paired_end=False)
    coverage = float(len(reads) * 50) / len(ref_genome)
    print "\tCoverage: " + str(coverage) + "x"
    if ref_name != reads_name:
        error_die("Reference genome id and reads genome id do not match")
    answer_key_name += (str(ref_name) + ".txt")

    ############################# START MAIN ALGORITHM ##############################

    # Algorithm Parameters
    seq_length = 10
    min_align_score = 0.4
    snp_thresh = 0.6
    local_alignment = True

    lookup_table = create_lookup_table(ref_genome, seq_length)
    read_map, inserts, deletes = indel.find_indels_read_map(
        ref_genome,
        reads,
        lookup_table,
        subseq_length=seq_length,
        min_score=min_align_score,
        coverage=coverage,
        local=local_alignment
    )

    ############################ WRITE RESULTS TO FILE ##############################
    with open(answer_file_name, "w") as answer_file:
        answer_file.write(">" + ref_name + "\n")

        # Write inserts
        answer_file.write(">INSERT:")
        for i in inserts:
            answer_file.write("\n1," + str(i[0]) + "," + str(i[1]))
        answer_file.write("\n")

        # Write deletes
        answer_file.write(">DELETE:")
        for d in deletes:
            answer_file.write("\n1," + str(d[0]) + "," + str(d[1]))
        answer_file.write("\n")

        # Find and write SNPs
        snp.find_and_write_snps(
            answer_file,
            ref_genome,
            read_map,
            thresh=snp_thresh
        )

    print "DONE\n"


def run_eval():
    try:
        with open(answer_file_name, "r") as student_ans:
            with open(answer_key_name, "r") as answer_key:
                print "EVALUATION"
                print "=========="
                eval.Eval(answer_key, student_ans)
    except IOError as e:
        sys.stderr.write("Couldn't open answer key \'" + answer_key_name + "\'\n")
        sys.stderr.write(e.message + "\n")


if __name__ == '__main__':
    main()
    run_eval()