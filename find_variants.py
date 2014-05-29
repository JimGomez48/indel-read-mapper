"""
The main module used to re-sequence and find variants in a donor genome
"""

__author__ = 'James Gomez'

import sys
import time
from argparse import ArgumentParser

from paired_end_read import PairedEndRead
import indel
import snp
import eval


allele_alphabet = "ACGT"
answer_file_name = "myanswers_"
answer_key_name = "ans_"


# def error_die(msg):
#     sys.stderr.write("ERROR: " + msg + "\n")
#     sys.exit(1)


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
                    # chr_name = line.translate(None, ">\r\n")
                    # print "\tLoading chromosome \'" + chr_name + "\'..."
                    continue
                line = line.rstrip("\r\n")
                genome += line

    return name, genome


def load_reads(filename, paired_end):
    """
    Loads the reads file into main memory as a either list of single reads or a
    list of paired-end reads, depending on the 'paired-end' flag. Returns the
    list and the name of the genome corresponding to the reads.

    :param filename: A string. The name of the file to load the reads from
    :param paired_end: A boolean. Set to True to load reads as paired-end reads,
    or False to load the reads as single reads
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
                # chr_name = line.translate(None, ">\r\n")
                # print "\tChromosome \'" + chr_name + "\'..."
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
    Creates a hash table to look up positions of allele sub-sequences of
    specified length within the genome. Returns a dictionary as the lookup
    hash-table, with string sub-sequences as keys and a list of positions as the
    values.
    Complexity: O(l/m) where l is length of the ref-genome, m is seq-length

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

    print "Created table of size: " + str(len(lookup_table)) + "keys"
    return lookup_table


def main():
    global answer_key_name
    global answer_file_name
    default_name = 'genome1'

    arg_parser = ArgumentParser(
        description="Performs resequencing of a donor genome with respect to a "
                    "reference genome and, optionally, identifies variations "
                    "within the donor genome. The resequenced genome will be "
                    "written to a file named \'donor_<genome_id>.txt\' and the "
                    "variants will be written to a file named "
                    "\'vars_<genome_id>.txt\'. The genome_id will be taken "
                    "from the first line of each file containing the '>' "
                    "delimiter (default genome_id=\'genome1\')."
    )
    arg_parser.add_argument(
        "ref_file",
        type=str,
        help="The file containing the reference genome")
    arg_parser.add_argument(
        "reads_file",
        type=str,
        help="The file containing the reads from the donor genome")
    arg_parser.add_argument(
        '-i', '--indels',
        help="Search for indels within the donor genome",
        action="store_true")
    arg_parser.add_argument(
        '-s', '--snps',
        help="Search for SNPs within the donor genome",
        action="store_true")
    args = arg_parser.parse_args()
    ref_file_name = args.ref_file
    reads_file_name = args.reads_file
    find_indels = args.indels
    find_snps = args.snps
    print "\nref file:    " + ref_file_name
    print "reads file:  " + reads_file_name
    print "Find Indels: " + str(find_indels)
    print "find SNPs:   " + str(find_snps)

    ref_name, ref_genome = load_genome(ref_file_name)
    if ref_name is None:
        ref_name = default_name
    reads_name, reads = load_reads(reads_file_name, paired_end=False)
    if reads_name is None:
        reads_name = default_name
    coverage = float(len(reads) * 50) / len(ref_genome)
    print "Read Coverage: " + str(coverage) + "x"
    if ref_name != reads_name:
        sys.stderr.write("WARNING: reference and donor genome_ids do not match")
    answer_key_name += (str(ref_name) + ".txt")
    answer_file_name += (str(ref_name) + ".txt")

    ########################## START MAIN ALGORITHM ############################
    # Algorithm Parameters
    seq_length = 10
    min_align_score = 0.4
    snp_thresh = 0.6
    local_alignment = True

    start_time = time.clock()
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

    ########################## WRITE RESULTS TO FILE ###########################
    with open(answer_file_name, "w") as answer_file:
        answer_file.write(">" + reads_name + "\n")

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
            reads,
            lookup_table,
            thresh=snp_thresh
        )
    total_time = time.clock() - start_time
    print "RUN-TIME"
    print "========"
    print "\tSeconds: " + str(total_time)
    print "\tMinutes: " + str(total_time / 60.0)

    run_eval()
    print "DONE\n"


def run_eval():
    try:
        with open(answer_file_name, "r") as student_ans:
            with open(answer_key_name, "r") as answer_key:
                print "EVALUATION"
                print "=========="
                grades = eval.eval(answer_key, student_ans)
                sum = 0
                keys = grades.keys()
                for k in keys:
                    grade = grades[k]
                    if grade > 0:
                        print "\t" + str(k) + ":\t " + str(grade)
                        sum += grades[k]
                print "\tOVERALL: " + str(sum / len(keys))
    except IOError as e:
        sys.stderr.write(
            "Couldn't open answer key \'" + answer_key_name + "\'\n")
        sys.stderr.write(e.message + "\n")


if __name__ == '__main__':
    main()