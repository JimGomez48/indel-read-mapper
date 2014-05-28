"""
Module for finding and writing SNPs
"""

__author__ = 'James Gomez'

__allele_alphabet__ = "ACGT"


def find_and_write_snps(answer_file, ref_genome, read_map, thresh=0.6):
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


def get_consensus_allele(alleles):
    """
    Returns the allele with the highest count.
    Complexity: O(len(alleles))
    :rtype : str
    """
    allele_counts = {}
    for a in alleles:
        assert a in __allele_alphabet__
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
