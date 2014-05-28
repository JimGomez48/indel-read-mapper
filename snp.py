"""
Module for finding and writing SNPs
"""

__author__ = 'James Gomez'

__allele_alphabet__ = "ACGT"


def find_and_write_snps(answer_file, ref_genome, reads, lookup_table, thresh=0.6):
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
    read_map = __map_reads__(ref_genome, reads, lookup_table)
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


def __map_reads__(ref_genome, reads, lookup_table):
    hash_key_length = 10
    thresh = 2
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

        read_position = get_best_read_position(ref_genome, read, positions, thresh)
        if not read_position is None:
            for i in range(0, len(read)):
                try:
                    read_map[read_position + i].append(str(read[i]))
                except KeyError:
                    read_map[read_position + i] = [str(read[i])]

    return read_map


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
