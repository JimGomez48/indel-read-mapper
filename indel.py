"""
Module for finding insertions and deletions in a donor genome with respect to a
reference genome
"""

__author__ = 'James'

import seq_align


def find_indels_read_map(ref_genome, reads, lookup_table, subseq_length, min_score):
    print "Finding indels and building read-map..."
    read_length = 50
    read_map = {}
    insert_map = {}
    delete_map = {}
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
                local=False
            )
            if current[SCORE] > best[SCORE] and current[SCORE] > min_score:
                best = current
                best_pos = p

        if not best_pos is None:
            # get indels from alignment
            ins, dels = __find_indels__(current[REF], current[TEST], best_pos)
            if not ins is None:   # Inserts
                for i in ins:
                    try:
                        insert_map[i[0]].append(i[1])
                    except KeyError:
                        insert_map[i[0]] = [i[1]]
            if not dels is None:   # Deletes
                for d in dels:
                    try:
                        delete_map[d[0]].append(d[1])
                    except KeyError:
                        delete_map[d[0]] = [d[1]]

            # Add aligned read to read map
            for i in range(0, len(best[TEST])):
                if best[TEST][i] == '-':  # skip gaps
                    continue
                try:
                    read_map[best_pos + i].append(str(best[TEST][i]))
                except KeyError:
                    read_map[best_pos + i] = [str(best[TEST][i])]

    inserts = __resolve_inserts__(insert_map)
    deletes = __resolve_deletes__(delete_map)

    return read_map, inserts, deletes


def __find_indels__(ref_seq, test_seq, start_pos):
    inserts = None
    deletes = None

    # FIND INSERTIONS
    insert_pos = position = start_pos
    inside = False
    insert_found = False
    current_insert = ""
    for i in range(len(ref_seq)):
        ref_c = ref_seq[i]
        test_c = test_seq[i]
        if not inside:
            if ref_c == '-':
                continue
            else:
                inside = True
        else:
            if ref_c == '-':
                if not insert_found:
                    insert_found = True
                    insert_pos = position
                current_insert += test_c
            else:
                if not current_insert == "":
                    if inserts is None:
                        inserts = []
                    inserts.append((current_insert, insert_pos))
                    current_insert = ""
                insert_found = False
                position += 1

    # FIND DELETIONS
    delete_pos = position = start_pos
    inside = False
    insert_found = False
    current_delete = ""
    for i in range(len(test_seq)):
        ref_c = ref_seq[i]
        test_c = test_seq[i]
        if not inside:
            if test_c == '-':
                continue
            else:
                inside = True
        else:
            if test_c == '-':
                if not insert_found:
                    insert_found = True
                    delete_pos = position
                current_delete += ref_c
            else:
                if not current_delete == "":
                    if deletes is None:
                        deletes = []
                    deletes.append((current_delete, delete_pos))
                    current_delete = ""
                insert_found = False
                position += 1

    return inserts, deletes


def __resolve_inserts__(insert_map):
    inserts = []
    return inserts


def __resolve_deletes__(delete_map):
    deletes = []
    return deletes