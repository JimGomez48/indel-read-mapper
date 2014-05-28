"""
Module used for alignment of two sequences. Can perform both global and local
alignment
"""

__author__ = 'James Gomez'

import sys


class MatrixEntry:
    DIAG = "DIAG"
    UP = "UP"
    LEFT = "LEFT"

    def __init__(self, score=0, parent=None):
        self.score = score
        self.parent = parent


def align(ref_seq, test_seq, local=False):
    """
    Aligns the two sequences with minimal error. If local=True, uses the
    Waterman-Smith local-alignment algorithm. Otherwise, uses the Needleman-Wunsch
    global-alignment algorithm

    both quadratic time  O(n^2) or O(n*m)
    """
    if local:
        matrix, score = __smith_waterman__(ref_seq, test_seq)
    else:
        matrix, score = __needleman_wunsch__(ref_seq, test_seq)

    ref_align, test_align = __align__(matrix, ref_seq=ref_seq, test_seq=test_seq)

    score = float(score) / len(ref_align)  # Normalize the score
    # if 30 < score < 50:
    # if 0.5 < score < 1.0:
        # print "ref-align:  " + ref_align
        # print "test-align: " + test_align
        # print "score: " + str(score)
        # print

    return ref_align, test_align, score


def __align__(matrix, ref_seq, test_seq):
    # perform alignment
    ref_align = ""
    test_align = ""
    i = len(matrix) - 1
    j = len(matrix[0]) - 1
    while i > 0 or j > 0:
        current = matrix[i][j]
        next = current.parent

        if next is None:
            if i == 0:
                ref_align += '-'
                test_align += test_seq[j - 1]
                j -= 1
            elif j == 0:
                ref_align += ref_seq[i - 1]
                test_align += '-'
                i -= 1
        elif i > 0 and j > 0 and next == MatrixEntry.DIAG:
            ref_align += ref_seq[i - 1]
            test_align += test_seq[j - 1]
            i -= 1
            j -= 1
        elif i > 0 and next == MatrixEntry.UP:
            # UP
            ref_align += ref_seq[i - 1]
            test_align += "-"
            i -= 1
        elif j > 0 and next == MatrixEntry.LEFT:
            # LEFT
            ref_align += "-"
            test_align += test_seq[j - 1]
            j -= 1

    ref_align = ref_align[::-1]
    test_align = test_align[::-1]

    return ref_align, test_align


def __needleman_wunsch__(ref_seq, test_seq):
    """
    Calculates the global alignment matrix

    Complexity: O(len(seq1*len(seq2))
    """
    match_score = 1
    mismatch_score = -1
    gap_score = -1

    rows = len(ref_seq) + 1
    cols = len(test_seq) + 1
    matrix = []

    # Create the matrix
    for row in range(0, rows):
        matrix.append(list())
        for col in range(cols):
            matrix[row].append(MatrixEntry())

    # Initialize first column
    for row in range(0, rows):
        matrix[row][0].score = row * gap_score

    # Initialize first row
    for col in range(0, cols):
        matrix[0][col].score = col * gap_score

    # Generate scores and parent pointers
    for row in range(1, rows):
        for col in range(1, cols):
            diag_score = matrix[row - 1][col - 1].score          # MATCH
            top_score = matrix[row - 1][col].score + gap_score   # DELETE
            left_score = matrix[row][col - 1].score + gap_score  # INSERT
            if ref_seq[row - 1] == test_seq[col - 1]:
                diag_score += match_score
            else:
                diag_score += mismatch_score

            best_score = max(top_score, diag_score, left_score)
            if best_score == diag_score:    # MATCH
                matrix[row][col].score = diag_score
                matrix[row][col].parent = MatrixEntry.DIAG
            elif best_score == left_score:  # INSERT
                matrix[row][col].score = left_score
                matrix[row][col].parent = MatrixEntry.LEFT
            else:                           # DELETE
                matrix[row][col].score = top_score
                matrix[row][col].parent = MatrixEntry.UP
    # print_matrix(matrix)
    align_score = matrix[rows - 1][cols - 1].score

    return matrix, align_score


def __smith_waterman__(ref_seq, test_seq):
    """
    Calculates the local alignment matrix

    Complexity: O(len(seq1*len(seq2))
    """
    match_score = 1
    mismatch_score = -1
    gap_score = -1

    rows = len(ref_seq) + 1
    cols = len(test_seq) + 1
    matrix = []

    # Create the matrix
    for row in range(0, rows):
        matrix.append(list())
        for col in range(cols):
            matrix[row].append(MatrixEntry())

    # Initialize first column
    for row in range(0, rows):
        matrix[row][0].score = 0

    # Initialize first row
    for col in range(0, cols):
        matrix[0][col].score = 0

    # Generate scores and parent pointers
    for row in range(1, rows):
        for col in range(1, cols):
            diag_score = matrix[row - 1][col - 1].score          # MATCH
            top_score = matrix[row - 1][col].score + gap_score   # DELETE
            left_score = matrix[row][col - 1].score + gap_score  # INSERT
            if ref_seq[row - 1] == test_seq[col - 1]:
                diag_score += match_score
            else:
                diag_score += mismatch_score

            best_score = max(top_score, diag_score, left_score)
            if best_score == diag_score:    # MATCH
                matrix[row][col].score = max(0, diag_score)
                matrix[row][col].parent = MatrixEntry.DIAG
            elif best_score == left_score:  # INSERT
                matrix[row][col].score = max(0, left_score)
                matrix[row][col].parent = MatrixEntry.LEFT
            else:                           # DELETE
                matrix[row][col].score = max(0, top_score)
                matrix[row][col].parent = MatrixEntry.UP
    # print_matrix(matrix)
    align_score = matrix[rows - 1][cols - 1].score

    return matrix, align_score


def print_matrix(matrix):
    for i in range(0, len(matrix)):
        sys.stdout.write("[ ")
        for j in range(0, len(matrix[0])):
            sys.stdout.write(str(matrix[i][j].score).zfill(3) + " ")
        sys.stdout.write("]\n")
    print