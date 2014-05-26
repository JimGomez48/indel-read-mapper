__author__ = 'james'

import sys


class SequenceAligner:

    class MatrixEntry:
        DIAG = "DIAG"
        UP = "UP"
        LEFT = "LEFT"

        def __init__(self, score=0, parent=None):
            self.score = score
            self.parent = parent

    def __init__(self):
        self.matrix = []

    def align(self, ref_seq, test_seq):
        """
        Cubic time algorithm
        """
        self.matrix = self.__needleman_wunsch_matrix__(ref_seq, test_seq)

        # perform alignment
        ref_align = ""
        test_align = ""
        i = len(self.matrix) - 1
        j = len(self.matrix[0]) - 1
        while i > 0 or j > 0:
            current = self.matrix[i][j]
            next = current.parent
            if i > 0 and j > 0 and next == self.MatrixEntry.DIAG:
                ref_align += ref_seq[i - 2]
                test_align += test_seq[j - 2]
                i -= 1
                j -= 1
            elif i > 0 and next == self.MatrixEntry.UP:
                # UP
                ref_align += ref_seq[i - 2]
                test_align += "-"
                i -= 1
            elif j > 0 and next == self.MatrixEntry.LEFT:
                # LEFT
                ref_align += "-"
                test_align += ref_seq[j - 2]
                j -= 1
            else:
                pass



        ref_align = ref_align[::-1]
        test_align = test_align[::-1]
        print ref_align
        print test_align

    def __needleman_wunsch_matrix__(self, ref_seq, test_seq):
        """
        Complexity: O(len(seq1*len(seq2))
        """
        gap_score = -1
        match_score = 1
        mismatch_score = -1

        rows = len(ref_seq) + 1
        cols = len(test_seq) + 1
        self.matrix = []

        # Create the matrix
        for row in range(0, rows):
            self.matrix.append(list())
            for col in range(cols):
                self.matrix[row].append(self.MatrixEntry())
        # self.__print_matrix__()

        # Initialize first column
        for row in range(0, rows):
            self.matrix[row][0].score = row * gap_score
        # self.__print_matrix__()

        # Initialize first row
        for col in range(0, cols):
            self.matrix[0][col].score = col * gap_score
        # self.__print_matrix__()

        # Generate scores and parent pointers
        for row in range(1, rows):
            for col in range(1, cols):
                diag_score = self.matrix[row - 1][col - 1].score          # MATCH
                top_score = self.matrix[row - 1][col].score + gap_score   # DELETE
                left_score = self.matrix[row][col - 1].score + gap_score  # INSERT
                if ref_seq[row - 1] == test_seq[col - 1]:
                    diag_score += match_score
                else:
                    diag_score += mismatch_score

                best_score = max(top_score, diag_score, left_score)
                if best_score == diag_score:    # MATCH
                    self.matrix[row][col].score = diag_score
                    # self.matrix[row][col].parent = [row - 1, col - 1]
                    self.matrix[row][col].parent = self.MatrixEntry.DIAG
                elif best_score == left_score:  # INSERT
                    self.matrix[row][col].score = left_score
                    # self.matrix[row][col].parent = [row, col - 1]
                    self.matrix[row][col].parent = self.MatrixEntry.LEFT
                else:                           # DELETE
                    self.matrix[row][col].score = top_score
                    # self.matrix[row][col].parent = [row - 1, col]
                    self.matrix[row][col].parent = self.MatrixEntry.UP

        # self.__print_matrix__()

        return self.matrix

    def __print_matrix__(self):
        for i in range(0, len(self.matrix)):
            sys.stdout.write("[ ")
            for j in range(0, len(self.matrix[0])):
                sys.stdout.write(str(self.matrix[i][j].score) + " ")
            sys.stdout.write("]\n")
        print