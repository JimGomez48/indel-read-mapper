__author__ = 'james'

import sys


class SequenceAligner:

    class MatrixEntry:
        def __init__(self, score=0, parents=[]):
            self.score = score
            self.parents = parents

    def __init__(self):
        self.matrix = []

    def align(self, seq1, seq2):
        self.matrix = self.__needleman_wunsch_matrix__(seq1, seq2)
        # TODO perform alignment

    def __needleman_wunsch_matrix__(self, seq1, seq2):
        """
        Complexity: O(len(seq1*len(seq2))
        """
        gap_score = -1
        match_score = 1
        mismatch_score = -1

        rows = len(seq1) + 1
        cols = len(seq2) + 1
        self.matrix = []

        # Create the matrix
        for row in range(0, rows):
            self.matrix.append(list())
            for col in range(cols):
                self.matrix[row].append(self.MatrixEntry())
        self.__print_matrix__()

        # Initialize first column
        for row in range(0, rows):
            self.matrix[row][0].score = row * gap_score
        self.__print_matrix__()

        # Initialize first row
        for col in range(0, cols):
            self.matrix[0][col].score = col * gap_score
        self.__print_matrix__()

        # Generate scores and parent pointers
        for row in range(1, rows):
            for col in range(1, cols):
                diag_score = self.matrix[row - 1][col - 1].score     # MATCH
                top_score = self.matrix[row - 1][col].score          # DELETE
                left_score = self.matrix[row][col - 1].score         # INSERT
                if seq1[row - 1] == seq2[col - 1]:
                    diag_score += match_score
                else:
                    diag_score += mismatch_score
                self.matrix[row][col].score = max(
                    top_score + gap_score,
                    diag_score,
                    left_score + gap_score
                )

        self.__print_matrix__()

        return self.matrix

    def __print_matrix__(self):
        for i in range(0, len(self.matrix)):
            sys.stdout.write("[ ")
            for j in range(0, len(self.matrix[0])):
                sys.stdout.write(str(self.matrix[i][j].score) + " ")
            sys.stdout.write("]\n")
        print