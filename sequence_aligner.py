__author__ = 'james'

import sys


class SequenceAligner:

    class MatrixEntry:
        def __init__(self, score=0, parent=None):
            self.score = score
            self.parent = parent

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
                diag_score = self.matrix[row - 1][col - 1].score          # MATCH
                top_score = self.matrix[row - 1][col].score + gap_score   # DELETE
                left_score = self.matrix[row][col - 1].score + gap_score  # INSERT
                if seq1[row - 1] == seq2[col - 1]:
                    diag_score += match_score
                else:
                    diag_score += mismatch_score

                best_score = max(top_score, diag_score, left_score)
                if best_score == top_score:
                    self.matrix[row][col].score = top_score
                    self.matrix[row][col].parent = self.matrix[row - 1][col]
                elif best_score == left_score:
                    self.matrix[row][col].score = left_score
                    self.matrix[row][col].parent = self.matrix[row][col - 1]
                else:
                    self.matrix[row][col].score = diag_score
                    self.matrix[row][col].parent = self.matrix[row - 1][col - 1]

        self.__print_matrix__()

        return self.matrix

    def __print_matrix__(self):
        for i in range(0, len(self.matrix)):
            sys.stdout.write("[ ")
            for j in range(0, len(self.matrix[0])):
                sys.stdout.write(str(self.matrix[i][j].score) + " ")
            sys.stdout.write("]\n")
        print