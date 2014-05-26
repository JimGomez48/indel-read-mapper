__author__ = 'james'


class PairedEndRead:
    def __init__(self, seq1=None, seq2=None):
        self.seq1 = seq1
        self.seq2 = seq2
        self.p1 = None
        self.p2 = None

    def set_positions(self, p1=0, p2=0):
        self.p1 = p1
        self.p2 = p2

    def gap_distance(self):
        return abs(self.p2 - self.p1)