class Solution:
    score = 0
    alignment = []

    def __init__(self, score, alignment):
        self.score = score
        self.alignment = alignment

    def get_immutable_aln(self):
        return self.alignment.copy()