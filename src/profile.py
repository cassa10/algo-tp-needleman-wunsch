import numpy as np

from src.msa_op import Op


def increment_value(sumValues, incValue):
    sumValues[incValue] = sumValues[incValue] + 1


def set_all_values(sumValues, value):
    for seq_key in sumValues:
        sumValues[seq_key] = value


class Profile:
    alignment = []
    values = {}
    total = 0
    gap_penalty = 0
    score = 0
    traceback_table = []

    def __init__(self, seq):
        self.alignment = [seq]
        self.total = len(self.alignment)
        self.values = {"-": [], "A": [], "C": [], "G": [], "T": []}
        sumValues = {"-": 0, "A": 0, "C": 0, "G": 0, "T": 0}
        for i in range(len(seq)):
            sumValues[seq[i]] = sumValues[seq[i]] + 1
            self.__store_values(sumValues)
            set_all_values(sumValues, 0)

    def __store_values(self, sumValues):
        for seq_key, seq_sum in sumValues.items():
            ls = self.values[seq_key]
            prob = seq_sum / self.total
            self.values[seq_key] = ls + [prob]

    # TODO: Refactor with store_values
    def __store_values_reverse(self, sumValues):
        for seq_key, seq_sum in sumValues.items():
            ls = self.values[seq_key]
            prob = seq_sum / self.total
            self.values[seq_key] = [prob] + ls

    def score_seq(self, score_mtx, pos, seq_param):
        score_pos = 0
        for seq_key in self.values:
            score_pos = score_pos + self.values[seq_key][pos] * score_mtx[seq_param + seq_key]
        return score_pos

    def length(self):
        return len(self.alignment[0])

    def init_traceback_table(self, seq):
        lenS = len(seq) + 1
        lenP = self.length() + 1
        self.traceback_table = np.repeat(None, lenS * lenP).reshape(lenS, lenP)

    def get_prob_values(self, pos_col):
        res = {"-": 0, "A": 0, "C": 0, "G": 0, "T": 0}
        for seq_key, prob_values in self.values.items():
            res[seq_key] = prob_values[pos_col]
        return res

    def do_traceback(self, lenS, lenP):
        # Traceback_table = [(Enum.Op, alignment: [col0, col1, ..., col_n])]
        self.total = self.total + 1
        new_alignment = [""] * self.total
        set_all_values(self.values, [])
        sumColValues = {"-": 0, "A": 0, "C": 0, "G": 0, "T": 0}
        while lenS > 0 or lenP > 0:
            last_op, aln_column = self.traceback_table[lenS][lenP]
            for i in range(len(new_alignment)):
                current_char = aln_column[i]
                new_alignment[i] = current_char + new_alignment[i]
                sumColValues[current_char] = sumColValues[current_char] + 1

            self.__store_values_reverse(sumColValues)
            set_all_values(sumColValues, 0)

            if last_op == Op.MA_MM:
                lenP = lenP - 1
                lenS = lenS - 1
            elif last_op == Op.GAP_S:
                lenP = lenP - 1
            else:
                lenS = lenS - 1

        self.alignment = new_alignment

    def get_chars_col(self, pos_col):
        return [self.alignment[i][pos_col] for i in range(self.total)]

    def get_score_and_alignment(self):
        return self.score, self.alignment
