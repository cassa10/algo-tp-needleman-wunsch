
class Profile:

    values = {
        "A": [],
        "C": [],
        "G": [],
        "T": [],
    }
    total = 0

    # alignment : [[Char]]
    def __init__(self, alignment):
        self.alignment = alignment
        self.total = len(alignment)
        sumValues = {"A": 0, "C": 0, "G": 0, "T": 0}

        for i in range(0, len(alignment)-1):
            for j in range(0, len(alignment)-1):
                self.increment_value(sumValues, alignment[i][j])
            self.store_values(sumValues)
            self.reset_sum_values(sumValues)

    def store_values(self, sumValues):
        for seq_key, seq_sum in sumValues:
            percentage_sum = seq_sum / self.total
            self.values[seq_key].append(percentage_sum)

    def increment_value(self, sumValues, incValue):
        sumValues[incValue] = sumValues[incValue] + 1

    def reset_sum_values(self, sumValues):
        for seq_key in sumValues:
            sumValues[seq_key] = 0

    def score(self, score_mtx, pos, seq_param):
        score_pos = 0
        for seq_key in self.values:
            score_pos = score_pos + self.values[seq_key][pos] * score_mtx[seq_param + seq_key]
        return score_pos

    def getValue(self, seq, pos):
        return self.values[seq][pos]

