class PossibilityResult:
    index_seq_snd = None
    index_seq_fst = None
    profile_result = None

    def __init__(self, profile_result, index_seq_snd, index_seq_fst=None):
        self.index_seq_fst = index_seq_fst
        self.index_seq_snd = index_seq_snd
        self.profile_result = profile_result

    def score(self):
        return self.profile_result.score

    def not_contain_index(self, index):
        return index != self.index_seq_snd and \
               (self.index_seq_fst is None or index != self.index_seq_fst)