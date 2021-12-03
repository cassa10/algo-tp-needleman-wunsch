import numpy as np
from enum import Enum

from src.profile import Profile
from src.score import Score


class Op(Enum):
    GAP_S = 1
    GAP_P = 2
    MA_MM = 3


def init(sequences, score_file_dir, gap_penalty=0):
    if len(sequences) <= 1:
        return sequences
    score_mtx = Score(score_file_dir).score
    seq_fst = sequences.pop(0)
    profile = Profile(seq_fst)
    # Mutates profile for each iteration
    for seq in sequences:
        nw_with_profile(profile, seq, score_mtx, gap_penalty)
    return profile.score, profile.alignment


def nw_with_profile(profile, seq, score_mtx, gap_penalty):
    lenS = len(seq) + 1
    lenP = profile.length() + 1
    table = np.repeat(0, lenS * lenP).reshape(lenS, lenP)
    profile.init_traceback_table(seq)

    for i in range(1, lenS):
        table[i][0] = get_score_probs(score_mtx, gap_penalty, seq[i - 1], {"-": [1]}) + table[i - 1][0]
        profile.traceback_table[i][0] = (Op.GAP_P, ["-"] * profile.total + [seq[i - 1]])

    for j in range(1, lenP):
        table[0][j] = get_score_probs(score_mtx, gap_penalty, "-", profile.get_prob_values(j-1)) + table[0][j-1]
        profile.traceback_table[0][j] = (Op.GAP_S, profile.get_chars_col(j-1) + ["-"])

    for i in range(1, lenS):
        for j in range(1, lenP):
            # match
            table[i][j] = get_score_probs(score_mtx, gap_penalty, seq[i - 1], profile.get_prob_values(j-1)) + \
                          table[i-1][j-1]
            profile.traceback_table[i][j] = (Op.MA_MM, profile.get_chars_col(j-1))

            # gap en profile
            score_gap_p = get_score_probs(score_mtx, gap_penalty, seq[i - 1], {"-": [1]}) + table[i-1][j]
            if score_gap_p > table[i][j]:
                table[i][j] = score_gap_p
                profile.traceback_table[i][j] = (Op.GAP_P, ["-"] * profile.total + [seq[i - 1]])

            # gap en seq
            score_gap_seq = get_score_probs(score_mtx, gap_penalty, "-", profile.get_prob_values(j - 1)) + table[i][j-1]
            if score_gap_seq > table[i][j]:
                table[i][j] = score_gap_seq
                profile.traceback_table[i][j] = (Op.GAP_S, profile.get_chars_col(j-1) + ["-"])

    profile.score = table[lenP][lenS]
    profile.do_traceback(lenP, lenS)


def get_score_probs(score_mtx, gap_penalty, char, probs):
    score = 0
    for char_prob, prob in probs:
        score += prob * get_score(score_mtx, gap_penalty, char, char_prob)
    return score


def get_score(score_mtx, gap_penalty, char_x, char_y):
    if char_x == "-" and char_y == "-":
        return gap_penalty
    return score_mtx[char_x + char_y]
