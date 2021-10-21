import numpy as np
from enum import Enum

class op(Enum):
    GAP_A = 1
    GAP_B = 2
    MA_MM = 3


def init(s1, s2):
    score_mtx = init_score_mtx(s1, s2)
    gap_score = 0
    return nw(s1, s2, score_mtx, gap_score)


def nw(seqA_str, seqB_str, score_mtx, gap_score=0):
    # _ = string vacio
    # seqA = "_CGTAA"
    # seqB = "_ATTCA"
    row = len(seqA_str) + 1
    col = len(seqB_str) + 1
    table = np.repeat(None, row * col).reshape(row, col)
    """ example:
       [
            _   A   T   C   T
          _[0,  0,  0,  0,  0],
          A[0,  0,  1,  0,  1],
          T[0,  1,  0,  0,  0],
          C[0,  1,  0,  0,  0]
       ]
    """

    table[0][0] = 0.0

    for i in range(1, row):
        table[i][0] = 0.0

    for j in range(1, col):
        table[0][j] = 0.0

    for i in range(1, row):
        for j in range(1, col):
            match_score = score_mtx[i-1][j-1]
            prev = table[i - 1][j - 1]
            if match_score > gap_score and prev is not None:
                table[i][j] = table[i-1][j-1] + match_score
            else: 
                table[i][j] = 0.0

    return table


def init_score_mtx(s1, s2):
    """ example:
    [   A   T   C    T
       A[1,  0,  0,  0],
       T[0,  1,  0,  1],
       A[1,  0,  0,  0]
    ]
    """
    row = len(s1)
    col = len(s2)
    matrix = np.repeat(0.0, row * col).reshape(row, col)

    for i in range(0, row):
        for j in range(0, col):
            matrix[i][j] = scoreOf(s1[i], s2[j])
    return matrix


def scoreOf(c1, c2):
    if c1 == c2:
        return 1
    return 0

if __name__ == '__main__':
    s1 = "ATC"
    s2 = "ATCT"
    #print(init_score_mtx(s1, s2))
    print(init(s1, s2))