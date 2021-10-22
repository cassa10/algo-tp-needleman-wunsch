import numpy as np

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
            if match_score > gap_score:
                table[i][j] = table[i-1][j-1] + match_score
            else:
                gap_en_a = table[i-1][j]
                gap_en_b = table[i][j-1]
                table[i][j] = max(gap_en_a, gap_en_b) + 0

    return [table[row-1, col-1], table]

def traceback(tabla_sol, seqA, seqB):
    i = len(seqA)
    j = len(seqB)

    aln_a = ""
    aln_b = ""

    while i > 0 and j > 0:
        ma_mm = tabla_sol[i][j]
        gap_a = tabla_sol[i-1][j]
        gap_b = tabla_sol[i][j-1]

        if (ma_mm > gap_a) and (ma_mm > gap_b):
            aln_a = seqA[i - 1] + aln_a
            aln_b = seqB[j - 1] + aln_b
        elif gap_a <= gap_b:
            aln_a = "-" + aln_a
            aln_b = seqB[j - 1] + aln_b
        else:
            aln_a = seqB[i - 1] + aln_a
            aln_b = "-" + aln_b
        i = i - 1
        j = j - 1

    return [aln_a, aln_b]

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