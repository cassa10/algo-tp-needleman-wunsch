import numpy as np
from enum import Enum

class Op(Enum):
    GAP_A = 1
    GAP_B = 2
    MA_MM = 3

def init(s1, s2):
    score_mtx = init_score_mtx()
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

    table[0][0] = gap_score

    for i in range(1, row):
        table[i][0] = table[i-1][0] + gap_score

    for j in range(1, col):
        table[0][j] = table[0][j-1] + gap_score

    for i in range(1, row):
        for j in range(1, col):
            tmp_match = seqA_str[i-1] + seqB_str[j-1]
            ma_mm_score = score_mtx[tmp_match]
            #TODO: Corregir y agregar otra tabla para hacer traceback
            if ma_mm_score > gap_score:
                table[i][j] = table[i-1][j-1] + ma_mm_score
            else:
                gap_en_a = table[i-1][j]
                gap_en_b = table[i][j-1]
                table[i][j] = max(gap_en_a, gap_en_b) + gap_score

    return [table[row-1, col-1], table]

def traceback(tabla_sol, seqA, seqB):
    #TODO: Utilizar el enum con otra tabla para hacer el recorrido inverso)
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
            i = i - 1
            j = j - 1
        elif gap_a <= gap_b:
            aln_a = "-" + aln_a
            aln_b = seqB[j - 1] + aln_b
            j = j-1
        else:
            aln_a = seqA[i - 1] + aln_a
            aln_b = "-" + aln_b
            i = i - 1

    return [aln_a, aln_b]

def init_score_mtx():
    """ example: MAAAAAL
    [   A   T   C    T
       A[1,  0,  0,  0],
       T[0,  1,  0,  1],
       A[1,  0,  0,  0]
    ]
    """
    #TODO:
    # Leer archivo de texto y llevarlo a un diccionario (MAP) en base a la tabla
    # Combinaciones de letras (AA=5 o AC=-1 o AG=-2, AT=-3, CC=4, CG=-1, etc en cada una)
    #values = ['A','C','G','T']
    score = {
        "AA": 1,
        "CC": 1,
        "GG": 1,
        "TT": 1,
        "AC": -2,
        "AG": -1,
        "AT": -4,
        "CA": -3,
        "CG": -5,
        "CT": 0,
        "GA": -1,
        "GC": -5,
        "GT": -3,
        "TA": -4,
        "TC": 0,
        "TG": -3,
    }

    return score