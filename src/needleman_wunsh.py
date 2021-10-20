import numpy as np

class op(Enum):
    GAP_A = 1
    GAP_B = 2
    MA_MM = 3

def init(s1, s2):
    score_mtx = init_score_mtx(s1, s2)
    print(score_mtx)
    return nw(s1, s2, score_mtx)

def nw(seqA_str, seqB_str, score_mtx, gap_score = 0):
    table = init_mem(seqA_str, seqB_str)
        
    #TODO: Fill matrix
    for i in range(1, len(seqA_str)):
        for j in range(1, len(seqB_str)):
            return 2


def init_mem(s1, s2):
    row = len(s1) + 2
    col = len(s2) + 2
    matrix = np.repeat(0.0, col * row).reshape(col, row)

    matrix[0][1] = ""
    matrix[1][1] = ""
    for i in range(2, row):
        matrix[0][i] = s1[i-2]
    
    for i in range(2, col):
        matrix[i][0] = s2[i-2]

    return matrix

def init_score_mtx(s1, s2):
    """ example:
    [
        [ ,  A,  T,  C,  T],
        [A,  1,  0,  0,  0],
        [T,  0,  1,  0,  1],
        [A,  1,  0,  0,  0]
    ]
    """
    row = len(s1) + 1
    col = len(s2) + 1
    matrix = np.repeat(0.0, col * row).reshape(col, row)

    for i in range(1, row):
        matrix[0][i] = s1[i-1]
    
    for i in range(1, col):
        matrix[i][0] = s2[i-1]
    
    for i in range(1, row): 
        for j in range(1, col): 
            matrix[i][j] = scoreOf(matrix[0][row], matrix[col][0])

    return matrix

def scoreOf(c1, c2):
    if c1 == c2:
        return 1 
    return 0