from Bio import pairwise2, SeqIO
import needleman_wunsh


def print_test_pass(n, isPassed):
    if isPassed:
        print(f"Test {n} - PASSED")
    else:
        print(f"Test {n} - ERROR")


def test_1():
    f_aln = SeqIO.parse("ab.fasta", 'fasta')
    fas_A = next(f_aln)
    fas_B = next(f_aln)

    seq_A = str(fas_A.seq)
    seq_B = str(fas_B.seq)

    aln = pairwise2.align.globalxx(seq_A, seq_B)
    aln_score = aln[0].score

    result = needleman_wunsh.init(seq_A, seq_B)
    score = result[0]

    return aln_score == score

def test_2():
    f_aln = SeqIO.parse("ab.fasta", 'fasta')
    fas_A = next(f_aln)
    fas_B = next(f_aln)

    seq_A = str(fas_A.seq)
    seq_B = str(fas_B.seq)

    aln = pairwise2.align.globalxx(seq_A, seq_B)
    aln_score = aln[0].score

    results = needleman_wunsh.init(seq_A, seq_B)
    score = results[0]
    table_result = results[1]

    str_results = needleman_wunsh.traceback(table_result, seq_A, seq_B)

    nw_A = str_results[0]
    nw_B = str_results[1]

    print("Test 2: ")
    print(f"seq_A: {seq_A}")
    print(f"seq_B: {seq_B}")
    print(f"nw(seq_A) => {nw_A}")
    print(f"nw(seq_B) => {nw_B}")

    return aln_score == score


def execute_tests():
    tests = [lambda : test_1(), lambda : test_2()]
    for i, test in enumerate(tests):
        print_test_pass(i + 1, test())
