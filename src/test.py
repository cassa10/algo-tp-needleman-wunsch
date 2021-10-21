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

    # aln_a, aln_b, \
    score = needleman_wunsh.init(seq_A, seq_B)

    return aln[0].score == score


def execute_tests():
    tests = [test_1()]
    for i, passed in enumerate(tests):
        print_test_pass(i + 1, passed)
