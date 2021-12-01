from Bio import pairwise2, SeqIO
import needleman_wunsh


def print_test(n):
    print(f"[Test {n}]")


def print_test_pass(n, is_passed):
    if is_passed:
        print(f"[Test {n}] - PASSED")
    else:
        print(f"[Test {n}] - ERROR")
    print(f"------------------------------")


def parse_fasta_and_validate(file_dir, score_mtx_file):
    f_aln = SeqIO.parse(file_dir, 'fasta')
    fas_A = next(f_aln)
    fas_B = next(f_aln)

    seq_A = str(fas_A.seq)
    seq_B = str(fas_B.seq)

    aln = pairwise2.align.globalxx(seq_A, seq_B)
    aln_score = aln[0].score

    results = needleman_wunsh.init(seq_A, seq_B, score_mtx_file)
    score = results[0]
    alignments = results[1]

    alignment_seq_A = alignments[0]
    alignment_seq_B = alignments[1]

    print("Test: ")
    print(f"seq_A: {seq_A}")
    print(f"seq_B: {seq_B}")
    print(f"nw(seq_A) => {alignment_seq_A}")
    print(f"nw(seq_B) => {alignment_seq_B}")

    return aln_score == score


def test_1():
    return parse_fasta_and_validate("../resources/ab.fasta", "../resources/test_nw_score_matrix")


def test_2():
    return parse_fasta_and_validate("../resources/sample.fasta", "../resources/test_nw_score_matrix")


def test_3():
    return parse_fasta_and_validate("../resources/ae.fasta", "../resources/test_nw_score_matrix")


def run():
    tests_nw = [lambda: test_1(), lambda: test_2(), lambda: test_3()]
    for i, test in enumerate(tests_nw):
        n = i + 1
        print_test(n)
        print_test_pass(n, test())
