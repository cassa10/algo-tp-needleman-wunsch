from Bio import pairwise2
import needleman_wunsh
from file_parser import fasta_first_two_seq, score_matrix
from test_util import execute_tests

RESOURCES_FOLDER = "./resources"
DEFAULT_SCORE_MTX = f"{RESOURCES_FOLDER}/test_nw_score_matrix"


def parse_fasta_and_validate(file_dir, score_mtx):
    fas_A, fas_B, seq_A, seq_B = fasta_first_two_seq(file_dir)

    aln = pairwise2.align.globalxx(seq_A, seq_B)
    aln_score = aln[0].score

    score, alignments = needleman_wunsh.init(seq_A, seq_B, score_mtx)

    alignment_seq_A = alignments[0]
    alignment_seq_B = alignments[1]

    print("Test: ")
    print(f"seq_A: {seq_A}")
    print(f"seq_B: {seq_B}")
    print(f"nw(seq_A) => {alignment_seq_A}")
    print(f"nw(seq_B) => {alignment_seq_B}")

    return aln_score == score


def test_1():
    return parse_fasta_and_validate(f"{RESOURCES_FOLDER}/ab.fasta", score_matrix(DEFAULT_SCORE_MTX))


def test_2():
    return parse_fasta_and_validate(f"{RESOURCES_FOLDER}/sample.fasta", score_matrix(DEFAULT_SCORE_MTX))


def test_3():
    return parse_fasta_and_validate(f"{RESOURCES_FOLDER}/ae.fasta", score_matrix(DEFAULT_SCORE_MTX))


def run_tests():
    passed = execute_tests(
        [
            lambda: test_1(),
            lambda: test_2(),
            lambda: test_3()
        ], "---- [ TESTS Needleman Wunsh ] ----")
    return passed
