from src import needleman_wunsh, msa
from src.needleman_wunsh_test import parse_fasta, score_matrix
from src.test_util import execute_tests


def parse_fasta_and_validate(file_dir, score_mtx):
    _, _, seq_A, seq_B = parse_fasta(file_dir)

    gap_penalty = 0
    nw_score, nw_alignments = needleman_wunsh.init(seq_A, seq_B, score_mtx, gap_penalty)
    nw_alignment_seq_A = nw_alignments.pop(0)
    nw_alignment_seq_B = nw_alignments.pop(0)

    msa_score, msa_alignments = msa.init([seq_A, seq_B], score_mtx, gap_penalty)
    msa_alignment_seq_A = msa_alignments.pop(0)
    msa_alignment_seq_B = msa_alignments.pop(0)

    print("Test: ")
    print(f"seq_A: {seq_A}")
    print(f"seq_B: {seq_B}")
    print(f" nw.score => {nw_score}")
    print(f" msa.score => {msa_score}")
    print(f" nw(seq_A) => {nw_alignment_seq_A}")
    print(f" nw(seq_B) => {nw_alignment_seq_B}")
    print(f"msa(seq_A) => {msa_alignment_seq_A}")
    print(f"msa(seq_B) => {msa_alignment_seq_B}")
    print("Asserting if score and alignment are identical in nw and msa with 2 seqs")

    return nw_score == msa_score and nw_alignment_seq_A == msa_alignment_seq_A and nw_alignment_seq_B == msa_alignment_seq_B


def test_1():
    return parse_fasta_and_validate("../resources/ab.fasta", score_matrix("../resources/test_nw_score_matrix"))


def test_2():
    return parse_fasta_and_validate("../resources/sample.fasta", score_matrix("../resources/test_nw_score_matrix"))


def test_3():
    return parse_fasta_and_validate("../resources/ae.fasta", score_matrix("../resources/test_nw_score_matrix"))


def run_tests():
    execute_tests(
        [
            lambda: test_2(),
            lambda: test_1(),

            lambda: test_3()
        ])
