from Bio import SeqIO

from src import needleman_wunsh, msa
from src.needleman_wunsh_test import score_matrix
from src.test_util import execute_tests

DEFAULT_SCORE_MTX_FILE = "NUC.4.2"


def print_test_results(seq_A, seq_B, nw_score, msa_score,
                       nw_alignment_seq_A, nw_alignment_seq_B,
                       msa_alignment_seq_A, msa_alignment_seq_B):
    print("Test: ")
    print(f" nw.score => {nw_score}")
    print(f" msa.score => {msa_score}")
    print(f"       seq_A: {seq_A}")
    print(f"       seq_B: {seq_B}")
    print(f" nw(seq_A) => {nw_alignment_seq_A}")
    print(f" nw(seq_B) => {nw_alignment_seq_B}")
    print(f"msa(seq_A) => {msa_alignment_seq_A}")
    print(f"msa(seq_B) => {msa_alignment_seq_B}")
    print("Asserting if score and alignment are identical in nw and msa with 2 seqs")


def print_test_results_msa(seqs, score, alignment):
    print("Test: ")
    for i, seq in enumerate(seqs):
        print(f"seq_{i}: {seq}")

    print(f"msa.score => {score}")

    for i, aln in enumerate(alignment):
        print(f"msa(seq_{i}) => {aln}")
    print("Asserting if msa do not raise exception")


def parse_fasta_multiple_seqs(file_dir):
    seqs = []
    for fas in SeqIO.parse(file_dir, 'fasta'):
        seqs.append(str(fas.seq))
    return seqs


def parse_fasta_and_validate(file_dir, score_mtx, gap_penalty=0):
    seqs = parse_fasta_multiple_seqs(file_dir)

    seq_A = seqs.pop(0)
    seq_B = seqs.pop(0)
    nw_score, nw_alignments = needleman_wunsh.init(seq_A, seq_B, score_mtx, gap_penalty)
    nw_alignment_seq_A = nw_alignments.pop(0)
    nw_alignment_seq_B = nw_alignments.pop(0)

    msa_score, msa_alignments = msa.init([seq_A, seq_B], score_mtx, gap_penalty)
    msa_alignment_seq_A = msa_alignments.pop(0)
    msa_alignment_seq_B = msa_alignments.pop(0)

    print_test_results(seq_A, seq_B, nw_score, msa_score,
                       nw_alignment_seq_A, nw_alignment_seq_B,
                       msa_alignment_seq_A, msa_alignment_seq_B)

    return nw_score == msa_score and \
           nw_alignment_seq_A == msa_alignment_seq_A and \
           nw_alignment_seq_B == msa_alignment_seq_B


def execute_msa(file_dir, score_mtx, gap_penalty=0):
    seqs = parse_fasta_multiple_seqs(file_dir)

    msa_score, msa_alignment = msa.init(seqs.copy(), score_mtx, gap_penalty)
    print_test_results_msa(seqs, msa_score, msa_alignment)
    return True


def test_1():
    return parse_fasta_and_validate("../resources/ab.fasta",
                                    score_matrix(f"../resources/{DEFAULT_SCORE_MTX_FILE}"))


def test_2():
    return parse_fasta_and_validate("../resources/sample.fasta",
                                    score_matrix(f"../resources/{DEFAULT_SCORE_MTX_FILE}"))


def test_3():
    return parse_fasta_and_validate("../resources/ae.fasta",
                                    score_matrix(f"../resources/{DEFAULT_SCORE_MTX_FILE}"))


def test_4():
    return parse_fasta_and_validate("../resources/ab.fasta",
                                    score_matrix(f"../resources/{DEFAULT_SCORE_MTX_FILE}"),
                                    -1)


def test_5():
    return parse_fasta_and_validate("../resources/sample.fasta",
                                    score_matrix(f"../resources/{DEFAULT_SCORE_MTX_FILE}"),
                                    -1)


def test_6():
    return parse_fasta_and_validate("../resources/ae.fasta",
                                    score_matrix(f"../resources/{DEFAULT_SCORE_MTX_FILE}"),
                                    -1)


def test_7_msa_with_multiple_sequences():
    return execute_msa("../resources/msa_test.fasta",
                       score_matrix(f"../resources/{DEFAULT_SCORE_MTX_FILE}"),
                       -1)


def test_8_msa_with_multiple_sequences():
    return execute_msa("../resources/10.fasta",
                       score_matrix(f"../resources/{DEFAULT_SCORE_MTX_FILE}"),
                       -1)


def run_tests():
    execute_tests(
        [
            lambda: test_1(),
            lambda: test_2(),
            lambda: test_3(),
            lambda: test_4(),
            lambda: test_5(),
            lambda: test_6(),
            lambda: test_7_msa_with_multiple_sequences(),
            lambda: test_8_msa_with_multiple_sequences()
        ])
