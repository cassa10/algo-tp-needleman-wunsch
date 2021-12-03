from Bio import pairwise2, SeqIO
import needleman_wunsh
from src.score import Score
from src.test_util import execute_tests


def parse_fasta(file_dir):
    f_aln = SeqIO.parse(file_dir, 'fasta')
    fas_A = next(f_aln)
    fas_B = next(f_aln)

    seq_A = str(fas_A.seq)
    seq_B = str(fas_B.seq)
    return fas_A, fas_B, seq_A, seq_B


def parse_fasta_and_validate(file_dir, score_mtx):
    fas_A, fas_B, seq_A, seq_B = parse_fasta(file_dir)

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


def score_matrix(file_dir):
    return Score(file_dir).score


def test_1():
    return parse_fasta_and_validate("../resources/ab.fasta", score_matrix("../resources/test_nw_score_matrix"))


def test_2():
    return parse_fasta_and_validate("../resources/sample.fasta", score_matrix("../resources/test_nw_score_matrix"))


def test_3():
    return parse_fasta_and_validate("../resources/ae.fasta", score_matrix("../resources/test_nw_score_matrix"))


def run_tests():
    execute_tests(
        [
            lambda: test_1(),
            lambda: test_2(),
            lambda: test_3()
        ])
