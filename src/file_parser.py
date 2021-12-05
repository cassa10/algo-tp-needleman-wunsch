from Bio import SeqIO

from src.score import Score


def fasta_multiple_seqs(file_dir):
    seqs = []
    for fas in SeqIO.parse(file_dir, 'fasta'):
        seqs.append(str(fas.seq))
    return seqs


def fasta_first_two_seq(file_dir):
    f_aln = SeqIO.parse(file_dir, 'fasta')
    fas_A = next(f_aln)
    fas_B = next(f_aln)

    seq_A = str(fas_A.seq)
    seq_B = str(fas_B.seq)
    return fas_A, fas_B, seq_A, seq_B


def score_matrix(file_dir):
    return Score(file_dir).score
