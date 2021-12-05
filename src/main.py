from src import tests, msa, file_parser

RESOURCES_DIR = "../resources/"


def init_msa_grasp(file_dir_fasta, file_dir_score_matrix, gap_penalty):
    seqs = file_parser.fasta_multiple_seqs(file_dir_fasta)
    score_mtx = file_parser.score_matrix(file_dir_score_matrix)
    return msa.init(seqs, score_mtx, gap_penalty)


def build_dir_file(file_name):
    return f"{RESOURCES_DIR}{file_name}"


def exec_tests(have_to_exec):
    if have_to_exec:
        if tests.run_all():
            print("All tests passed :D\n")
        else:
            print("ERROR - All tests have no passed\n")


if __name__ == '__main__':
    have_run_tests = False
    exec_tests(have_run_tests)

    score_matrix_file = build_dir_file("NUC.4.2")
    fasta_file = build_dir_file("10.fasta")
    gap_penalty = -1
    score, aln = init_msa_grasp(fasta_file, score_matrix_file, gap_penalty)
    print(score)
    for e in aln:
        print(e)
