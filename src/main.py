from src import tests, file_parser, grasp, msa

RESOURCES_DIR = "../resources/"


def init_msa_grasp(file_dir_fasta, file_dir_score_matrix, _gap_penalty):
    seqs = file_parser.fasta_multiple_seqs(file_dir_fasta)
    score_mtx = file_parser.score_matrix(file_dir_score_matrix)

    # TODO: Delete dummy score, only for debug
    dummy_score = msa.init(seqs.copy(), score_mtx, gap_penalty)
    print(f"msa dummy score: {dummy_score.score}")

    results = grasp.init(seqs.copy(), score_mtx, _gap_penalty)
    # TODO: Imprimir grafico de results
    # TODO: Crear archivo que tenga los alineamientos del resultado final.


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

    # TODO: Make it configurable with script params
    score_matrix_file = "NUC.4.2"
    fasta_file = "10.fasta"
    gap_penalty = -1
    init_msa_grasp(build_dir_file(fasta_file), build_dir_file(score_matrix_file), gap_penalty)
