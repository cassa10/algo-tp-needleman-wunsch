import datetime

from src import tests, file_parser, grasp, msa
from src.chart import bar_chart

RESOURCES_DIR = "../resources"
OUTPUT_DIR = "../output"


def print_aln_as_output(alignment, fullpath_file):
    with open(fullpath_file, 'w') as f:
        for seq in alignment:
            f.write(f"{seq}\n")


def init_msa_grasp(file_dir_fasta, file_dir_score_matrix, _gap_penalty):
    seqs = file_parser.fasta_multiple_seqs(file_dir_fasta)
    score_mtx = file_parser.score_matrix(file_dir_score_matrix)

    # TODO: Delete dummy score, only for debug
    dummy_score = msa.init(seqs.copy(), score_mtx, gap_penalty)
    print(f"msa dummy score: {dummy_score.score}")

    results = grasp.init(seqs.copy(), score_mtx, _gap_penalty)

    bar_chart("GRASP results", "iterations", "scores",
              map(lambda profile: profile.score, results),
              f"{OUTPUT_DIR}/chart_results-{time}.png",
              save_chart_out_file)
    # TODO: Crear archivo que tenga los alineamientos del resultado final.
    if save_aln_out_file:
        print_aln_as_output(results, f"{OUTPUT_DIR}/alignment-{time}.txt")


def build_dir_file(file_name):
    return f"{RESOURCES_DIR}/{file_name}"


def exec_tests():
    if have_run_tests:
        if tests.run_all():
            print("All tests passed :D\n")
        else:
            print("ERROR - All tests have no passed\n")


if __name__ == '__main__':
    time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    # TODO: Make vars configurable with script params
    have_run_tests = False
    save_chart_out_file = True
    save_aln_out_file = True
    score_matrix_file = "NUC.4.2"
    fasta_file = "10.fasta"
    gap_penalty = -1

    exec_tests()

    # print_aln_as_output(["AGT","-G-","--T","-GT","A-T"], f"{OUTPUT_DIR}/alignment-{time}.txt") bar_chart("GRASP
    # results", "iterations", "scores", [390, 430, 483, 655], f"{OUTPUT_DIR}/chart_results-{time}.png",
    # save_chart_out_file)

    init_msa_grasp(build_dir_file(fasta_file), build_dir_file(score_matrix_file), gap_penalty)
