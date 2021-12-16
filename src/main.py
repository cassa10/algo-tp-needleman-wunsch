import datetime
import time
import file_parser
import grasp
import tests
from output import bar_chart, print_aln

RESOURCES_DIR = "./resources"
OUTPUT_DIR = "./output"


def init_msa_grasp(file_dir_fasta, file_dir_score_matrix, _gap_penalty):
    seqs = file_parser.fasta_multiple_seqs(file_dir_fasta)
    score_mtx = file_parser.score_matrix(file_dir_score_matrix)

    results = grasp.init(seqs.copy(), score_mtx, _gap_penalty)
    make_output_files(results)


def make_output_files(results):
    bar_chart("GRASP results", "# de iterations of local search", "score",
              results,
              build_output_dir_file(chart_out_file),
              save_chart_out_file)
    if save_aln_out_file:
        print_aln(get_best_solution(results), build_output_dir_file(aln_out_file))


def get_best_solution(results):
    best_sol = results[0][-1]
    for result in results:
        current_sol = result[-1]
        if best_sol.score < current_sol.score:
            best_sol = current_sol
    return best_sol


def build_resource_dir_file(file_name):
    return f"{RESOURCES_DIR}/{file_name}"


def build_output_dir_file(file_name):
    return f"{OUTPUT_DIR}/{file_name}"


def exec_tests():
    if have_run_tests:
        if tests.run_all():
            print("All tests passed :D\n")
        else:
            print("ERROR - All tests have no passed\n")


if __name__ == '__main__':
    cur_datetime = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    aln_out_file = f"alignment-{cur_datetime}.txt"
    chart_out_file = f"chart_results-{cur_datetime}.png"
    # TODO: Make vars configurable with script params
    save_chart_out_file = True
    save_aln_out_file = True
    score_matrix_file = "NUC.4.2"
    fasta_file = "10.fasta"
    gap_penalty = -1

    have_run_tests = False
    exec_tests()

    start_time = time.time()
    init_msa_grasp(build_resource_dir_file(fasta_file), build_resource_dir_file(score_matrix_file), gap_penalty)
    print(f"Finish MSA GRASP in {time.time() - start_time} seconds")
