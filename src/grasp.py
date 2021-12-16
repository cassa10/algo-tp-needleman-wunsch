from random import randint
from numpy import ceil, log2
import msa
from possibility_result import PossibilityResult
from profile import Profile
from solution import Solution


def init(seqs, score_mtx, gap_penalty):
    if len(seqs) <= 1:
        return [Solution(0, seqs)]

    n_iteration = get_iterations(seqs)
    select_best_in_greedy = 5
    results = []
    iteration_result = []
    while n_iteration > 0:
        # greedy solution
        greedy_sol = do_greedy_heuristic(seqs, score_mtx, gap_penalty, select_best_in_greedy)
        print(f"DEBUG - #{n_iteration}-greedy-score: {greedy_sol.score}")
        iteration_result.append(greedy_sol)
        # local solution
        termination_criteria = 100
        local_sol = do_local_heuristic(iteration_result, greedy_sol, score_mtx, gap_penalty, termination_criteria)
        print(f"DEBUG - #{n_iteration}-best-local-score: {local_sol.score}")
        results.append(iteration_result)
        iteration_result = []
        n_iteration -= 1

    return results


def do_greedy_heuristic(seqs, score_mtx, gap_penalty, select_best):
    possibilities = []
    for i1, seq_a in enumerate(seqs):
        for i2, seq_b in enumerate(seqs):
            if i1 != i2:
                possibilities.append(
                    PossibilityResult(msa.nw_with_profile(Profile(seq_a), seq_b, score_mtx, gap_penalty), i2, i1)
                )

    seqs, profile_result = select_possibility(seqs, possibilities, select_best)

    while len(seqs) > 0:
        possibilities = []
        for i, seq in enumerate(seqs):
            possibilities.append(
                PossibilityResult(msa.nw_with_profile(profile_result.copy(), seq, score_mtx, gap_penalty), i)
            )
        seqs, profile_result = select_possibility(seqs, possibilities, select_best)

    return profile_result.map_to_solution()


def select_possibility(seqs, possibilities, select_best):
    possibilities.sort(key=lambda possibility_res: possibility_res.score(), reverse=True)
    select_top = select_best
    if len(possibilities) < select_best:
        select_top = len(possibilities)

    top_poss = possibilities[:select_top]
    # Random
    selector = randint(0, select_top - 1)
    selected_poss_res = top_poss[selector]

    seqs = [seq for i, seq in enumerate(seqs) if selected_poss_res.not_contain_index(i)]

    return seqs, selected_poss_res.profile_result


def get_iterations(seqs):
    return int(ceil(log2(len(seqs)) * 1.33))


def do_local_heuristic(iteration_res, solution, score_mtx, gap_penalty, termination_criteria):
    best_solution = solution
    n_iterations = termination_criteria
    while 0 < n_iterations:
        index_cols, aln = get_neighbor_aln(best_solution)
        current_solution = get_neighbor_solution(best_solution, index_cols, aln, score_mtx, gap_penalty)

        if best_solution.score < current_solution.score:
            best_solution = current_solution

        iteration_res.append(best_solution)
        n_iterations -= 1
    return best_solution


# Neighbor solution - Switch nucleotide with gap, eg: (-X) => (X-) where X in [A,C,G,T]
def get_neighbor_aln(solution):
    neighbor_aln = solution.get_immutable_aln()
    # Use random because always return same solution in case of return a bad solution
    seq_selector = randint(0, len(neighbor_aln) - 1)

    seq_selected = neighbor_aln[seq_selector]

    prev_nuc = False
    prev_nuc_col = 0
    gap_col = 0
    nuc = ""

    select_start_col = randint(0, int(len(seq_selected) * 0.75))
    for i, c in enumerate(seq_selected[select_start_col:]):
        if c != "-":
            prev_nuc = True
            prev_nuc_col = select_start_col + i

        if prev_nuc and c == "-":
            gap_col = select_start_col + i
            nuc = c
            break

    # Not found swap
    if prev_nuc_col == 0 and gap_col == 0:
        return [], neighbor_aln

    seq = neighbor_aln[seq_selector]
    new_seq = seq[:prev_nuc_col] + "-" + nuc + seq[gap_col+1:]
    neighbor_aln[seq_selector] = new_seq
    return [prev_nuc_col, gap_col], neighbor_aln


def get_neighbor_solution(solution, index_cols, aln, score_mtx, gap_penalty):
    new_score = solution.score
    for index_col in index_cols:
        new_score -= calculate_col_score(solution.alignment, index_col, score_mtx, gap_penalty)
        new_score += calculate_col_score(aln, index_col, score_mtx, gap_penalty)

    return Solution(new_score, aln)


def calculate_col_score(aln, index_col, score_mtx, gap_penalty):
    acc = 0
    for i1, s1 in enumerate(aln):
        for i2, s2 in enumerate(aln):
            if i1 != i2:
                acc += msa.get_score(score_mtx, gap_penalty, s1[index_col], s2[index_col])
    return acc
