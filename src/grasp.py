from random import randint
from numpy import ceil, log2
from src import msa
from src.possibility_result import PossibilityResult
from src.profile import Profile
from src.solution import Solution


def init(seqs, score_mtx, gap_penalty):
    if len(seqs) <= 1:
        return seqs

    results = []
    n_iteration = get_iterations(seqs)
    select_best_in_greedy = 5
    while n_iteration > 0:
        # greedy solution
        greedy_sol = do_greedy_heuristic(seqs, score_mtx, gap_penalty, select_best_in_greedy)
        print(f"DEBUG - #{n_iteration}-greedy-score: {greedy_sol.score}")
        # local solution
        result = do_local_heuristic(greedy_sol, score_mtx, gap_penalty)
        results.append(result)
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
    return int(ceil(log2(len(seqs)) * 2.33))


def do_local_heuristic(solution, score_mtx, gap_penalty):
    last_score = 0
    current_solution = solution
    while last_score < current_solution.score:
        last_score = current_solution.score
        index_cols, aln = get_neighbor_aln(current_solution)
        current_solution = get_neighbor_solution(current_solution, index_cols, aln, score_mtx, gap_penalty)
    return current_solution


def get_neighbor_aln(solution):
    index_cols = []
    neighbor_aln = solution.get_immutable_aln()

    # TODO: make neighbor solution

    return index_cols, neighbor_aln


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
