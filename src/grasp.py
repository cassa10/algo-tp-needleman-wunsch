from random import randint
from numpy import ceil, log2
from src import msa
from src.possibility_result import PossibilityResult
from src.profile import Profile


def init(seqs, score_mtx, gap_penalty):
    if len(seqs) <= 1:
        return seqs

    results = []
    n_iteration = get_iterations(seqs)
    select_best_in_greedy = 5
    while n_iteration > 0:
        # greedy solution
        greedy_sol = do_greedy_heuristic(seqs, score_mtx, gap_penalty, select_best_in_greedy)
        print(f"{n_iteration}-greedy-score: {greedy_sol.score}")
        #result = do_local_heuristic(greedy_sol, score_mtx, gap_penalty)
        #results.append(result)
        results.append(greedy_sol)
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

    return profile_result


def select_possibility(seqs, possibilities, select_best):

    possibilities.sort(key=lambda possibility_res: possibility_res.score(), reverse=True)
    select_top = select_best
    if len(possibilities) < select_best:
        select_top = len(possibilities)

    top_poss = possibilities[:select_top]
    # Random
    selector = randint(0, select_top-1)
    selected_poss_res = top_poss[selector]

    seqs = [seq for i, seq in enumerate(seqs) if selected_poss_res.not_contain_index(i)]

    return seqs, selected_poss_res.profile_result


def get_iterations(seqs):
    return int(ceil(log2(len(seqs)) * 2.33))


def do_local_heuristic(solution, score_mtx, gap_penalty):
    # while siga obteniendo mejor scoring sigo y sino corto y retorno lo mejor que encontre
    return None


