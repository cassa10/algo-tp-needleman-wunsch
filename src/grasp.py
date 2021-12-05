from random import randint
from numpy import ceil
from src import msa
from src.profile import Profile


def init(seqs, score_mtx, gap_penalty):
    if len(seqs) <= 1:
        return seqs
    n_iteration = get_iterations(seqs)
    results = []
    while n_iteration > 0:
        # greedy solution
        greedy_sol = do_greedy_heuristic(seqs, score_mtx, gap_penalty)

        result = do_local_heuristic(greedy_sol, score_mtx, gap_penalty)
        results.append(result)
        n_iteration -= 1

    return results


def do_greedy_heuristic(seqs, score_mtx, gap_penalty):
    select_best = 3
    possibilities = []

    for i1, seq_a in enumerate(seqs):
        seqs.pop(i1)
        for i2, seq_b in enumerate(seqs):
            possibilities.append((i1, i2, msa.nw_with_profile(Profile(seq_a), seq_b, score_mtx, gap_penalty)))
        seqs.insert(i1, seq_a)

    possibilities.sort(key=lambda i_seq_fst, i_seq_snd, prof_res: prof_res.score)
    select_top = select_best
    if len(possibilities) < select_best:
        select_top = len(possibilities)

    top_poss = possibilities[:select_top]
    # Random greedy
    selector = randint(0, select_top - 1)
    index_seq_fst, index_seq_snd, profile_result = top_poss[selector]

    seqs = [seq for i, seq in enumerate(seqs) if i != index_seq_fst and i != index_seq_snd]

    if len(seqs) == 0:
        return profile_result
    if len(seqs) == 1:
        return msa.nw_with_profile(profile_result, seqs.pop(0), score_mtx, gap_penalty)

    # TODO: VOLVER A ITERAR

    result = None
    return result


def get_iterations(seqs):
    return ceil(len(seqs) * 2 / 1.6)


def do_local_heuristic(solution, score_mtx, gap_penalty):
    # while siga obteniendo mejor scoring sigo y sino corto y retorno lo mejor que encontre
    return None


