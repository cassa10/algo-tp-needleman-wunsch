SCORE_MTX_FILE = "../resources/NUC.4.2"


class Score:
    score = {
        "AA": 1, "CC": 1, "GG": 1, "TT": 1,
        "AC": 0, "AG": 0, "AT": 0,
        "CA": 0, "CG": 0, "CT": 0,
        "GA": 0, "GC": 0, "GT": 0,
        "TA": 0, "TC": 0, "TG": 0,
        # TODO: Ver que score poner en este caso
        "--": 1
    }

    def __init__(self, score_mtx_file_dir):
        score_col_order = []
        col_readed = False
        for line in open(score_mtx_file_dir, "r"):
            li = line.strip()
            if not li.startswith("#") and not li.__eq__(""):
                line.rstrip()
                if not col_readed:
                    score_col_order = list(filter(lambda x: x != "", li.split(" ")))
                    col_readed = True
                else:
                    file_scores = li.split(" ")
                    current_l = file_scores.pop(0)
                    file_scores = enumerate(list(filter(lambda x: x != "", file_scores)))
                    for i, scr_str in file_scores:
                        self.score[current_l + score_col_order[i]] = int(scr_str)
