import matplotlib.pyplot as plt


def print_aln(solution, fullpath_file):
    with open(fullpath_file, 'w') as f:
        f.write(f"SCORE = {solution.score}\n\n")
        for seq in solution.alignment:
            f.write(f"{seq}\n")


# results :: [[Solution]]
def bar_chart(title, x_label, y_label, results, fullpath_file, save_file=False):
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    for i, result in enumerate(results):
        scores = [s.score for s in result]
        plt.plot(list(range(0, len(result))), scores, label=f"GRASP {i}")

    plt.xlim(1, max([len(result) for result in results]) + 10)
    plt.ylim(0, max([s.score for result in results for s in result]) + 50)
    plt.legend()
    if save_file:
        plt.savefig(fullpath_file)
    plt.show()
