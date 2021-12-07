import matplotlib.pyplot as plt


def print_aln(alignment, fullpath_file):
    with open(fullpath_file, 'w') as f:
        for seq in alignment:
            f.write(f"{seq}\n")


# data : [(n_iteration : Int, score : Int )]
def bar_chart(title, x_label, y_label, scores, fullpath_file, save_file=False):
    bar_width = 0.35
    fig, ax = plt.subplots()

    ind = range(1, len(scores) + 1)
    rects = ax.bar(ind, scores, bar_width, label="GRASP")
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xticks(ind, labels=range(1, len(scores) + 1))
    ax.legend()

    ax.bar_label(rects)
    fig.tight_layout()
    if save_file:
        plt.savefig(fullpath_file)

    plt.show()
