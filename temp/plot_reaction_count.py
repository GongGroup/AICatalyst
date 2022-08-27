from functools import wraps

from matplotlib import pyplot as plt


def plot_wrapper(func):
    @wraps(func)
    def wrapper(*args, **kargs):
        bwidth = 2
        fontsize = 22
        plt.rc('font', family="Arial", weight="regular")  # set the font globally
        plt.rcParams['mathtext.default'] = 'regular'  # set the math-font globally
        plt.rcParams['lines.linewidth'] = 2  # set line-width
        plt.rcParams['lines.markersize'] = 9.0
        ax = plt.gca()
        ax.spines['bottom'].set_linewidth(bwidth)  # set border
        ax.spines['left'].set_linewidth(bwidth)
        ax.spines['top'].set_linewidth(bwidth)
        ax.spines['right'].set_linewidth(bwidth)
        # plt.xlim() if self.xlim is None else plt.xlim(self.xlim)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        # plt.xlabel(fontsize=fontsize + 1)
        # plt.ylabel(fontsize=fontsize + 1)
        # plt.title(self.title, fontsize=self.fontsize + 2)
        func(*args, **kargs)

    return wrapper


@plot_wrapper
def main():
    count = [14, 8, 10, 20, 75, 326]
    label = ['>5', '5', '4', '3', '2', '1']
    plt.bar(label, count)
    for x, y in zip(label, count):
        plt.text(x, y, y, ha='center', va='bottom', fontdict={"size": 20})

    plt.ylim([0, 360])
    plt.xlabel("Sample Range", fontsize=24)
    plt.ylabel("Sample Count", fontsize=24)
    plt.show()


if __name__ == '__main__':
    main()
