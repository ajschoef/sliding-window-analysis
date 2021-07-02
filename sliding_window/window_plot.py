import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


class WindowPlot:

    def __init__(self, sliding_window):
        self.sliding_window = sliding_window
        self.filtered_data = sliding_window.filtered_data
        self.n_subset = sliding_window.n_subset
        self.window_range = sliding_window.window_range
        self.colors = sliding_window.colors
        self.target = sliding_window.target
        self.window_size = sliding_window.window_size

    def plot(self):
        plt.rcParams['figure.figsize'] = [15, 6]
        sns.scatterplot(data=self.filtered_data, x='x_ticks', y='frequency',
                        hue='subset', palette=self.colors)
        # get min and max frequency for each window
        min_freqs = self.filtered_data.groupby('position').min()[
            'frequency']
        max_freqs = self.filtered_data.groupby('position').max()[
            'frequency']
        axes = plt.gca()
        _, y_max = axes.get_ylim()
        plt.ylim((min_freqs.min(), y_max))

        # plot lines (with alternating linestyle) extending down to x axis from maximum frequency
        plt.vlines(x=self.window_range[0::2], ymin=-0.1,
                   ymax=max_freqs[0::2], color='grey', zorder=0, lw=0.6)
        plt.vlines(x=self.window_range[1::2], ymin=-0.1,
                   ymax=max_freqs[1::2], color='gray', zorder=0, lw=0.6, linestyle='dotted')
        axes.tick_params(axis='x', which='major',
                         labelsize=6, labelrotation=-90)
        axes.set_xticks(self.window_range)
        axes.set_xticklabels(self.filtered_data['position']
                             [::self.n_subset])
        plt.title(
            f"Target{'s' if self.target.size > 1 else ''}: {', '.join(self.target)}")
        plt.xlabel(
            f"Window position (window size = {self.window_size})")
        plt.ylabel('Mean window frequency')
        plt.tight_layout()
        plt.show()
        plt.savefig('results/window_comparisons.png')
        plt.savefig('results/window_comparisons.pdf')
