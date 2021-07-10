import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


class WindowPlot:

    def __init__(self, sliding_window):
        self.sw = sliding_window
        #### helper variables ####
        self.plots_path = 'results/plots/'
        # get min and max frequency for each window
        self.min_freqs = self.sw.filtered_data.groupby('position').min()[
            'window_mean']
        self.max_freqs = self.sw.filtered_data.groupby('position').max()[
            'window_mean']
        self.window_range = self.calc_window_range()
        # add x-axis tick positions to filtered data
        self.add_x_ticks()

    def calc_window_range(self):
        return np.arange(self.sw.filtered_data['position'].nunique())

    def add_x_ticks(self):
        self.sw.filtered_data['x_ticks'] = np.repeat(
            self.window_range, self.sw.n_subset)

    def save_plot(self, file_name):
        plt.savefig(f'{self.plots_path}{file_name}.png')
        plt.savefig(f'{self.plots_path}{file_name}.pdf')

    def ecdf_plot(self):
        plot = sns.displot(
            data=self.sw.window_data,
            x="window_mean",
            hue="subset",
            kind="ecdf",
            palette=self.sw.colors)
        plot._legend.set_title('Subset')
        plt.xlabel('Mean window frequency')
        plt.tight_layout()
        self.save_plot('ecdf_plot')
        plt.show()

    def kde_plot(self):
        max_freq = self.sw.window_data['window_mean'].max()
        plot = sns.displot(
            data=self.sw.window_data,
            x="window_mean",
            hue="subset",
            kind="kde",
            common_norm=False,
            palette=self.sw.colors).set(xlim=(0, max_freq))
        plot._legend.set_title('Subset')
        plt.xlabel('Mean window frequency')
        plt.tight_layout()
        self.save_plot('kde_plot')
        plt.show()

    def alternate_vlines(self):
        # plot lines (with alternating linestyle) extending down to x axis from maximum frequency
        plt.vlines(x=self.window_range[0::2], ymin=-0.1,
                   ymax=self.max_freqs[0::2], color='grey', zorder=0, lw=0.6)
        plt.vlines(x=self.window_range[1::2], ymin=-0.1,
                   ymax=self.max_freqs[1::2], color='gray', zorder=0, lw=0.6, linestyle='dotted')

    def dot_plot(self):
        plt.rcParams['figure.figsize'] = [15, 6]
        plot = sns.scatterplot(
            data=self.sw.filtered_data,
            x='x_ticks',
            y='window_mean',
            hue='subset',
            palette=self.sw.colors)
        plot.legend(title='Subset')
        axes = plt.gca()
        _, y_max = axes.get_ylim()
        plt.ylim((self.min_freqs.min(), y_max))
        self.alternate_vlines()
        axes.tick_params(axis='x', which='major',
                         labelsize=6, labelrotation=-90)
        axes.set_xticks(self.window_range)
        axes.set_xticklabels(self.sw.filtered_data['position']
                             [::self.sw.n_subset])
        plt.title(
            f"Target{'s' if self.sw.target.size > 1 else ''}: {', '.join(self.sw.target)}")
        plt.xlabel(
            f"Window position (window size = {self.sw.window_size})")
        plt.ylabel('Mean window frequency')
        plt.tight_layout()
        self.save_plot('dot_plot')
        plt.show()

    def make_plots(self):
        self.dot_plot()
        self.kde_plot()
        self.ecdf_plot()
