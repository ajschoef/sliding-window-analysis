import matplotlib.pyplot as plt
import seaborn as sns


class WindowPlot:

    def __init__(self, sliding_window):
        self.sliding_window = sliding_window
        self.data = sliding_window.freq_means_long
        self.filtered_data = sliding_window.filtered_data
        self.n_subset = sliding_window.n_subset
        self.window_range = sliding_window.window_range
        self.colors = sliding_window.colors
        self.target = sliding_window.target
        self.window_size = sliding_window.window_size
        #### helper variables ####
        # get min and max frequency for each window
        self.min_freqs = self.filtered_data.groupby('position').min()[
            'frequency']
        self.max_freqs = self.filtered_data.groupby('position').max()[
            'frequency']

    def save_plot(self, file_name):
        plt.savefig(f'results/{file_name}.png')
        plt.savefig(f'results/{file_name}.pdf')

    def ecdf_plot(self):
        plot = sns.displot(data=self.filtered_data, x="frequency",
                           hue="subset", kind="ecdf", palette=self.colors)
        plot._legend.set_title('Subset')
        plt.xlabel('Frequency')
        plt.tight_layout()
        self.save_plot('ecdf_plot')
        plt.show()

    def kde_plot(self):
        plot = sns.displot(data=self.data,
                           x="frequency", hue="subset", kind="kde", palette=self.colors)
        plot._legend.set_title('Subset')
        plt.xlabel('Frequency')
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
        plot = sns.scatterplot(data=self.filtered_data, x='x_ticks', y='frequency',
                               hue='subset', palette=self.colors)
        plot.legend(title='Subset')
        axes = plt.gca()
        _, y_max = axes.get_ylim()
        plt.ylim((self.min_freqs.min(), y_max))
        self.alternate_vlines()
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
        self.save_plot('dot_plot')
        plt.show()
