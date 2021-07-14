import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.widgets import Slider


class WindowPlot:

    def __init__(self, sliding_window, colors=None):
        self.sw = sliding_window
        self.colors = colors
        #### helper variables ####
        self.plots_path = 'results/plots/'
        self.filtered_range = self.calc_window_range(self.sw.filtered_data)
        self.unfiltered_range = self.calc_window_range(self.sw.window_data)
        # add x-axis tick positions to data
        self.add_x_ticks(self.sw.filtered_data, self.filtered_range)
        self.add_x_ticks(self.sw.window_data, self.unfiltered_range)

    def num_windows(self, df):
        return df.shape[0] / self.sw.n_subset

    def calc_window_range(self, df):
        return np.arange(self.num_windows(df))

    def add_x_ticks(self, df, window_range):
        df['x_ticks'] = np.repeat(window_range, self.sw.n_subset)

    def get_min(self, df):
        return df.groupby('window').min()['window_mean']

    def get_max(self, df):
        return df.groupby('window').max()['window_mean']

    def save_plot(self, file_name):
        plt.savefig(f'{self.plots_path}{file_name}.png')
        plt.savefig(f'{self.plots_path}{file_name}.pdf')

    def ecdf_plot(self):
        plot = sns.displot(
            data=self.sw.window_data,
            x="window_mean",
            hue="subset",
            kind="ecdf",
            palette=self.colors)
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
            palette=self.colors).set(xlim=(0, max_freq))
        plot._legend.set_title('Subset')
        plt.xlabel('Mean window frequency')
        plt.tight_layout()
        self.save_plot('kde_plot')
        plt.show()

    def alternating_vlines(self, window_range, max_values):
        # plot lines (with alternating linestyle) extending down to x axis from maximum frequency
        plt.vlines(x=window_range[0::2], ymin=-0.1,
                   ymax=max_values[0::2], color='grey', zorder=0, lw=0.6)
        plt.vlines(x=window_range[1::2], ymin=-0.1,
                   ymax=max_values[1::2], color='gray', zorder=0, lw=0.6, linestyle='dotted')

    def dot_plot(self, df, window_range, axes=None, save_plot=True, show=True):
        plt.rcParams['figure.figsize'] = [15, 6]
        label = "Window"
        if not axes:
            label = f"Filtered window"
        axes = sns.scatterplot(
            data=df,
            x='x_ticks',
            y='window_mean',
            hue='subset',
            palette=self.colors,
            ax=axes)
        axes.legend(title='Subset')
        min_values = self.get_min(df)
        max_values = self.get_max(df)
        _, y_max = axes.get_ylim()
        plt.ylim((min_values.min(), y_max))
        self.alternating_vlines(window_range, max_values)
        axes.tick_params(
            axis='x',
            which='major',
            labelsize=6,
            labelrotation=-90)
        axes.set_xticks(window_range)
        axes.set_xticklabels(df['window'][::self.sw.n_subset])
        title = f"Target{'s' if self.sw.target.size > 1 else ''}: {', '.join(self.sw.target)}"
        plt.title(title)
        x_label = f'{label} position (window size = {self.sw.window_size})'
        plt.xlabel(x_label, labelpad=10)
        plt.ylabel('Mean window frequency')
        plt.tight_layout()
        if save_plot:
            self.save_plot('dot_plot')
        if show:
            plt.show()

    def line_plot(self, axes):
        sns.lineplot(
            data=self.sw.window_data,
            x='x_ticks',
            y='window_mean',
            hue='subset',
            sort=False,
            palette=self.colors,
            ax=axes)
        axes.legend(title='Subset')
        axes.set(xlabel=None)
        axes.set(ylabel=None)
        plt.tight_layout()

    def sliding_plot(self, plot_scroll_size):
        fig, axes = plt.subplots(2)
        self.dot_plot(self.sw.window_data, self.unfiltered_range,
                      axes=axes[1], save_plot=False, show=False)
        self.line_plot(axes[0])
        axes[1].set(title=None)
        axes[1].set(ylabel=None)
        axes[0].set(
            title=f"Target{'s' if self.sw.target.size > 1 else ''}: {', '.join(self.sw.target)}")
        plt.subplots_adjust(left=0.15, right=0.9, bottom=0.2, top=0.90)
        # add common y label
        fig.text(0.05, 0.5, 'Mean window frequency', va='center',
                 ha='center', rotation='vertical')
        # set the axis and slider position in the plot
        axis_pos = plt.axes([0.15, 0.01, 0.75, 0.03], facecolor='White')
        slider_pos = Slider(axis_pos, 'Window', 0.0,
                            self.num_windows(self.sw.window_data))
        # update() function to change the graph when the slider is in use
        y_min, y_max = axes[1].get_ylim()

        def update(val):
            pos = slider_pos.val
            axes[1].axis([pos - 1, pos + plot_scroll_size, y_min, y_max])
            fig.canvas.draw_idle()
        update(slider_pos.val)
        # update function called using on_changed() function
        slider_pos.on_changed(update)

        plt.show()

    def make_plots(self):
        self.sliding_plot(plot_scroll_size=50)
        self.dot_plot(self.sw.filtered_data, self.filtered_range)
        self.kde_plot()
        self.ecdf_plot()
