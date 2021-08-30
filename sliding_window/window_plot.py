import numpy as np
from scipy.stats.kde import gaussian_kde
from bokeh.io import show
from bokeh.layouts import row, gridplot
from bokeh.models import RangeTool, FuncTickFormatter, Div
from bokeh.plotting import figure, output_file
from bokeh.palettes import Category20


class WindowPlot:

    def __init__(self, sliding_window, colors=None):
        self.sw = sliding_window
        self.colors = colors
        # helper variables
        self.plots_path = 'results/plots/'
        self.groupby_subset = self.sw.window_data.groupby('subset')
        self.groupby_subset_filtered = self.sw.window_data_filtered.groupby(
            'subset')
        self.parameter_text = self.make_parameter_text()
        # append x-ticks for dot plots
        self.add_x_ticks(self.sw.window_data)
        self.add_x_ticks(self.sw.window_data_filtered)

    def make_parameter_text(self):
        target_text = ', '.join(self.sw.target)
        window_size_text = f"Window size: {self.sw.window_size}"
        is_plural = 's' if self.sw.target.size > 1 else ''
        stride_text = f"Window stride: {self.sw.stride}"
        return f"Target{is_plural}: {target_text}</br>{window_size_text}</br>{stride_text}"

    def num_windows(self, df):
        return df.shape[0] / self.sw.n_subset

    def calc_window_range(self, df):
        return np.arange(1, self.num_windows(df)+1)

    def add_x_ticks(self, df):
        df['x_ticks'] = np.repeat(self.calc_window_range(
            df), self.sw.n_subset)

    def ecdf_plot(self):
        p = figure(x_axis_label='Window frequency mean',
                   y_axis_label='Fraction of data')
        for name, group in self.groupby_subset:
            x = np.sort(group['window_mean'])
            y = np.arange(len(x))/float(len(x))
            color = group['color'].values[0]
            p.step(x, y, line_width=2, mode="center",
                   color=color, muted_color=color,
                   muted_alpha=0.2, legend_label=name)
        p.legend.click_policy = "hide"
        p.legend.location = "bottom_right"
        return p

    def kde_plot(self):
        maxs = self.groupby_subset['window_mean'].max()
        mins = self.groupby_subset['window_mean'].min()
        p = figure(x_axis_label='Window frequency mean',
                   y_axis_label='Density', x_range=(np.min(mins), np.max(maxs)))
        for name, group in self.groupby_subset:
            frequencies = group['window_mean']
            x = np.linspace(mins.loc[name], maxs.loc[name], 1000)
            pdf = gaussian_kde(frequencies)
            y = pdf(x)
            color = group['color'].values[0]
            p.line(x, y, legend_label=name, line_width=2,
                   line_color=color, muted_color=color,
                   muted_alpha=0.2)
            p.legend.click_policy = "hide"
        return p

    def distribution_plots(self):
        output_file("results/plots/distributions.html",
                    title="Sliding window plot")
        p = self.ecdf_plot()
        q = self.kde_plot()
        show(row(p, q))

    def make_color_palette(self):
        colors = Category20[20]
        color_map = {name: color for name,
                     color in zip(self.sw.subset_names, colors)}
        return self.sw.window_data['subset'].map(color_map)

    def dot_plot(self, groupby_subset, tools, x_range):
        tooltips = [
            ("Subset", "@subset"),
            ("Window position", "@window_start"),
            ("Mean frequency", "@window_mean"),
        ]
        p = figure(
            plot_height=400,
            plot_width=1000,
            tools=tools,
            toolbar_location="right",
            x_axis_location="below",
            background_fill_color="#efefef",
            x_range=x_range,
            tooltips=tooltips)

        for name, group in groupby_subset:
            p.segment('x_ticks', 0, 'x_ticks', 'window_mean', line_color='#bdbdbd',
                      line_width=1, legend_label=name, source=group)
        circles = []
        for name, group in groupby_subset:
            circles.append(p.circle('x_ticks', 'window_mean',
                                    size=7,
                                    legend_label=name,
                                    color='color',
                                    muted_color='color',
                                    muted_alpha=0.2,
                                    source=group))
        p.legend.click_policy = "hide"
        p.hover.renderers = circles
        p.yaxis.axis_label = 'Window frequency mean'
        return p

    def filtered_dot_plot(self, df):
        output_file("results/plots/filtered_dot_plot.html",
                    title="Sliding window plot")
        if self.sw.n_largest > 50:
            x_range = (0, 50)
        else:
            x_range = (0, self.sw.n_largest+1)
        p = self.dot_plot(self.groupby_subset_filtered,
                          ['xpan', 'save'], x_range)
        p.xaxis.axis_label = 'Window position'
        p.xaxis.ticker = df['x_ticks']
        self.format_x_axis(p, df, 'x_ticks', 'window_start')
        p.xaxis.major_label_orientation = "vertical"
        div = Div(text=self.parameter_text, width=200, height=100)
        show(row(p, div))

    def format_x_axis(self, p, df, ticks, tick_labels):
        x_label = {key: value for key, value in zip(
            df[ticks], df[tick_labels])}
        p.xaxis.formatter = FuncTickFormatter(code="""
            var labels = %s;
            return labels[tick] || tick;
            """ % x_label)

    def sliding_panel(self, p, df, groupby_subset):
        self.format_x_axis(p, df, 'x_ticks', 'window_start')
        select = figure(
            plot_height=150,
            plot_width=1000,
            y_range=p.y_range,
            y_axis_type=None,
            tools="",
            toolbar_location=None,
            background_fill_color="#efefef"
        )
        select.xaxis.axis_label = 'Window position'
        self.format_x_axis(select, df, 'x_ticks', 'window_start')
        range_tool = RangeTool(x_range=p.x_range)
        range_tool.overlay.fill_color = "navy"
        range_tool.overlay.fill_alpha = 0.2
        for name, group in groupby_subset:
            select.circle('x_ticks', 'window_mean',
                          color=group['color'].values[0], alpha=0.6, source=group)
        select.ygrid.grid_line_color = None
        select.add_tools(range_tool)
        select.toolbar.active_multi = range_tool
        return select

    def sliding_window_plot(self, df):
        output_file("results/plots/dot_plot.html",
                    title="Sliding window plot")
        initial_panel_size = int(
            self.sw.summary_stats['alignment_length'][0] * 0.1)
        p = self.dot_plot(self.groupby_subset, [
                          "xpan", "save"], (0, initial_panel_size))
        select = self.sliding_panel(p, df, self.groupby_subset)
        div = Div(text=self.parameter_text, width=200, height=100)
        grid = gridplot([[p, div], [select, None]])
        show(grid)

    def make_plots(self):
        self.sw.window_data['color'] = self.make_color_palette()
        self.sw.window_data_filtered['color'] = self.make_color_palette()
        self.sliding_window_plot(self.sw.window_data)
        self.filtered_dot_plot(self.sw.window_data_filtered)
        self.distribution_plots()
