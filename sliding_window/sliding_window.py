import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pickle
from sklearn.linear_model import MultiTaskElasticNetCV
from pathlib import Path


class SlidingWindow:

    def __init__(self, path_to_data, target, window_size, cutoff, stride=None, colors=None, fit=None):
        # parameters
        self.path_to_data = Path(path_to_data)
        self.target = np.asarray(list(target))
        self.window_size = window_size
        self.cutoff = cutoff
        self.stride = stride
        self.colors = colors
        self.fit = fit
        # data
        self.frequencies = None
        self.frequency_means = None
        self.summary = []
        self.filtered_indices = None
        self.fasta_extensions = (
            '.fasta', '.fa', '.fna', '.ffn', '.faa', '.frn')
        with self.path_to_data as entries:
            self.subset_names = [f.name.split('.')[0]
                                 for f in entries.iterdir()
                                 if f.name.lower().endswith(self.fasta_extensions)]

    def seq_to_np_array(self, sequences):
        # convert sequence strings to 2d array
        return np.asarray(list(map(list, sequences)))

    def to_boolean(self, sequence_arrays):
        # convert sequence arrays to boolean values indicating target presence
        return np.where(np.isin(sequence_arrays, self.target), 1, 0)

    def window_counts(self, sequence):
        # convolution weights
        weights = np.ones(self.window_size, dtype=int)
        # calculate sum of boolean values for each window in a sequence
        return np.convolve(sequence, weights, mode='valid')

    def sequence_window_counts(self, is_target):
        return np.asarray([self.window_counts(seq) for seq in is_target])

    def window_means(self, window_stat):
        return window_stat.mean(axis=0)

    def window_var(self, window_stat):
        return window_stat.var(axis=0)

    def filter_stride(self, sequences):
        self.filtered_indices = np.arange(
            start=1, stop=sequences.shape[1]+1, step=self.stride)
        return sequences[:, 0::self.stride]

    def window_frequencies(self, sequences):
        sequences_arrays = self.seq_to_np_array(sequences)
        is_target = self.to_boolean(sequences_arrays)
        self.summary.append(self.summary_stats(is_target))
        window_counts = self.sequence_window_counts(is_target)
        window_counts = self.filter_stride(window_counts)
        window_frequencies = window_counts / self.window_size
        return window_frequencies

    def read_fasta(self, infile):
        # read in raw FASTA file
        with open(infile, "r") as f:
            return f.read()

    def process_fasta(self, file_text):
        # extract protein sequence data from FASTA file
        return [''.join(sequence.split('\n')[1:]) for sequence in file_text.split('>')][1:]

    def subset_sequences(self):
        with self.path_to_data as entries:
            fastas = [self.read_fasta(f) for f in entries.iterdir()
                      if f.name.lower().endswith(self.fasta_extensions)]
        return [self.process_fasta(lines) for lines in fastas]

    def subset_frequencies(self):
        # calculate window frequencies for each subset
        self.frequencies = [pd.DataFrame(self.window_frequencies(subset))
                            for subset in self.subset_sequences()]
        # concatenate subset frequencies and add subset column to dataframe
        self.frequencies = pd.concat(
            self.frequencies, keys=self.subset_names).reset_index()
        self.frequencies.drop('level_1', axis=1, inplace=True)
        self.frequencies.rename(columns={'level_0': 'subset'}, inplace=True)
        # rename columns so they match the original positions in their sequence before filtering
        col_names = {col: idx for idx, col in zip(
            self.filtered_indices, self.frequencies.columns[1:])}
        self.frequencies.rename(columns=col_names, inplace=True)

    def freq_window_means(self):
        self.frequency_means = self.frequencies.groupby('subset').mean()

    def summary_stats(self, is_target):
        alignment_length = is_target.shape[1]
        seq_count = is_target.shape[0]
        target_count = np.sum(is_target)
        seq_sums = np.sum(is_target, axis=0)
        print(seq_sums)
        mean = np.mean(seq_sums)
        variance = np.var(seq_sums)
        std = np.sqrt(variance)
        return pd.DataFrame([[alignment_length, seq_count, target_count, mean, variance, std]],
                            columns=['alignment_length', 'sequence_count', 'target_count',
                                     'target_mean', 'target_variance', 'target_std'])

    def fit_elastic_net(self):
        elastic_net = MultiTaskElasticNetCV()
        y = pd.get_dummies(self.frequencies['subset'], drop_first=True)
        elastic_net.fit(self.frequencies.drop('subset', 1), y)
        with open('results/multitask_elastic_net.pickle', 'wb') as f:
            pickle.dump(elastic_net, f)

    def run_data_pipeline(self):
        self.subset_sequences()
        self.subset_frequencies()
        self.frequencies.to_csv('results/window_frequencies.csv', index=False)
        self.freq_window_means()
        self.frequency_means.to_csv('results/window_frequency_means.csv')
        self.summary = pd.concat(self.summary)
        self.summary.index = self.subset_names
        self.summary.index.name = 'subset'
        self.summary.reset_index(inplace=True)
        print(self.summary)
        self.summary.to_csv('results/summary_stats.csv', index=False)
        if self.fit:
            self.fit_elastic_net()

    def plot(self, df):
        plt.rcParams['figure.figsize'] = [15, 6]
        # filter out data below user specified cutoff
        df = df[df > self.cutoff].dropna(axis=1).T
        # add index for plotting
        window_range = np.arange(len(df))
        df = df.stack().reset_index()
        n_subset = df['subset'].nunique()
        df['plot_index'] = np.repeat(window_range, n_subset)
        df.columns = ['Position', 'Subset', 'Frequency', 'Index']
        sns.scatterplot(data=df, x='Index', y='Frequency',
                        hue='Subset', palette=self.colors)
        # get min and max frequency for each window
        min_freqs = df.groupby('Position').min()['Frequency']
        max_freqs = df.groupby('Position').max()['Frequency']
        axes = plt.gca()
        _, y_max = axes.get_ylim()
        plt.ylim((min_freqs.min(), y_max))
        # plot lines extending down to x axis from maximum frequency
        plt.vlines(x=window_range[0::2], ymin=-0.1,
                   ymax=max_freqs[0::2], color='grey', zorder=0, lw=0.6)
        plt.vlines(x=window_range[1::2], ymin=-0.1,
                   ymax=max_freqs[1::2], color='gray', zorder=0, lw=0.6, linestyle='dotted')
        axes.tick_params(axis='x', which='major',
                         labelsize=6, labelrotation=-90)
        axes.set_xticks(window_range)
        axes.set_xticklabels(df['Position'][::n_subset])
        plt.title(
            f"Target{'s' if self.target.size > 1 else ''}: {', '.join(self.target)}")
        plt.xlabel(f"Window position (window size = {self.window_size})")
        plt.ylabel('Mean window frequency')
        plt.tight_layout()
        plt.savefig('results/window_comparisons.png')
        plt.savefig('results/window_comparisons.pdf')
