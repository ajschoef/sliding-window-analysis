import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pickle
import os
from sklearn.linear_model import MultiTaskElasticNetCV
from pathlib import Path


class SlidingWindow:

    def __init__(self, directory_path, amino_acids, window_size, cutoff, stride=None):
        np.warnings.filterwarnings(
            'error', category=np.VisibleDeprecationWarning)
        self.directory_path = Path(directory_path)
        self.amino_acids = np.asarray(list(amino_acids))
        self.window_size = window_size
        self.cutoff = cutoff
        self.stride = stride
        self.frequencies = None
        self.frequency_means = None
        self.filtered_indices = None
        with self.directory_path as entries:
            self.taxa_names = [f.name.split('.')[0] for f in entries.iterdir()]
        self.fasta_extensions = (
            '.fasta', '.fa', '.fna', '.ffn', '.faa', '.frn')

    def process_fasta(self, lines):
        # extract protein sequence data from FASTA file
        return [line[:-1] for line in lines if line[0] not in ('\n', '>')]

    def read_fasta(self, infile):
        # read in raw FASTA file
        with open(infile, "r") as f:
            return f.readlines()

    def seq_to_np_array(self, sequences):
        # convert sequence strings to 2d array
        return np.asarray(list(map(list, sequences)))

    def to_boolean(self, sequence_arrays):
        # convert sequence arrays to boolean values indicating amino_acid presence
        return np.where(np.isin(sequence_arrays, self.amino_acids), 1, 0)

    def window_counts(self, sequence):
        # convolution weights
        weights = np.ones(self.window_size, dtype=int)
        # calculate sum of boolean values for each window in a sequence
        return np.convolve(sequence, weights, mode='valid')

    def sequence_window_counts(self, is_amino_acid):
        return np.asarray([self.window_counts(seq) for seq in is_amino_acid])

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
        is_amino_acid = self.to_boolean(sequences_arrays)
        window_counts = self.sequence_window_counts(is_amino_acid)
        window_counts = self.filter_stride(window_counts)
        window_frequencies = window_counts / self.window_size
        return window_frequencies

    def taxa_sequences(self):
        with self.directory_path as entries:
            fastas = [self.read_fasta(f) for f in entries.iterdir()
                      if f.name.lower().endswith(self.fasta_extensions)]
        return [self.process_fasta(lines) for lines in fastas]

    def taxa_frequencies(self):
        # calculate window frequencies for each taxa
        self.frequencies = [pd.DataFrame(self.window_frequencies(taxa))
                            for taxa in self.taxa_sequences()]
        # concatenate taxa frequencies and add taxa column to dataframe
        self.frequencies = pd.concat(
            self.frequencies, keys=self.taxa_names).reset_index()
        self.frequencies.drop('level_1', axis=1, inplace=True)
        self.frequencies.rename(columns={'level_0': 'taxa'}, inplace=True)
        # rename columns so they match the original positions in their sequence before filtering
        col_names = {col: idx for idx, col in zip(
            self.filtered_indices, self.frequencies.columns[1:])}
        self.frequencies.rename(columns=col_names, inplace=True)

    def freq_window_means(self):
        return self.frequencies.groupby('taxa').mean()

    def fit_elastic_net(self):
        elastic_net = MultiTaskElasticNetCV()
        y = pd.get_dummies(self.frequencies['taxa'], drop_first=True)
        elastic_net.fit(self.frequencies.drop('taxa', 1), y)
        if not os.path.exists('results'):
            os.makedirs('results')
        with open('results/multitask_elastic_net.pickle', 'wb') as f:
            pickle.dump(elastic_net, f)

    def run_data_pipeline(self):
        self.taxa_sequences()
        self.taxa_frequencies()
        return self.freq_window_means()

    def summary_stats(self, is_amino_acid):
        target_count = np.sum(is_amino_acid)
        seq_count = is_amino_acid.shape[0]
        seq_sums = np.sum(is_amino_acid, axis=0)
        mean = np.mean(seq_sums)
        variance = np.var(seq_sums)
        std = np.sqrt(variance)
        arr = np.array([target_count, seq_count, mean, variance, std])
        stat_names = ['target_count', 'seq_count', 'mean', 'variance', 'std']
        stat_names = ', '.join(stat_names)
        if not os.path.exists('results'):
            os.makedirs('results')
        np.savetxt('results/summary_stats.csv', arr, delimiter=', ',
                   header=stat_names, comments='', fmt='%.6f')

    def plot(self, df):
        plt.rcParams['figure.figsize'] = [15, 6]
        # filter out data below user specified cutoff
        df = df[df > self.cutoff].dropna(axis=1).T
        # add index for plotting
        window_range = np.arange(len(df))
        df = df.stack().reset_index()
        df['plot_index'] = np.repeat(window_range, 3)
        df.columns = ['Position', 'Taxa', 'Frequency', 'Index']
        sns.scatterplot(data=df, x='Index', y='Frequency', hue='Taxa')
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
        axes.set_xticklabels(df['Position'][::3])
        plt.title(
            f"Amino acid{'s' if self.amino_acids.size > 1 else ''}: {', '.join(self.amino_acids)}")
        plt.xlabel(f"Window position (window size = {self.window_size})")
        plt.ylabel('Mean window frequency')
        plt.tight_layout()
        plt.show()
