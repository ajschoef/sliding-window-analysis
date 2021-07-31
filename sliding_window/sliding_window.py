import numpy as np
import pandas as pd
from pathlib import Path


class SlidingWindow:

    def __init__(self, path_to_data, target, window_size, n_largest, stride=1):
        # parameters
        self.path_to_data = Path(path_to_data)
        self.target = np.array(list(target))
        self.window_size = window_size
        self.n_largest = n_largest
        self.stride = stride
        # data
        self.frequencies = None
        self.window_data = None
        self.filtered_data = None
        self.summary = []
        # helper variables
        self.processed_data_path = 'results/processed_data/'
        self.fasta_extensions = (
            '.fasta', '.fa', '.fna', '.ffn', '.faa', '.frn')
        # generate subset labels from filenames
        with self.path_to_data as entries:
            self.subset_names = [f.name.split('.')[0]
                                 for f in entries.iterdir()
                                 if self.is_fasta(f)]
        self.n_subset = len(self.subset_names)

    def is_fasta(self, file):
        return file.name.lower().endswith(self.fasta_extensions)

    def seq_to_np_array(self, sequences):
        # convert sequence strings to 2d array
        return np.asarray(list(map(list, sequences)))

    def to_boolean(self, sequence_arrays):
        # convert sequence arrays to boolean values indicating target presence
        return np.where(np.isin(sequence_arrays, self.target), 1, 0)

    def window_counts(self, sequence):
        # convolution weights; vector of ones with length of window size
        weights = np.ones(self.window_size, dtype=int)
        # calculate sum of boolean values for each window in a sequence
        return np.convolve(sequence, weights, mode='valid')

    def sequence_window_counts(self, is_target):
        return np.asarray([self.window_counts(seq) for seq in is_target])

    def filter_stride(self, sequences):
        return sequences[:, 0::self.stride]

    # read and return raw FASTA file text
    def read_fasta(self, infile):
        with open(infile, "r") as f:
            return f.read()

    # returns a tuble containing the sequence header and the sequence itself
    def separate_header(self, sequence):
        return sequence[0], ''.join(sequence[1:])

    # returns list tuples containing the headers and sequences from a FASTA file
    def process_fasta(self, file_text):
        sequences = file_text.split('>')[1:]
        sequences = [sequence.split('\n') for sequence in sequences]
        return [self.separate_header(sequence) for sequence in sequences]

    def subset_sequences(self):
        with self.path_to_data as entries:
            fastas = [self.read_fasta(f) for f in entries.iterdir()
                      if self.is_fasta(f)]
        return [self.process_fasta(fasta) for fasta in fastas]

    def add_header(self, header, window_frequencies):
        return pd.concat([pd.Series(header, name='header'),
                          pd.DataFrame(window_frequencies)], axis=1)

    def window_frequencies(self, sequences):
        header, sequences = zip(*sequences)
        sequences_arrays = self.seq_to_np_array(sequences)
        is_target = self.to_boolean(sequences_arrays)
        self.summary.append(self.subset_summary_stats(is_target))
        window_counts = self.sequence_window_counts(is_target)
        window_counts = self.filter_stride(window_counts)
        window_frequencies = window_counts / self.window_size
        return self.add_header(header, window_frequencies)

    # calculate window frequencies for each subset
    def make_frequencies(self):
        return [self.window_frequencies(subset) for subset in self.subset_sequences()]

    # concatenate subset frequencies and add subset column to dataframe
    def subset_frequencies(self):
        self.frequencies = (
            pd.concat(self.make_frequencies(), keys=self.subset_names)
            .reset_index()
            .drop('level_1', axis=1)
            .rename(columns={'level_0': 'subset'})
        )

    # returns a symmetric matrix of the differences between group means for a window
    def delta_matrix(self, group):
        return group['window_mean'].values - group['window_mean'].values[:, None]

    def make_delta_names(self):
        delta_names = self.subset_names.copy()
        delta_names.sort()
        return [subset.lower() + '_delta' for subset in delta_names]

    def make_deltas(self, groupby_window):
        return [self.delta_matrix(group) for _, group in groupby_window]

    def append_deltas(self, group, diff_matrix, names):
        delta_matrix = pd.DataFrame(diff_matrix, columns=names)
        return pd.concat([group.reset_index(), delta_matrix], axis=1)

    def concat_deltas(self, deltas, groupby_window):
        delta_names = self.make_delta_names()
        self.window_data = [self.append_deltas(group[1], diff_matrix, delta_names)
                            for diff_matrix, group in zip(deltas, groupby_window)]
        self.window_data = (
            pd.concat(self.window_data)
            .reset_index(drop=True)
            .drop(labels='index', axis=1)
        )

    # returns the lower half of symmetric matrix (upper half is redundent) as a flattened array
    def lower_tri_values(self, diff_matrix):
        return diff_matrix[np.triu_indices(diff_matrix.shape[0], k=1)]

    # returns the largest difference between the group means for a window
    def pairwise_abs_max(self, diff_matrix):
        diff_array = self.lower_tri_values(diff_matrix)
        return np.max(abs(diff_array))

    # filter windows with the n (user specified) largest differences between group means
    def make_filtered_data(self, deltas):
        maxs = [self.pairwise_abs_max(delta) for delta in deltas]
        indices = np.argpartition(maxs, -self.n_largest)[-self.n_largest:]
        mask = self.window_data['window'].isin(indices)
        self.filtered_data = (
            self.window_data[mask]
            .copy()
            .reset_index(drop=True)
        )

    def make_window_data(self):
        groupby_subset = self.frequencies.groupby('subset')
        self.window_data = groupby_subset.mean().T.stack().reset_index()
        variance = groupby_subset.var().T.stack().reset_index(drop=True)
        std = np.sqrt(variance)
        self.window_data['variance'] = variance
        self.window_data['std'] = std
        self.window_data.columns = [
            'window', 'subset', 'window_mean', 'window_variance', 'window_std']
        groupby_window = self.window_data.groupby('window')
        deltas = self.make_deltas(groupby_window)
        self.concat_deltas(deltas, groupby_window)
        self.make_filtered_data(deltas)

    def subset_summary_stats(self, is_target):
        alignment_length = is_target.shape[1]
        seq_count = is_target.shape[0]
        target_count = np.sum(is_target)
        proportion = np.mean(is_target)
        return pd.DataFrame([[alignment_length, seq_count, target_count, proportion]],
                            columns=['alignment_length', 'sequence_count', 'target_count',
                                     'target_proportion'])

    def make_summary_stats(self):
        self.summary = pd.concat(self.summary)
        self.summary.index = self.subset_names
        self.summary.index.name = 'subset'
        subsets = [self.frequencies.loc[self.frequencies['subset'] == subset]
                   for subset in self.subset_names]
        self.summary['frequency_mean'] = [
            np.mean(subset.iloc[:, 2:].values) for subset in subsets]
        self.summary['frequency_variance'] = [
            np.var(subset.iloc[:, 2:].values) for subset in subsets]
        self.summary['frequency_std'] = np.sqrt(
            self.summary['frequency_variance'])
        self.summary.reset_index(inplace=True)

    def save_csv(self, df, filename):
        df.to_csv(f'{self.processed_data_path}{filename}.csv', index=False)

    def save_data(self):
        self.save_csv(self.frequencies, 'window_frequencies')
        self.save_csv(self.window_data, 'window_data')
        self.save_csv(self.filtered_data, 'window_data_filtered')
        self.save_csv(self.summary, 'summary_stats')

    def run_pipeline(self):
        self.subset_sequences()
        self.subset_frequencies()
        self.make_window_data()
        self.make_summary_stats()
        self.save_data()
