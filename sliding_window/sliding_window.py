import numpy as np
import pandas as pd
from pathlib import Path


class SlidingWindow:

    def __init__(self, path_to_data, target, window_size, n_largest, stride=1):
        # parameters
        self.path_to_data = Path(path_to_data)
        self.target = np.asarray(list(target))
        self.window_size = window_size
        self.n_largest = n_largest
        self.stride = stride
        # data
        self.frequencies = None
        self.groupby_subset = None
        self.window_data = None
        self.filtered_data = None
        self.stride_indices = None
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
        # convolution weights
        weights = np.ones(self.window_size, dtype=int)
        # calculate sum of boolean values for each window in a sequence
        return np.convolve(sequence, weights, mode='valid')

    def sequence_window_counts(self, is_target):
        return np.asarray([self.window_counts(seq) for seq in is_target])

    def filter_stride(self, sequences):
        self.stride_indices = np.arange(
            start=1, stop=sequences.shape[1]+1, step=self.stride)
        return sequences[:, 0::self.stride]

    def read_fasta(self, infile):
        # read in raw FASTA file
        with open(infile, "r") as f:
            return f.read()

    def process_fasta(self, file_text):
        # extract protein sequence data from FASTA file
        sequences = file_text.split('>')[1:]
        sequences = [sequence.split('\n') for sequence in sequences]
        return [(sequence[0], ''.join(sequence[1:]))
                for sequence in sequences]
        # headers, sequences = zip(*sequences)

    def subset_sequences(self):
        with self.path_to_data as entries:
            fastas = [self.read_fasta(f) for f in entries.iterdir()
                      if self.is_fasta(f)]
        return [self.process_fasta(fasta) for fasta in fastas]

    def window_frequencies(self, sequences):
        header, sequences = zip(*sequences)
        sequences_arrays = self.seq_to_np_array(sequences)
        is_target = self.to_boolean(sequences_arrays)
        self.summary.append(self.subset_summary_stats(is_target))
        window_counts = self.sequence_window_counts(is_target)
        window_counts = self.filter_stride(window_counts)
        window_frequencies = window_counts / self.window_size
        return pd.concat([pd.Series(header, name='header'),
                          pd.DataFrame(window_frequencies)], axis=1)

    def subset_frequencies(self):
        # calculate window frequencies for each subset
        self.frequencies = [self.window_frequencies(subset)
                            for subset in self.subset_sequences()]
        # concatenate subset frequencies and add subset column to dataframe
        self.frequencies = (
            pd.concat(self.frequencies, keys=self.subset_names)
            .reset_index()
            .drop('level_1', axis=1)
            .rename(columns={'level_0': 'subset'})
        )
        # rename columns so they match the original positions in their sequence before filtering by stride
        col_names = {col: idx for idx, col in zip(
            self.stride_indices, self.frequencies.columns[2:])}
        self.frequencies.rename(columns=col_names, inplace=True)
        self.frequencies.to_csv(
            f'{self.processed_data_path}window_frequencies.csv', index=False)
        self.groupby_subset = self.frequencies.groupby('subset')

    def make_window_data(self):
        self.window_data = self.groupby_subset.mean().T.stack().reset_index()
        variance = self.groupby_subset.var().to_numpy().flatten()
        std = np.sqrt(variance)
        self.window_data['variance'] = variance
        self.window_data['std'] = std
        self.window_data.columns = [
            'window', 'subset', 'window_mean', 'window_variance', 'window_std']
        self.window_data.to_csv(
            f'{self.processed_data_path}window_data.csv', index=False)

    # returns a symmetric matrix of the differences between group means for a window
    def difference_matrix(self, group):
        return abs(group['window_mean'].values - group['window_mean'].values[:, None])

    # returns lower half of symmetric matrix (upper half is redundent) as a flattened array
    def lower_tri_values(self, diff_matrix):
        return diff_matrix[np.triu_indices(diff_matrix.shape[0], k=1)]

    # returns the largest difference between the group means for a window
    def pairwise_max(self, group):
        diff_matrix = self.difference_matrix(group)
        diff_array = self.lower_tri_values(diff_matrix)
        return np.max(diff_array)

    def filter_data(self):
        groupby_window = self.window_data.groupby('window')
        max_diff = [self.pairwise_max(group) for _, group in groupby_window]
        # get windows with the n (user specified) largest differences between group means
        indices = np.argpartition(max_diff, -self.n_largest)[-self.n_largest:]
        self.filtered_data = (
            self.window_data[self.window_data['window'].isin(indices+1)]
            .copy()
            .reset_index(drop=True)
        )
        self.filtered_data.to_csv(
            f'{self.processed_data_path}window_data_filtered.csv', index=False)

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
        self.summary = (
            self.summary
            .reset_index()
            .to_csv(f'{self.processed_data_path}summary_stats.csv', index=False)
        )

    def run_pipeline(self):
        self.subset_sequences()
        self.subset_frequencies()
        self.make_window_data()
        self.filter_data()
        self.make_summary_stats()
