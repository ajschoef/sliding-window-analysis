import numpy as np
import pandas as pd
from pathlib import Path


class SlidingWindow:
    # FIXME epsilon is not currently in program arguments
    # FIXME make max_iter, epsilon, and minimum_regularization all parameters for --fit
    def __init__(self, path_to_data, target, window_size, cutoff, stride=1, colors=None, fit=None, epsilon=0.05):
        #### parameters ####
        self.path_to_data = Path(path_to_data)
        self.target = np.asarray(list(target))
        self.window_size = window_size
        self.cutoff = cutoff
        self.stride = stride
        self.colors = colors
        self.fit = fit
        self.epsilon = epsilon
        #### data ####
        self.frequencies = None
        self.groupby_subset = None
        self.window_data = None
        self.filtered_data = None
        self.stride_indices = None
        self.summary = []
        #### helper variables ####
        self.processed_data_path = 'results/processed_data/'
        self.fasta_extensions = (
            '.fasta', '.fa', '.fna', '.ffn', '.faa', '.frn')
        # generate subset labels from filenames
        with self.path_to_data as entries:
            self.subset_names = [f.name.split('.')[0]
                                 for f in entries.iterdir()
                                 if f.name.lower().endswith(self.fasta_extensions)]
        self.n_subset = len(self.subset_names)

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
                      if f.name.lower().endswith(self.fasta_extensions)]
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
        self.frequencies = pd.concat(
            self.frequencies, keys=self.subset_names).reset_index()
        self.frequencies.drop('level_1', axis=1, inplace=True)
        self.frequencies.rename(columns={'level_0': 'subset'}, inplace=True)
        # rename columns so they match the original positions in their sequence before filtering by stride
        col_names = {col: idx for idx, col in zip(
            self.stride_indices, self.frequencies.columns[2:])}
        self.frequencies.rename(columns=col_names, inplace=True)
        self.frequencies.to_csv(
            f'{self.processed_data_path}window_frequencies.csv', index=False)
        self.groupby_subset = self.frequencies.groupby('subset')

    def freq_window_means(self):
        self.window_data = self.groupby_subset.mean().T.stack().reset_index()

    def freq_window_variance(self):
        return self.groupby_subset.var().T.stack().reset_index()[0]

    def freq_window_std(self):
        return self.groupby_subset.std().T.stack().reset_index()[0]

    def make_dataframe(self):
        self.freq_window_means()
        self.window_data['variance'] = self.freq_window_variance()
        self.window_data['std'] = self.freq_window_std()
        self.window_data.columns = [
            'position', 'subset', 'window_mean', 'variance', 'std']
        self.window_data.to_csv(
            f'{self.processed_data_path}window_data.csv', index=False)

    def greater_than_cutoff(self, groupby_position):
        is_cutoff = []
        deltas = []
        # calculate pairwise differences between means of subsets for each window
        for _, group in groupby_position:
            diff_matrix = abs(
                group['window_mean'].values - group['window_mean'].values[:, None])
            deltas.append(diff_matrix)
            # if any of the differences in the means of subsets for a window is less than the cutoff, filter it out
            is_cutoff.append(np.any(np.max(diff_matrix, axis=0) > self.cutoff))
        return is_cutoff, deltas

    def filter_by_delta_cutoff(self):
        groupby_position = self.window_data.groupby('position')
        is_cutoff, deltas = self.greater_than_cutoff(groupby_position)
        # if all values are below the user specified cutoff, then filter values below the 3rd quartile of subset differences instead
        if (~np.array(is_cutoff)).all():
            self.cutoff = np.percentile(
                np.array(deltas).reshape(-1), 75)
            is_cutoff, _ = self.greater_than_cutoff(groupby_position)
        is_cutoff = pd.Series(np.repeat(is_cutoff, self.n_subset))
        # remove values below delta cutoff
        self.filtered_data = self.window_data[is_cutoff].copy()
        self.filtered_data.reset_index(inplace=True, drop=True)
        self.filtered_data.to_csv(
            f'{self.processed_data_path}window_data_filtered.csv', index=False)

    def subset_summary_stats(self, is_target):
        alignment_length = is_target.shape[1]
        seq_count = is_target.shape[0]
        target_count = np.sum(is_target)
        seq_sums = np.sum(is_target, axis=0)
        mean = np.mean(seq_sums)
        variance = np.var(seq_sums)
        std = np.sqrt(variance)
        return pd.DataFrame([[alignment_length, seq_count, target_count, mean, variance, std]],
                            columns=['alignment_length', 'sequence_count', 'target_count',
                                     'target_mean', 'target_variance', 'target_std'])

    def make_summary_stats(self):
        self.summary = pd.concat(self.summary)
        self.summary.index = self.subset_names
        self.summary.index.name = 'subset'
        self.summary.reset_index(inplace=True)
        self.summary.to_csv(
            f'{self.processed_data_path}summary_stats.csv', index=False)

    def run_pipeline(self):
        self.subset_sequences()
        self.subset_frequencies()
        self.make_dataframe()
        self.filter_by_delta_cutoff()
        self.make_summary_stats()
