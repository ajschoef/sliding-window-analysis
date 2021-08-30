import numpy as np
import pandas as pd
from pathlib import Path


class SlidingWindow:

    def __init__(self, path_to_data, target, window_size, n_largest, stride):
        # parameters
        self.path_to_data = Path(path_to_data)
        self.target = np.array(list(target))
        self.window_size = window_size
        self.n_largest = n_largest
        self.stride = stride
        # data
        self.frequencies = None
        self.window_data = None
        self.window_data_filtered = None
        self.summary_stats = []
        # helper variables
        self.groupby_window = None
        self.processed_data_path = 'results/processed_data/'
        self.fasta_extensions = (
            '.fasta', '.fa', '.fna', '.ffn', '.faa', '.frn')
        # generate subset labels from filenames
        self.subset_names = self.apply_if_fasta(self.make_subset_name)
        self.subset_names.sort()
        self.n_subset = len(self.subset_names)

    def make_subset_name(self, file):
        return file.name.split('.')[0]

    # returns a boolean indicating if a file's extension matches a FASTA format
    def is_fasta(self, file):
        return file.name.lower().endswith(self.fasta_extensions)

    # applies a function to each file in the data directory if it is a FASTA file
    def apply_if_fasta(self, f):
        with self.path_to_data as entries:
            return [f(file) for file in entries.iterdir() if self.is_fasta(file)]

    # converts multiple alignment strings to a 2d numpy array
    def seq_to_np_array(self, sequences):
        return np.asarray(list(map(list, sequences)))

    # convert alignment arrays to boolean values indicating target presence
    def to_boolean(self, sequence_arrays):
        return np.where(np.isin(sequence_arrays, self.target), 1, 0)

    # returns an array of target counts for each window in a sequence
    def window_counts(self, sequence):
        # convolution weights; vector of ones with length of window size
        weights = np.ones(self.window_size, dtype=int)
        # calculate sum of boolean values for each window in a sequence
        return np.convolve(sequence, weights, mode='valid')

    # generates summary statistics aggregated across all alighments for a single subset
    def subset_summary_stats(self, is_target):
        alignment_length = is_target.shape[1]
        seq_count = is_target.shape[0]
        target_count = np.sum(is_target)
        proportion = np.mean(is_target)
        return pd.DataFrame([[alignment_length, seq_count, target_count, proportion]],
                            columns=['alignment_length', 'sequence_count', 'target_count',
                                     'target_proportion'])

    # returns a 2d array of the target counts in each window (column) for each sequence (row)
    def sequence_window_counts(self, is_target):
        return np.array([self.window_counts(seq) for seq in is_target])

    # returns sequences at each window stride step
    def filter_stride(self, sequences):
        return sequences[:, 0::self.stride]

    # returns a dataframe of sequence headers appended to window frequencies dataframe
    def add_header(self, header, window_frequencies):
        header = pd.Series(header, name='header')
        # name_range = range(1, window_frequencies.shape[1] + 1)
        # [f'window_{i}' for i in name_range]
        columns = range(1, window_frequencies.shape[1] + 1)
        window_frequencies = pd.DataFrame(
            window_frequencies, columns=columns)
        return pd.concat([header, window_frequencies], axis=1)

    # returns a data frame of window frequencies for each sequence in a subset
    def window_frequencies(self, sequences):
        header, sequences = zip(*sequences)
        sequences_arrays = self.seq_to_np_array(sequences)
        is_target = self.to_boolean(sequences_arrays)
        self.summary_stats.append(self.subset_summary_stats(is_target))
        window_counts = self.sequence_window_counts(is_target)
        window_counts = self.filter_stride(window_counts)
        window_frequencies = window_counts / self.window_size
        return self.add_header(header, window_frequencies)

    # read and return raw FASTA file text
    def read_fasta(self, infile):
        with open(infile, "r") as f:
            return f.read()

    # returns a tuble containing the sequence header and the sequence itself
    def separate_header(self, sequence):
        return sequence[0], ''.join(sequence[1:])

    # returns a list of tuples containing a header and sequence for each sequence in a FASTA file
    def process_fasta(self, file_text):
        sequences = file_text.split('>')[1:]
        sequences = [sequence.split('\n') for sequence in sequences]
        return [self.separate_header(sequence) for sequence in sequences]

    # process each FASTA file
    def subset_sequences(self):
        fastas = self.apply_if_fasta(self.read_fasta)
        return [self.process_fasta(fasta) for fasta in fastas]

    # calculate window frequencies for each subset
    def subset_frequencies(self):
        return [self.window_frequencies(subset) for subset in self.subset_sequences()]

    # concatenate subset frequencies and add subset column to dataframe
    def make_frequencies(self):
        self.frequencies = (
            pd.concat(self.subset_frequencies(), keys=self.subset_names)
            .reset_index()
            .drop('level_1', axis=1)
            .rename(columns={'level_0': 'subset'})
        )

    # generate aggregate summary statistics for each subset
    def make_summary_stats(self):
        self.summary_stats = pd.concat(self.summary_stats)
        self.summary_stats.index = self.subset_names
        self.summary_stats.index.name = 'subset'
        subsets = [self.frequencies.loc[self.frequencies['subset'] == subset]
                   for subset in self.subset_names]
        self.summary_stats['frequency_mean'] = [
            np.mean(subset.iloc[:, 2:].values) for subset in subsets]
        self.summary_stats['frequency_variance'] = [
            np.var(subset.iloc[:, 2:].values) for subset in subsets]
        self.summary_stats['frequency_std'] = np.sqrt(
            self.summary_stats['frequency_variance'])
        self.summary_stats.reset_index(inplace=True)

    def make_column_names(self, postfix):
        return [subset.lower() + postfix for subset in self.subset_names]

    # returns a symmetric matrix of the differences between group means for a window
    def delta_matrix(self, group):
        return group['window_mean'].values - group['window_mean'].values[:, None]

    def calculate_deltas(self):
        return [self.delta_matrix(group) for _, group in self.groupby_window]

    def concat_dataframes(self, deltas, effect_sizes):
        self.window_data = pd.concat(
            [self.window_data, deltas, effect_sizes], axis=1)

    # scale variance of each subset relative to their sample size
    def scale_variance(self, sample_size, variance):
        return (sample_size - 1) * variance

    # returns a 2d array of all pairwise sums of a vector's elements
    def pairwise_sum(self, vector):
        return vector.values + vector.values[:, None]

    # check for division by zero and replace with 0 if encountered
    def safe_zero_divide(self, dividend, divisor):
        dividend = dividend.astype(float)
        divisor = divisor.astype(float)
        out = np.zeros_like(divisor)
        return np.divide(dividend, divisor, out=out, where=divisor != 0)

    # calculates pooled standard deviation for each pairwise combination of subset frequency variances
    def pool_std(self, variance, sample_size_sums):
        variance_sums = self.pairwise_sum(variance)
        return np.sqrt(self.safe_zero_divide(variance_sums, sample_size_sums - 2))

    def calculate_cohens_d(self, deltas, pooled_std):
        return self.safe_zero_divide(deltas, pooled_std)

    # returns the bias corrected Cohen's d effect size
    def correct_bias(self, effect_size, n_sums):
        divisor = (4 * n_sums) - 9
        correction_factor = 1 - self.safe_zero_divide(np.array([3.0]), divisor)
        return correction_factor * effect_size

    # calculates Hedge's g effect size for groups with unequal sample sizes
    def calculate_hedges_g(self, delta, group, sample_size, sample_size_sums):
        variance = group['window_variance'].values
        scaled_variance = self.scale_variance(sample_size, variance)
        pooled_std = self.pool_std(scaled_variance, sample_size_sums)
        cohens_d = self.calculate_cohens_d(delta, pooled_std)
        return self.correct_bias(cohens_d, sample_size_sums)

    # calculate Hedge's g effect size for each pairwise combination of subset frequencies
    def calculate_effect_sizes(self, deltas):
        sample_size = self.summary_stats['alignment_length']
        sample_size_sums = self.pairwise_sum(sample_size)
        return [self.calculate_hedges_g(delta, group[1], sample_size, sample_size_sums)
                for delta, group in zip(deltas, self.groupby_window)]

    # calculate and append deltas and effect sizes to window dataframe
    def make_effect_sizes(self):
        deltas = self.calculate_deltas()
        deltas_temp = deltas.copy()
        effect_sizes = self.calculate_effect_sizes(
            deltas)
        deltas = pd.DataFrame(np.concatenate(
            deltas), columns=self.make_column_names('_delta'))
        effect_sizes = pd.DataFrame(np.concatenate(
            effect_sizes), columns=self.make_column_names('_hedges_g'))
        self.concat_dataframes(deltas, effect_sizes)
        return deltas_temp

    # returns the lower half of a symmetric matrix (upper half is redundent) as a flattened array
    def lower_tri_values(self, diff_matrix):
        return diff_matrix[np.triu_indices(diff_matrix.shape[0], k=1)]

    # returns the largest difference between the group means for a window
    def pairwise_abs_max(self, diff_matrix):
        diff_array = self.lower_tri_values(diff_matrix)
        return np.max(abs(diff_array))

    # filter windows with the n (user specified) largest differences between group means
    def make_filtered_data(self, deltas):
        max_deltas = [self.pairwise_abs_max(delta) for delta in deltas]
        indices = np.argpartition(
            max_deltas, -self.n_largest)[-self.n_largest:] + 1
        mask = self.window_data['window'].isin(indices)
        self.window_data_filtered = (
            self.window_data[mask]
            .copy()
            .reset_index(drop=True)
        )

    # convert a scalar from one range to another
    def interpolate_range(self, value, old_min, old_max, new_min, new_max):
        return int(np.interp(value, [old_min, old_max], [new_min, new_max]))

    def make_window_start(self):
        alighment_length = self.summary_stats['alignment_length'][0]
        old_min = self.window_data["window"].values[0]
        old_max = self.window_data["window"].values[-1]
        return [self.interpolate_range(value, old_min, old_max, 1, alighment_length)
                for value in self.window_data["window"]]

    def make_window_data(self):
        groupby_subset = self.frequencies.groupby('subset')
        self.window_data = groupby_subset.mean().T.stack().reset_index()
        variance = groupby_subset.var().T.stack().reset_index(drop=True)
        self.window_data['variance'] = variance
        std = np.sqrt(variance)
        self.window_data['std'] = std
        self.window_data.columns = [
            'window', 'subset', 'window_mean', 'window_variance', 'window_std']
        self.window_data.insert(1, "window_start", self.make_window_start())
        self.groupby_window = self.window_data.groupby('window')
        deltas = self.make_effect_sizes()
        self.make_filtered_data(deltas)

    def save_csv(self, df, filename):
        df.to_csv(f'{self.processed_data_path}{filename}.csv', index=False)

    def save_data(self):
        self.save_csv(self.frequencies, 'window_frequencies')
        self.save_csv(self.window_data, 'window_data')
        self.save_csv(self.window_data_filtered, 'window_data_filtered')
        self.save_csv(self.summary_stats, 'summary_stats')

    def run_pipeline(self):
        self.make_frequencies()
        self.make_summary_stats()
        self.make_window_data()
        self.save_data()
