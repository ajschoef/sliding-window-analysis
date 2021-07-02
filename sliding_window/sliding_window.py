import numpy as np
import pandas as pd
import pickle
from numpy import linalg as LA
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
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
        self.freq_means_wide = None
        self.freq_means_long = None
        self.filtered_data = None
        self.stride_indices = None
        self.summary = []
        #### helper variables ####
        self.fasta_extensions = (
            '.fasta', '.fa', '.fna', '.ffn', '.faa', '.frn')
        # generate subset labels from filenames
        with self.path_to_data as entries:
            self.subset_names = [f.name.split('.')[0]
                                 for f in entries.iterdir()
                                 if f.name.lower().endswith(self.fasta_extensions)]
        self.n_subset = len(self.subset_names)
        self.window_range = None

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
        self.stride_indices = np.arange(
            start=1, stop=sequences.shape[1]+1, step=self.stride)
        return sequences[:, 0::self.stride]

    def window_frequencies(self, sequences):
        sequences_arrays = self.seq_to_np_array(sequences)
        is_target = self.to_boolean(sequences_arrays)
        self.summary.append(self.subset_summary_stats(is_target))
        window_counts = self.sequence_window_counts(is_target)
        window_counts = self.filter_stride(window_counts)
        return window_counts / self.window_size

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
            self.stride_indices, self.frequencies.columns[1:])}
        self.frequencies.rename(columns=col_names, inplace=True)
        self.frequencies.to_csv('results/window_frequencies.csv', index=False)

    def freq_window_means(self):
        self.freq_means_wide = self.frequencies.groupby('subset').mean()

    def to_long_format(self):
        self.freq_means_long = self.freq_means_wide.T
        self.freq_means_long = self.freq_means_long.stack().reset_index()
        self.freq_means_long.columns = [
            'position', 'subset', 'frequency']
        self.freq_means_long.to_csv('results/window_frequency_means.csv')

    def filter_by_delta_cutoff(self):
        grouped_df = self.freq_means_long.groupby('position')
        greater_than_cutoff = []
        for _, group in grouped_df:
            difference_matrix = abs(
                group['frequency'].values - group['frequency'].values[:, None])
            greater_than_cutoff.append(
                np.any(np.max(difference_matrix, axis=0) > self.cutoff))
        # FIXME breaks when all are False. Use try/except and set cutoff to median in except
        greater_than_cutoff = pd.Series(
            np.repeat(greater_than_cutoff, self.n_subset))
        self.filtered_data = self.freq_means_long[greater_than_cutoff].copy()
        self.filtered_data.reset_index(inplace=True, drop=True)
        self.filtered_data.to_csv(
            'results/filtered_window_frequency_means.csv')
        self.window_range = np.arange(
            self.filtered_data['position'].nunique())
        self.filtered_data['x_ticks'] = np.repeat(
            self.window_range, self.n_subset)

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
        self.summary.to_csv('results/summary_stats.csv', index=False)

    # FIXME save how many collinear neighbors were filtered for display later
    def filter_collinear_neighbors(self, df):

        def unit_vector(vector):
            return vector / LA.norm(vector)

        def angle_between(v0, v1):
            return np.arccos(np.clip(np.dot(unit_vector(v0), unit_vector(v1)), -1.0, 1.0))

        mask = [angle_between(df[i], df[i+1]) >
                self.epsilon for i in df.iloc[:, :-1]]
        mask = [True] + mask
        indices = np.where(mask)
        return indices, df.iloc[:, mask]

    # FIXME make modeling its own file/class
    def fit_elastic_net(self):
        # set design matrix and outcome
        y = self.frequencies['subset']
        X = self.frequencies.drop('subset', 1)
        # FIXME indices are currently not used; will be used to find neighbors
        indices, X = self.filter_collinear_neighbors(X)
        # scale data to 0 mean and unit variance
        scaler = StandardScaler()
        scaler.fit(X)
        X = scaler.transform(X)
        # grid search ranges for regularization parameters
        Cs = 1 / np.power(10, np.arange(3, 7, step=0.1))
        l1_ratio = [0.1, 0.25, 0.5, 0.75, 0.9]
        """ tune hyperparameters for weighted multinomial logistic regression with
            elastic net penalty via stratified k-fold cross validation"""
        elastic_net = LogisticRegressionCV(multi_class='multinomial',
                                           penalty='elasticnet',
                                           class_weight='balanced',
                                           solver='saga',
                                           max_iter=1000,
                                           Cs=Cs,
                                           l1_ratios=l1_ratio)
        elastic_net.fit(X, y)
        with open('results/logistic_regression_elastic_net.pickle', 'wb') as f:
            pickle.dump(elastic_net, f)

    def run_pipeline(self):
        self.subset_sequences()
        self.subset_frequencies()
        self.freq_window_means()
        self.to_long_format()
        self.filter_by_delta_cutoff()
        self.make_summary_stats()
        if self.fit:
            self.fit_elastic_net()
