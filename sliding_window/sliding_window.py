import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt


class SlidingWindow:

    def __init__(self, directory_path, codon, window_size):
        self.directory_path = Path(directory_path)
        self.codon = codon
        self.window_size = window_size

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

    def to_boolean(self, sequences_arrays):
        # convert sequence arrays to boolean values indicating codon presence
        return np.where(sequences_arrays == self.codon, 1, 0)

    def window_counts(self, is_codon):
        # convolution weights
        weights = np.ones(self.window_size, dtype=int)
        # calculate sum of boolean values for each window in each sequence
        return np.asarray([np.convolve(row, weights, mode='valid')
                           for row in is_codon])

    def windows_freqs(self, window_counts):
        return window_counts / self.window_size

    def window_means(self, window_counts):
        return window_counts.mean(axis=0)

    def window_var(self, window_counts):
        return window_counts.var(axis=0)

    def window_freq_means(self, sequences):
        # calculate the mean frequency in each window of each sequence
        sequences_arrays = self.seq_to_np_array(sequences)
        is_codon = self.to_boolean(sequences_arrays)
        window_counts = self.window_counts(is_codon)
        window_frequencies = self.windows_freqs(window_counts)
        window_freq_means = self.window_means(window_frequencies)
        return window_freq_means

    def run_data_pipeline(self):
        with self.directory_path as entries:
            fastas = [self.read_fasta(f) for f in entries.iterdir()]
            taxa_sequences = [self.process_fasta(lines) for lines in fastas]
            return [self.window_freq_means(taxa) for taxa in taxa_sequences]

    def plot_freqs(self, taxa_window_data):
        plt.rcParams['figure.figsize'] = [15, 5]
        plt.stackplot(range(taxa_window_data[0].size), *taxa_window_data)
        plt.title(f"Amino acid: {self.codon}")
        plt.xlabel(f"Window position (window size = {self.window_size})")
        plt.ylabel('Stacked mean window frequency')
        with self.directory_path as entries:
            taxa_names = [f.name.split('.')[0] for f in entries.iterdir()]
        plt.legend(taxa_names)
        plt.show()
