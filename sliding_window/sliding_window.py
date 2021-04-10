import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt


class SlidingWindow:

    def __init__(self, directory_path, amino_acid, window_size):
        self.directory_path = Path(directory_path)
        self.amino_acid = amino_acid
        self.window_size = window_size
        self.num_fasta_files = 0

    # create new FASTA file for edited protein sequence data
    def process_fasta(self, lines, new_file):
        for line in lines:
            if line[0] != '>':
                if line[0] != '\n':
                    line = line[:-1]
                    new_file.write(line)
            else:
                if line != lines[0]:
                    new_file.write('\n')
        new_file.write('\n')

    def write_processed_fasta(self):
        with self.directory_path as entries:
            for file_idx, infile in enumerate(entries.iterdir()):
                if infile.is_file():
                    self.num_fasta_files += 1
                # read in raw FASTA file
                input_filename = f'{self.directory_path}/{infile.name}'
                with open(input_filename, "r") as f:
                    lines = f.readlines()
                # write processed FASTA file
                output_filename = f'data/processed/processed_data_{file_idx}.fa'
                with open(output_filename, "w") as f:
                    self.process_fasta(lines, f)

    def calc_freqs_in_windows(self, sequences):
        # convert sequence strings to 2d array
        sequences = np.asarray(list(map(list, sequences)))
        # convert sequence arrays to boolean values indicating amino acid presence
        is_amino = np.where(sequences == self.amino_acid, 1, 0)
        # convolution weights
        weights = np.ones(self.window_size, dtype=int)
        # calculate sum of boolean values for each window in each sequence
        window_sums = [np.convolve(a, weights, mode='valid') for a in is_amino]
        window_sums = np.asarray(window_sums)
        # calculate frequencies
        windows_freqs = window_sums / self.window_size
        return windows_freqs

    def calculate_taxa_freqs(self):
        taxa_window_freqs = np.full((self.num_fasta_files), None)
        with Path('data/processed') as entries:
            for file_idx, taxa_file in enumerate(entries.iterdir()):
                # read in processed fasta file(s)
                with open(taxa_file, "r") as f:
                    sequences = f.readlines()
                # calculate the frequencies in each window of each sequence
                freqs_in_windows = self.calc_freqs_in_windows(sequences)
                taxa_window_freqs[file_idx] = freqs_in_windows
        return taxa_window_freqs

    def window_means(self, taxa_window_freqs):
        return np.asarray([taxa.mean(axis=0) for taxa in taxa_window_freqs])

    def window_var(self, taxa_window_freqs):
        return np.asarray([taxa.var(axis=0) for taxa in taxa_window_freqs])

    def plot_freqs(self, names, freq_means):
        plt.rcParams['figure.figsize'] = [15, 5]
        plt.stackplot(range(freq_means[0].size), *freq_means)
        plt.title(f"Amino acid: {self.amino_acid}")
        plt.xlabel(f"Window position (window size = {self.window_size})")
        plt.ylabel('Stacked mean window frequency')
        plt.legend(names)
        plt.show()
