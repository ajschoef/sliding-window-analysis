import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt


class SlidingWindow:

    def __init__(self, directory_path, codon, window_size):
        self.directory_path = Path(directory_path)
        self.codon = codon
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

    def taxa_freq_means(self):
        taxa_window_freqs = np.full((self.num_fasta_files), None)
        with Path('data/processed') as entries:
            for file_idx, taxa_file in enumerate(entries.iterdir()):
                # read in processed fasta file(s)
                with open(taxa_file, "r") as f:
                    sequences = f.readlines()
                # calculate the frequencies in each window of each sequence
                sequences_arrays = self.seq_to_np_array(sequences)
                is_codon = self.to_boolean(sequences_arrays)
                window_counts = self.window_counts(is_codon)
                window_frequencies = self.windows_freqs(window_counts)
                taxa_freq_means = self.window_means(
                    window_frequencies)
                taxa_window_freqs[file_idx] = taxa_freq_means
        return taxa_window_freqs

    def plot_freqs(self, names, freq_means):
        plt.rcParams['figure.figsize'] = [15, 5]
        plt.stackplot(range(freq_means[0].size), *freq_means)
        plt.title(f"Amino acid: {self.codon}")
        plt.xlabel(f"Window position (window size = {self.window_size})")
        plt.ylabel('Stacked mean window frequency')
        plt.legend(names)
        plt.show()
