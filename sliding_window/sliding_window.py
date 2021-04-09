import pandas as pd
import numpy as np
from pathlib import Path


class SlidingWindow:
    def __init__(self, directory_path, amino_acid, window_size):
        self.directory_path = Path(directory_path)
        self.amino_acid = amino_acid
        self.window_size = window_size
        self.num_fasta_files = 0

    # create new FASTA file for edited protein sequence data
    def write_processed_fasta(self, lines, file_idx):
        new_file = open(f'data/processed/processed_data_{file_idx}.fa', "w")
        for line in lines:
            if line[0] != '>':
                if line[0] != '\n':
                    line = line[:-1]
                    new_file.write(line)
            else:
                if line != lines[0]:
                    new_file.write('\n')
        new_file.write('\n')
        new_file.close()

    def process_fasta(self):
        with self.directory_path as entries:
            for file_idx, infile in enumerate(entries.iterdir()):
                if infile.is_file():
                    self.num_fasta_files += 1
                # read in raw FASTA file
                content = open(f'{self.directory_path}/{infile.name}', "r")
                lines = content.readlines()
                content.close()
                self.write_processed_fasta(lines, file_idx)

    def calc_freqs_in_windows(self, sequences):
        seq_window_freqs = np.full(len(sequences), None)
        for seq_idx, seq in enumerate(sequences):
            # convert sequence to a series of booleans indicating presence of specified amino acid
            s = pd.Series(np.where(np.array(list(seq))
                          == self.amino_acid, 1, 0))
            # calculate frequencies of each window in the sequence
            sum_in_windows = s.rolling(self.window_size).sum()
            # remove NAs and reset index to 0
            sum_in_windows = sum_in_windows.dropna().reset_index(drop=True)
            # calculate frequencies
            freqs_in_windows = sum_in_windows / self.window_size
            seq_window_freqs[seq_idx] = freqs_in_windows
        return seq_window_freqs

    def calculate_taxa_freqs(self):
        with Path('data/processed') as entries:
            taxa_window_freqs = np.full(self.num_fasta_files, None)
            for file_idx, taxa_file in enumerate(entries.iterdir()):
                # read in processed fasta file(s)
                new_file = open(taxa_file, "r")
                sequences = new_file.readlines()
                new_file.close()
                # calculate the frequencies in each window of each sequence
                freqs_in_windows = self.calc_freqs_in_windows(sequences)
                taxa_window_freqs[file_idx] = freqs_in_windows
        return taxa_window_freqs

    def create_df_from_freqs(self, taxa_freqs):
        return pd.DataFrame.from_dict(
            {f"seq_{i}": freq for i, freq in enumerate(taxa_freqs, 1)})

    def create_taxa_dict(self, names, taxa_window_freqs):
        return {name: self.create_df_from_freqs(taxa_freqs) for name, taxa_freqs in zip(names, taxa_window_freqs)}

    def window_means(self, taxa_dict):
        return pd.DataFrame({name: df.mean(axis=1) for (name, df) in taxa_dict.items()})

    def window_variances(self, taxa_dict):
        return pd.DataFrame({name: df.var(axis=1) for (name, df) in taxa_dict.items()})
