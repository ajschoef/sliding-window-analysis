from sliding_window import sliding_window
import matplotlib.pyplot as plt


def main():

    directory_path = 'data/raw'
    window_size = 10
    amino_acid = "G"
    names = ['psychro', 'meso', 'thermo']

    # create sliding window instance
    sw = sliding_window.SlidingWindow(directory_path, amino_acid, window_size)
    # processed and write FASTA files to local
    sw.write_processed_fasta()
    # calculate window frequencies for each taxa
    taxa_window_freqs = sw.calculate_taxa_freqs()
    # calculate window means for each taxa
    window_freq_means = sw.window_means(taxa_window_freqs)
    # plot window frequencies of all taxa
    sw.plot_freqs(names, window_freq_means)


if __name__ == "__main__":
    main()
