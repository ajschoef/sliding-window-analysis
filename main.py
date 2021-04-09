from sliding_window import sliding_window
import matplotlib.pyplot as plt


def main():

    directory_path = 'data/raw'
    window_size = 10
    amino_acid = "G"

    # create sliding window instance
    sw = sliding_window.SlidingWindow(directory_path, amino_acid, window_size)
    # processed and write FASTA files to local
    sw.write_processed_fasta()
    # calculate window frequencies for each taxa
    taxa_window_freqs = sw.calculate_taxa_freqs()
    # taxa names
    names = ['psychro', 'meso', 'thermo']
    # calculate window means for each taca
    taxa_dict = sw.create_taxa_dict(names, taxa_window_freqs)
    window_freq_means = sw.window_means(taxa_dict)
    # plot window frequencies of all taxa
    sw.plot_freqs(window_freq_means)


if __name__ == "__main__":
    main()
