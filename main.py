from sliding_window import sliding_window
import matplotlib.pyplot as plt


def main():

    directory_path = 'data/raw'
    window_size = 10
    amino_acid = "G"

    sw = sliding_window.SlidingWindow(directory_path, amino_acid, window_size)
    sw.process_fasta()
    taxa_window_freqs = sw.calculate_taxa_freqs()

    names = ['psychro', 'meso', 'thermo']
    taxa_dict = sw.create_taxa_dict(names, taxa_window_freqs)

    window_freq_means = sw.window_means(taxa_dict)

    # Plot mean window frequency for each taxa
    plt.rcParams['figure.figsize'] = [15, 5]
    cols = [window_freq_means[col] for col in window_freq_means]
    plt.stackplot(window_freq_means.index, *cols,
                  labels=window_freq_means.columns)
    plt.title(f"Amino acid: {amino_acid}")
    plt.xlabel(f"Window position (window size = {window_size})")
    plt.ylabel('Stacked mean window frequency')
    plt.legend(["Psychro", "Meso", "Thermo"])
    plt.show()


if __name__ == "__main__":
    main()
