from sliding_window import sliding_window
import sys


def main():

    try:
        directory_path = sys.argv[1]
        window_size = int(sys.argv[2])
        amino_acid = sys.argv[3]
    except IndexError:
        raise SystemExit(
            f"Usage: {sys.argv[0]} <path_to_fasta_files> <codon> <window_size>")

    # create sliding window instance
    sw = sliding_window.SlidingWindow(directory_path, amino_acid, window_size)
    # calculate window frequencies for each taxa
    taxa_window_freq_means = sw.run_data_pipeline()
    # plot window frequencies of all taxa
    sw.plot_freqs(taxa_window_freq_means)


if __name__ == "__main__":
    main()
