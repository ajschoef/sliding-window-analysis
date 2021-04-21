from sliding_window import sliding_window
import argparse


def main():

    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "directory_path",
        help="the path to the directory containing FASTA files")
    parser.add_argument(
        "amino_acids",
        help="the amino acid(s) to calculate statistics for in each window")
    parser.add_argument(
        "window_size",
        help="the number of positions to include in each window", type=int)
    args = parser.parse_args()

    # create sliding window instance
    sw = sliding_window.SlidingWindow(
        args.directory_path, args.codon, args.window_size)
    # calculate window frequencies for each taxa
    taxa_window_freq_means = sw.run_data_pipeline()
    # plot window frequencies of all taxa
    sw.plot_freqs(taxa_window_freq_means)


if __name__ == "__main__":
    main()
