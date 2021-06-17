from sliding_window import sliding_window
import argparse


def main():

    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "directory_path",
        type=str,
        help="The path to the directory containing FASTA files")
    parser.add_argument(
        "amino_acids",
        type=str,
        help="The amino acid(s) to calculate statistics for in each window")
    parser.add_argument(
        "window_size",
        type=int,
        help="The number of sequence positions (>= 1) to include in each window")
    parser.add_argument(
        "cutoff",
        type=float,
        help="Average window frequencies below the cutoff value [0-1] are filtered out")
    parser.add_argument(
        "--stride",
        type=int,
        help="The number of sequence positions (>= 1) between adjacent windows (default = 1)")
    parser.add_argument(
        "--fit",
        help="A multitask elastic net model is fit and saved to: results/multitask_elastic_net.pickle")

    # parse arguments and validate them
    args = parser.parse_args()
    if args.window_size < 1:
        parser.error("Minimum allowed window size is 1")
    if 0 > args.cutoff > 1:
        parser.error("The cutoff must be between 0 and 1 (inclusive)")
    # create sliding window instance with stride if stride argument is provided
    if args.stride:
        if args.stride < 1:
            parser.error("Minimum allowed is 1")
        sw = sliding_window.SlidingWindow(
            args.directory_path, args.amino_acids, args.window_size, args.cutoff, args.stride)
    else:  # otherwise, create instance with stride = 1
        sw = sliding_window.SlidingWindow(
            args.directory_path, args.amino_acids, args.window_size, args.cutoff)

    # calculate window frequencies for each taxa
    taxa_window_freq_means = sw.run_data_pipeline()
    # fit multitask elastic net model and saves it to: results/multitask_elastic_net.pickle
    if args.fit:
        sw.fit_elastic_net()
    # plot window frequencies of all taxa
    sw.plot(taxa_window_freq_means)


if __name__ == "__main__":
    main()
