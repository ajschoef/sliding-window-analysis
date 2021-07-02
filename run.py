from sliding_window.sliding_window import SlidingWindow
from sliding_window.window_plot import WindowPlot
import argparse
import os


def main():

    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "path_to_data",
        type=str,
        help="The path to the directory containing FASTA files")
    parser.add_argument(
        "target",
        type=str,
        help="The target(s) to calculate statistics for in each window")
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
        action='store_true',
        help="A multitask elastic net model is fit and saved to: results/multitask_elastic_net.pickle")
    parser.add_argument(
        "--colors",
        nargs='*',
        help="The colors of each sequence in the comparision plots. Can be any color or code accepted by Matplotlib")

    # parse arguments and validate them
    args = parser.parse_args()
    if args.window_size < 1:
        parser.error("'window_size' must be greater than 0")
    if 0 > args.cutoff > 1:
        parser.error("'cutoff' must be between 0 and 1 (inclusive)")
    if args.stride and args.stride < 1:
        parser.error("'stride' must be greater than 0")

    # make results directory if it doesn not already exist
    if not os.path.exists('results'):
        os.makedirs('results')
    # create sliding window instance and run data pipeline
    sw = SlidingWindow(**vars(args))
    sw.run_pipeline()
    # plot window frequencies of all subsets
    window_plot = WindowPlot(sw)
    window_plot.dot_plot()
    window_plot.ecdf_plot()
    window_plot.kde_plot()


if __name__ == "__main__":
    main()
