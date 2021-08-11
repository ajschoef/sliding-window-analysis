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
        help="The number of sequence positions (integer >= 1) to include in each window")
    parser.add_argument(
        "n_largest",
        type=int,
        help="Display n (integer >= 1) windows with the largest mean frequency difference between subsets")
    parser.add_argument(
        "--stride",
        type=int,
        help="The number of sequence positions (integer >= 1) between adjacent windows (default = 1)")
    parser.add_argument(
        "--colors",
        nargs='*',
        help="The colors of each subset in the comparision plots. Can be any color or code accepted by Matplotlib")

    # parse arguments and validate them
    args = vars(parser.parse_args())
    colors = args.pop('colors')

    pos_int_text = "must be a positive integer"
    if args['window_size'] < 1:
        parser.error(f"'window_size' {pos_int_text}")
    if args['n_largest'] < 1:
        parser.error(f"'n_largest' {pos_int_text}")
    if args['stride'] and args['stride'] < 1:
        parser.error(f"'stride' {pos_int_text}")
    if not args['stride']:
        args['stride'] = 1

    # makes a directory if it doesn't already exist
    def make_dir(dir_name):
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

    # make results, plots, processed data directories
    directories = ['', 'plots', 'processed_data']
    for dir_name in directories:
        make_dir(f'results/{dir_name}')

    # create sliding window instance
    sw = SlidingWindow(**args)
    # if colors are passed, check if the number of colors matches the number of subsets
    if colors and len(colors) != sw.n_subset:
        parser.error(
            f'the number of colors must equal to the number of subsets ({sw.n_subset})')
    # run data processing pipeline
    print("Processing data...")
    sw.run_pipeline()
    print(f"Processed data is in {sw.processed_data_path}")
    # generate all plots
    wp = WindowPlot(sw, colors)
    print("Rendering plots...")
    wp.make_plots()
    print(f"Plots are in {wp.plots_path}")


if __name__ == "__main__":
    main()
