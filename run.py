from sliding_window.sliding_window import SlidingWindow
from sliding_window.window_plot import WindowPlot
from sliding_window.model import ElasticNet
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
        help="A multinomial logistic regression with elastic net penalty model is fit and saved to: results/elastic_net.pickle")
    parser.add_argument(
        "--colors",
        nargs='*',
        help="The colors of each sequence in the comparision plots. Can be any color or code accepted by Matplotlib")

    # parse arguments and validate them
    args = vars(parser.parse_args())
    colors = args.pop('colors')
    fit = args.pop('fit')

    if args['window_size'] < 1:
        parser.error("'window_size' must be greater than 0")
    if 0 > args['cutoff'] > 1:
        parser.error("'cutoff' must be between 0 and 1 (inclusive)")
    if args['stride'] and args['stride'] < 1:
        parser.error("'stride' must be greater than 0")

    def make_dir(dir_name):
        # makes a directory if it doesn't already exist
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

    # make results, plots, processed data, and model output directories
    make_dir('results')
    make_dir('results/plots')
    make_dir('results/processed_data')
    make_dir('results/model_output')

    # create sliding window instance
    sw = SlidingWindow(**args)
    # if colors are passed, check if the number of colors matches the number of subsets
    if colors and len(colors) != sw.n_subset:
        parser.error(
            f'the number of colors must equal to the number of subsets ({sw.n_subset})')
    # run data processing pipeline
    print("Processing data...")
    sw.run_pipeline()
    # fit weighted multinomial logistic regression with elastic net penalty
    if fit:
        print("Fitting model...")
        elastic_net = ElasticNet(sw)
        elastic_net.fit_elastic_net()
    # plot window frequencies of all subsets
    window_plot = WindowPlot(sw, colors)
    print("Rendering plots...")
    window_plot.make_plots()


if __name__ == "__main__":
    main()
