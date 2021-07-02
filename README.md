# sliding-window-analysis

![Example plot](/example_plot.png)

### Setup and run instuctions

Download or clone this repository.

Download and install [Docker](https://www.docker.com/).

Open a command line interface and change the working directory to the sliding-window-analysis directory. For example, in Linux/Unix

    cd example/path/to/directory/sliding-window-analysis

Build the docker image and specify the path to the directory containing your [FASTA files](https://en.wikipedia.org/wiki/FASTA_format)

    docker-compose build --build-arg path_to_data=example/path/

For documentation on the program's arguments, run

    docker-compose run app -h

To excute the program, pass the arguments listed below.

    docker-compose run app path_to_data target window_size cutoff --stride --fit --colors

The path_to_data argument should be the same as the one provided in the docker-compose build command. After running this command, the program's output should now be in the generated "results" subdirectory.

### Examples

With a single target and no optional arguments

    docker-compose run app example/path/ G 20 0.15

With multiple targets and optional window `--stride` and `--colors` arguments

    docker-compose run app example/path/ YRV 10 0.0 --stride 10 --colors blue green yellow

If the optional `--fit` argument is passed, a weighted multinomial logistic regression with an [elastic net](https://en.wikipedia.org/wiki/Elastic_net_regularization) penalty is fit to the data and the resulting model is saved to the results folder. Elastic net hyperparameters are chosen via [cross-validation](<https://en.wikipedia.org/wiki/Cross-validation_(statistics)>)

    docker-compose run app example/path/ R 15 0.2 --fit
