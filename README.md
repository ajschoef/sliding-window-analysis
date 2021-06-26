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

    docker-compose run app path_to_data amino_acid window_size cutoff --stride --fit

For example,

    docker-compose run app example/path/ G 10 0.15

The path_to_data argument should be the same as the one provided in the docker-compose build command. After running this command, the program's output should now be in the generated "results" subdirectory.
