# CB_Mini_Project
Comp 388 Mini Project

## Getting Started

This project uses Kallisto, Sleuth, Bowtie2 and BLAST via BioPython to
function so all of these programs are required to run the pipeline.

### Installing Kallisto

Download the Kallisto program for your system from the
[Patcher Lab Website](http://pachterlab.github.io/kallisto/download) and put
it in your bin directory. To check if it is successfully installing open
your terminal and run `kallisto --help`. You should get a standard help
screen.

### Installing Sleuth
Sleuth is run as an R package. Follow the instructions provided by the
[Pachter Lab](https://pachterlab.github.io/sleuth/download).

### Installing Bowtie2
Use the command
```
sudo apt-get install bowtie2
```
to install bowtie2. You can also visit the [Bowtie Website](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
for other installation modes.

### Installing BioPython
Assuming you have a version of python3 already installing use
`pip install biopython` to install this package. We will use it for running
BLAST searches over the internet. To test it is working run `import Bio`
in python. If you get an error try `pip3 install biopython`.


## Running the Pipline

To run the pipeline in its entirety your command will look like this.
```
python3 main.py -o [Path\to\your\output\directory]
```
This will download all files and build any required indices which can take
quite some time. So if you would like to just see the functionality
test files can be downloaded from LINK COMING SOON and pass in these files
with program arguements. See the args section for a more detailed look at what
each arguement does. You can also run `python3 main.py --help` for a similar
output.