# CB_Mini_Project
Comp 388 Mini Project

## Getting Started

This project uses a number of differnet programs to function and include...
- [Kallisto](#Installing-Kallisto)
- [Bowtie2](#Installing-Bowtie2)
- [Samtools](#Installing-Samtools)
- [fastq-dump](#Installing-fastq-dump)
- [Sleuth (R package)](#Installing-Sleuth)
- [optparse (R package)](#Installing-Sleuth)
- [Biopython (Python package)](#Installing-Biopython)

If you are missing any of these programs click the respective link for info on how to download.

## Running with Test Data

The most straightforward way to run the pipeline is using the command below.
```
python3 main.py -test 1 -local 1 
```
Setting `test 1` tells the pipeline to use the included test data and `-local 1` tells the pipeline to run BLAST locally. This will write all results to a new directory called miniProject_Ethan_Holleman in your current working directory. Additionally if `main.py` is not in your current working directory please specify the path to it.  
To run the pipeline with test data and BLAST over the internet use the command below.
```
python3 main.py -test 1
```
In both cases a log file named `miniProject.log` will be written to miniProject_Ethan Holleman.

## Running the Pipeline in Full

To run the pipeline in its entirety your command will look like this.
```
python3 main.py -o [Path/to/your/output/directory -l [Path/to/your/output/directory]
```
Where -o is the path to output directory and -l is path where the log file will be written. 
 
This will download all files and build any required indices which can take
quite some time. So if you would like to just see the functionality you can run the program using the included test data. More details on that below.


# Install Guide

### Installing Kallisto

Download the Kallisto program for your system from the
[Patcher Lab Website](http://pachterlab.github.io/kallisto/download) and put
it in your bin directory. To check if it is successfully installing open
your terminal and run `kallisto --help`. You should get a standard help
screen.

### Installing Sleuth and/or Optparse
Sleuth is run as an R package. Follow the instructions provided by the
[Pachter Lab](https://pachterlab.github.io/sleuth/download).
In addition make sure you have the R package optparse installed. If not you can install it using the command `install.packages('optparse')`. This will allow the pipeline to pass file paths from python to the R script where sleuth is run from.

### Installing Bowtie2
Use the command
```
sudo apt-get install bowtie2
```
to install bowtie2. You can also visit the [Bowtie Website](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
for other installation modes.

### Installing Biopython
Assuming you have a version of python3 already installing use
`pip install biopython` to install this package. We will use it for running
BLAST searches over the internet. To test it is working run `import Bio`
in python. If you get an error try `pip3 install biopython`.

### Installing Samtools
Run the command `sudo apt-get install samtools ` or [follow this guide](https://www.biostars.org/p/328831/)

### Installing fastq-dump
Fastq-dump is a command available from SRA toolkit you can access the downloads page [here](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software) or install using the command `sudo apt-get install sra-toolkit`
