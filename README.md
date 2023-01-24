# CNAsim

CNAsim is a cancer evolution simulator for generating single-cell data. There are two types of data that can be generated: realistic copy number data mimicing single-cell CNA detection algorithms, and sequencing reads generated from different single-cell sequencing platforms.

More information can be found in our paper, located here: [link to paper].

# Table of Contents

[Installation](#Installation)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [Python packages](#python-packages)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [External packages](#external-packages)

[Usage](#Usage)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [I/O and Utilities](#io-and-utilities)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [Stage 1: simulating cell lineage tree](#stage-1)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [Stage 2: simulating genomes and mutations](#stage-2)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [Stage 3: generating single-cell data](#stage-3)

[Examples](#Examples)

# Installation

CNAsim is written in python and can run reliably on version 3.7 or later. You can download the source code by cloning this repository:

```
git clone 
```

### Python packages

CNAsim requires the installation of the following python packages:
* [Numpy](https://numpy.org/)
* [Scipy](https://scipy.org/)
* [Pyfaidx](https://github.com/mdshw5/pyfaidx)

If not already installed, it is highly recommended that you install python packages with either `pip` or `conda`. If you use `pip`, you can install the packages with

```
pip install numpy scipy pyfaidx
```

If you use conda, the pyfaidx needs to be installed form the `bioconda` channel. For best practices, create a new environment before installing. You can install the packages with
```
conda create -n CNAsim
conda activate CNAsim

conda config --env --add channels conda-forge 
conda config --env --add channels bioconda

conda install numpy scipy
conda install -c bioconda pyfaidx
```

### External packages

Additionally, CNAsim requires that the following binaries are installed and configured on the environment.
* [samtools](http://www.htslib.org/download/) 
* [ms](https://home.uchicago.edu/~rhudson1/source/mksamples.html)
* [dwgsim](https://github.com/nh13/DWGSIM)

You may follow the installation guides with the links provided. For convienience, we also describe below how to install all the necessary packages.
```
Not implemented.
```
# Usage

CNAsim is run through the command line by simply calling 

```
python main.py 
```

CNAsim consists of three phases: 1) Generate the cell lineage tree, 2) simulate genomes and mutations, and 3) generate single-cell data. The following documentation details the relevant parameters associated with each phase.

### I/O and Utilities

* Commands for specifying the output directory, saving simulation details, and an option to parse parameters from a file. The parameter file passed must follow a specific format to be parsed correctly. An example is found in the *parameters* file, located in this repository. One can simply edit this file, or follow the format in another file.

   `-o, --out-path` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Path to output directory where data will be | |saved to. Will create one if no directory exists. Default: current directory

   `-s, --summary` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Saves the details of the simulation to a file. Default: False

   `-F, --param-file` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Parses the file given by the provided path for all simulation parameter. See the *parameters* file for the format. Takes input parameters directly from input if None is provided. Default: None

### Stage 1

* Commands for generating the cell lineage tree. Note that, invoking *ms* requires the directory of the installed binary to be part of the user's PATH variable. The user can specify any combination of tumor, normal, and pseudonormal cells using the -n, -n1, and -n2 parameters. The number of tumor cells will equal the value passed with -n if both normal and pseudonormal fractions are set to 0. Otherwise, the number of tumor cells will be the remaining fraction. For example, if -n is set to 100, -n1 is set to 0.1, and -n2 is set to 0.05, there will be 85 tumor cells, 10 normal cells, and 5 pseudonormal cells. 

   `-t, --tree-type` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Method to generate topology of main tumor lineage. Option 0 uses *ms* to generate a tree under neutral coalescence. Option 1 generates a random tree topology. Option 2 takes an existing tree from the user in newick format (see the following parameter). Default: 0 

   `-T, --tree-path` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Path to input tree if user selection option 2 as the tree type.

   `-g, --growth-rate` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The exponential growth rate parameter used to generate the tree under neutral coalescence. Default: 15.1403

   `-n, --num-cells` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Number of observed cells in the experiment. Includes all types of cells (tumor, normal, and pseudonormal). Default: 250

   `-n1, --normal-fraction` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Proportion of the observed cells that are normal diploid cells. Default: 0

   `-n2, --pseudonormal-fraction` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Proportion of observed cells that are pseudonormal cells. Default: 0

   `-c, --num-clones` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Number of ancestral nodes in the tumor lineage to be selected as clonal founders. Default: 0

### Stage 2

* Commands for simulating genomes. The minimum copy number length is considered the minimum resolution of the genome used in the simulator. Genomes are subsequently represented as a sequence of regions with length equal to the value given by -k. If a reference genome is given, the number of chromosomes and chromosome-lengths are derived from the reference. Without a given reference genome, chromosome lengths derived from hg38 can be used by toggling -U. If this parameter is not toggled, the number of chromosomes will be equal to the value given by -N with set lengths given by -L. Chromosome arm ratios can either be fixed with the -A parameter, which details the ratio of the long arm, or by using the chromosome arm ratios from hg38 again by toggling -U.

   `-k, --min-cn-length` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The minimum length of a copy number event. Default: 1000 bp

   `-r, --reference` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The path to the input reference genome.

   `-U, --use-hg38-static` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use chromosome lengths and chromosome arm ratios derived from the hg38 reference genome. Default: False

   `-N, --num-chromosomes` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A set number of chromosomes to represent the genome. Default: 22

   `-L, --chrom-length` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A fixed chromosome length if set number of chromosomes used. Default: 100,000,000 bp

   `-A, --chrom-arm-ratio` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A fixed chromosome-arm ratio if arm-ratios are not derived from the hg38 reference. Default: 0.5

* Commands for simulating focal copy number aberrations. The number of focal CNAs chosen along each edge is controlled with the -p1 and -p2 parameters. If option 1 is chosen for the placement type, this means that edges directly above leaves in the tree have a mean number of events equal to the value given by -p2. Under neutral coalescence, this likely means more events will be assigned to edges higher up in the tree. The length of focal CNAs cannot drop below the minimum given by -k, and cannot exceed the length of the chromosome. The -s parameter controls the rate at which an event is chosen to be an amplifciation. If this value is given by *x*, then the rate of deletion is simply *1-x*. 

   `-p1, --placement-type` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The method for choosing the number of events along each edge. Option 0 draws from a Poisson distribution with a fixed mean. Option 1 draws from a Poisson with a mean proportional to edge length. Option 2 is a fixed number per edge. Default: 0

   `-p2, --placement-param` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The value to be used in conjunction with the placement method chosen. Default: 2

   `-l, --cn-length-mean` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The average length of a focal CNA in number of base pairs. Default: 5,000,000 bp

   `-s, --cn-event-rate` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The probability that an event is an amplification. Default: 0.5

   `-a, --cn-copy-param` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Parameter used in a geometric distribution to choose the number of additional copies to add for amplification events. Default: 0.5

   `-j, --founder-event-mult` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A scalar multiplier to increase the number of focal CNAs into the founder genome. Default: 10

* Commands for simulating large scale copy number aberrations. Whole-genome duplications (WGD) and chromosomal CNAs can be included in the simulation by toggling -w and -v, respectively. Chromosomal CNAs can either be chromosome-arm level events, with probability given by -q, or whole-chromosome level events, with probability 1 minus that given by -q. Chromosomal CNAs can occur in two locations: edges into and out of the founder cell, and edges into clonal ancestors. The rate at which chromosomal CNAs occur in these locations are controlled by separate parameters in -i1 and -i2, respectively. A chromosomal CNA is either a duplication, with probability given by -u, or a deletion, with probability 1 minus the value given by -u.

   `-w, --WGD` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Include a WGD event. Default: False

   `-v, --chrom-level-event` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Include chromosomal CNAs. Default: False

   `-q, --chrom-arm-rate` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The probability that a given chromosomal CNA is a chromosome-arm event. Default: 0.75

   `-i1, --chrom-rate-founder` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The mean number of chromosomal CNAs along the edges into and out of the founder cell. Default: 2

   `-i2, --chrom-rate-clone` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The mean number of chromosomal CNAs along edges into ancestral nodes representing clonal founders. Default: 1.5

   `-u, --chrom-event-type` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The probability that a chromosomal event is a duplication. Default: 0.5

### Stage 3

* CNAsim can either generate copy number data or sequencing data. This is controlled with the *mode* parameter. Option 0 is for the sequencing mode, and option 1 is for the copy number mode.

   `-m, --mode` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Data generation mode. To generate sequencing data, use option 0. To generate copy number data, use option 1. Default: 0

* Commands generating copy number data. Regions from the starting diploid genome are grouped into contiguous fixed size bins, the size of which is given by -B. The copy number of a bin from an observed cell is the average number of copies of all regions in that bin. CNAsim also implements an error model meant to mimic the effects of noise from CNA detection methods on real sequencing data. The level of noise is controlled by two parameters, -E1 and -E2. The -E1 parameter controls the boundary model, which shortens or lengthens continuous copy number segments. The -E2 parameter controls the jitter model, which increases or decreases the value of each copy number independently. 

   `-B, --bin-length` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The bin length to use in the copy number profiles. Default: 1,000,000 bp

   `-E1, --error-rate-1` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The error rate used for adding noise at the boundaries of copy number segments. Default: 0

   `-E2, --error-rate-2` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The error rate used for adding random jitter in the copy number profiles. Default: 0

* Commands for generating sequencing read data. When generating reads for an observed cell's genome, the genome is divided into non-overlapping windows with length given by -W. CNAsim then makes a system call to *dwgsim* to generate reads for each window independently, allowing for varying read counts. The average read depth across the entire genome is given by -C, and is used to compute the expected read count of each window. The -I parameter controls at what interval a window draws an independent read depth from the coverage distribution. The read depths of windows inbetween each interval are smoothed out. Basically, smaller intervals results in more peaks and troughs, while larger intervals results in fewer peaks and troughs. The degree of read depth non-uniformity is given by a point on the lorenz curve, which can measure a wide variety of single-cell sequencing technologies. Generating whole genome read data can require a lot of computational resources, so there is the option to parallelize this step. This is recommended for larger datasets. 

   `-C, --coverage` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sequencing coverage across the entire genome. Default: 0.1

   `-W, --window-size` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The size of genomic segments to generate reads for at each iteration. Default: 100000

   `-I, --interval` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Draws from the coverage distribution every *interval* number of windows.

   `-X, --lorenz-x` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The value of the x-coordinate for the point on the lorenz curve. Default: 0.5

   `-Y, --lorenz-y` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The value of the y-coordinate for the point on the lorenz curve. Default: 0.4

   `-R, --read-length` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The length of paired-end short reads. Default: 35

   `-P, --num-processors` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The number of processors to use to generate sequencing reads in parallel. Default: 1

For a complete list of parameters, use the *--help* parameter.
```
python main.py --help
```

# A few demonstrations

## 1. Quick experiment for obtaining copy number profiles of a single chromosome for 50 cells

```
python main.py -m 1 -n 50 -N 1  
```

## 2. Mimicing whole-genome copy number profiles of 100 tumor cells with resolution and noise expected of Ginkgo/HMMcopy

```
python main.py -m 1 -n 100 -U -B 500000 -E1 0.04 -E2 0.1
```

## 3. Generating whole-genome sequencing reads of 100 tumor cells using MALBAC

```
python main.py -m 0 -n 100 -r hg38.fa -U -C 25 -X 0.5 -Y 0.27 
```

## 4. Mimicing large-scale whole-genome 10x genomics breast cancer dataset with ultra-low coverage

```
python main.py -m 0 -n 10000 -n1 0.4 -n2 0.05 -c 7 -p2 1 -r hg38.fa -U -w -v -u 0.2 -C 0.03
```


## 
