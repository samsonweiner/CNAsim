# CNAsim

CNAsim is a software package for improved simulation of single-cell copy number alteration (CNA) data from tumors. CNAsim can be used to generate copy number profiles with noise patterns that mimic those of single-cell CNA detection algorithms, and to generate DNA-seq data for sampled cells. It offers significantly improved scalability, a high degree of customizability, and improved biological realism of simulated data. CNAsim can be cited as follows:


CNAsim: Improved simulation of single-cell copy number profiles and DNA-seq data from tumors

Samson Weiner and Mukul S. Bansal

Under review

<!--
More information can be found in our paper, located here: [link to paper].
-->

# Table of Contents

[Installation](#Installation)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [Python packages](#python-packages)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [External packages](#external-packages)

[Usage](#Usage)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [I/O and Utilities](#io-and-utilities)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [Stage 1: simulating cell lineage tree and tumor population](#stage-1)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [Stage 2: simulating genomes and mutations](#stage-2)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [Stage 3: generating single-cell data](#stage-3)

[Examples](#Example-commands)

# Installation

CNAsim is written in python and can run reliably on versions 3.7 or later. You can download the source code by cloning this repository:

```
git clone https://github.com/samsonweiner/CNAsim.git
```

### Python packages

CNAsim requires the installation of the following python packages:
* [Numpy](https://numpy.org/)
* [Scipy](https://scipy.org/)
* [msprime](https://tskit.dev/msprime/docs/stable/intro.html)
* [Biopython](https://biopython.org/)
* [Pyfaidx](https://github.com/mdshw5/pyfaidx)

Python packages can be easily installed with a package manager such as `pip` or `conda`. If using `pip`, you can install the packages by running:

```
pip install numpy scipy msprime biopython pyfaidx
```

If using `conda`, it is highly recommended that you create a fresh environment with a compatible python version before installing the packages. You can do so with the following.
```
conda create -n CNAsim python=3.7
conda activate CNAsim

conda config --env --add channels conda-forge 
conda config --env --add channels bioconda
```

Now you can install the packages by running:
```
conda install numpy scipy msprime biopython pyfaidx
```

### External packages

Additionally, CNAsim requires that the following packages are installed and configured on the environment.
* [samtools](http://www.htslib.org/download/) 
* [dwgsim](https://github.com/nh13/DWGSIM)

You may follow the installation guides with the links provided. Note that both packages can also be installed with conda. If the former option is used, the downloaded binaries must be compiled and the directories containing each binary must be added to your `$PATH` variable. For example,
```
export PATH=/path/to/msdir:$PATH
```
You may also wish to add this line to your *~/.bashrc* or */.bash_profile* configuration file to avoid having to retype this command on login. 

# Usage

CNAsim roughly follows three stages: simulate a cell lineage tree and subclonal structure, simulate genomes and evolution, and generate single-cell data. The following documentation details the relevant parameters organized according to these three stages.

To run CNAsim, call the main.py file with the desired parameters in the command line. Users must specify the simulator mode (-m, --mode) of which there are three choices. Pass **0** to generate copy number profiles (CNP mode), **1** to generate sequencing data (seq mode), or **2** to generate both.

```
python3 main.py -m mode [options]
```
Alternatively, you can modify the provided *parameters* file and toggle the *-F* option.
```
python3 main.py -F
```

### I/O and Utilities

* Commands for controlling the simulator mode and all inputs/outputs. The mode parameter -m is required and takes the integer values 0,1, or 2. All output files will be deposited in the directory located by the path given by -o. If the directory does not exist, one will be created. The remaining parameters are only circumstainstally required or entirely optional. If run in CNP mode using standard features, no inputs are required.

   `-m, --mode` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Simulator mode for generating data. 0: CNP data, 1: seq data, 2: both. Required.

   `-o, --out-path` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Path to output directory where data will be saved to. Will create one if no directory exists. Default: CNA_output/

   `-T, --tree-path` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Users may use a precomputed cell lineage tree by selecting option 2 as the tree type (see below). Specifies the path to the tree in newick format.

   `-r, --reference` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The path to the input reference genome. Required if run in seq mode.

   `-P, --num-processors` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The number of processors to use for generating sequencing reads in parallel. Does not affect the generation of CNPs. Default: 1

   `-d, --disable-info` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Toggle to omit the output simulation log, cell types, and ground truth events. Default: False

   `-F, --param-file` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Parses the provided `parameters` file in the source directory for all simulation parameters. Default: False

### Simulate a cell lineage tree and subclonal structure

* Commands for generating the cell lineage tree and subclonal populations. The user can specify any combination of tumor, normal, and pseudonormal cells using the -n, -n1, and -n2 parameters. The number of tumor cells will equal the value given by -n if both normal and pseudonormal fractions are set to 0. Otherwise, the number of tumor cells will be the remaining fraction. For example, if -n is set to 100, -n1 is set to 0.1, and -n2 is set to 0.05, there will be 85 tumor cells, 10 normal cells, and 5 pseudonormal cells. Ancestral nodes are selected to represent diverging subclonal populations, the number of which is given by -c. The probability of selecting ancestral nodes are can either be proportional to the size of the induced subtree, or the node's branch length. If the former, users can further control clone sizes by defining a normal distribution (parameters -c2 and -c3). Otherwise, nodes are selected uniformly with respect to tree depth. 

   `-t, --tree-type` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Method to generate topology of main tumor lineage. Option 0 uses *ms* to generate a tree under neutral coalescence. Option 1 generates a random tree topology. Option 2 takes an existing tree from the user in newick format (see the following parameter). Default: 0 

   `-g, --growth-rate` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The exponential growth rate parameter used to generate the tree under neutral coalescence. Default: 15.1403

   `-n, --num-cells` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Number of observed cells in the experiment. Includes all types of cells (tumor, normal, and pseudonormal). Default: 250

   `-n1, --normal-fraction` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Proportion of the observed cells that are normal diploid cells. Default: 0

   `-n2, --pseudonormal-fraction` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Proportion of observed cells that are pseudonormal cells. Default: 0

   `-c, --num-clones` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Number of ancestral nodes in the tumor lineage to be selected as clonal founders. Default: 0

   `-c1, --clone-criteria` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Criteria to choose clonal ancesters. 0: proportional to number of leaves in subtree. 1: proportional to edge length. Default: 0

   `-c2, --clone-mu` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  Mean number of leaves in a subclone. Must select 0 for clone-criteria. Default: 0
   
   `-c3, --clone-sd` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  SD in subclone size. Must select 0 for clone-criteria. Default: 0.25*mu

### Simulate genomes and evolution

* Commands for simulating genomes. The minimum copy number length is considered the minimum resolution of the genome used in the simulator. Genomes are subsequently represented as a sequence of regions with length equal to the value given by -k. If a reference genome is provided, the number of chromosomes and chromosome-lengths are derived from the reference. Without a given reference genome, pre-computed chromosome lengths derived from hg38 can be used by toggling -U. If this parameter is not toggled, the number of chromosomes will be equal to the value given by -N with set lengths given by -L. Chromosome arm ratios can either be fixed with the -A parameter, which details the ratio of the short arm, or by using the chromosome arm ratios from hg38 again by toggling -U. Even if a reference genome is given, it is still recommended to toggle -U for more accurate chromosome-arm ratios.

   `-k, --region-length` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Region length in bp. CNA events occur in region-sized units. Essentually controls the resolution of the genome. Default: 1000 bp

   `-U, --use-hg38-static` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use chromosome lengths and chromosome arm ratios derived from the hg38 reference genome. Default: False

   `-N, --num-chromosomes` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A set number of chromosomes to represent the genome. Default: 22

   `-L, --chrom-length` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A fixed chromosome length if set number of chromosomes used. Default: 100 Mbp

   `-A, --chrom-arm-ratio` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A fixed chromosome-arm ratio if arm-ratios are not derived from the hg38 reference. Default: 0.5

* Commands for simulating focal copy number aberrations. The number of focal CNAs chosen along each edge is controlled with the -p1 and -p2 parameters. If option 1 is chosen for the placement type, this means that edges directly above leaves in the tree have a mean number of events equal to the value given by -p2. Under neutral coalescence, this likely means more events will be assigned to edges higher up in the tree. The length of focal CNAs cannot drop below the minimum given by -k, and cannot exceed the length of the chromosome. The -s parameter controls the rate at which an event is chosen to be an amplifciation. If this value is given by *x*, then the rate of deletion is simply *1-x*. 

   `-p1, --placement-type` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The method for choosing the number of events along each edge. Option 0 draws from a Poisson distribution with a fixed mean. Option 1 draws from a Poisson with a mean proportional to edge length. Option 2 is a fixed number per edge. Default: 0

   `-p2, --placement-param` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The value to be used in conjunction with the placement method chosen. Default: 2

   `-l, --cn-length-mean` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The average length of a focal CNA in number of base pairs. Default: 5 Mbp

   `-l2, --min-cn-length` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Minimum copy number event length in bp. Should be at minimum the region length and no greater than the mean. Default: 1000 bp

   `-s, --cn-event-rate` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The probability that an event is an amplification. Default: 0.5

   `-a, --cn-copy-param` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Parameter used in a geometric distribution to choose the number of additional copies to add for amplification events. Default: 0.5

   `-j, --founder-event-mult` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A scalar multiplier to increase the number of focal CNAs into the founder genome. Default: 10

* Commands for simulating large scale copy number aberrations. Whole-genome duplications (WGD) and chromosomal CNAs can be included in the simulation by toggling -w and -v, respectively. Note that if the number of clones is > 0, -v will be toggled regardless. Chromosomal CNAs can either be chromosome-arm level events, with probability given by -q, or whole-chromosome level events, with probability 1 minus that given by -q. Chromosomal CNAs can occur in three locations: the edge into the founder cell, the two edges out of the founder cell, and edges into clonal ancestors. The rate at which chromosomal CNAs occur in these locations are controlled by separate parameters in -i1, -i2, and -i3, respectively. Adding chromosomal CNAs along the two edges out of the founder cell are intended to induce a super-clonal structure across the population, and should be used accordingly. A chromosomal CNA is either a duplication, with probability given by -u, or a deletion, with probability 1 minus the value given by -u.

   `-w, --WGD` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Include a WGD event. Default: False

   `-v, --chrom-level-event` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Include chromosomal CNAs. Default: False

   `-q, --chrom-arm-rate` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The probability that a given chromosomal CNA is a chromosome-arm event. Default: 0.75

   `-i1, --chrom-rate-founder` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The mean number of chromosomal CNAs along the edge into the founder cell. Default: 2

   `-i2, --chrom-rate-super-clone` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The mean number of chromosomal CNAs along the two edges out of the founder cell. Default: 1

   `-i3, --chrom-rate-clone` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The mean number of chromosomal CNAs along edges into ancestral nodes representing clonal founders. Default: 1

   `-u, --chrom-event-type` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The probability that a chromosomal event is a duplication. Default: 0.5

### Generate single-cell data

* Commands generating copy number data. Regions from the starting diploid genome are grouped into contiguous fixed size bins, the size of which is given by -B. The copy number of a bin from an observed cell is the average number of copies of all regions in that bin. CNAsim also implements an error model meant to mimic the effects of noise from CNA detection methods on real sequencing data. The level of noise is controlled by two parameters, -E1 and -E2. The -E1 parameter controls the boundary model, which shortens or lengthens continuous copy number segments. The -E2 parameter controls the jitter model, which increases or decreases the value of each copy number independently. 

   `-B, --bin-length` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The bin length to use in the copy number profiles. Default: 1,000,000 bp

   `-E1, --error-rate-1` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The error rate used for adding noise at the boundaries of copy number segments. Default: 0

   `-E2, --error-rate-2` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The error rate used for adding random jitter in the copy number profiles. Default: 0

     `-O, --output-clean-CNP` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; If run in CNP mode with the error model, outputs the clean CNPs in addition to the noisy ones. Default: False

* Commands for generating sequencing reads. For each observed cell, CNAsim builds a reference based on its evolved genome and makes system calls to *dwgsim* to generate the reads themselves. The average read depth (coverage) across the entire genome is given by -C, and is used to compute the expected read count of a given region. To use uniform coverage, toggle the -M parameter. In practice, this is an oversimplification of single-cell sequencing technologies, however doing so will reduce the time it takes to generate reads by roughly 1/4th. 

   `-C, --coverage` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sequencing coverage across the entire genome. Default: 0.1

   `-M, --use-uniform-coverage` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Assumes uniform coverage across the genome. Default: False

* Commands for controlling read counts and coverage. If non-uniform coverage is used, the genome is divided into non-overlapping windows with length given by -W. Reads are generated for each window independently thereby creating region-specific variation in read counts. The degree of coverage non-uniformity is defined by a point on the lorenz curve, which can measure a wide variety of single-cell sequencing technologies (see [this paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008012)). This is used to create a coverage distribution across the windows. The -I parameter controls at what interval a window draws an independent read depth from the coverage distribution. The read depths of windows inbetween each interval are connected via bezier curves to create smooth transitions. Basically, using smaller intervals results in more peaks and troughs, while using larger intervals results in fewer peaks and troughs. Generating whole genome sequencing reads can require a lot of computational resources, so there is the option to parallelize this step. This is recommended for larger datasets.

   `-W, --window-size` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The size of genomic segments to generate reads for at each iteration. Default: 1 Mbp

   `-I, --interval` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Draws from the coverage distribution every *interval* number of windows.

   `-X, --lorenz-x` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The value of the x-coordinate for the point on the lorenz curve. Default: 0.5

   `-Y, --lorenz-y` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The value of the y-coordinate for the point on the lorenz curve. Default: 0.4

   `-R, --read-length` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The length of paired-end short reads. Default: 150


For a complete list of parameters, use the *--help* parameter.
```
python main.py --help
```

# Example commands

## 1. Quick experiment for obtaining copy number profiles of a single chromosome for 50 cells

```
python main.py -m 0 -n 50 -N 1  
```

## 2. Whole-genome copy number profiles of 100 tumor cells and 3 subclones with resolution and noise as may result from using Ginkgo/HMMcopy

```
python main.py -m 0 -n 100 -c 3 -v -U -B 500000 -E1 0.04 -E2 0.1
```

## 3. Generating whole-genome sequencing reads of 100 tumor cells assuming read depth nonuniformity of MALBAC

```
python main.py -m 1 -n 100 -r hg38.fa -U -C 25 -X 0.5 -Y 0.27 
```

## 4. Mimicing large-scale whole-genome 10x genomics breast cancer dataset with ultra-low coverage with synthetic reads and ground truth copy number profiles leveraging 8 cores.

```
python main.py -m 2 -n 10000 -n1 0.4 -n2 0.05 -c 7 -c2 100 -c3 50 -r path/to/hg38.fa -U -w -v -u 0.2 -C 0.02 -k 10000 -B 5000000 -P 8
```
It should be noted that even with very low coverage, the above command will take a substantial amount of time (estimated at ~100 hours). When generating sequencing data for a large number of cells, it is *strongly* recommended to make use of parallelization.

## 
