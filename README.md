# CNAsim

CNAsim is a software package for improved simulation of single-cell copy number alteration (CNA) data from tumors. CNAsim can be used to generate copy number profiles with noise patterns that mimic those of single-cell CNA detection algorithms, readcount data with ground truth CNPs, and DNA sequencing reads for sampled cells with ground truth CNPs. It offers significantly improved scalability, a high degree of customizability, and improved biological realism of simulated data.

CNAsim can be cited as follows:

<a href="https://doi.org/10.1093/bioinformatics/btad434">CNAsim: Improved simulation of single-cell copy number profiles and DNA-seq data from tumors</a><br>
Samson Weiner and Mukul S. Bansal<br>
Bioinformatics, Volume 39, Issue 7, July 2023, btad434.

Current Version: 1.3.1

### New Features

* Sequencing reads can be sampled from two distinct haploid reference genomes, enabling the analysis of allelic frequency. To make sure of this feature, run CNAsim in mode **2** and specify the paths to the reference sequences using the `-r1` and `-r2` parameters.
* CNAsim can now be used to simulate readcounts directly without having to generate synthetic sequencing reads. This can be achieved by using mode **1**. Passing a reference genome is optional.

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

[Using the error model](#Using-the-error-model)

# Installation

CNAsim is written in python and runs reliably on versions 3.7 - 3.10. For Mac and Linux users, a standalone executable of CNAsim is available and can be downloaded from [here](https://compbio.engr.uconn.edu/software/CNAsim/) (OUTDATED). 

Otherwise, users with a python interpreter can use CNAsim by downloading the source code and setting up an environment with all the necessary packages. The remainder of this section provides instructions on how to install CNAsim this way.

You can download the source code by cloning this repository:

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
conda create -n CNAsim
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

Both packages may be installed with `conda`.
```
conda install -c bioconda samtools dwgsim
```
You can also install the packages individually by following the instructions found on the package homepage. If this option is used, the downloaded binaries must be compiled and the directories containing each binary must be added to your `$PATH` variable. For example,
```
export PATH=/path/to/bin:$PATH
```
You may also wish to add this line to your *~/.bashrc* or */.bash_profile* configuration file to avoid having to retype this command on login. 

# Usage

CNAsim roughly follows three stages: simulate a cell lineage tree and subclonal structure, simulate genomes and evolution, and generate single-cell data. The following documentation details the relevant parameters organized according to these three stages.

To run CNAsim, run the executable with the desired parameters in the command line. Users must specify the simulator mode (-m, --mode) of which there are three choices. Pass **0** to generate copy number profiles only (CNP mode), **1** to generate readcounts and CNPs (count mode), or **2** to generate synthetic DNA sequencing reads and CNPs (seq mode).

```
python main.py -m mode [options]
```

Users running CNAsim using an executable should instead run:

```
./CNAsim -m mode [options]
```

Instead of inputing each parameter in the command line, you can modify the provided *parameters* file and toggle the *-F* option.
```
python main.py -F
./CNAsim -F
```

### A note on the input reference
If using CNAsim to generate sequence reads, please ensure that the provided reference genome(s) contain only the desired chromosomes (remove any small alternate/unorganized sequences that you do not want to be represented as chromosomes).

For convenience, if using human reference genome hg38 with chromosomes named chr1, chr2, ect..., simply toggle the --use-hg38-static parameter to bypass having to filter the fasta file.

### I/O and Utilities

* Commands for controlling the simulator mode and all inputs/outputs. The mode parameter -m is required and takes the integer values 0,1, or 2. All output files will be deposited in the directory located by the path given by -o. If the directory does not exist, one will be created. The remaining parameters are only circumstainstally required or entirely optional. If run in CNP mode using standard features, no inputs are required.

   `-m, --mode` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Simulator mode for generating data. 0: CNP data, 1: count data, 2: seq data. Required.

   `-o, --out-path` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Path to output directory where data will be saved to. Will create one if no directory exists. Default: CNA_output/

   `-T, --tree-path` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Users may use a precomputed cell lineage tree by selecting option 2 as the tree type (see below). Specifies the path to the tree in newick format.

   `-r1, --reference` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Path to input reference genome as the primary haplotype in fasta format. Will be duplicated as both haplotypes if an alternate is not provided. Required if run in seq mode.

   `-r2, --alt-reference` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Path to an alternate reference genome to be used as a secondary haplotype, also in fasta format.

   `-P, --num-processors` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The number of processors to use for generating sequencing reads in parallel. Does not affect the generation of CNPs. Default: 1

   `-d, --disable-info` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Toggle to omit the output simulation log, cell types, and ground truth events. Default: False

   `-F, --param-file` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Parses the provided `parameters` file in the source directory for all simulation parameters. Default: False

### Simulate a cell lineage tree and subclonal structure

* Commands for generating the cell lineage tree and subclonal populations. The user can specify any combination of tumor, normal, and pseudonormal cells using the -n, -n1, and -n2 parameters. The number of tumor cells will equal the value given by -n if both normal and pseudonormal fractions are set to 0. Otherwise, the number of tumor cells will be the remaining fraction. For example, if -n is set to 100, -n1 is set to 0.1, and -n2 is set to 0.05, there will be 85 tumor cells, 10 normal cells, and 5 pseudonormal cells. To generate the main tumor lineage, CNAsim simulates an exponentially growing population under neutral coalescence with growth rate given by the -g parameter. Selective sweeps can be introduced at random points during the tumor lifespan using the -s parameter. Once the cell-lineage tree is generated, ancestral nodes can be selected to represent diverging subclonal populations, the number of which is given by -c. The probability of selecting ancestral nodes are can either be proportional to the size of the induced subtree, or the node's branch length. If the former, users can further control clone sizes by defining a normal distribution (parameters -c2 and -c3) over the number of cells belonging to the clone. Otherwise, nodes are selected uniformly with respect to tree depth. 

   `-t, --tree-type` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Method to generate topology of main tumor lineage. Option 0 uses *ms* to generate a tree under neutral coalescence. Option 1 generates a random tree topology. Option 2 takes an existing tree from the user in newick format (see the following parameter). Default: 0 

   `-g, --growth-rate` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The exponential growth rate parameter used to generate the tree under neutral coalescence. Default: 0.003785

   `-s, --num-sweep` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Number of selective sweeps to introduce throughout the tumor lineage. Default: 0

   `-s1, --selective-strength` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Relative selective strength of the sweeps. Recommended values are 0.01 for weak sweeps, and 0.25 for strong sweeps. For more information, see the msprime documentation. Default: 0.01

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

   `-U, --use-hg38-static` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use chromosome lengths and chromosome arm ratios derived from the hg38 reference genome, excluding sex chromosomes. Default: False

   `-N, --num-chromosomes` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A set number of chromosomes to represent the genome. Default: 22

   `-L, --chrom-length` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A fixed chromosome length if set number of chromosomes used. Default: 100 Mbp

   `-A, --chrom-arm-ratio` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A fixed chromosome-arm ratio if arm-ratios are not derived from the hg38 reference. Default: 0.5

* Commands for simulating focal copy number aberrations. The number of focal CNAs chosen along each edge is controlled with the -p1 and -p2 parameters. If option 1 is chosen for the placement type, this means that edges directly above leaves in the tree have a mean number of events equal to the value given by -p2. Under neutral coalescence, this likely means more events will be assigned to edges higher up in the tree. The length of focal CNAs cannot drop below the minimum given by -k, and cannot exceed the length of the chromosome. The -s parameter controls the rate at which an event is chosen to be an amplifciation. If this value is given by *x*, then the rate of deletion is simply *1-x*. 

   `-p, --placement-type` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The method for choosing the number of events along each edge. Option 0 draws from a Poisson distribution with a fixed mean. Option 1 draws from a Poisson with a mean proportional to edge length. Option 2 is a fixed number per edge. Default: 0

   `-p1, --placement-param` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The value to be used in conjunction with the placement method chosen. Default: 2

   `-l, --cn-length-mean` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The average length of a focal CNA in number of base pairs. Default: 5 Mbp

   `-l1, --min-cn-length` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Minimum copy number event length in bp. Should be at minimum the region length and no greater than the mean. Default: 1000 bp

   `-b, --cn-event-rate` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The probability that an event is an amplification. Default: 0.5

   `-a, --cn-copy-param` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Parameter used in a geometric distribution to choose the number of additional copies to add for amplification events. Default: 0.5

   `-j, --founder-event-mult` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A scalar multiplier to increase the number of focal CNAs into the founder genome. Default: 10

* Commands for simulating large scale copy number aberrations. Whole-genome duplications (WGD) and chromosomal CNAs can be included in the simulation by toggling -w and -v, respectively. Note that if the number of clones is > 0, -v will be toggled regardless. Chromosomal CNAs can either be chromosome-arm level events, with probability given by -q, or whole-chromosome level events, with probability 1 minus that given by -q. Chromosomal CNAs can occur in three locations: the edge into the founder cell, the two edges out of the founder cell, and edges into clonal ancestors. The rate at which chromosomal CNAs occur in these locations are controlled by separate parameters in -i1, -i2, and -i3, respectively. Adding chromosomal CNAs along the two edges out of the founder cell are intended to induce a super-clonal structure across the population, and should be used accordingly. A chromosomal CNA is either a duplication, with probability given by -u, or a deletion, with probability 1 minus the value given by -u.

   `-w, --WGD` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Include a WGD event. Default: False

   `-v, --chrom-level-event` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Include chromosomal CNAs. Default: False

   `-q, --chrom-arm-rate` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The probability that a given chromosomal CNA is a chromosome-arm event. Default: 0.75

   `-i, --chrom-rate-founder` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The mean number of chromosomal CNAs along the edge into the founder cell. Default: 2

   `-i1, --chrom-rate-super-clone` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The mean number of chromosomal CNAs along the two edges out of the founder cell. Default: 1

   `-i2, --chrom-rate-clone` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The mean number of chromosomal CNAs along edges into ancestral nodes representing clonal founders. Default: 1

   `-u, --chrom-event-type` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The probability that a chromosomal event is a duplication. Default: 0.5

### Generate single-cell data

* Commands generating copy number data. Regions from the starting diploid genome are grouped into contiguous fixed size bins, the size of which is given by -B. The copy number of a bin from an observed cell is the average number of copies of all regions in that bin. CNAsim also implements an error model meant to mimic the effects of noise from CNA detection methods on real sequencing data. The level of noise is controlled by two parameters, -E1 and -E2. The -E1 parameter controls the boundary model, which shortens or lengthens continuous copy number segments. The -E2 parameter controls the jitter model, which increases or decreases the value of each copy number independently. 

   `-B, --bin-length` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The bin length to use in the copy number profiles. Default: 1,000,000 bp

   `-E1, --error-rate-1` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The error rate used for adding noise at the boundaries of copy number segments. Default: 0

   `-E2, --error-rate-2` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The error rate used for adding random jitter in the copy number profiles. Default: 0

* Commands for generating sequencing reads. For each observed cell, CNAsim builds a reference based on its evolved genome and makes system calls to *dwgsim* to generate the reads themselves. The average read depth (coverage) across the entire genome is given by -C, and is used to compute the expected read count of a given region. To use uniform coverage, toggle the -M parameter. In practice, this is an oversimplification of single-cell sequencing technologies, however doing so will reduce the time it takes to generate reads by roughly 1/4th. 

   `-C, --coverage` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sequencing coverage across the entire genome. Default: 0.1

   `-M, --use-uniform-coverage` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Assumes uniform coverage across the genome. Default: False

   `-R, --read-length` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The length of paired-end short reads. Default: 150

   `-S, --seq-error` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The rate of sequencing errors. Default: 0.02

* Commands for controlling coverage non-uniformity. If non-uniform coverage is used, the genome is divided into non-overlapping windows with length given by -W. Reads are generated for each window independently thereby creating region-specific variation in read counts. The degree of coverage non-uniformity is defined by a point on the lorenz curve, which can measure a wide variety of single-cell sequencing technologies (see [this paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008012)). This is used to create a coverage distribution across the windows. The -I parameter controls at what interval a window draws an independent read depth from the coverage distribution. The read depths of windows inbetween each interval are connected via bezier curves to create smooth transitions. Basically, using smaller intervals results in more peaks and troughs, while using larger intervals results in fewer peaks and troughs. Generating whole genome sequencing reads can require a lot of computational resources, so there is the option to parallelize this step. This is recommended for larger datasets.

   `-W, --window-size` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The size of genomic segments to generate reads for at each iteration. Default: 1 Mbp

   `-I, --interval` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Draws from the coverage distribution every *interval* number of windows.

   `-X, --lorenz-x` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The value of the x-coordinate for the point on the lorenz curve. Default: 0.5

   `-Y, --lorenz-y` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The value of the y-coordinate for the point on the lorenz curve. Default: 0.4



For a complete list of parameters, use the *--help* parameter.
```
./CNAsim --help
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

## 3. Generating whole-genome sequencing reads of 100 tumor cells assuming highly nonuniform coverage

```
python main.py -m 2 -n 100 -r1 hg38.fa -U -C 25 -X 0.5 -Y 0.4 
```

## 4. Mimicing large-scale whole-genome 10x genomics breast cancer dataset with ultra-low coverage with synthetic reads and ground truth copy number profiles leveraging 8 cores.

```
python main.py -m 2 -n 10000 -n1 0.4 -n2 0.05 -c 7 -c2 100 -c3 50 -r1 path/to/hg38.fa -U -w -v -u 0.2 -C 0.02 -k 10000 -B 5000000 -P 8
```
It should be noted that even with very low coverage, the above command will take a substantial amount of time (estimated at ~100 hours). When generating sequencing data for a large number of cells, it is *strongly* recommended to make use of parallelization.

## 5. Generate whole-genome read counts of 100 tumor cells directly without generating any reads

```
python main.py -m 1 -n 100 -U -C 0.1 -M -W 1000000 -B 5000000
```

# Using the error model

CNAsim can bypass the explicit generation of reads and simulate copy number profiles directly with error patterns of existing CNA detection algorithms. Assuming the precision and recall values of a particular variant caller are known, the same precision and recall values can be achieved in the noisy CNPs outputted by CNAsim with the correct combination of boundary error rate (-E1) and jitter error rate (-E2).

Included in the source code is the script `eval_error.py` which can be used to obtain the precision/recall of a particular set of simulation parameters. To use the script, run the simulator with the -O flag to output the clean CNPs in addition to the noisy ones. Then, use the script as follows:
```
python eval_error.py -p /path/to/data/directory -t tolerance
```
Here, *tolerance* represents the distance from the predicted breakpoint to the groundtruth breakpoint in terms of number of bins. If the predicted breakpoint is within the distance threshold, it is counted as correct. By default, the tolerance is set to 0, meaning the predicted breakpoint must match up directly with the ground truth.

Finding the correct combination of boundary and jitter error rates to yield a particular precision/recall can be tedious. The following table can be used as a starting point. In general, increasing the jitter error decreases precision, and increasing the boundary error decreases recall. Note that these precision/recall rates were derived using CNAsim with default parameter values and 100 cells.

## Sample error rate combinations
| E1   | E2   | Precision | Recall |         | E1   | E2   | Precision | Recall |
|------|------|-----------|--------|---------|------|------|-----------|--------|
| 0.01 | 0.05 |   0.963   |  0.956 |         | 0.02 | 0.05 |   0.905   |  0.860 |
|      |  0.1 |   0.850   |  0.940 |         |      |  0.1 |   0.774   |  0.856 |
|      | 0.15 |   0.616   |  0.934 |         |      | 0.15 |   0.619   |  0.858 |
| 0.04 | 0.05 |   0.841   |  0.762 |         | 0.06 | 0.05 |   0.830   |  0.727 |
|      |  0.1 |   0.738   |  0.784 |         |      |  0.1 |   0.687   |  0.670 |
|      | 0.15 |   0.527   |  0.757 |         |      | 0.15 |   0.568   |  0.725 |
| 0.08 | 0.05 |   0.775   |  0.662 |         |  0.1 | 0.05 |   0.744   |  0.598 |
|      |  0.1 |   0.643   |  0.654 |         |      |  0.1 |   0.645   |  0.619 |
|      | 0.15 |   0.512   |  0.662 |         |      | 0.15 |   0.478   |  0.624 |
| 0.12 | 0.05 |   0.721   |  0.580 |         | 0.14 | 0.05 |   0.709   |  0.550 |
|      |  0.1 |   0.604   |  0.587 |         |      |  0.1 |   0.624   |  0.577 |
|      | 0.15 |   0.437   |  0.583 |         |      | 0.15 |   0.464   |  0.568 |

In practice, the precision and recall of an existing method may not be known. While less efficient than generating CNPs directly, CNAsim can be used to evaluate the precision/recall of existing CNA detection algorithms using sequencing data. First,run CNAsim under mode set to 2 in order to output both ground truth copy number profiles and sequencing reads. Second, apply the CNA detection algorithm to the sequencing reads to obtain 'noisy' copy number profiles. Third, apply the `eval_error.py` script as described above. This will yield the precision and recall of the algorithm. The benefit of doing this process is that for any subsequent experiments that require CNP data, the user can generate the CNPs directly with error patterns representative of the existing CNA detection algorithm.
