# GeFaST

Clustering tool primarily designed for processing nucleotide sequences obtained from high-throughput sequencing.
**GeFaST** (**Ge**neralised **Fa**stidious **S**warming **T**ool) picks up the two-phased, iterative clustering strategy of 
[Swarm](https://github.com/torognes/swarm) and generalises resp. extends both the clustering and the refinement phase.
**GeFaST** considers the iterative strategy as an abstract scheme, which can be coupled with different notions of distance.
The refinement phase has been generalised by broadening the applicability of the original fastidious approach
and introducing further refinement methods.


## Required software
 * C++ (GCC 4.9.2 or higher)
 * [StatsLib C++ library](https://www.kthohr.com/statslib.html)
 * [GCE-Math C++ library](URLhttps://www.kthohr.com/gcem.html)


## Installation
Get the source code from GitHub and compile GeFaST:

```sh
git clone https://github.com/romueller/gefast.git
cd gefast

# install StatsLib / GCE-Math library
git clone -b master --single-branch https://github.com/kthohr/stats libs/stats
git clone https://github.com/kthohr/gcem libs/gcem
mkdir libs/gcem/build
cd libs/gcem/build
cmake .. -DCMAKE_INSTALL_PREFIX=../../stats
make install
cd ../../..

make
```


## Usage

**GeFaST** is controlled via a simple command-line interface:

```sh
GeFaST <mode> <input> <config> [option]...
```

The first three arguments are mandatory and fixed in their order:
 * ```<mode>```: The abbreviation of the mode (see below).
 * ```<input>```: By default, **GeFaST** expects a comma-separated list of file paths as its second argument.
 * ```<config>```: The file path of the configuration file containing the (basic) configuration for this execution of **GeFaST**.

These arguments are followed by an arbitrary list of options described in the next section.
When the same option is provided several times, the value of the last occurrence is used.
When an option is provided in the configuration file and via command line, 
the value specified on the command line has priority.
    
### Available modes

**GeFaST** offers the following modes, which group related clustering techniques and notions of distance or provide specific functionalities (e.g. dereplication).
The chosen mode influences the fundamental behaviour of **GeFaST** as well as the available parameters.
 * Levenshtein (```lev```): Cluster amplicons based on the number of edit operations in optimal pairwise alignments.
 * Alignment score (```as```): Cluster amplicons based on the score of optimal pairwise alignments.
 * Quality Levenshtein (```qlev```): Cluster amplicons based on the number of edit operations in optimal pairwise alignments
    considering the quality scores associated with the sequences.
 * Quality alignment score (```qas```): Cluster amplicons based on the score of optimal pairwise alignments
    considering the quality scores associated with the sequences.
 * Consistency (```cons```): Cluster amplicons using a notion of consistency considering the quality and abundance of amplicons.
 * Dereplication (```derep```): Group amplicons based on exact sequence equality.

A full description of the available options and more details on the different modes are provided in the manual. 
Moreover, examples of how to run GeFaST in the different modes are shown in the `examples/` subdirectory.

## Citation

To cite **GeFaST**, please refer to:

Müller, R., & Nebel, M. (2018). GeFaST: An improved method for OTU assignment by generalising Swarm’s fastidious clustering approach. 
		*BMC Bioinformatics*, 19(1), 321. doi:[10.1186/s12859-018-2349-1](https://doi.org/10.1186/s12859-018-2349-1)

## Version history

### version 2.0.1 ###
**GeFaST** v2.0.1 provides additional options for how to combine quality scores while dereplicating FASTQ sequences 
and changes the default scaling factor of the quality-weighted Frith cost function to a value that corresponds to the default scoring function.
It also includes some smaller performance improvements and adds examples of how to run **GeFaST** in the different modes.


### version 2.0.0 ###
**GeFaST** v2.0.0 further generalises the iterative clustering strategy by considering additional notions of distance
and quality-aware methods for clustering and refinement.
The previously considered quality-unaware, alignment-based notions of distance are still accessible in Levenshtein (and dereplication) mode.
The structure of **GeFaST** has been extensively revised to incorporate the new (and future) clustering and refinement methods.


### version 1.0.0 ###
**GeFaST** v1.0.0 improves the runtime and memory consumption through a range of optimisations and
contains some modifications that reflect changes to Swarm since version 2.1.13 like an upper-case normalisation,
changed amplicon sorting criteria and modified outputs.
It also fixes small bugs that could lead to slightly wrong results on undereplicated data.
Moreover, the dereplication procedure, the on-screen logging and the document are improved.


### version 0.9.0 ###
**GeFaST** v0.9.0 incorporates several small optimisations and adds (a first version of) alternative code for the clustering steps that makes use of succinct data structures.


### version 0.8.0 ###
**GeFaST** v0.8.0 reduces the memory footprint of **GeFaST**.
It also adds the possibility to use an additional q-gram filter similar to Swarm.


### version 0.7.0 ###
**GeFaST** v0.7.0 is the first publicly available version of **GeFaST**.

**GeFaST** implements Swarm's clustering capabilities and also expands the fastidious clustering mechanism to arbitrary input thresholds,
makes the threshold of the fastidious clustering step adjustable, and introduces an alternative edit-distance mode.
