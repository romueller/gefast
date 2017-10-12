# GeFaSt

Clustering tool using [Swarm](https://github.com/torognes/swarm)'s clustering strategy and [Pass-Join](http://dx.doi.org/10.1145/2487259.2487261)'s segment filter. 
It extends and generalises the fastidious clustering capabilities - hence the name **GeFaST** (**Ge**neralised **Fa**stidious **S**warming **T**ool).


## Installation
Get the source code from GitHub and compile GeFaST:

```sh
git clone https://github.com/romueller/gefast.git
cd gefast
make
```


## Required software
 * C++ (GCC 4.9.2 or higher)
 * make


## Version history

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
