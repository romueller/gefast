# GeFaSt

Clustering tool using [Swarm](https://github.com/torognes/swarm)'s clustering strategy and [Pass-Join](http://dx.doi.org/10.1145/2487259.2487261)'s segment filter. 
It extends and generalises the fastidious clustering capabilities - hence the name **GeFaST** (**Ge**neralised **Fa**stidious **S**warming **T**ool).


## Installation
Get the source code from GitHub and compile GeFaST:

```sh
git clone https://github.com/romueller/gefast.git
cd gefast
cmake .
make
```


## Required software
 * C++ (developed with GCC 4.9.2)
 * cmake, make


## Version history

### version 0.7.0 ###
**GeFaST** v0.7.0 is the first publicly available version of **GeFaST**. 

**GeFaST** implements Swarm's clustering capabilities and also expands the fastidious clustering mechanism to arbitrary input thresholds, makes the threshold of the fastidious clustering step adjustable, and introduces an alternative edit-distance mode.
