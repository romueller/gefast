## How to run GeFaST in different modes

This directory exemplarily describes how to run GeFaST in different modes and with different parameters on a toy data set.
For each mode, a possible scenario is described, followed by a configuration and a command suitable to perform the desired clustering.
Here, the full configuration of GeFaST is provided through the configuration file (if possible).
The described configuration files together with the exemplary outputs are provided in the respective subdirectories.

The toy data set comprises 57 sequences of varying abundance.
The data set is provided in three versions to show different input processing options: 
1) single file in FASTA format
2) single file in FASTQ format (Illumina 1.8 quality encoding)
3) two files in FASTQ format (Illumina 1.8 quality encoding)

See the manual for a full list of available configuration options.


### Levenshtein (lev) mode

*Scenario:*  
Sequences (provided in FASTA format) should be clustered based on their Levenshtein (edit) distance.
Two sequences are considered similar when they differ by at most 2 edit operations.
The clusters should be refined (using the default fastidious refiner) applying a refinement threshold of 4.
All main output files should be generated afterwards.

The corresponding configuration contains the following entries:
```sh
threshold=2
refinement_threshold=4

output_internal=lev/internal.txt
output_otus=lev/otus.txt
output_statistics=lev/statistics.txt
output_seeds=lev/seeds.fasta
```

When the sequences and the configuration are stored in `reads.fasta` and `lev/lev.conf`, respectively,
the following command performs the clustering described above:
```sh
GeFaST lev reads.fasta lev/lev.conf
```

The corresponding output files are stored in the `lev/` subdirectory.


### Alignment-score (as) mode

*Scenario:*  
Sequences (provided in FASTQ format) are to be clustered using the score of an optimal pairwise alignment as the distance between two sequences.
In order to be similar, this score must not be larger than 60.
The quality scores in the FASTQ file use the Illumina 1.8 encoding.
All main output files should be generated afterwards.

The corresponding configuration contains the following entries:
```sh
preprocessor=fastq
quality_encoding=illumina1.8
threshold=60

output_internal=as/internal.txt
output_otus=as/otus.txt
output_statistics=as/statistics.txt
output_seeds=as/seeds.fasta
```

When the sequences and the configuration are stored in `reads.fastq` and `as/as.conf`, respectively,
the following command performs the clustering described above:
```sh
GeFaST as reads.fastq as/as.conf
```

The corresponding output files are stored in the `as/` subdirectory.


### Alignment-free (af) mode

*Scenario:*  
Sequences (provided in FASTA format) should be clustered based on the word-composition vectors representing them.
The feature vectors are built using a word length of 4 and two sequences are considered similar when their Euclidean distance is at most 3.2.
All main output files should be generated afterwards.

The corresponding configuration contains the following entries:
```sh
threshold=3.2
distance=euclidean
representation=wcv
misc=wcv_word_length:4

output_internal=af/internal.txt
output_otus=af/otus.txt
output_statistics=af/statistics.txt
output_seeds=af/seeds.fasta
```

When the sequences and the configuration are stored in `reads.fasta` and `af/af.conf`, respectively,
the following command performs the clustering described above:
```sh
GeFaST af reads.fasta af/af.conf
```

The corresponding output files are stored in the `af/` subdirectory.


### Quality Levenshtein (qlev) mode

*Scenario:*  
Sequences (provided in FASTQ format) should be clustered based on the number of edit operations in an optimal pairwise alignment
using a quality-weighted scoring function (namely the Clement version).
Two sequences are considered similar when they differ by at most 3 edit operations.
The quality scores in the FASTQ file again use the Illumina 1.8 encoding.
The breaking mechanism should be deactivated and all main output files should be generated afterwards.

The corresponding configuration contains the following entries:
```sh
preprocessor=fastq
quality_encoding=illumina1.8
break_swarms=0
distance=clement
threshold=3

output_internal=qlev/internal.txt
output_otus=qlev/otus.txt
output_statistics=qlev/statistics.txt
output_seeds=qlev/seeds.fasta
```

When the sequences and the configuration are stored in `reads.fastq` and `qlev/qlev.conf`, respectively,
the following command performs the clustering described above:
```sh
GeFaST qlev reads.fastq qlev/qlev.conf
```

The corresponding output files are stored in the `qlev/` subdirectory.


### Quality alignment-score (qas) mode

*Scenario:*  
Sequences (provided in two FASTQ files) have to be clustered using the score of an optimal pairwise alignment as the distance between two sequences,
while using a quality-weighted scoring function (more precisely, the Malde-A version).
In addition, use the root boosting function with degree 3.
Two sequences are considered similar only when the alignment score is at most 30. 
The quality scores in the FASTQ file use the Illumina 1.8 encoding.
All main output files should be generated afterwards.

The corresponding configuration contains the following entries:
```sh
preprocessor=fastq
quality_encoding=illumina1.8
distance=score-malde-a
threshold=30
misc=boosting_method:root$root_degree:3

output_internal=qas/internal.txt
output_otus=qas/otus.txt
output_statistics=qas/statistics.txt
output_seeds=qas/seeds.fasta
```

When the sequences and the configuration are stored in `reads_pt1.fastq`,  `reads_pt2.fastq` and `qas/qas.conf`, respectively,
the following command performs the clustering described above:
```sh
GeFaST qas reads_pt1.fastq,reads_pt2.fastq qas/qas.conf
```

Alternatively, multiple input files can be listed in an extra file, from which the file paths are then extracted.
This is indicated by the `--list_file` option.
When the extra file `parts.txt` contains the two lines
```sh
reads_pt1.fastq
reads_pt2.fastq
```
and when the sequences and the configuration are still stored in `reads_pt1.fastq`,  `reads_pt2.fastq` and `qas/qas.conf`, respectively,
the following command performs the same clustering as the previous command:
```sh
GeFaST qas parts.txt qas/qas.conf --list_file
```

The corresponding output files are stored in the `qas/` subdirectory.


### Consistency (cons) mode

*Scenario:*  
Sequences (provided in FASTQ format) should be clustered based on the notion of consistency alone (no distance threshold needed).
Keep only sequences with an abundance of at least 2 (after dereplication). 
The quality scores in the FASTQ file once more use the Illumina 1.8 encoding.
All main output files should be generated afterwards.

The corresponding configuration contains the following entries:
```sh
preprocessor=fastq
quality_encoding=illumina1.8
min_abundance=2

output_internal=cons/internal.txt
output_otus=cons/otus.txt
output_statistics=cons/statistics.txt
output_seeds=cons/seeds.fasta
```

When the sequences and the configuration are stored in `reads.fastq` and `cons/cons.conf`, respectively,
the following command performs the clustering described above:
```sh
GeFaST cons reads.fastq cons/cons.conf
```

The corresponding output files are stored in the `cons/` subdirectory.

### Dereplication (derep) mode

*Scenario:*  
Dereplicate the sequences based on their sequence (requiring exact identity) 
and output the dereplicated sequences (with aggregated abundances).

The corresponding configuration contains the following entry:
```sh
output_seeds=derep/dereplicated.fasta
```

When the sequences and the configuration are stored in `reads.fasta` and `derep/derep.conf`, respectively,
the following command performs the clustering described above:
```sh
GeFaST derep reads.fasta derep/derep.conf
```

The corresponding output file is stored in the `derep/` subdirectory.
