# Dashing [![Build Status](https://travis-ci.com/dnbaker/dashing.svg?branch=master)](https://travis-ci.com/dnbaker/dashing)

dashing sketches and computes distances between fasta and fastq data.

# Build
Clone this repository recursively, and use make.

```bash
git clone https://github.com/dnbaker/dashing
cd dashing && make update dashing
```

Dashing is written in C++14, which means that it requires a relatively new compiler.
Dashing is tested under gcc`{5.4-9}`, but fails for gcc4, which is installed by default on many machines.
For OSX, we recommend using Homebrew to install gcc-8.
On Linux, we recommend package managers. (For instance, our Travis-CI Ubuntu example upgrades to a sufficiently new GCC using `sudo update-alternatives`.

# Usage

To see all usage options, use `./dashing <subcommand>`, for subcommand in `[sketch, dist, hll, union, printmat]`.
Of most interest is probably the dist command, which can take either genomes or pre-built sketches as arguments.

## dist
For the simplest case of unspaced, unminimized kmers for a set of genomes with `k = 31` and 13 threads:

```
dashing dist -k31 -p13 -Odistance_matrix.txt -osize_estimates.txt genome1.fna.gz genome2.fna genome3.fasta <...>
```

The genomes can be omitted as positional arguments if `-F genome_paths.txt` is provided, where `genome_paths.txt` is a file containing a path to a genome per line.
This can avoid system limits on the number of arguments in a shell command.

These can be cached with `-c`, which saves the sketches for later use. These sketch filenames are based on spacing, kmer size, and sketch size, so there is no risk of overwriting each other.

### dist (asymmetric mode)

`dashing dist` performs all pairwise jaccard index estimates by default. By providing the `-Q` flag, dashing performs a core
comparison operation between all queries and all references, where references are provided by `-F`.

This is necessary to provide containment.

For example:

```
dashing dist --containment-index -k21 -Odistmat.txt -ofsizes.txt -Q query_paths.txt -F ref_paths.txt
```

To generate a full, asymmetric distance matrix, provide the same path to -F and -Q.



## sketch
The sketch command largely mirrors dist, except that only sketches are computed.

```
dashing sketch -k31 -p13 -F genome_paths.txt
```

## hll
The hll command simply estimates the number of unique elements in a set of files. This can be useful for estimating downstream database sizes based on spacing schemes.

## union
The union command takes a set of pre-sketched HLLs and performs unions between them. Currently, the sketches must be of the same size.
We may modify this in future releases to allow a merger of different sizes by flattening larger sketches to the smallest sketch size discovered.
This would involve a loss of precision from the larger models.
This currently doesn't support data structures besides HLLs, but we plan to make this change at a later date.


## Alterative Data Structures

Dashing supports comparisons with a variety of data structures, which have speed and accuracy tradeoffs for given situations.
By default, HyperLogLog sketches are used, while b-bit minhashing, bottom-k minhashing, bloom filters, and hash sets are supported. 
Using hash sets provides a ground truth at the expense of greatly increased runtime costs.

```
b-bit minhashing:             --use-bb-minhash
bottom-k minhashing:          --use-range-minhash
weighted bottom-k minhashing: --use-counting-range-minhash
SuperMinHash:                 --use-super-minhash
hash sets:                    --use-full-khash-sets
bloom filters:                --use-bloom-filter
```

References:
[SuperMinHash](https://arxiv.org/abs/1706.05698), modified. (Use 32-bit register instead of float between 0 and 1 to make use of more information.)
[Bloom Filter Jaccard Index](https://www.ncbi.nlm.nih.gov/pubmed/17444629)
