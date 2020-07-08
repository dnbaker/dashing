# Dashing 🕺 [![Build Status](https://travis-ci.com/dnbaker/dashing.svg?branch=main)](https://travis-ci.com/dnbaker/dashing) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/dashing/README.html)

dashing sketches and computes distances between fasta and fastq data.

Our paper is available [here](https://www.biorxiv.org/content/10.1101/501726v2) as a preprint and [here](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1875-0) at Genome Biology.

# Use

The easiest way to use dashing is to grab a binary release. These are located in `dashing/release/{osx,linux}/dashing_s{128,256,512}`, where `dashing_s128`, `dashing_s256`, and `dashing_s512`
work, respectively, on systems supporting SSE2, AVX2, and AVX512BW. If these don't work, you'll need to build from source.

# Build
Clone this repository recursively, and use make.

```bash
git clone --recursive https://github.com/dnbaker/dashing
cd dashing && make dashing
```

If you clone without submodules, this can be corrected by `git submodule update --init --recursive`.

Dashing is written in C++14, which means that it requires a relatively new compiler.
Dashing is tested under gcc`{5.4-9}`, but fails for gcc4, which is installed by default on many machines.
For OSX, we recommend using Homebrew to install gcc-8.
On Linux, we recommend package managers. (For instance, our Travis-CI Ubuntu example upgrades to a sufficiently new GCC using `sudo update-alternatives`.)

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


## Features

Transparently consumes uncompressed, zlib- or zstd-compressed files.

Caching of sketches to disk (in compressed form)

Calculation of a variety of (dis)similarity measures:
1. Jaccard Similarity
2. Mash distance
3. Containment index 
4. Containment distance (log transformed containment index)
5. Symmetric Containment Index (`\frac{|A \bigcap B|}{\min{|A|,|B|}}`) (The maximum of each containment index)
6. Symmetric Continament Distance (log transformed SCI)
7. Intersection size

Additionally, supports all the above under the weighted/multiset Jaccard index via labeled w-shingling. (See Broder, 1997 "On the Resemblance and Containment of Documents" for more details.)

#### Filtering
Filtering of of rare k-mer events via count-min sketch point query estimates. This is primarily desirable for raw sequencing datasets rather than genome assemblies. This is enabled with `-y/--countmin`, and the number of hashes (`--nhashes`), sketch size (`--cm-sketch-size`) and min count `--min-count` can all be controlled by command-line parameters.

#### Encoding options
1. Exact k-mer encoding (`k <= 32`)
2. Rolling hashing encoding for any k
3. Spaced seed encoding (Hamming weight <= 32)
4. Windowed/minimized k-mers

See the [bonsai](https://github.com/dnbaker/bonsai) for more details on encoding.

#### Output formats

Dashing defaults to upper triangular TSV matrix emission, but it also suppurts upper triangular PHYLIP format, packed binary encoding, and top k-nearest-neighbor emission formats.

This is supported for all symmetric measures (Mash distance, Jaccard, Intersection size, and their multiset equivalents), whereas asymmetric measures and nearest neighbor forms (all variations of containment) have two emission options: tabular and binary.


## Alternative Data Structures

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

Wide HLL:                     --use-wide-hll
```

References:
[SuperMinHash](https://arxiv.org/abs/1706.05698), modified. (Use 32-bit register instead of float between 0 and 1 to make use of more information.)
[Bloom Filter Jaccard Index](https://www.ncbi.nlm.nih.gov/pubmed/17444629)


## To Cite:

```tex
@Article{pmid31801633,
   Author="Baker, D. N.  and Langmead, B. ",
   Title="{{D}ashing: fast and accurate genomic distances with {H}yper{L}og{L}og}",
   Journal="Genome Biol.",
   Year="2019",
   Volume="20",
   Number="1",
   Pages="265",
   Month="12"
}
```
