#include "sketch_and_cmp.h"

using namespace sketch;
using hll::hll_t;


namespace bns {
GlobalArgs gargs;
uint64_t global_hash_seed = 137;
// sketch_core forward declaration
extern template void sketch_core<mh::RangeMinHash<uint64_t>>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);
extern template void sketch_core<mh::CountingRangeMinHash<uint64_t>>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);
extern template void sketch_core<SuperMinHashType>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);
extern template void sketch_core<hll::hll_t>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);
extern template void sketch_core<HLLH>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);
extern template void sketch_core< bf::bf_t>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);
extern template void sketch_core<khset64_t>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);
extern template void sketch_core<mh::BBitMinHasher<uint64_t>>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);

extern template void sketch_core<wj::WeightedSketcher<mh::RangeMinHash<uint64_t>>>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);
extern template void sketch_core<wj::WeightedSketcher<mh::CountingRangeMinHash<uint64_t>>>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);
extern template void sketch_core<wj::WeightedSketcher<SuperMinHashType>>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);
extern template void sketch_core<wj::WeightedSketcher<hll::hll_t>>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);
extern template void sketch_core<wj::WeightedSketcher<HLLH>>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);
extern template void sketch_core<wj::WeightedSketcher<bf::bf_t>>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);
extern template void sketch_core<wj::WeightedSketcher<khset64_t>>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);
extern template void sketch_core<wj::WeightedSketcher<mh::BBitMinHasher<uint64_t>, wj::ExactCountingAdapter>>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);

extern template void sketch_core<wj::WeightedSketcher<mh::RangeMinHash<uint64_t>, wj::ExactCountingAdapter>>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);
extern template void sketch_core<wj::WeightedSketcher<mh::CountingRangeMinHash<uint64_t>, wj::ExactCountingAdapter>>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);
extern template void sketch_core<wj::WeightedSketcher<SuperMinHashType>>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);
extern template void sketch_core<wj::WeightedSketcher<hll::hll_t, wj::ExactCountingAdapter>>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);
extern template void sketch_core<wj::WeightedSketcher<HLLH, wj::ExactCountingAdapter>>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);
extern template void sketch_core<wj::WeightedSketcher<bf::bf_t, wj::ExactCountingAdapter>>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);
extern template void sketch_core<wj::WeightedSketcher<khset64_t, wj::ExactCountingAdapter>>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);
extern template void sketch_core<wj::WeightedSketcher<mh::BBitMinHasher<uint64_t>, wj::ExactCountingAdapter>>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct, std::string);
// sketch_by_seq_core forward declaration
#if 0
extern template void sketch_by_seq_core<mh::RangeMinHash<uint64_t>>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct);
extern template void sketch_by_seq_core<mh::CountingRangeMinHash<uint64_t>>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct);
extern template void sketch_by_seq_core<SuperMinHashType>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct);
extern template void sketch_by_seq_core<hll::hll_t>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct);
extern template void sketch_by_seq_core< bf::bf_t>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct);
extern template void sketch_by_seq_core<khset64_t>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct);
extern template void sketch_by_seq_core<mh::BBitMinHasher<uint64_t>>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, int sketch_flags, uint32_t mincount, EncodingType enct);
#endif


void main_usage(char **argv) {
    std::fprintf(stderr, "Usage: %s <subcommand> [options...]. Use %s <subcommand> for more options. [Subcommands: sketch, cmp, hll, mkdist, union, view, flatten, printmat.]\[cmp was formerly dist]\n",
                 *argv, *argv);
    std::exit(EXIT_FAILURE);
}

size_t posix_fsize(const char *path) {
    struct stat st;
    stat(path, &st);
    return st.st_size;
}

void dist_usage(const char *arg) {
    std::fprintf(stderr, "Usage: %s <opts> [genome1 genome2 seq.fq [...] if not provided from a file with -F]\n"
                         "Flags:\n"
                         "-h/-?, --help\tUsage\n\n\n"
                         "===Encoding Options===\n\n"
                         "-k, --kmer-length\tSet kmer size [31], max 32\n"
                         "-s, --spacing\tadd a spacer of the format <int>x<int>,<int>x<int>,"
                         "..., where the first integer corresponds to the space "
                         "-w, --window-size\tSet window size [max(size of spaced kmer, [parameter])]\n"
                         "-S, --sketch-size\tSet sketch size [10, for 2**10 bytes each]\n"
                         "--use-nthash\tUse nthash for encoding. (not reversible, but fast, rolling, and specialized for DNA). Allows k to be unbounded\n"
                         "--use-cyclic-hash\tUses a cyclic hash for encoding. Not reversible, but efficient. Allows k to be unbounded\n"
                         "-C, --no-canon\tDo not canonicalize. [Default: canonicalize]\n\n\n"
                         "===Output Files===\n\n"
                         "-o, --out-sizes\tOutput for genome size estimates [stdout]\n"
                         "-O, --out-dists\tOutput for genome distance matrix [stdout]\n\n\n"
                         "===Filtering Options===\n\n"
                         "-y, --countmin\tFilter all input data by count-min sketch.\n"
                         "--sketch-by-fname\tAutodetect fastq or fasta data by filename (.fq or .fastq within filename).\n"
                         " When filtering with count-min sketches by either -y or -N, set minimum count:"
                         "-c, --min-count\tSet minimum count for kmers to pass count-min filtering.\n"
                         "-q, --nhashes\tSet count-min number of hashes. Default: [1]\n"
                         "-t, --cm-sketch-size\tSet count-min sketch size (log2). Default: 20\n"
                         "-R, --seed\tSet seed for seeds for count-min sketches\n\n\n"
                         "===Runtime Options\n\n"
                         "-F, --paths\tGet paths to genomes from file rather than positional arguments\n"
                         "-W, --cache-sketches\tCache sketches/use cached sketches\n"
                         "-p, --nthreads\tSet number of threads [1]\n"
                         "-Q, --query-paths\tSets query paths to use for performing query against reference paths\n"
                         "                 \tParticularly for asymmetric distances or panel queries.\n"
                         "For use with option -F.\n"
                         "--presketched\tTreat provided paths as pre-made sketches.\n"
                         "-P, --prefix\tSet prefix for sketch file locations [empty]\n"
                         "-x, --suffix\tSet suffix in sketch file names [empty]\n"
                         "--avoid-sorting\tAvoid sorting files by genome sizes. This avoids a computational step, but can result in degraded load-balancing.\n\n\n"
                         "===Emission Formats===\n\n"
                         "-b, --emit-binary\tEmit distances in binary (default: human-readable, upper-triangular)\n"
                         "-U, --phylip\tEmit distances in PHYLIP upper triangular format(default: human-readable, upper-triangular)\n"
                         "between bases repeated the second integer number of times\n"
                         "-T, --full-tsv\tpostprocess binary format to human-readable TSV (not upper triangular) [Square matrix]\n\n\n"
                         "===Emission Details===\n\n"
                         "-e, --emit-scientific\tEmit in scientific notation\n\n\n"
                         "===Data Structures===\n\n"
                         "Default: HyperLogLog. Alternatives:\n"
                         "--use-hyperminhash\tUse HyperMinHash. Defaults to 16-bit registers. Use --bbits/-B bits to modify remainder size to 8, 16, 32, or 64.\n"
                         "--use-bb-minhash/-8\tCreate b-bit minhash sketches\n"
                         "--use-bloom-filter\tCreate bloom filter sketches\n"
                         "--use-range-minhash\tCreate range minhash sketches\n"
                         "--use-super-minhash\tCreate b-bit superminhash sketches\n"
                         "--use-counting-range-minhash\tCreate range minhash sketches\n"
                         "--use-full-khash-sets\tUse full khash sets for comparisons, rather than sketches. This can take a lot of memory and time!\n\n\n"
                         "===Sketch-specific Options===\n\n"
                         "-I, --improved      \tUse Ertl's Improved Estimator for HLL\n"
                         "-E, --original      \tUse Ertl's Original Estimator for HLL\n"
                         "-J, --ertl-joint-mle\tUse Ertl's JMLE Estimator for HLL[default:Uses Ertl-MLE]\n\n\n"
                         "===b-bit Minhashing Options (apply for b-bit minhash and b-bit superminhash) ===\n\n"
                         "--bbits,-B\tSet `b` for b-bit minwise hashing to <int>. Default: 16. For HyperMinHash, this sets the full register size.\n\n\n"
                         "===Distance Emission Types===\n\n"
                         "Default: Jaccard Index\n"
                         "Alternatives:\n"
                         "-M, --mash-dist    \tEmit Mash distance [ji ? (-log(2. * ji / (1. + ji)) / k) : 1.]\n"
                         "--full-mash-dist   \tEmit full (not approximate) Mash distance. [1. - (2.*ji/(1. + ji))^(1/k)]\n"
                         "--sizes            \tEmit intersection sizes (default: jaccard index)\n"
                         "--containment-index\tEmit Containment Index (|A & B| / |A|)\n"
                         "--containment-dist \tEmit distance metric using containment index. [Let C = (|A & B| / |A|). C ? -log(C) / k : 1.] \n"
                         "--symmetric-containment-dist\tEmit symmetric containment index symcon(A, B) = max(C(A, B), C(B, A))\n"
                         "--symmetric-containment-index\ttEmit distance metric using maximum containment index. symdist(A, B) = min(cdist(A,B), cdist(B, A))\n"
                         "--full-containment-dist \tEmit distance metric using containment index, without log approximation. [Let C = (|A & B| / |A|). C ? 1. - C^(1/k) : 1.] \n"
                         "\n\n"
                         "===Streaming Weighted Jaccard===\n"
                         "--wj               \tEnable weighted jaccard adapter using the count-min sketch\n"
                         "--wj-cm-sketch-size\tSet count-min sketch size for count-min streaming weighted jaccard [16]\n"
                         "--wj-cm-nhashes    \tSet count-min sketch number of hashes for count-min streaming weighted jaccard [8]\n"
                         "===Miscellaneous===\n"
                         "--defer-hll        \tMaintain k-partition MinHash and produce an HLL at the end. May be faster (fewer instructions) or slower (more memory).\n"
                , arg);
    std::exit(EXIT_FAILURE);
}


// Usage, utilities
void sketch_usage(const char *arg) {
    std::fprintf(stderr, "Usage: %s <opts> [genomes if not provided from a file with -F]\n"
                         "Flags:\n"
                         "-h/-?:\tEmit usage\n"
                         "\n\n"
                         "Sketch options --\n\n"
                         "--kmer-length/-k\tSet kmer size [31], max 32\n"
                         "--spacing/-s\tadd a spacer of the format <int>x<int>,<int>x<int>,"
                         "..., where the first integer corresponds to the space "
                         "between bases repeated the second integer number of times\n"
                         "--window-size/-w\tSet window size [max(size of spaced kmer, [parameter])]\n"
                         "--sketch-size/-S\tSet log2 sketch size in bytes [10, for 2**10 bytes each]\n"
                         "--no-canon/-C\tDo not canonicalize. [Default: canonicalize]\n"
                         "--bbits/-B\tSet `b` for b-bit minwise hashing to <int>. Default: 16\n\n\n"
                         "Run options --\n\n"
                         "--nthreads/-p\tSet number of threads [1]\n"
                         "--prefix/-P\tSet prefix for sketch file locations [empty]\n"
                         "--suffix/-x\tSet suffix in sketch file names [empty]\n"
                         "--paths/-F\tGet paths to genomes from file rather than positional arguments\n"
                         "--skip-cached/-c\tSkip alreday produced/cached sketches (save sketches to disk in directory of the file [default] or in folder specified by -P\n"
                         "--avoid-sorting\tAvoid sorting files by genome sizes. This avoids a computational step, but can result in degraded load-balancing.\n\n\n"
                         "\n\n"
                         "Estimation methods --\n\n"
                         "--original/-E\tUse Flajolet with inclusion/exclusion quantitation method for hll. [Default: Ertl MLE]\n"
                         "--improved/-I\tUse Ertl Improved estimator [Default: Ertl MLE]\n"
                         "--ertl-jmle/-J\tUse Ertl JMLE\n\n\n"
                         "Filtering Options --\n\n"
                         "Default: consume all kmers. Alternate options: \n"
                         "--sketch-by-fname\tAutodetect fastq or fasta data by filename (.fq or .fastq within filename).\n"
                         "--countmin/-b\tFilter all input data by count-min sketch.\n\n\n"
                         "Options for count-min filtering --\n\n"
                         "--nhashes/-H\tSet count-min number of hashes. Default: [1]\n"
                         "--cm-sketch-size/-q\tSet count-min sketch size (log2). Default: 20\n"
                         "--min-count/-n\tProvide minimum expected count for fastq data. If unspecified, all kmers are passed.\n"
                         "--seed/-R\tSet seed for seeds for count-min sketches\n\n\n"
                         "Sketch Type Options --\n\n"
                         "--use-hyperminhash\tUse HyperMinHash. Defaults to 16-bit registers. Use --bbits/-B bits to modify remainder size to 8, 16, 32, or 64.\n"
                         "--use-bb-minhash/-8\tCreate b-bit minhash sketches\n"
                         "--use-bloom-filter\tCreate bloom filter sketches\n"
                         "--use-range-minhash\tCreate range minhash sketches\n"
                         "--use-super-minhash\tCreate b-bit super minhash sketches\n"
                         "--use-counting-range-minhash\tCreate range minhash sketches\n"
                         "--use-full-khash-sets\tUse full khash sets for comparisons, rather than sketches. This can take a lot of memory and time!\n"
                         "\n\n"
                         "===Streaming Weighted Jaccard===\n"
                         "--wj               \tEnable weighted jaccard adapter using the count-min sketch\n"
                         "--wj-cm-sketch-size\tSet count-min sketch size for count-min streaming weighted jaccard [16]\n"
                         "--wj-cm-nhashes    \tSet count-min sketch number of hashes for count-min streaming weighted jaccard [8]\n"
                         "--wj-exact         \tEnable exact weighted jaccard using a hash map\n"
                         "===Miscellaneous===\n"
                         "--defer-hll        \tMaintain k-partition MinHash and produce an HLL at the end. May be faster (fewer instructions) or slower (more memory).\n"
                , arg);
    std::exit(EXIT_FAILURE);
}


void sketch_by_seq_usage(const char *arg) {
    std::fprintf(stderr, "Usage: %s <opts> [genomes if not provided from a file with -F]\n"
                         "Creates a set of sketches, one per sequence, for a given file.\n"
                         "\n"
                         "Flags:\n"
                         "-h/-?:\tEmit usage\n"
                         "\n\n"
                         "Output options --\n\n"
                         "-o\tWrite sketches to [file] and names to [file].names rather than /dev/stdout and stdout.name\n"
                         "Sketch options --\n\n"
                         "--kmer-length/-k\tSet kmer size [31], max 32\n"
                         "--spacing/-s\tadd a spacer of the format <int>x<int>,<int>x<int>,"
                         "..., where the first integer corresponds to the space "
                         "between bases repeated the second integer number of times\n"
                         "--window-size/-w\tSet window size [max(size of spaced kmer, [parameter])]\n"
                         "--sketch-size/-S\tSet log2 sketch size in bytes [10, for 2**10 bytes each]\n"
                         "--no-canon/-C\tDo not canonicalize. [Default: canonicalize]\n"
                         "--bbits/-B\tSet `b` for b-bit minwise hashing to <int>. Default: 16\n\n\n"
                         "Run options --\n\n"
                         "--nthreads/-p\tSet number of threads [1]\n"
                         "\n\n"
                         "Estimation methods --\n\n"
                         "--original/-E\tUse Flajolet with inclusion/exclusion quantitation method for hll. [Default: Ertl MLE]\n"
                         "--improved/-I\tUse Ertl Improved estimator [Default: Ertl MLE]\n"
                         "--ertl-jmle/-J\tUse Ertl JMLE\n\n\n"
                         "Filtering Options --\n\n"
                         "Default: consume all kmers. Alternate options: \n"
                         "--sketch-by-fname\tAutodetect fastq or fasta data by filename (.fq or .fastq within filename).\n"
                         "--countmin/-b\tFilter all input data by count-min sketch.\n\n\n"
                         "Options for count-min filtering --\n\n"
                         "--nhashes/-H\tSet count-min number of hashes. Default: [1]\n"
                         "--cm-sketch-size/-q\tSet count-min sketch size (log2). Default: 20\n"
                         "--min-count/-n\tProvide minimum expected count for fastq data. If unspecified, all kmers are passed.\n"
                         "--seed/-R\tSet seed for seeds for count-min sketches\n\n\n"
                         "Sketch Type Options --\n\n"
                         "--use-hyperminhash\tUse HyperMinHash. Defaults to 16-bit registers. Use --bbits/-B bits to modify remainder size to 8, 16, 32, or 64.\n"
                         "--use-bb-minhash/-8\tCreate b-bit minhash sketches\n"
                         "--use-bloom-filter\tCreate bloom filter sketches\n"
                         "--use-range-minhash\tCreate range minhash sketches\n"
                         "--use-super-minhash\tCreate b-bit super minhash sketches\n"
                         "--use-counting-range-minhash\tCreate range minhash sketches\n"
                         "--use-full-khash-sets\tUse full khash sets for comparisons, rather than sketches. This can take a lot of memory and time!\n"
                         "\n\n"
                         "===Count-min-based Streaming Weighted Jaccard===\n"
                         "--wj               \tEnable weighted jaccard adapter\n"
                         "--wj-cm-sketch-size\tSet count-min sketch size for count-min streaming weighted jaccard [16]\n"
                         "--wj-cm-nhashes    \tSet count-min sketch number of hashes for count-min streaming weighted jaccard [8]\n"
                , arg);
    std::exit(EXIT_FAILURE);
}

bool fname_is_fq(const std::string &path) {
    static const std::string fq1 = ".fastq", fq2 = ".fq";
    return path.find(fq1) != std::string::npos || path.find(fq2) != std::string::npos;
}




#define SKETCH_LONG_OPTS \
static option_struct sketch_long_options[] = {\
    LO_FLAG("countmin", 'b', sm, CBF)\
    LO_FLAG("sketch-by-fname", 'f', sm, BY_FNAME)\
    LO_FLAG("no-canon", 'C', canon, false)\
    LO_FLAG("skip-cached", 'c', skip_cached, true)\
    LO_FLAG("by-entropy", 'e', entropy_minimization, true) \
    LO_FLAG("use-bb-minhash", '8', sketch_type, BB_MINHASH)\
    LO_ARG("bbits", 'B')\
    LO_ARG("paths", 'F')\
    LO_ARG("prefix", 'P')\
    LO_ARG("nhashes", 'H')\
    LO_ARG("original", 'E')\
    LO_ARG("improved", 'I')\
    LO_ARG("ertl-joint-mle", 'J')\
    LO_ARG("seed", 'R')\
    LO_ARG("sketch-size", 'S')\
    LO_ARG("kmer-length", 'k')\
    LO_ARG("min-count", 'n')\
    LO_ARG("nthreads", 'p')\
    LO_ARG("cm-sketch-size", 'q')\
    LO_ARG("spacing", 's')\
    LO_ARG("window-size", 'w')\
    LO_ARG("suffix", 'x')\
    LO_ARG("wj-cm-sketch-size", 136)\
    LO_ARG("wj-cm-nhashes", 137)\
    LO_ARG("suffix", 'x')\
\
    LO_FLAG("use-range-minhash", 128, sketch_type, RANGE_MINHASH)\
    LO_FLAG("use-counting-range-minhash", 129, sketch_type, COUNTING_RANGE_MINHASH)\
    LO_FLAG("use-full-khash-sets", 130, sketch_type, FULL_KHASH_SET)\
    LO_FLAG("use-bloom-filter", 131, sketch_type, BLOOM_FILTER)\
    LO_FLAG("use-super-minhash", 132, sketch_type, BB_SUPERMINHASH)\
    LO_FLAG("use-nthash", 133, enct, NTHASH)\
    LO_FLAG("use-cyclic-hash", 134, enct, NTHASH)\
    LO_FLAG("avoid-sorting", 135, avoid_fsorting, true)\
    LO_FLAG("wj", 138, weighted_jaccard, true)\
    SHARED_OPTS \
    {0,0,0,0}\
};

// Main functions
int sketch_main(int argc, char *argv[]) {
    int wsz(0), k(31), sketch_size(10), skip_cached(false), co, nthreads(1), mincount(1), nhashes(1), cmsketchsize(-1);
    int canon(true);
    int entropy_minimization = false, avoid_fsorting = false, weighted_jaccard = false;
    hll::EstimationMethod estim = hll::EstimationMethod::ERTL_MLE;
    hll::JointEstimationMethod jestim = static_cast<hll::JointEstimationMethod>(hll::EstimationMethod::ERTL_MLE);
    std::string spacing, paths_file, suffix, prefix, output_file;
    sketching_method sm = EXACT;
    Sketch sketch_type = HLL;
    EncodingType enct = BONSAI;
    uint64_t seedseedseed = 1337u;
    int option_index = 0;
    SKETCH_LONG_OPTS
    while((co = getopt_long(argc, argv, "n:P:F:o:p:x:R:s:S:k:w:H:q:B:8JbfjEIcCeh?", sketch_long_options, &option_index)) >= 0) {
        switch(co) {
            case 'B': gargs.bbnbits = std::atoi(optarg); break;
            case 'F': paths_file = optarg; break;
            case 'H': nhashes = std::atoi(optarg); break;
            case 'E': jestim = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ORIGINAL); break;
            case 'I': jestim = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ERTL_IMPROVED); break;
            case 'J': jestim = hll::JointEstimationMethod::ERTL_JOINT_MLE; break;
            case 'P': prefix = optarg; break;
            case 'R': seedseedseed = std::strtoull(optarg, nullptr, 10); break;
            case 'S': sketch_size = std::atoi(optarg); break;
            case 'k': k = std::atoi(optarg); break;
            case 'o': output_file = optarg; break;
            case '8': sketch_type = BB_MINHASH; break;
            case 'b': sm = CBF; break;
            case 136:
                gargs.weighted_jaccard_cmsize  = std::atoi(optarg); weighted_jaccard = true; break;
            case 137:
                gargs.weighted_jaccard_nhashes = std::atoi(optarg); weighted_jaccard = true; break;
            case 'n':
                      mincount = std::atoi(optarg);
                      std::fprintf(stderr, "mincount: %d\n", mincount);
                      break;
            case 'p': nthreads = std::atoi(optarg); break;
            case 'q': cmsketchsize = std::atoi(optarg); break;
            case 's': spacing = optarg; break;
            case 'w': wsz = std::atoi(optarg); break;
            case 'x': suffix = optarg; break;
            case 'h': case '?': sketch_usage(*argv); break;
        }
    }
    if(k > 32 && enct == BONSAI)
        UNRECOVERABLE_ERROR("k must be <= 32 for non-rolling hashes.");
    if(k > 32 && spacing.size())
        UNRECOVERABLE_ERROR("kmers must be unspaced for k > 32");
    nthreads = std::max(nthreads, 1);
    omp_set_num_threads(nthreads);
    std::fprintf(stderr, "Using %d threads\n", nthreads);
    Spacer sp(k, wsz, parse_spacing(spacing.data(), k));
    std::vector<bool> use_filter;
    std::vector<CountingSketch> cms;
    std::vector<std::string> inpaths(paths_file.size() && isfile(paths_file)
        ? get_paths(paths_file.data())
        : std::vector<std::string>(argv + optind, argv + argc));
    LOG_INFO("Sketching genomes with sketch: %d/%s\n", sketch_type, sketch_names[sketch_type]);
    if(inpaths.empty()) {
        std::fprintf(stderr, "No paths. See usage.\n");
        sketch_usage(*argv);
    }
    if(!avoid_fsorting)
        detail::sort_paths_by_fsize(inpaths);
    if(sm != EXACT) {
        if(cmsketchsize < 0) {
            cmsketchsize = 20;
            LOG_WARNING("Note: count-min sketch size not set. Defaulting to 20 for log2(sketch_size).\n");
        }
        if(sm == CBF)
            use_filter = std::vector<bool>(inpaths.size(), true);
        else // BY_FNAME
            for(const auto &path: inpaths) use_filter.emplace_back(fname_is_fq(path));
        while(cms.size() < unsigned(nthreads))
#if DASHING_USE_HK
            cms.emplace_back(cmsketchsize, nhashes, 1.05, (cms.size() ^ seedseedseed) * 1337uL);
#else
            cms.emplace_back(16, cmsketchsize, nhashes, (cms.size() ^ seedseedseed) * 1337uL);
#endif
    }
    KSeqBufferHolder kseqs(nthreads);
    if(wsz < (int)sp.c_) wsz = sp.c_;
    const int sketch_flags = skip_cached | (int(canon) << 1) | (int(entropy_minimization) << 2);
#define SKETCH_CORE(type) \
    do {\
        if(gargs.exact_weighted) {\
            sketch_core<wj::WeightedSketcher<type, wj::ExactCountingAdapter>>(sketch_size, nthreads, wsz, k, sp, inpaths,\
                                                                              suffix, prefix, cms, estim, jestim,\
                                                                              kseqs, use_filter, spacing, sketch_flags, mincount, enct, output_file);\
        } else if(weighted_jaccard) {\
            sketch_core<wj::WeightedSketcher<type>>(sketch_size, nthreads, wsz, k, sp, inpaths,\
                                                    suffix, prefix, cms, estim, jestim,\
                                                    kseqs, use_filter, spacing, sketch_flags, mincount, enct, output_file);\
        } else {\
        sketch_core<type>(sketch_size, nthreads, wsz, k, sp, inpaths,\
                          suffix, prefix, cms, estim, jestim,\
                          kseqs, use_filter, spacing, sketch_flags, mincount, enct, output_file);\
        } \
    } while(0)
    switch(sketch_type) {
        case HLL: if(gargs.defer_hll_creation) SKETCH_CORE(HLLH); else SKETCH_CORE(hll::hll_t); break;
        case WIDE_HLL: SKETCH_CORE(sketch::WideHyperLogLogHasher<>); break;
        case BLOOM_FILTER: SKETCH_CORE(bf::bf_t); break;
        case RANGE_MINHASH: SKETCH_CORE(mh::RangeMinHash<uint64_t>); break;
        case COUNTING_RANGE_MINHASH: SKETCH_CORE(mh::CountingRangeMinHash<uint64_t>); break;
        case BB_MINHASH: SKETCH_CORE(mh::BBitMinHasher<uint64_t>); break;
        case BB_SUPERMINHASH: SKETCH_CORE(SuperMinHashType); break;
        case FULL_KHASH_SET: SKETCH_CORE(khset64_t); break;
        case HYPERMINHASH: SKETCH_CORE(sketch::HyperMinHash); break;
        default: {
            char buf[128];
            std::sprintf(buf, "Sketch %s not yet supported.\n", (size_t(sketch_type) >= (sizeof(sketch_names) / sizeof(char *)) ? "Not such sketch": sketch_names[sketch_type]));
            UNRECOVERABLE_ERROR(buf);
        }
    }
#undef SKETCH_CORE
    LOG_INFO("Successfully finished sketching from %zu files\n", inpaths.size());
    return EXIT_SUCCESS;
}




namespace {
enum CompReading: unsigned {
    UNCOMPRESSED,
    GZ,
    AUTODETECT
};
}




int print_binary_main(int argc, char *argv[]) {
    int c;
    bool use_scientific = false;
    std::string outpath;
    for(char **p(argv); *p; ++p) if(std::strcmp(*p, "-h") && std::strcmp(*p, "--help") == 0) goto usage;
    if(argc == 1) {
        usage:
        std::fprintf(stderr, "%s printmat <path to binary file> [- to read from stdin]\n"
                             "-o\tSpecify output file (default: stdout)\n"
                             "-s\tEmit in scientific notation\n",
                     argv ? static_cast<const char *>(*argv): "dashing");
        std::exit(EXIT_FAILURE);
    }
    while((c = getopt(argc, argv, ":o:sh?")) >= 0) {
        switch(c) {
            case 'o': outpath = optarg; break;
            case 's': use_scientific = true; break;
            case 'h': case '?': goto usage;
        }
    }
    std::FILE *fp;
    if(outpath.empty()) outpath = "/dev/stdout";
    dm::DistanceMatrix<float> mat(argv[optind]);
    if((fp = std::fopen(outpath.data(), "wb")) == nullptr) UNRECOVERABLE_ERROR(ks::sprintf("[print_binary_main] Could not open file at %s", outpath.data()).data());
    mat.printf(fp, use_scientific);
    std::fclose(fp);
    return EXIT_SUCCESS;
}


#if 0
void flatten_usage() {
    std::fprintf(stderr, "Usage: dashing flatten <output.bin> [in1.bin in2.bin...]\n");
    std::exit(1);
}

int flatten_main(int argc, char *argv[]) {
    if(argc < 3 || std::find_if(argv, argv + argc, [](auto x){return std::strcmp(x, "-h") == 0;}) != argv + argc) flatten_usage();
    std::vector<std::string> fpaths(argv + 2, argv + argc);
    std::vector<unsigned> ks(fpaths.size(), -1);
    omp_set_num_threads(std::thread::hardware_concurrency());
    return flatten_all(fpaths, fpaths.size(), argv[1], ks);
}
#endif

int setdist_main(int argc, char *argv[]) {
    UNRECOVERABLE_ERROR("setdist_main was deprecated and has ben removed. Instead, call `dashing dist` with --use-full-khash-sets to use hash sets instead of sketches.\n");
    return 1;
}

void union_usage [[noreturn]] (char *ex) {
    std::fprintf(stderr, "Usage: %s genome1 <genome2>...\n"
                         "Flags:\n"
                         "-o: Write union sketch to file [/dev/stdout]\n"
                         "-z: Emit compressed sketch\n"
                         "-Z: Set gzip compression level\n"
                         "-r: RangeMinHash sketches\n"
                         "-H: Full Khash Sets\n"
                         "-b: Bloom Filters\n"
                ,
                 ex);
    std::exit(1);
}

int sketch_by_seq_main(int argc, char *argv[]) {
    int wsz(0), k(31), sketch_size(10), skip_cached(false), co, nthreads(1), mincount(1), nhashes(1), cmsketchsize(-1);
    int canon(true);
    int entropy_minimization = false, avoid_fsorting = false, weighted_jaccard = false;
    hll::EstimationMethod estim = hll::EstimationMethod::ERTL_MLE;
    hll::JointEstimationMethod jestim = static_cast<hll::JointEstimationMethod>(hll::EstimationMethod::ERTL_MLE);
    std::string spacing, paths_file, suffix, prefix;
    sketching_method sm = EXACT;
    Sketch sketch_type = HLL;
    EncodingType enct = BONSAI;
    uint64_t seedseedseed = 1337u;
    int option_index = 0;
    std::string outpath = "/dev/stdout";
    SKETCH_LONG_OPTS
    while((co = getopt_long(argc, argv, "o:n:P:p:x:R:s:S:k:w:H:q:B:8JbfjEIcCeh?", sketch_long_options, &option_index)) >= 0) {
        switch(co) {
            case 'B': gargs.bbnbits = std::atoi(optarg); break;
            case 'H': nhashes = std::atoi(optarg); break;
            case 'E': jestim = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ORIGINAL); break;
            case 'I': jestim = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ERTL_IMPROVED); break;
            case 'J': jestim = hll::JointEstimationMethod::ERTL_JOINT_MLE; break;
            case 'P': prefix = optarg; break;
            case 'R': seedseedseed = std::strtoull(optarg, nullptr, 10); break;
            case 'S': sketch_size = std::atoi(optarg); break;
            case 'k': k = std::atoi(optarg); break;
            case '8': sketch_type = BB_MINHASH; break;
            case 'b': sm = CBF; break;
            case 'o': outpath = optarg; break;
            case 136:
                gargs.weighted_jaccard_cmsize  = std::atoi(optarg); weighted_jaccard = true; break;
            case 137:
                gargs.weighted_jaccard_nhashes = std::atoi(optarg); weighted_jaccard = true; break;
            case 'n':
                      mincount = std::atoi(optarg);
                      std::fprintf(stderr, "mincount: %d\n", mincount);
                      break;
            case 'p': nthreads = std::atoi(optarg); break;
            case 'q': cmsketchsize = std::atoi(optarg); break;
            case 's': spacing = optarg; break;
            case 'w': wsz = std::atoi(optarg); break;
            case 'x': suffix = optarg; break;
            case 'h': case '?': sketch_usage(*argv); break;
        }
    }
    if(k > 32 && enct == BONSAI)
        UNRECOVERABLE_ERROR("k must be <= 32 for non-rolling hashes.");
    if(k > 32 && spacing.size())
        UNRECOVERABLE_ERROR("kmers must be unspaced for k > 32");
    if(argc != optind + 1) sketch_by_seq_usage(*argv);
    if(nthreads > 1)
        std::fprintf(stderr, "note: sketch_by_seq isn't parallelized.");
    nthreads = std::max(nthreads, 1);
    omp_set_num_threads(nthreads);
    auto leftover = argc - optind;
    std::string inpath;
    switch(leftover) {
        case 0: inpath = "/dev/stdin"; break;
        case 1: inpath = argv[optind]; break;
        default: throw std::runtime_error("sketch_by_seq_main only takes one path");
    }
    const bool use_filter = sm == CBF ? true: sm == BY_FNAME ? fname_is_fq(inpath): false;
    std::unique_ptr<CountingSketch> cs(use_filter ? new CountingSketch(cmsketchsize, nhashes, 1.08, seedseedseed): nullptr);
    //if(use_filter)
    //    cs = new CountingSketch(cmsketchsize, nhashes, 1.08, seedseedseed);
    const int sketch_flags = skip_cached | (int(canon) << 1) | (int(entropy_minimization) << 2);
    const Spacer sp(k, wsz, parse_spacing(spacing.data(), k));
#define SKETCH_BY_SEQ_CORE(type) \
    sketch_by_seq_core<type>(sketch_size, nthreads, sp, inpath, outpath,\
                             cs.get(), estim, jestim,\
                             use_filter, sketch_flags, mincount, enct)
    switch(sketch_type) {
        case HLL: if(gargs.defer_hll_creation) SKETCH_BY_SEQ_CORE(hll::hll_t);
                  else                         SKETCH_BY_SEQ_CORE(HLLH);
        break;
        case WIDE_HLL: SKETCH_BY_SEQ_CORE(sketch::WideHyperLogLogHasher<>); break;
        case BLOOM_FILTER: SKETCH_BY_SEQ_CORE(bf::bf_t); break;
        case RANGE_MINHASH: SKETCH_BY_SEQ_CORE(mh::RangeMinHash<uint64_t>); break;
        case COUNTING_RANGE_MINHASH: SKETCH_BY_SEQ_CORE(mh::CountingRangeMinHash<uint64_t>); break;
        case BB_MINHASH: SKETCH_BY_SEQ_CORE(mh::BBitMinHasher<uint64_t>); break;
        case BB_SUPERMINHASH: SKETCH_BY_SEQ_CORE(SuperMinHashType); break;
        case FULL_KHASH_SET: SKETCH_BY_SEQ_CORE(khset64_t); break;
        case HYPERMINHASH: SKETCH_BY_SEQ_CORE(sketch::HyperMinHash); break;
        default: {
            char buf[128];
            std::sprintf(buf, "Sketch %s not yet supported.\n", (size_t(sketch_type) >= (sizeof(sketch_names) / sizeof(char *)) ? "Not such sketch": sketch_names[sketch_type]));
            UNRECOVERABLE_ERROR(buf);
        }
    }
    return EXIT_SUCCESS;
}

int view_main(int argc, char *argv[]) {
    if(argc < 2) UNRECOVERABLE_ERROR("Usage: dashing view f1.hll [f2.hll ...]. Only HLLs currently supported.");
    for(int i = 1; i < argc; hll::hll_t(argv[i++]).printf(stdout));
    return 0;
}

} // namespace bns

