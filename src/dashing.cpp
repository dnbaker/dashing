#include "tinythreadpp/source/fast_mutex.h"
#include <fstream>
#include <omp.h>
#include "bonsai/bonsai/include/util.h"
#include "bonsai/bonsai/include/database.h"
#include "bonsai/bonsai/include/bitmap.h"
#include "bonsai/bonsai/include/setcmp.h"
#include "hll/bbmh.h"
#include "hll/mh.h"
#include "hll/ccm.h"
#include "khset/khset.h"
#include "distmat/distmat.h"
#include <sstream>
#include "getopt.h"
#include <sys/stat.h>

using namespace sketch;
using circ::roundup;
using hll::hll_t;
using sketch::common::NotImplementedError;

#ifndef BUFFER_FLUSH_SIZE
#define BUFFER_FLUSH_SIZE (1u << 18)
#endif

#if ZWRAP_USE_ZSTD
#define COMPRESSED_FILE_SUFFIX ".zst"
#else
#define COMPRESSED_FILE_SUFFIX ".gz"
#endif

using option_struct = struct option;
namespace bns {
using sketch::common::WangHash;
static const char *executable = nullptr;
enum EmissionType {
    MASH_DIST = 0,
    JI        = 1,
    SIZES     = 2,
    FULL_MASH_DIST = 3
};

enum EmissionFormat: unsigned {
    UT_TSV = 0,
    BINARY   = 1,
    UPPER_TRIANGULAR = 2,
    PHYLIP_UPPER_TRIANGULAR = 2,
    FULL_TSV = 3,
};

enum Sketch: int {
    HLL,
    BLOOM_FILTER,
    RANGE_MINHASH,
    FULL_KHASH_SET,
    COUNTING_RANGE_MINHASH,
    BB_MINHASH,
    COUNTING_BB_MINHASH, // TODO make this work.
};
static const char *sketch_names [] {
    "HLL/HyperLogLog",
    "BF/BloomFilter",
    "RMH/Range Min-Hash/KMV",
    "FHS/Full Hash Set",
    "CRHM/Counting Range Minhash",
    "BB/B-bit Minhash",
    "CBB/Counting B-bit Minhash",
};
using CBBMinHashType = mh::CountingBBitMinHasher<uint64_t, uint16_t>; // Is counting to 65536 enough for a transcriptome?
template<typename T> struct SketchEnum;

template<> struct SketchEnum<hll::hll_t> {static constexpr Sketch value = HLL;};
template<> struct SketchEnum<bf::bf_t> {static constexpr Sketch value = BLOOM_FILTER;};
template<> struct SketchEnum<mh::RangeMinHash<uint64_t>> {static constexpr Sketch value = RANGE_MINHASH;};
template<> struct SketchEnum<mh::CountingRangeMinHash<uint64_t>> {static constexpr Sketch value = COUNTING_RANGE_MINHASH;};
template<> struct SketchEnum<mh::BBitMinHasher<uint64_t>> {static constexpr Sketch value = BB_MINHASH;};
template<> struct SketchEnum<CBBMinHashType> {static constexpr Sketch value = COUNTING_BB_MINHASH;};

static uint32_t bbnbits = 16;

template<typename T>
double cardinality_estimate(T &x);
template<>
double cardinality_estimate(hll::hll_t &x) {return x.report();}
template<>
double cardinality_estimate(mh::FinalBBitMinHash &x) {return x.est_cardinality_;}
template<>
double cardinality_estimate(mh::BBitMinHasher<uint64_t> &x) {return x.cardinality_estimate();}

static size_t bytesl2_to_arg(int nblog2, Sketch sketch) {
    switch(sketch) {
        case HLL: return nblog2;
        case BLOOM_FILTER: return nblog2 + 3; // 8 bits per byte
        case RANGE_MINHASH: return size_t(1) << (nblog2 - 3); // 8 bytes per minimizer
        case COUNTING_RANGE_MINHASH: return (size_t(1) << (nblog2)) / (sizeof(uint64_t) + sizeof(uint32_t));
        case BB_MINHASH: return nblog2 - std::ceil(std::log2(bbnbits / 8));
        case FULL_KHASH_SET: return 16; // Reserve hash set size a bit. Mostly meaningless, resizing as necessary.
        default: {
            char buf[128];
            std::sprintf(buf, "Sketch %s not yet supported.\n", (size_t(sketch) >= (sizeof(sketch_names) / sizeof(char *)) ? "Not such sketch": sketch_names[sketch]));
            RUNTIME_ERROR(buf);
            return -1337;
        }
    }

}

struct khset64_t: public kh::khset64_t {
    void addh(uint64_t v) {this->insert(v);}
    khset64_t(): kh::khset64_t() {}
    khset64_t(size_t reservesz): kh::khset64_t(reservesz) {}
    double jaccard_index(const khset64_t &other) const {
        auto p1 = this, p2 = &other;
        if(size() > other.size())
            std::swap(p1, p2);
        size_t olap = 0;
        p1->for_each([&](auto v) {olap += p2->contains(v);});
        return static_cast<double>(olap) / (p1->size() + p2->size() - olap);
    }
};

template<typename SketchType>
struct FinalSketch {
    using final_type = SketchType;
};
#define FINAL_OVERLOAD(type) \
template<> struct FinalSketch<type> { \
    using final_type = typename type::final_type;}
FINAL_OVERLOAD(mh::CountingRangeMinHash<uint64_t>);
FINAL_OVERLOAD(mh::RangeMinHash<uint64_t>);
FINAL_OVERLOAD(mh::BBitMinHasher<uint64_t>);
FINAL_OVERLOAD(CBBMinHashType);
template<typename T>struct SketchFileSuffix {static constexpr const char *suffix = ".sketch";};
#define SSS(type, suf) template<> struct SketchFileSuffix<type> {static constexpr const char *suffix = suf;}
SSS(mh::CountingRangeMinHash<uint64_t>, ".crmh");
SSS(mh::RangeMinHash<uint64_t>, ".rmh");
SSS(khset64_t, ".khs");
SSS(bf::bf_t, ".bf");
SSS(mh::BBitMinHasher<uint64_t>, ".bmh");
SSS(CBBMinHashType, ".cbmh");
SSS(mh::HyperMinHash<uint64_t>, ".hmh");
SSS(hll::hll_t, ".hll");

using CRMFinal = mh::FinalCRMinHash<uint64_t, std::greater<uint64_t>, uint32_t>;
template<typename T> INLINE double similarity(const T &a, const T &b) {
    return jaccard_index(a, b);
}

template<> INLINE double similarity<CRMFinal>(const CRMFinal &a, const CRMFinal &b) {
    return a.histogram_intersection(b);
}

namespace us {
template<typename T> INLINE double union_size(const T &a, const T &b) {
    throw NotImplementedError(std::string("union_size not available for type ") + __PRETTY_FUNCTION__);
    [[unreachable]];
    return 0.;
}

template<> INLINE double union_size<hll::hllbase_t<>> (const hll::hllbase_t<> &a, const hll::hllbase_t<> &b) {
    return a.union_size(b);
}
template<> INLINE double union_size<mh::FinalBBitMinHash> (const mh::FinalBBitMinHash &a, const mh::FinalBBitMinHash &b) {
    return (a.est_cardinality_ + b.est_cardinality_ ) / (1. + a.jaccard_index(b));
}
}


void main_usage(char **argv) {
    std::fprintf(stderr, "Usage: %s <subcommand> [options...]. Use %s <subcommand> for more options. [Subcommands: sketch, dist, setdist, hll, printmat.]\n",
                 *argv, *argv);
    std::exit(EXIT_FAILURE);
}

size_t posix_fsize(const char *path) {
    struct stat st;
    stat(path, &st);
    return st.st_size;
}

namespace detail {
struct path_size {
    friend void swap(path_size&, path_size&);
    std::string path;
    size_t size;
    path_size(std::string &&p, size_t sz): path(std::move(p)), size(sz) {}
    path_size(const std::string &p, size_t sz): path(p), size(sz) {}
    path_size(path_size &&o): path(std::move(o.path)), size(o.size) {}
    path_size(): size(0) {}
    path_size &operator=(path_size &&o) {
        std::swap(o.path, path);
        std::swap(o.size, size);
        return *this;
    }
};

inline void swap(path_size &a, path_size &b) {
    std::swap(a.path, b.path);
    std::swap(a.size, b.size);
}

void sort_paths_by_fsize(std::vector<std::string> &paths) {
    uint32_t *fsizes = static_cast<uint32_t *>(std::malloc(paths.size() * sizeof(uint32_t)));
    if(!fsizes) throw std::bad_alloc();
    #pragma omp parallel for
    for(size_t i = 0; i < paths.size(); ++i)
        fsizes[i] = posix_fsize(paths[i].data());
    std::vector<path_size> ps(paths.size());
    #pragma omp parallel for
    for(size_t i = 0; i < paths.size(); ++i)
        ps[i] = path_size(paths[i], fsizes[i]);
    std::free(fsizes);
    std::sort(ps.begin(), ps.end(), [](const auto &x, const auto &y) {return x.size > y.size;});
    paths.clear();
    for(const auto &p: ps) paths.emplace_back(std::move(p.path));
}
} // namespace detail

void dist_usage(const char *arg) {
    std::fprintf(stderr, "Usage: %s <opts> [genomes if not provided from a file with -F]\n"
                         "Flags:\n"
                         "-h/-?\tUsage\n"
                         "-k\tSet kmer size [31]\n"
                         "-W\tCache sketches/use cached sketches\n"
                         "-p\tSet number of threads [1]\n"
                         "-b\tEmit distances in binary (default: human-readable, upper-triangular)\n"
                         "-U\tEmit distances in PHYLIP upper triangular format(default: human-readable, upper-triangular)\n"
                         "-s\tadd a spacer of the format <int>x<int>,<int>x<int>,"
                         "..., where the first integer corresponds to the space "
                         "between bases repeated the second integer number of times\n"
                         "-w\tSet window size [max(size of spaced kmer, [parameter])]\n"
                         "-S\tSet sketch size [10, for 2**10 bytes each]\n"
                         "-H\tTreat provided paths as pre-made sketches.\n"
                         "-C\tDo not canonicalize. [Default: canonicalize]\n"
                         "-P\tSet prefix for sketch file locations [empty]\n"
                         "-x\tSet suffix in sketch file names [empty]\n"
                         "-o\tOutput for genome size estimates [stdout]\n"
                         "-I\tUse Ertl's Improved Estimator\n"
                         "-E\tUse Ertl's Original Estimator\n"
                         "-J\tUse Ertl's JMLE Estimator [default:Uses Ertl-MLE]\n"
                         "-O\tOutput for genome distance matrix [stdout]\n"
                         "-e\tEmit in scientific notation\n"
                         "-F\tGet paths to genomes from file rather than positional arguments\n"
                         "-M\tEmit Mash distance (default: jaccard index)\n"
                         "-l\tEmit full (not approximate) Mash distance. default: jaccard index\n"
                         "-T\tpostprocess binary format to human-readable TSV (not upper triangular)\n"
                         "-Z\tEmit genome sizes (default: jaccard index)\n"
                         "-N\tAutodetect fastq or fasta data by filename (.fq or .fastq within filename).\n"
                         "-n\tAvoid sorting files by genome sizes. This avoids a computational step, but can result in degraded load-balancing.\n"
                         "-y\tFilter all input data by count-min sketch.\n"
                         "-q\tSet count-min number of hashes. Default: [4]\n"
                         "-c\tSet minimum count for kmers to pass count-min filtering.\n"
                         "-t\tSet count-min sketch size (log2). Default: ceil(log2(max_filesize)) + 2\n"
                         "-R\tSet seed for seeds for count-min sketches\n"
                         "-b\tSet `b` for b-bit minwise hashing to <int>. Default: 16\n"
                , arg);
    std::exit(EXIT_FAILURE);
}


// Usage, utilities
void sketch_usage(const char *arg) {
    std::fprintf(stderr, "Usage: %s <opts> [genomes if not provided from a file with -F]\n"
                         "Flags:\n"
                         "-h/-?:\tEmit usage\n"
                         "\n"
                         "Sketch options --\n"
                         "--kmer-length/-k\tSet kmer size [31]\n"
                         "--spacing/-s\tadd a spacer of the format <int>x<int>,<int>x<int>,"
                         "..., where the first integer corresponds to the space "
                         "between bases repeated the second integer number of times\n"
                         "--window-size/-w\tSet window size [max(size of spaced kmer, [parameter])]\n"
                         "--sketch-size/-S\tSet log2 sketch size in bytes [10, for 2**10 bytes each]\n"
                         "--no-canon/-C\tDo not canonicalize. [Default: canonicalize]\n"
                         "--bbits/-B\tSet `b` for b-bit minwise hashing to <int>. Default: 16\n"
                         "Run options --\n"
                         "--nthreads/-p\tSet number of threads [1]\n"
                         "--prefix/-P\tSet prefix for sketch file locations [empty]\n"
                         "--suffix/-x\tSet suffix in sketch file names [empty]\n"
                         "--paths/-F\tGet paths to genomes from file rather than positional arguments\n"
                         "--skip-cached/-c\tSkip alreday produced/cached sketches (save sketches to disk in directory of the file [default] or in folder specified by -P\n"
                         "Estimation methods --\n"
                         "--original/-E\tUse Flajolet with inclusion/exclusion quantitation method for hll. [Default: Ertl MLE]\n"
                         "--improved/-I\tUse Ertl Improved estimator [Default: Ertl MLE]\n"
                         "--ertl-jmle/-J\tUse Ertl JMLE\n"
                         "Filtering Options --\n"
                         "Default: consume all kmers. Alternate options: \n"
                         "--sketch-by-fname/-f\tAutodetect fastq or fasta data by filename (.fq or .fastq within filename).\n"
                         "--count-min/-b\tFilter all input data by count-min sketch.\n"
                         "Options for count-min filtering --\n"
                         "--nhashes/-H\tSet count-min number of hashes. Default: [4]\n"
                         "--cm-sketch-size/-q\tSet count-min sketch size (log2). Default: ceil(log2(max_filesize)) + 2\n"
                         "--min-count/-n\tProvide minimum expected count for fastq data. If unspecified, all kmers are passed.\n"
                         "--seed/-R\tSet seed for seeds for count-min sketches\n"
                         "Sketch Type Options --\n"
                         "--use-bb-minhash/-8\tCreate bbit minhash sketches\n"
                         "--use-range-minhash\tCreate range minhash sketches\n"
                         "--use-counting-range-minhash\tCreate range minhash sketches\n"
                         "----\n"
                , arg);
    std::exit(EXIT_FAILURE);
}

bool fname_is_fq(const std::string &path) {
    static const std::string fq1 = ".fastq", fq2 = ".fq";
    return path.find(fq1) != std::string::npos || path.find(fq2) != std::string::npos;
}

template<typename SketchType>
std::string make_fname(const char *path, size_t sketch_p, int wsz, int k, int csz, const std::string &spacing, const std::string &suffix="", const std::string &prefix="") {
    std::string ret(prefix);
    {
        const char *p;
        if(ret.size() && (p = strrchr(get_cstr(path), '/')))
            ret += std::string(p);
        else
            ret += get_cstr(path);
    }
    ret += ".w";
    ret + std::to_string(std::max(csz, wsz));
    ret += ".";
    ret += std::to_string(k);
    ret += ".spacing";
    ret += spacing;
    ret += '.';
    if(suffix.size()) {
        ret += "suf";
        ret += suffix;
        ret += '.';
    }
    ret += std::to_string(sketch_p);
    ret += SketchFileSuffix<SketchType>::suffix;
    return ret;
}

enum sketching_method: int {
    EXACT = 0,
    CBF   = 1,
    BY_FNAME = 2
};


size_t fsz2countcm(uint64_t fsz, double factor=1.) {
    return std::log2(roundup(size_t(fsz * factor))) + 2; // plus 2 to account for the fact that the file is likely compressed
}

size_t fsz2count(uint64_t fsz) {
    static constexpr size_t mul = 4; // Take these estimates and multiply by 4 just to be safe
    // This should be adapted according to the error rate of the pcbf
    if(fsz < 300ull << 20) return 1 * mul;
    if(fsz < 500ull << 20) return 3 * mul;
    if(fsz < 1ull << 30)  return 10 * mul;
    if(fsz < 3ull << 30)  return 20 * mul;
    return static_cast<size_t>(std::pow(2., std::log((fsz >> 30))) * 30.) * mul;
    // This likely does not account for compression.
}

/*
 *
  enc.for_each([&](u64 kmer){h.addh(kmer);}, inpaths[i].data(), &kseqs[tid]);\
 */

template<typename T>
INLINE void set_estim_and_jestim(T &x, hll::EstimationMethod estim, hll::JointEstimationMethod jestim) {}

template<typename Hashstruct>
INLINE void set_estim_and_jestim(hll::hllbase_t<Hashstruct> &h, hll::EstimationMethod estim, hll::JointEstimationMethod jestim) {
    h.set_estim(estim);
    h.set_jestim(jestim);
}
using hll::EstimationMethod;
using hll::JointEstimationMethod;
template<typename T>
T construct(size_t ssarg) {
    return T(ssarg);
}
template<> mh::BBitMinHasher<uint64_t> construct<mh::BBitMinHasher<uint64_t>>(size_t p) {return mh::BBitMinHasher<uint64_t>(p, bbnbits);}

template<typename SketchType>
void sketch_core(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<cm::ccm_t> &cms, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, bool skip_cached, bool canon, uint32_t mincount, bool entropy_minimization) {
    std::vector<SketchType> sketches;
    uint32_t sketch_size = bytesl2_to_arg(ssarg, SketchEnum<SketchType>::value);
    while(sketches.size() < (u32)nthreads) sketches.push_back(construct<SketchType>(sketch_size)), set_estim_and_jestim(sketches.back(), estim, jestim);
    std::vector<std::string> fnames(nthreads);

#define MAIN_SKETCH_LOOP(MinType)\
    for(size_t i = 0; i < inpaths.size(); ++i) {\
        const int tid = omp_get_thread_num();\
        std::string &fname = fnames[tid];\
        fname = make_fname<SketchType>(inpaths[i].data(), sketch_size, wsz, k, sp.c_, spacing, suffix, prefix);\
        if(skip_cached && isfile(fname)) continue;\
        Encoder<MinType> enc(nullptr, 0, sp, nullptr, canon);\
        auto &h = sketches[tid];\
        if(use_filter.size() && use_filter[i]) {\
            auto &cm = cms[tid];\
            enc.for_each([&](u64 kmer){if(cm.addh(kmer) >= mincount) h.addh(kmer);}, inpaths[i].data(), &kseqs[tid]);\
            cm.clear();  \
        } else {\
            enc.for_each([&](u64 kmer){h.addh(kmer);}, inpaths[i].data(), &kseqs[tid]);\
        }\
        h.write(fname.data());\
        h.clear();\
    }
    if(entropy_minimization) {
        #pragma omp parallel for schedule(dynamic)
        MAIN_SKETCH_LOOP(bns::score::Entropy)
    } else {
        #pragma omp parallel for schedule(dynamic)
        MAIN_SKETCH_LOOP(bns::score::Lex)
    }
#undef MAIN_SKETCH_LOOP
}

#if 0
#ifndef no_argument
#define no_argument 0
#endif
#ifndef required_argument
#define required_argument 1
#endif
#ifndef optional_argument
#define optional_argument 2
#endif
#endif

#define LO_ARG(LONG, SHORT) {LONG, required_argument, 0, SHORT},

#define LO_NO(LONG, SHORT) {LONG, no_argument, 0, SHORT},

#define LO_FLAG(LONG, SHORT, VAR, VAL) {LONG, no_argument, (int *)&VAR, VAL},

static_assert(sizeof(option) >= sizeof(void *), "must be bigger");

#define SKETCH_LONG_OPTS \
static option_struct sketch_long_options[] = {\
    LO_FLAG("countmin", 'b', sm, CBF)\
    LO_FLAG("no-canon", 'C', canon, false)\
    LO_FLAG("sketch-by-fname", 'f', sm, BY_FNAME)\
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
\
    LO_FLAG("use-range-minhash", 128, sketch_type, RANGE_MINHASH)\
    LO_FLAG("use-counting-range-minhash", 129, sketch_type, COUNTING_RANGE_MINHASH)\
};

// Main functions
int sketch_main(int argc, char *argv[]) {
    int wsz(0), k(31), sketch_size(10), skip_cached(false), co, nthreads(1), mincount(1), nhashes(1), cmsketchsize(-1);
    int canon(true);
    int entropy_minimization = false;
    hll::EstimationMethod estim = hll::EstimationMethod::ERTL_MLE;
    hll::JointEstimationMethod jestim = static_cast<hll::JointEstimationMethod>(hll::EstimationMethod::ERTL_MLE);
    std::string spacing, paths_file, suffix, prefix;
    sketching_method sm = EXACT;
    Sketch sketch_type = HLL;
    uint64_t seedseedseed = 1337u;
    int option_index = 0;
    SKETCH_LONG_OPTS
    while((co = getopt_long(argc, argv, "n:P:F:p:x:R:s:S:k:w:H:q:B:8JbfjEIcCeh?", sketch_long_options, &option_index)) >= 0) {
        switch(co) {
            case 'B': bbnbits = std::atoi(optarg); break;
            case 'E': jestim = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ORIGINAL); break;
            case 'F': paths_file = optarg; break;
            case 'H': nhashes = std::atoi(optarg); break;
            case 'I': jestim = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ERTL_IMPROVED); break;
            case 'P': prefix = optarg; break;
            case 'R': seedseedseed = std::strtoull(optarg, nullptr, 10); break;
            case 'S': sketch_size = std::atoi(optarg); break;
            case 'J': jestim = hll::JointEstimationMethod::ERTL_JOINT_MLE; break;
            case 'k': k = std::atoi(optarg); break;
            case '8': sketch_type = BB_MINHASH; break;
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
    nthreads = std::max(nthreads, 1);
    omp_set_num_threads(nthreads);
    Spacer sp(k, wsz, parse_spacing(spacing.data(), k));
    std::vector<bool> use_filter;
    std::vector<cm::ccm_t> cms;
    std::vector<std::string> inpaths(paths_file.size() ? get_paths(paths_file.data())
                                                       : std::vector<std::string>(argv + optind, argv + argc));
    LOG_INFO("Sketching genomes with sketch: %d/%s\n", sketch_type, sketch_names[sketch_type]);
    if(inpaths.empty()) {
        std::fprintf(stderr, "No paths. See usage.\n");
        sketch_usage(*argv);
    }
    detail::sort_paths_by_fsize(inpaths);
    if(sm != EXACT) {
        if(cmsketchsize < 0) {
            cmsketchsize = fsz2countcm(
                std::accumulate(inpaths.begin(), inpaths.end(), 0u,
                                [](unsigned x, const auto &y) ->unsigned {return std::max(x, (unsigned)posix_fsize(y.data()));})
            );
        }
        if(sm == CBF)
            use_filter = std::vector<bool>(inpaths.size(), true);
        else // BY_FNAME
            for(const auto &path: inpaths) use_filter.emplace_back(fname_is_fq(path));
        auto nbits = std::log2(mincount) + 1;
        while(cms.size() < unsigned(nthreads))
            cms.emplace_back(nbits, cmsketchsize, nhashes, (cms.size() ^ seedseedseed) * 1337uL);
    }
    KSeqBufferHolder kseqs(nthreads);
    if(wsz < sp.c_) wsz = sp.c_;
#define SKETCH_CORE(type) \
    sketch_core<type>(sketch_size, nthreads, wsz, k, sp, inpaths,\
                            suffix, prefix, cms, estim, jestim,\
                            kseqs, use_filter, spacing, skip_cached, canon, mincount, entropy_minimization)
    switch(sketch_type) {
        case HLL: SKETCH_CORE(hll::hll_t); break;
        case BLOOM_FILTER: SKETCH_CORE(bf::bf_t); break;
        case RANGE_MINHASH: SKETCH_CORE(mh::RangeMinHash<uint64_t>); break;
        case COUNTING_RANGE_MINHASH: SKETCH_CORE(mh::CountingRangeMinHash<uint64_t>); break;
        case BB_MINHASH: SKETCH_CORE(mh::BBitMinHasher<uint64_t>); break;
        default: {
            char buf[128];
            std::sprintf(buf, "Sketch %s not yet supported.\n", (size_t(sketch_type) >= (sizeof(sketch_names) / sizeof(char *)) ? "Not such sketch": sketch_names[sketch_type]));
            RUNTIME_ERROR(buf);
        }
    }
#undef SKETCH_CORE
    LOG_INFO("Successfully finished sketching from %zu files\n", inpaths.size());
    return EXIT_SUCCESS;
}


template<typename FType, typename=typename std::enable_if<std::is_floating_point<FType>::value>::type>
size_t submit_emit_dists(int pairfi, const FType *ptr, u64 hs, size_t index, ks::string &str, const std::vector<std::string> &inpaths, EmissionFormat emit_fmt, bool use_scientific, const size_t buffer_flush_size=BUFFER_FLUSH_SIZE) {
    if(emit_fmt & BINARY) {
        const ssize_t nbytes = sizeof(FType) * (hs - index - 1);
        LOG_DEBUG("Writing %zd bytes for %zu items\n", nbytes, (hs - index - 1));
        ssize_t i = ::write(pairfi, ptr, nbytes);
        if(i != nbytes) {
            std::fprintf(stderr, "written %zd bytes instead of expected %zd\n", i, nbytes);
        }
    } else {
        auto &strref = inpaths[index];
        str += strref;
        if(emit_fmt == UT_TSV) {
            const char *fmt = use_scientific ? "\t%e": "\t%f";
            {
                u64 k;
                for(k = 0; k < index + 1;  ++k, kputsn_("\t-", 2, reinterpret_cast<kstring_t *>(&str)));
                for(k = 0; k < hs - index - 1; str.sprintf(fmt, ptr[k++]));
            }
        } else { // emit_fmt == UPPER_TRIANGULAR
            const char *fmt = use_scientific ? " %e": " %f";
            if(strref.size() < 9)
                str.append(9 - strref.size(), ' ');
            for(u64 k = 0; k < hs - index - 1; str.sprintf(fmt, ptr[k++]));
        }
        str.putc_('\n');
        str.flush(pairfi);
    }
    return index;
}

template<typename FType1, typename FType2,
         typename=typename std::enable_if<
            std::is_floating_point<FType1>::value && std::is_floating_point<FType2>::value
          >::type
         >
typename std::common_type<FType1, FType2>::type dist_index(FType1 ji, FType2 ksinv) {
    // Adapter from Mash https://github.com/Marbl/Mash
    return ji ? -std::log(2. * ji / (1. + ji)) * ksinv: 1.;
}

template<typename FType1, typename FType2,
         typename=typename std::enable_if<
            std::is_floating_point<FType1>::value && std::is_floating_point<FType2>::value
          >::type
         >
typename std::common_type<FType1, FType2>::type containment_index(FType1 containment, FType2 ksinv) {
    // Adapter from Mash https://github.com/Marbl/Mash
    return containment ? -std::log(containment) * ksinv: 1.;
}

template<typename FType1, typename FType2,
         typename=typename std::enable_if<
            std::is_floating_point<FType1>::value && std::is_floating_point<FType2>::value
          >::type
         >
typename std::common_type<FType1, FType2>::type full_dist_index(FType1 ji, FType2 ksinv) {
    return 1. - std::pow(2.*ji/(1. + ji), ksinv);
}
template<typename FType1, typename FType2,
         typename=typename std::enable_if<
            std::is_floating_point<FType1>::value && std::is_floating_point<FType2>::value
          >::type
         >
typename std::common_type<FType1, FType2>::type full_containment_index(FType1 containment, FType2 ksinv) {
    return 1. - std::pow(containment, ksinv);
}

template<typename SketchType, typename T, typename Func>
INLINE void perform_core_op(T &dists, size_t nhlls, SketchType *hlls, const Func &func, size_t i) {
    auto &h1 = hlls[i];
    #pragma omp parallel for
    for(size_t j = i + 1; j < nhlls; ++j)
        dists[j - i - 1] = func(hlls[j], h1);
    h1.free();
}

#if defined(ENABLE_COMPUTED_GOTO) && !defined(__clang__)
#define CORE_ITER(zomg) do {{\
        static constexpr void *TYPES[] {&&mash##zomg, &&ji##zomg, &&sizes##zomg, &&full_mash##zomg};\
        goto *TYPES[result_type];\
        mash##zomg:\
            perform_core_op(dists, nsketches, hlls, [ksinv](const auto &x, const auto &y) {return dist_index(similarity<const SketchType>(x, y), ksinv);}, i);\
            goto next_step##zomg;\
        ji##zomg:\
            perform_core_op(dists, nsketches, hlls, similarity<const SketchType>, i);\
           goto next_step##zomg;\
        sizes##zomg:\
            perform_core_op(dists, nsketches, hlls, us::union_size<SketchType>, i);\
           goto next_step##zomg;\
        full_mash##zomg:\
            perform_core_op(dists, nsketches, hlls, [ksinv](const auto &x, const auto &y) {return full_dist_index(similarity<const SketchType>(x, y), ksinv);}, i);\
           goto next_step##zomg;\
        next_step##zomg: ;\
    } } while(0);
#else
#define CORE_ITER(zomg) do {\
        switch(result_type) {\
            case MASH_DIST: {\
                perform_core_op(dists, nsketches, hlls, [ksinv](const auto &x, const auto &y) {return dist_index(similarity<const SketchType>(x, y), ksinv);}, i);\
                break;\
            }\
            case JI: {\
            perform_core_op(dists, nsketches, hlls, similarity<const SketchType>, i);\
                break;\
            }\
            case SIZES: {\
            perform_core_op(dists, nsketches, hlls, us::union_size<SketchType>, i);\
                break;\
            }\
            case FULL_MASH_DIST:\
                perform_core_op(dists, nsketches, hlls, [ksinv](const auto &x, const auto &y) {return full_dist_index(similarity<const SketchType>(x, y), ksinv);}, i);\
                break;\
            default:\
                __builtin_unreachable();\
        } } while(0)
#endif

template<typename SketchType>
void dist_loop(std::FILE *ofp, SketchType *hlls, const std::vector<std::string> &inpaths, const bool use_scientific, const unsigned k, const EmissionType result_type, EmissionFormat emit_fmt, int nthreads, const size_t buffer_flush_size=BUFFER_FLUSH_SIZE) {
    const float ksinv = 1./ k;
    const int pairfi = fileno(ofp);
    omp_set_num_threads(nthreads);
    const size_t nsketches = inpaths.size();
    if((emit_fmt & BINARY) == 0) {
        std::future<size_t> submitter;
        std::array<std::vector<float>, 2> dps;
        dps[0].resize(nsketches - 1);
        dps[1].resize(nsketches - 2);
        ks::string str;
        for(size_t i = 0; i < nsketches; ++i) {
            std::vector<float> &dists = dps[i & 1];
            CORE_ITER(_a);
            LOG_DEBUG("Finished chunk %zu of %zu\n", i + 1, nsketches);
            if(i) submitter.get();
            submitter = std::async(std::launch::async, submit_emit_dists<float>,
                                   pairfi, dists.data(), nsketches, i,
                                   std::ref(str), std::ref(inpaths), emit_fmt, use_scientific, buffer_flush_size);
        }
        submitter.get();
    } else {
        dm::DistanceMatrix<float> dm(nsketches);
        for(size_t i = 0; i < nsketches; ++i) {
            auto span = dm.row_span(i);
            auto &dists = span.first;
            CORE_ITER(_b);
        }
        if(emit_fmt == FULL_TSV) dm.printf(ofp, use_scientific, &inpaths);
        else {
            assert(emit_fmt == BINARY);
            dm.write(ofp);
        }
    }
}

namespace {
enum CompReading: unsigned {
    UNCOMPRESSED,
    GZ,
    AUTODETECT
};
}

//dist_sketch_and_cmp<sketchtype>(inpaths, cms, kseqs, ofp, pairofp, sp, sketch_size, mincount, estim, jestim, cache_sketch, result_type, emit_fmt, presketched_only, nthreads,
//use_scientific, suffix, prefix, canon, entropy_minimization, spacing);
template<typename SketchType>
void dist_sketch_and_cmp(const std::vector<std::string> &inpaths, std::vector<sketch::cm::ccm_t> &cms, KSeqBufferHolder &kseqs, std::FILE *ofp, std::FILE *pairofp,
                         Spacer sp,
                         unsigned ssarg, unsigned mincount, EstimationMethod estim, JointEstimationMethod jestim, bool cache_sketch, EmissionType result_type, EmissionFormat emit_fmt,
                         bool presketched_only, unsigned nthreads, bool use_scientific, std::string suffix, std::string prefix, bool canon, bool entropy_minimization, std::string spacing) {
    using final_type = typename FinalSketch<SketchType>::final_type;
    std::vector<SketchType> sketches;
    sketches.reserve(inpaths.size());
    uint32_t sketch_size = bytesl2_to_arg(ssarg, SketchEnum<SketchType>::value);
    while(sketches.size() < inpaths.size()) {
        sketches.emplace_back(construct<SketchType>(sketch_size));
        set_estim_and_jestim(sketches.back(), estim, jestim);
    }
    static constexpr bool samesketch = std::is_same<SketchType, final_type>::value;
    final_type *final_sketches =
        samesketch ? reinterpret_cast<final_type *>(sketches.data())
                   : static_cast<final_type *>(std::malloc(sizeof(*final_sketches) * inpaths.size()));
    CONST_IF(samesketch) {
        if(final_sketches == nullptr) throw std::bad_alloc();
    }
#if 0
    for(auto &&pair: map) {
        *bcp++ = pair.first;
        new(pp++) FinalType(std::move(pair.second));
    }
    {decltype(map) tmap(std::move(map));} // Free map now that it's not needed
#endif

    std::atomic<uint32_t> ncomplete;
    ncomplete.store(0);
    const unsigned k = sp.k_;
    const unsigned wsz = sp.w_;
    #pragma omp parallel for schedule(dynamic)
    for(size_t i = 0; i < sketches.size(); ++i) {
        const std::string &path(inpaths[i]);
        auto &sketch = sketches[i];
        if(presketched_only)  {
            CONST_IF(samesketch) {
                sketch.read(path);
                set_estim_and_jestim(sketch, estim, jestim); // HLL is the only type that needs this, and it's the same
            } else new(final_sketches + i) final_type(path.data()); // Read from path
        } else {
            const std::string fpath(make_fname<SketchType>(path.data(), sketch_size, wsz, k, sp.c_, spacing, suffix, prefix));
            const bool isf = isfile(fpath);
            if(cache_sketch && isf) {
                LOG_DEBUG("Sketch found at %s with size %zu, %u\n", fpath.data(), size_t(1ull << sketch_size), sketch_size);
                CONST_IF(samesketch) {
                    sketch.read(fpath);
                    set_estim_and_jestim(sketch, estim, jestim);
                } else {
                    new(final_sketches + i) final_type(fpath);
                }
            } else {
                const int tid = omp_get_thread_num();
#define FILL_SKETCH_MIN(MinType) \
        Encoder<MinType> enc(nullptr, 0, sp, nullptr, canon);\
        if(cms.empty()) {\
            auto &h = sketch;\
            enc.for_each([&](u64 kmer){h.addh(kmer);}, inpaths[i].data(), &kseqs[tid]);\
        } else {\
            sketch::cm::ccm_t &cm = cms[tid];\
            enc.for_each([&,mincount](u64 kmer){if(cm.addh(kmer) >= mincount) sketch.addh(kmer);}, inpaths[i].data(), &kseqs[tid]);\
            cm.clear();\
        }\
        CONST_IF(!samesketch) new(final_sketches + i) final_type(std::move(sketch));
                if(entropy_minimization) {
                    FILL_SKETCH_MIN(score::Entropy);
                } else {
                    FILL_SKETCH_MIN(score::Lex);
                }
#undef FILL_SKETCH_MIN
                CONST_IF(samesketch) {
                    if(cache_sketch && !isf) sketch.write(fpath);
                } else if(cache_sketch) final_sketches[i].write(fpath);
            }
        }
        ++ncomplete; // Atomic
    }
    kseqs.free();
    ks::string str("#Path\tSize (est.)\n");
    assert(str == "#Path\tSize (est.)\n");
    str.resize(BUFFER_FLUSH_SIZE);
    {
        const int fn(fileno(ofp));
        for(size_t i(0); i < sketches.size(); ++i) {
            double card;
            CONST_IF(samesketch) {
                card = cardinality_estimate(sketches[i]);
            } else {
                card = cardinality_estimate(final_sketches[i]);
            }
            str.sprintf("%s\t%lf\n", inpaths[i].data(), card);
            if(str.size() >= BUFFER_FLUSH_SIZE) str.flush(fn);
        }
        str.flush(fn);
    }
    if(ofp != stdout) std::fclose(ofp);
    str.clear();
    if(emit_fmt == UT_TSV) {
        str.sprintf("##Names\t");
        for(const auto &path: inpaths) str.sprintf("%s\t", path.data());
        str.back() = '\n';
        str.write(fileno(pairofp)); str.free();
    } else if(emit_fmt == UPPER_TRIANGULAR) { // emit_fmt == UPPER_TRIANGULAR
        std::fprintf(pairofp, "%zu\n", inpaths.size());
        std::fflush(pairofp);
    }
    // binary formats don't have headers we handle here
    static constexpr uint32_t buffer_flush_size = BUFFER_FLUSH_SIZE;
    dist_loop<final_type>(pairofp, final_sketches, inpaths, use_scientific, k, result_type, emit_fmt, nthreads, buffer_flush_size);
    CONST_IF(!samesketch) {
#if __cplusplus >= 201703L
        std::destroy_n(
#  ifdef USE_PAR_EX
            std::execution::par_unseq,
#  endif
            final_sketches, inpaths.size());
#else
        std::for_each(final_sketches, final_sketches + inpaths.size(), [](auto &sketch) {
            using destructor_type = typename std::decay<decltype(sketch)>::type;
            sketch.~destructor_type();
        });
#endif
        std::free(final_sketches);
    }
}

#define DIST_LONG_OPTS \
static option_struct sketch_long_options[] = {\
    LO_FLAG("full-tsv", 'T', emit_fmt, FULL_TSV)\
    LO_FLAG("emit-binary", 'b', emit_fmt, BINARY)\
    LO_FLAG("no-canon", 'C', canon, false)\
    LO_FLAG("by-entropy", 'g', entropy_minimization, true) \
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
\
    LO_FLAG("use-range-minhash", 128, sketch_type, RANGE_MINHASH)\
    LO_FLAG("use-counting-range-minhash", 129, sketch_type, COUNTING_RANGE_MINHASH)\
};

#if 0
    while((co = getopt(argc, argv, "n:Q:P:x:F:c:p:o:s:w:O:S:k:=:t:R:8TgDazlICbMEeHJhZBNyUmqW?")) >= 0) {
        switch(co) {
            case '8': use_bbmh = true;                 break;
            case 'T': emit_fmt = FULL_TSV;             break; // This also sets the emit_fmt bit for BINARY
            case 'B': bbnbits = std::atoi(optarg);     break;
            case 'b': emit_fmt = BINARY;               break;
            case 'C': canon = false;                   break;
            case 'E': jestim = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ORIGINAL); break;
            case 'F': paths_file = optarg;             break;
            case 'H': presketched_only = true;         break;
            case 'I': jestim = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ERTL_IMPROVED); break;
            case 'J': jestim = hll::JointEstimationMethod::ERTL_JOINT_MLE; break;
            case 'M': result_type = MASH_DIST;         break;
            case 'N': sm = BY_FNAME;                   break;
            case 'P': prefix = optarg;                 break;
            case 'Q': querypaths.emplace_back(optarg); LOG_EXIT("Error: querypaths method temporarily removed before integration into a separate subcommand."); break;
            case 'R': seedseedseed = std::strtoull(optarg, nullptr, 10); break;
            case 'S': sketch_size = std::atoi(optarg); break;
            case 'U': emit_fmt = UPPER_TRIANGULAR;     break;
            case 'W': cache_sketch = true;             break;
            case 'Z': result_type = SIZES;             break;
            case 'c': mincount = std::atoi(optarg);    break;
            case 'e': use_scientific = true;           break;
            case 'g': entropy_minimization = true; LOG_WARNING("Entropy-based minimization is probably theoretically ill-founded, but it might be of practical value.\n"); break;
            case 'k': k = std::atoi(optarg);           break;
            case 'l': result_type = FULL_MASH_DIST;    break;
            case 'm': jestim = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ERTL_MLE); break;
            case 'n': avoid_fsorting = true;           break;
            case 'o': if((ofp = fopen(optarg, "w")) == nullptr) LOG_EXIT("Could not open file at %s for writing.\n", optarg); break;
            case 'p': nthreads = std::atoi(optarg);    break;
            case 'q': nhashes = std::atoi(optarg);     break;
            case 't': cmsketchsize = std::atoi(optarg); break;
            case 's': spacing = optarg;                break;
            case 'w': wsz = std::atoi(optarg);         break;
            case 'x': suffix = optarg;                 break;
            case 'y': sm = CBF;                        break;
            //case 'z': reading_type = GZ;               break;
            case 'O': if((pairofp = fopen(optarg, "wb")) == nullptr)
                          LOG_EXIT("Could not open file at %s for writing.\n", optarg);
                      pairofp_labels = std::string(optarg) + ".labels";
                      pairofp_path = optarg;
                      break;
            case 'h': case '?': dist_usage(*argv);
#endif

int dist_main(int argc, char *argv[]) {
    int wsz(0), k(31), sketch_size(10), use_scientific(false), co, cache_sketch(false),
        nthreads(1), mincount(30), nhashes(4), cmsketchsize(-1);
    int canon(true), presketched_only(false), entropy_minimization(false),
         avoid_fsorting(false);
    Sketch sketch_type = HLL;
         // bool sketch_query_by_seq(true);
    EmissionFormat emit_fmt = UT_TSV;
    double factor = 1.;
    EmissionType result_type(JI);
    hll::EstimationMethod estim = hll::EstimationMethod::ERTL_MLE;
    hll::JointEstimationMethod jestim = static_cast<hll::JointEstimationMethod>(hll::EstimationMethod::ERTL_MLE);
    std::string spacing, paths_file, suffix, prefix, pairofp_labels, pairofp_path;
    FILE *ofp(stdout), *pairofp(stdout);
    sketching_method sm = EXACT;
    std::vector<std::string> querypaths;
    uint64_t seedseedseed = 1337u;
    if(argc == 1) dist_usage(*argv);
    while((co = getopt(argc, argv, "n:Q:P:x:F:c:p:o:s:w:O:S:k:=:t:R:8TgDazlICbMEeHJhZBNyUmqW?")) >= 0) {
        switch(co) {
            case 'B': bbnbits = std::atoi(optarg);     break;
            case 'b': emit_fmt = BINARY;               break;
            case 'C': canon = false;                   break;
            case 'E': jestim = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ORIGINAL); break;
            case 'F': paths_file = optarg;             break;
            case 'H': presketched_only = true;         break;
            case 'I': jestim = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ERTL_IMPROVED); break;
            case 'J': jestim = hll::JointEstimationMethod::ERTL_JOINT_MLE; break;
            case 'M': result_type = MASH_DIST;         break;
            case 'N': sm = BY_FNAME;                   break;
            case 'P': prefix = optarg;                 break;
            case 'Q': querypaths.emplace_back(optarg); LOG_EXIT("Error: querypaths method temporarily removed before integration into a separate subcommand."); break;
            case 'R': seedseedseed = std::strtoull(optarg, nullptr, 10); break;
            case 'S': sketch_size = std::atoi(optarg); break;
            case 'U': emit_fmt = UPPER_TRIANGULAR;     break;
            case 'W': cache_sketch = true;             break;
            case 'Z': result_type = SIZES;             break;
            case 'c': mincount = std::atoi(optarg);    break;
            case 'e': use_scientific = true;           break;
            case 'g': entropy_minimization = true; LOG_WARNING("Entropy-based minimization is probably theoretically ill-founded, but it might be of practical value.\n"); break;
            case 'k': k = std::atoi(optarg);           break;
            case 'l': result_type = FULL_MASH_DIST;    break;
            case 'm': jestim = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ERTL_MLE); break;
            case 'n': avoid_fsorting = true;           break;
            case 'o': if((ofp = fopen(optarg, "w")) == nullptr) LOG_EXIT("Could not open file at %s for writing.\n", optarg); break;
            case 'p': nthreads = std::atoi(optarg);    break;
            case 'q': nhashes = std::atoi(optarg);     break;
            case 't': cmsketchsize = std::atoi(optarg); break;
            case 's': spacing = optarg;                break;
            case 'w': wsz = std::atoi(optarg);         break;
            case 'x': suffix = optarg;                 break;
            case 'y': sm = CBF;                        break;
            //case 'z': reading_type = GZ;               break;
            case 'O': if((pairofp = fopen(optarg, "wb")) == nullptr)
                          LOG_EXIT("Could not open file at %s for writing.\n", optarg);
                      pairofp_labels = std::string(optarg) + ".labels";
                      pairofp_path = optarg;
                      break;
            case 'h': case '?': dist_usage(*argv);
        }
    }
    if(nthreads < 0) nthreads = 1;
    std::vector<std::string> inpaths(paths_file.size() ? get_paths(paths_file.data())
                                                       : std::vector<std::string>(argv + optind, argv + argc));
    if(inpaths.empty())
        std::fprintf(stderr, "No paths. See usage.\n"), dist_usage(*argv);
    omp_set_num_threads(nthreads);
    Spacer sp(k, wsz, parse_spacing(spacing.data(), k));
    if(!presketched_only && !avoid_fsorting)
        detail::sort_paths_by_fsize(inpaths);
    std::vector<sketch::cm::ccm_t> cms;
    KSeqBufferHolder kseqs(nthreads);
    switch(sm) {
        case CBF: case BY_FNAME: {
            if(cmsketchsize < 0) {
                cmsketchsize = fsz2countcm(
                std::accumulate(inpaths.begin(), inpaths.end(), 0u,
                                [](unsigned x, const auto &y) ->unsigned {return std::max(x, (unsigned)posix_fsize(y.data()));}), factor);
            }
            unsigned nbits = std::log2(mincount) + 1;
            while(cms.size() < static_cast<unsigned>(nthreads))
                cms.emplace_back(nbits, cmsketchsize, nhashes, (cms.size() ^ seedseedseed) * 1337uL);
            break;
        }
        case EXACT: default: break;
    }
#define CALL_DIST(sketchtype) \
        dist_sketch_and_cmp<sketchtype>(inpaths, cms, kseqs, ofp, pairofp, sp, sketch_size, mincount, estim, jestim, cache_sketch, result_type, emit_fmt, presketched_only, nthreads, use_scientific, suffix, prefix, canon, entropy_minimization, spacing);
    switch(sketch_type) {
    case BB_MINHASH:
        CALL_DIST(mh::BBitMinHasher<uint64_t>); break;
    case HLL:
        CALL_DIST(hll::hll_t); break;
    case RANGE_MINHASH:
        CALL_DIST(mh::RangeMinHash<uint64_t>); break;
    case COUNTING_RANGE_MINHASH:
        CALL_DIST(mh::CountingRangeMinHash<uint64_t>); break;
    case BLOOM_FILTER:
        CALL_DIST(bf::bf_t); break;
    default: {
            char buf[128];
            std::sprintf(buf, "Sketch %s not yet supported.\n", (size_t(sketch_type) >= (sizeof(sketch_names) / sizeof(char *)) ? "Not such sketch": sketch_names[sketch_type]));
            RUNTIME_ERROR(buf);
    }
    }

    std::future<void> label_future;
    if(emit_fmt == BINARY) {
        if(pairofp_labels.empty()) pairofp_labels = "unspecified";
        label_future = std::async(std::launch::async, [&inpaths](const std::string &labels) {
            std::FILE *fp = std::fopen(labels.data(), "wb");
            if(fp == nullptr) RUNTIME_ERROR(std::string("Could not open file at ") + labels);
            for(const auto &path: inpaths) std::fwrite(path.data(), path.size(), 1, fp), std::fputc('\n', fp);
            std::fclose(fp);
        }, pairofp_labels);
    }
    if(pairofp != stdout) std::fclose(pairofp);
    if(label_future.valid()) label_future.get();
    return EXIT_SUCCESS;
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
#define PRINTMAT_INNER(type) \
        dm::DistanceMatrix<type> mat(argv[optind]);\
        LOG_DEBUG("Name of found: %s\n", dm::DistanceMatrix<type>::magic_string());\
        if((fp = std::fopen(outpath.data(), "wb")) == nullptr) RUNTIME_ERROR(ks::sprintf("Could not open file at %s", outpath.data()).data());\
        mat.printf(fp, use_scientific);
    try {
        PRINTMAT_INNER(float);
    } catch(const std::runtime_error &re) {
        PRINTMAT_INNER(double);
    }
    std::fclose(fp);
    return EXIT_SUCCESS;
}

int setdist_main(int argc, char *argv[]) {
    int wsz(0), k(31), use_scientific(false), co;
    bool canon(true);
    EmissionFormat emit_fmt = UT_TSV;
    EmissionType emit_type = JI;
    unsigned bufsize(BUFFER_FLUSH_SIZE);
    int nt(1);
    std::string spacing, paths_file;
    FILE *ofp(stdout), *pairofp(stdout);
    while((co = getopt(argc, argv, "F:c:p:o:O:S:B:k:s:lTfCMeZh?")) >= 0) {
        switch(co) {
            case 'B': std::stringstream(optarg) << bufsize; break;
            case 'k': k = std::atoi(optarg);                break;
            case 'p': nt = std::atoi(optarg);       break;
            case 's': spacing = optarg;             break;
            case 'C': canon = false;                break;
            case 'w': wsz = std::atoi(optarg);      break;
            case 'F': paths_file = optarg;          break;
            case 'o': ofp = fopen(optarg, "w");     break;
            case 'O': pairofp = fopen(optarg, "w"); break;
            case 'e': use_scientific = true;        break;
            case 'l': emit_type = FULL_MASH_DIST;   break;
            case 'M': emit_type = MASH_DIST;        break;
            case 'Z': emit_type = SIZES;            break;
            case 'U': emit_fmt = UPPER_TRIANGULAR;  break;
            case 'b': emit_fmt = BINARY;            break;
            case 'T': emit_fmt = FULL_TSV;          break;
            case 'h': case '?': dist_usage(*argv);
        }
    }
    omp_set_num_threads(nt);
    std::vector<char> rdbuf(bufsize);
    spvec_t sv(parse_spacing(spacing.data(), k));
    Spacer sp(k, wsz, sv);
    std::vector<std::string> inpaths(paths_file.size() ? get_paths(paths_file.data())
                                                       : std::vector<std::string>(argv + optind, argv + argc));
    KSeqBufferHolder h(nt);
    std::vector<khash_t(all)> hashes;
    while(hashes.size() < inpaths.size()) hashes.emplace_back(khash_t(all){0, 0, 0, 0, 0, 0, 0});
    const size_t nhashes(hashes.size());
    if(wsz < sp.c_) wsz = sp.c_;
    if(inpaths.size() == 0) {
        std::fprintf(stderr, "No paths. See usage.\n");
        dist_usage(*argv);
    }
    #pragma omp parallel for
    for(size_t i = 0; i < nhashes; ++i) {
        fill_set_genome<score::Lex>(inpaths[i].data(), sp, &hashes[i], i, nullptr, canon, h.data() + omp_get_thread_num());
    }
    LOG_DEBUG("Filled genomes. Now analyzing data.\n");
    ks::string str;
    str.sprintf("#Path\tCardinality (exact)\n");
    {
        const int fn(fileno(ofp));
        for(size_t i(0); i < nhashes; ++i) {
            str.sprintf("%s\t%zu\n", inpaths[i].data(), kh_size(&hashes[i]));
            if(str.size() > BUFFER_FLUSH_SIZE) str.flush(fn);
        }
        str.flush(fn);
    }
    // TODO: Emit overlaps and symmetric differences.
    if(ofp != stdout) std::fclose(ofp);
    std::vector<float> dists(nhashes - 1);
    str.clear();
    const double ksinv = 1./static_cast<double>(k);
    if((emit_fmt & BINARY) == 0) {
        if(emit_fmt == UPPER_TRIANGULAR) throw NotImplementedError("Not Implemented: upper triangular phylip for setdist.");
        str.sprintf("##Names\t");
        for(auto &path: inpaths) str.sprintf("%s\t", path.data());
        str.back() = '\n';
            str.write(fileno(pairofp)); str.free();
        setvbuf(pairofp, rdbuf.data(), _IOLBF, rdbuf.size());
        for(size_t i = 0; i < nhashes; ++i) {
            const khash_t(all) *h1(&hashes[i]);
            size_t j;
#define DO_LOOP(val) for(j = i + 1; j < nhashes; ++j) dists[j - i - 1] = (val)
            if(emit_type == JI) {
                #pragma omp parallel for
                DO_LOOP(jaccard_index(&hashes[j], h1));
            } else if(emit_type == MASH_DIST) {
                #pragma omp parallel for
                DO_LOOP(dist_index(jaccard_index(&hashes[j], h1), ksinv));
            } else if(emit_type == FULL_MASH_DIST) {
                #pragma omp parallel for
                DO_LOOP(full_dist_index(jaccard_index(&hashes[j], h1), ksinv));
            } else {
                #pragma omp parallel for
                DO_LOOP(union_size(&hashes[j], h1));
            }
#undef DO_LOOP
            submit_emit_dists<float>(fileno(pairofp), dists.data(), hashes.size(), i, str, inpaths, emit_fmt, use_scientific);
            std::free(h1->keys); std::free(h1->vals); std::free(h1->flags);
        }
    } else {
        dm::DistanceMatrix<float> dm(nhashes);
        for(size_t i = 0; i < nhashes; ++i) {
            auto span = dm.row_span(i);
            auto &dists = span.first;
            auto h1 = &hashes[i];
            size_t j;
#define DO_LOOP(val) for(j = i + 1; j < nhashes; ++j) dists[j - i - 1] = (val)
            if(emit_type == JI) {
                #pragma omp parallel for
                DO_LOOP(jaccard_index(&hashes[j], h1));
            } else if(emit_type == MASH_DIST) {
                #pragma omp parallel for
                DO_LOOP(dist_index(jaccard_index(&hashes[j], h1), ksinv));
            } else if(emit_type == FULL_MASH_DIST) {
                #pragma omp parallel for
                DO_LOOP(full_dist_index(jaccard_index(&hashes[j], h1), ksinv));
            } else {
                #pragma omp parallel for
                DO_LOOP(union_size(&hashes[j], h1));
            }
#undef DO_LOOP
            std::free(h1->keys); std::free(h1->vals); std::free(h1->flags);
        }
        if(emit_fmt == FULL_TSV) dm.printf(ofp, use_scientific, &inpaths);
        else {
            assert(emit_fmt == BINARY);
            dm.write(ofp);
        }
    }
    return EXIT_SUCCESS;
}

int hll_main(int argc, char *argv[]) {
    int c, wsz(0), k(31), num_threads(-1), sketch_size(24);
    bool canon(true);
    std::string spacing, paths_file;
    if(argc < 2) {
        usage: LOG_EXIT("Usage: %s <opts> <paths>\nFlags:\n"
                        "-k:\tkmer length (Default: 31. Max: 31)\n"
                        "-w:\twindow size (Default: -1)  Must be -1 (ignored) or >= kmer length.\n"
                        "-s:\tspacing (default: none). format: <value>x<times>,<value>x<times>,...\n"
                        "   \tOmitting x<times> indicates 1 occurrence of spacing <value>\n"
                        "-S:\tsketch size (default: 24). (Allocates 2 << [param] bytes of memory per HyperLogLog.\n"
                        "-p:\tnumber of threads.\n"
                        "-F:\tPath to file which contains one path per line\n"
                        , argv[0]);
    }
    while((c = getopt(argc, argv, "Cw:s:S:p:k:tfh?")) >= 0) {
        switch(c) {
            case 'C': canon = false; break;
            case 'h': case '?': goto usage;
            case 'k': k = std::atoi(optarg); break;
            case 'p': num_threads = std::atoi(optarg); break;
            case 's': spacing = optarg; break;
            case 'S': sketch_size = std::atoi(optarg); break;
            case 'w': wsz = std::atoi(optarg); break;
            case 'F': paths_file = optarg; break;
        }
    }
    if(wsz < k) wsz = k;
    std::vector<std::string> inpaths(paths_file.size() ? get_paths(paths_file.data())
                                                       : std::vector<std::string>(argv + optind, argv + argc));
    spvec_t sv(parse_spacing(spacing.data(), k));
    LOG_INFO("Processing %zu paths with %i threads\n", inpaths.size(), num_threads);
    const double est(estimate_cardinality<bns::score::Lex>(inpaths, k, wsz, sv, canon, nullptr, num_threads, sketch_size));
    std::fprintf(stdout, "Estimated number of unique exact matches: %lf\n", est);
    return EXIT_SUCCESS;
}

void union_usage [[noreturn]] (char *ex) {
    std::fprintf(stderr, "Usage: %s genome1 <genome2>...\n"
                         "Flags:\n"
                         "-o: Write union sketch to file [/dev/stdout]\n"
                         "-z: Emit compressed sketch\n",
                 ex);
    std::exit(1);
}

int union_main(int argc, char *argv[]) {
    if(std::find_if(argv, argc + argv, [](const char *s) {return std::strcmp(s, "--help") == 0;}) != argc + argv)
        union_usage(*argv);
    bool compress = false;
    int compression_level = 6;
    const char *opath = "/dev/stdout";
    std::vector<std::string> paths;
    for(int c;(c = getopt(argc, argv, "o:F:zZ:h?")) >= 0;) {
        switch(c) {
            case 'h': union_usage(*argv);
            case 'Z': compression_level = std::atoi(optarg); [[fallthrough]];
            case 'z': compress = true; break;
            case 'o': opath = optarg; break;
            case 'F': paths = get_paths(optarg); break;
        }
    }
    if(argc == optind && paths.empty()) union_usage(*argv);
    std::for_each(argv + optind, argv + argc, [&](const char *s){paths.emplace_back(s);});
    hll::hll_t hll(paths.back());
    paths.pop_back();
    for(auto &path: paths) {
        hll += hll_t(path);
    }
    char mode[6];
    if(compress && compression_level)
        std::sprintf(mode, "wb%d", compression_level % 23);
    else
        std::sprintf(mode, "wT");
    gzFile ofp = gzopen(opath, mode);
    if(!ofp) throw std::runtime_error(std::string("Could not open file at ") + opath);
    hll.write(ofp);
    gzclose(ofp);
    return 0;
}

int partdist_usage() {
    std::fprintf(stderr, "%s partdist <opts> <args>...\n[This command has not been implemented.]\n", bns::executable);
    return EXIT_FAILURE;
}

int partdist_main(int argc, char *argv[]) {
    if(argc == 0) return partdist_usage();
    int c;
    while((c = getopt(argc, argv, "h?")) >= 0) {
        switch(c) {
            // WHEEEEEEEEE
        }
    }
    std::vector<std::string> qpaths = get_lines(argv[optind]), refpaths = get_lines(argv[optind + 1]);
    throw NotImplementedError("partdist_main has not been implemented.");
}

} // namespace bns

using namespace bns;


int main(int argc, char *argv[]) {
    bns::executable = argv[0];
    if(argc == 1) main_usage(argv);
    if(std::strcmp(argv[1], "sketch") == 0) return sketch_main(argc - 1, argv + 1);
    else if(std::strcmp(argv[1], "dist") == 0) return dist_main(argc - 1, argv + 1);
    else if(std::strcmp(argv[1], "union") == 0) return union_main(argc - 1, argv + 1);
    else if(std::strcmp(argv[1], "setdist") == 0) return setdist_main(argc - 1, argv + 1);
    else if(std::strcmp(argv[1], "hll") == 0) return hll_main(argc - 1, argv + 1);
    else if(std::strcmp(argv[1], "printmat") == 0) return print_binary_main(argc - 1, argv + 1);
    else if(std::strcmp(argv[1], "partdist") == 0) return partdist_main(argc - 1, argv + 1);
    else {
        for(const char *const *p(argv + 1); *p; ++p)
            if(std::string(*p) == "-h" || std::string(*p) == "--help") main_usage(argv);
        std::fprintf(stderr, "Usage: %s <subcommand> [options...]. Use %s <subcommand> for more options. [Subcommands: sketch, dist, setdist, hll, union, printmat.]\n",
                     *argv, *argv);
        RUNTIME_ERROR(std::string("Invalid subcommand ") + argv[1] + " provided.");
    }
}
