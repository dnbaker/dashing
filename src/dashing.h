#ifndef DASHING_H__
#define DASHING_H__
#include <omp.h>
#include "sketch/bbmh.h"
#include "sketch/mh.h"
#include "sketch/mult.h"
#include "sketch/hk.h"
#include "khset/khset.h"
#include "bonsai/util.h"
#include "bonsai/database.h"
#include "bonsai/bitmap.h"
#include "bonsai/setcmp.h"
#include "khset/khset.h"
#include "distmat/distmat.h"
#include <sstream>
#include "getopt.h"
#include <sys/stat.h>
#include "substrs.h"
#include "khset64.h"
#include "enums.h"

#if __cplusplus >= 201703L && __cpp_lib_execution
#include <execution>
#endif
#ifndef BUFFER_FLUSH_SIZE
#define BUFFER_FLUSH_SIZE (1u << 18)
#endif

#define LO_ARG(LONG, SHORT) {LONG, required_argument, 0, SHORT},
#define LO_NO(LONG, SHORT) {LONG, no_argument, 0, SHORT},
#define LO_FLAG(LONG, SHORT, VAR, VAL) {LONG, no_argument, (int *)&VAR, VAL},
#define DIST_LONG_OPTS \
static option_struct dist_long_options[] = {\
    LO_FLAG("avoid-sorting", 'n', avoid_fsorting, true)\
    LO_FLAG("by-entropy", 'g', entropy_minimization, true) \
    LO_FLAG("cache-sketches", 'W', cache_sketch, true)\
    LO_FLAG("countmin", 'y', sm, CBF)\
    LO_FLAG("emit-binary", 'b', emit_fmt, BINARY)\
    LO_FLAG("full-mash-dist", 'l', result_type, FULL_MASH_DIST)\
    LO_FLAG("full-tsv", 'T', emit_fmt, FULL_TSV)\
    LO_FLAG("no-canon", 'C', canon, false)\
    LO_FLAG("phylip", 'U', emit_fmt, UPPER_TRIANGULAR)\
    LO_FLAG("presketched", 'H', presketched_only, true)\
    LO_FLAG("sizes", 'Z', result_type, SIZES)\
    LO_FLAG("sketch-by-fname", 'N', sm, BY_FNAME)\
    LO_FLAG("use-bb-minhash", '8', sketch_type, BB_MINHASH)\
    LO_FLAG("use-scientific", 'e', use_scientific, true)\
    LO_ARG("bbits", 'B')\
    LO_ARG("cm-sketch-size", 't')\
    LO_ARG("ertl-joint-mle", 'J')\
    LO_ARG("ertl-mle", 'm')\
    LO_ARG("improved", 'I')\
    LO_ARG("kmer-length", 'k')\
    LO_ARG("min-count", 'c')\
    LO_ARG("nhashes", 'q')\
    LO_ARG("nthreads", 'p')\
    LO_ARG("original", 'E')\
    LO_ARG("out-dists", 'O') \
    LO_ARG("out-sizes", 'o') \
    LO_ARG("paths", 'F')\
    LO_ARG("prefix", 'P')\
    LO_ARG("query-paths", 'Q') \
    LO_ARG("seed", 'R')\
    LO_ARG("sketch-size", 'S')\
    LO_ARG("spacing", 's')\
    LO_ARG("suffix", 'x')\
    LO_ARG("window-size", 'w')\
    LO_ARG("help", 'h')\
    LO_ARG("mkdist", 1337)\
    LO_FLAG("use-range-minhash", 128, sketch_type, RANGE_MINHASH)\
    LO_FLAG("use-counting-range-minhash", 129, sketch_type, COUNTING_RANGE_MINHASH)\
    LO_FLAG("use-full-khash-sets", 130, sketch_type, FULL_KHASH_SET)\
    LO_FLAG("full-containment-dist", 133, result_type, FULL_CONTAINMENT_DIST) \
    LO_FLAG("use-bloom-filter", 134, sketch_type, BLOOM_FILTER)\
    LO_FLAG("use-super-minhash", 135, sketch_type, BB_SUPERMINHASH)\
    LO_FLAG("use-nthash", 136, enct, NTHASH)\
    LO_FLAG("containment-index", 131, result_type, CONTAINMENT_INDEX) \
    LO_FLAG("containment-dist", 132, result_type, CONTAINMENT_DIST) \
    LO_FLAG("mash-dist", 'M', result_type, MASH_DIST)\
    LO_FLAG("symmetric-containment-index", 137, result_type, SYMMETRIC_CONTAINMENT_INDEX) \
    LO_FLAG("symmetric-containment-dist", 138, result_type, SYMMETRIC_CONTAINMENT_DIST) \
    LO_FLAG("use-cyclic-hash", 139, enct, NTHASH)\
    LO_ARG("wj-cm-sketch-size", 140)\
    LO_ARG("wj-cm-nhashes", 141)\
    LO_FLAG("wj", 142, weighted_jaccard, true)\
    LO_ARG("nearest-neighbors", 143)\
    LO_FLAG("use-wide-hll", 144, sketch_type, WIDE_HLL) \
    {0,0,0,0}\
};



namespace bns {
int flatten_all(const std::vector<std::string> &fpaths, const std::string outpath, std::vector<unsigned> &k_values);
namespace detail {void sort_paths_by_fsize(std::vector<std::string> &paths);}
size_t posix_fsizes(const std::string &path, const char sep=FNAME_SEP);
using namespace sketch;
using namespace hll;
using option_struct = struct option;
using sketch::common::WangHash;
using CRMFinal = mh::FinalCRMinHash<uint64_t, uint32_t>;
using RMFinal = mh::FinalRMinHash<uint64_t, sketch::common::Allocator<uint64_t>>;
template<typename BaseHash>
struct SeededHash {
    BaseHash wh_;
    const uint64_t seed_;
    SeededHash(uint64_t seed): seed_(seed) {}
    uint64_t operator()(uint64_t x) const {return wh_(x ^ seed_);}
};

using CountingSketch = sketch::hk::HeavyKeeper<6, 10, SeededHash<sketch::common::WangHash>>;
template<typename T> inline double similarity(const T &a, const T &b) {
    return a.jaccard_index(b);
}

template<> inline double similarity<CRMFinal>(const CRMFinal &a, const CRMFinal &b) {
    return a.histogram_intersection(b);
}

template<typename T>
inline void sketch_finalize(T &x) {}
template<> inline void sketch_finalize<khset64_t>(khset64_t &x) {x.cvt2shs();}
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
typename std::common_type<FType1, FType2>::type containment_dist(FType1 containment, FType2 ksinv) {
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
typename std::common_type<FType1, FType2>::type full_containment_dist(FType1 containment, FType2 ksinv) {
    return 1. - std::pow(containment, ksinv);
}

template<typename T>
inline double containment_index(const T &a, const T &b) {
    return a.containment_index(b);
}
template<typename T>
inline std::array<double, 3> set_triple(const T &a, const T &b) {
    return a.full_set_comparison(b);
}
template<typename T>
inline std::array<double, 3> set_triple(const wj::WeightedSketcher<T> &a, const wj::WeightedSketcher<T> &b) {
    RUNTIME_ERROR("This should only be called on the finalized sketches");
}

#define CONTAIN_OVERLOAD_FAIL(x)\
template<>\
inline double containment_index<x>(const x &b, const x &a) {\
    RUNTIME_ERROR(std::string("Containment index not implemented for ") + __PRETTY_FUNCTION__);\
}\
template<>\
inline std::array<double, 3> set_triple<x>(const x &b, const x &a) {\
    RUNTIME_ERROR(std::string("set_triple not implemented for ") + __PRETTY_FUNCTION__);\
}
CONTAIN_OVERLOAD_FAIL(RMFinal)
CONTAIN_OVERLOAD_FAIL(bf::bf_t)
CONTAIN_OVERLOAD_FAIL(wj::WeightedSketcher<RMFinal>)
CONTAIN_OVERLOAD_FAIL(wj::WeightedSketcher<bf::bf_t>)
CONTAIN_OVERLOAD_FAIL(CRMFinal)
#undef CONTAIN_OVERLOAD_FAIL

using CBBMinHashType = mh::CountingBBitMinHasher<uint64_t, uint16_t>; // Is counting to 65536 enough for a transcriptome?
using SuperMinHashType = mh::SuperMinHash<>;



enum Sketch: int {
    HLL,
    BLOOM_FILTER,
    RANGE_MINHASH,
    FULL_KHASH_SET,
    COUNTING_RANGE_MINHASH,
    BB_MINHASH,
    BB_SUPERMINHASH,
    COUNTING_BB_MINHASH, // TODO make this work.
    WIDE_HLL
};
static constexpr const char *const sketch_names [] {
    "HLL/HyperLogLog",
    "BF/BloomFilter",
    "RMH/Range Min-Hash/KMV",
    "FHS/Full Hash Set",
    "CRHM/Counting Range Minhash",
    "BB/B-bit Minhash",
    "BBS/B-bit SuperMinHash",
    "CBB/Counting B-bit Minhash",
    "WHLL/Wide HLL",
};



struct GlobalArgs {
    size_t weighted_jaccard_cmsize = 22;
    size_t weighted_jaccard_nhashes = 8;
    uint32_t bbnbits = 16;
    uint32_t number_neighbors = 0;
    void show() const {
        std::fprintf(stderr, "Global Arguments: %zu wjcm, %zu wjnh, %u bbits %u nn\n", weighted_jaccard_cmsize, weighted_jaccard_nhashes, bbnbits, number_neighbors);
    }
};

extern GlobalArgs gargs;
extern uint64_t global_hash_seed;


INLINE static constexpr
NNType emt2nntype(EmissionType result_type) {
    switch(result_type) {
        case MASH_DIST: case FULL_MASH_DIST: case CONTAINMENT_DIST:
        case FULL_CONTAINMENT_DIST: case SYMMETRIC_CONTAINMENT_DIST:
            return DIST_MEASURE;
        case JI: case SIZES: case CONTAINMENT_INDEX: case SYMMETRIC_CONTAINMENT_INDEX:
        default:
            return SIMILARITY_MEASURE;
    }
    __builtin_unreachable();
    return SIMILARITY_MEASURE;
}

static const char *emt2str(EmissionType result_type) {
    switch(result_type) {
        case MASH_DIST: return "MASH_DIST";
        case JI: return "JI";
        case SIZES: return "SIZES";
        case FULL_MASH_DIST: return "FULL_MASH_DIST";
        case FULL_CONTAINMENT_DIST: return "FULL_CONTAINMENT_DIST";
        case CONTAINMENT_INDEX: return "CONTAINMENT_INDEX";
        case CONTAINMENT_DIST: return "CONTAINMENT_DIST";
        case SYMMETRIC_CONTAINMENT_INDEX: return "SYMMETRIC_CONTAINMENT_INDEX";
        case SYMMETRIC_CONTAINMENT_DIST: return "SYMMETRIC_CONTAINMENT_DIST";
        default: break;
    }
    return "ILLEGAL_EMISSION_FMT";
}

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
FINAL_OVERLOAD(WideHyperLogLogHasher<>);
FINAL_OVERLOAD(SuperMinHashType);
FINAL_OVERLOAD(CBBMinHashType);
FINAL_OVERLOAD(wj::WeightedSketcher<mh::CountingRangeMinHash<uint64_t>>);
FINAL_OVERLOAD(wj::WeightedSketcher<mh::RangeMinHash<uint64_t>>);
FINAL_OVERLOAD(wj::WeightedSketcher<mh::BBitMinHasher<uint64_t>>);
FINAL_OVERLOAD(wj::WeightedSketcher<SuperMinHashType>);
FINAL_OVERLOAD(wj::WeightedSketcher<WideHyperLogLogHasher<>>);
FINAL_OVERLOAD(wj::WeightedSketcher<CBBMinHashType>);
template<typename T>struct SketchFileSuffix {static constexpr const char *suffix = ".sketch";};
#define SSS(type, suf) template<> struct SketchFileSuffix<type> {static constexpr const char *suffix = suf;}
SSS(mh::CountingRangeMinHash<uint64_t>, ".crmh");
SSS(mh::RangeMinHash<uint64_t>, ".rmh");
SSS(khset64_t, ".khs");
SSS(bf::bf_t, ".bf");
SSS(mh::BBitMinHasher<uint64_t>, ".bmh");
SSS(SuperMinHashType, ".bbs");
SSS(CBBMinHashType, ".cbmh");
SSS(mh::HyperMinHash<uint64_t>, ".hmh");
SSS(hll::hll_t, ".hll");
#undef SSS
#undef FINAL_OVERLOAD


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

} // namespace detail

static constexpr bool is_symmetric(EmissionType result_type) {
    switch(result_type) {
        case MASH_DIST: case JI: case SIZES:
        case FULL_MASH_DIST: case SYMMETRIC_CONTAINMENT_INDEX: case SYMMETRIC_CONTAINMENT_DIST:
            return true;

        case CONTAINMENT_INDEX: case CONTAINMENT_DIST: case FULL_CONTAINMENT_DIST:
        default: break;
    }
    return false;
}



template<typename T> struct SketchEnum;
template<> struct SketchEnum<hll::hll_t> {static constexpr Sketch value = HLL;};
template<> struct SketchEnum<bf::bf_t> {static constexpr Sketch value = BLOOM_FILTER;};
template<> struct SketchEnum<mh::RangeMinHash<uint64_t>> {static constexpr Sketch value = RANGE_MINHASH;};
template<> struct SketchEnum<mh::CountingRangeMinHash<uint64_t>> {static constexpr Sketch value = COUNTING_RANGE_MINHASH;};
template<> struct SketchEnum<mh::BBitMinHasher<uint64_t>> {static constexpr Sketch value = BB_MINHASH;};
template<> struct SketchEnum<CBBMinHashType> {static constexpr Sketch value = COUNTING_BB_MINHASH;};
template<> struct SketchEnum<khset64_t> {static constexpr Sketch value = FULL_KHASH_SET;};
template<> struct SketchEnum<SuperMinHashType> {static constexpr Sketch value = BB_SUPERMINHASH;};
template<> struct SketchEnum<WideHyperLogLogHasher<>> {static constexpr Sketch value = WIDE_HLL;};
template<> struct SketchEnum<wj::WeightedSketcher<hll::hll_t>> {static constexpr Sketch value = HLL;};
template<> struct SketchEnum<wj::WeightedSketcher<bf::bf_t>> {static constexpr Sketch value = BLOOM_FILTER;};
template<> struct SketchEnum<wj::WeightedSketcher<mh::RangeMinHash<uint64_t>>> {static constexpr Sketch value = RANGE_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<mh::CountingRangeMinHash<uint64_t>>> {static constexpr Sketch value = COUNTING_RANGE_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<mh::BBitMinHasher<uint64_t>>> {static constexpr Sketch value = BB_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<CBBMinHashType>> {static constexpr Sketch value = COUNTING_BB_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<khset64_t>> {static constexpr Sketch value = FULL_KHASH_SET;};
template<> struct SketchEnum<wj::WeightedSketcher<SuperMinHashType>> {static constexpr Sketch value = BB_SUPERMINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<WideHyperLogLogHasher<>>> {static constexpr Sketch value = WIDE_HLL;};

template<typename T>
inline void set_estim_and_jestim(T &x, hll::EstimationMethod estim, hll::JointEstimationMethod jestim) {}

template<typename Hashstruct>
inline void set_estim_and_jestim(hll::hllbase_t<Hashstruct> &h, hll::EstimationMethod estim, hll::JointEstimationMethod jestim) {
    h.set_estim(estim);
    h.set_jestim(jestim);
}
template<typename T> inline T construct(size_t ssarg);
template<typename T, bool is_weighted>
struct Constructor;
template<typename T> struct Constructor<T, false> {
    static auto create(size_t ssarg) {
        return T(ssarg);
    }
};
template<typename T> struct Constructor<T, true> {
    static auto create(size_t ssarg) {
        using base_type = typename T::base_type;
        using cm_type = typename T::cm_type;
        return T(cm_type(16, gargs.weighted_jaccard_cmsize, gargs.weighted_jaccard_nhashes), construct<base_type>(ssarg));
    }
};

template<typename T>
inline T construct(size_t ssarg) {
    Constructor<T, wj::is_weighted_sketch<T>::value> constructor;
    return constructor.create(ssarg);
}


template<typename T>
inline double cardinality_estimate(T &x) {
    return x.cardinality_estimate();
}
template<> inline double cardinality_estimate(hll::hll_t &x) {return x.report();}
template<> inline double cardinality_estimate(mh::FinalBBitMinHash &x) {return x.est_cardinality_;}
template<> inline double cardinality_estimate(mh::FinalDivBBitMinHash &x) {return x.est_cardinality_;}
//extern template double cardinality_estimate(hll::hll_t &x);
//extern template double cardinality_estimate(mh::FinalBBitMinHash &x);
//extern template double cardinality_estimate(mh::FinalDivBBitMinHash &x);

template<typename SketchType>
static inline std::string make_fname(const char *path, size_t sketch_p, int wsz, int k, int csz, const std::string &spacing,
                       const std::string &suffix="", const std::string &prefix="",
                       EncodingType enct=BONSAI) {
    std::string ret(prefix);
    if(ret.size()) ret += '/';
    {
        const char *p, *p2;
		p = (p = std::strchr(path, FNAME_SEP)) ? p + 1: path;
        if(ret.size() && (p2 = strrchr(p, '/'))) ret += std::string(p2 + 1);
        else                                     ret += p;
    }
    ret += ".w";
    ret + std::to_string(std::max(csz, wsz));
    ret += ".";
    ret += std::to_string(k);
    ret += ".spacing";
    ret += spacing;
    ret += '.';
    ret += enct == BONSAI ? "": enct == NTHASH ? "nt.": "cyclic.";
    // default: no note. cyclic or nthash otherwise.
    if(suffix.size()) {
        ret += "suf";
        ret += suffix;
        ret += '.';
    }
    ret += std::to_string(sketch_p);
    ret += SketchFileSuffix<SketchType>::suffix;
    return ret;
}

namespace us {
template<typename T> inline double union_size(const T &a, const T &b) {
    throw NotImplementedError(std::string("union_size not available for type ") + __PRETTY_FUNCTION__);
}
template<typename T> inline double intersection_size(const T &a, const T &b) {
    throw NotImplementedError(std::string("intersection_size not available for type ") + __PRETTY_FUNCTION__);
}


#define US_DEC(type) \
template<> inline double union_size<type> (const type &a, const type &b) { \
    return a.union_size(b); \
} \
template<> inline double intersection_size<type> (const type &a, const type &b) { \
    return a.intersection_size(b); \
}

US_DEC(RMFinal)
US_DEC(CRMFinal)
US_DEC(khset64_t)
#undef US_DEC
template<> inline double union_size<hll::hllbase_t<>>(const hll::hllbase_t<> &h1, const hll::hllbase_t<> &h2) {
    return h1.union_size(h2);
}
template<> inline double intersection_size<hll::hllbase_t<>>(const hll::hllbase_t<> &h1, const hll::hllbase_t<> &h2) {
    return std::max(0., h1.creport() + h2.creport() - h1.union_size(h2));
}
template<> inline double union_size<mh::FinalBBitMinHash> (const mh::FinalBBitMinHash &a, const mh::FinalBBitMinHash &b) {
    return (a.est_cardinality_ + b.est_cardinality_ ) / (1. + a.jaccard_index(b));
}
template<> inline double intersection_size<mh::FinalBBitMinHash> (const mh::FinalBBitMinHash &a, const mh::FinalBBitMinHash &b) {
    const double ji = a.jaccard_index(b);
    return ji * (a.est_cardinality_ + b.est_cardinality_ ) / (1. + ji);
}
} // namespace us

template<typename T>
inline auto symmetric_containment_func(const T &x, const T &y) {
    auto tmp = set_triple(x, y);
    return tmp[2] / (std::min(tmp[0], tmp[1]) + tmp[2]);
}

template<typename SketchType>
void partdist_loop(std::FILE *ofp, SketchType *hlls, const std::vector<std::string> &inpaths, const bool use_scientific, const unsigned k, const EmissionType result_type, EmissionFormat emit_fmt, const size_t buffer_flush_size,
                   size_t nq)
{
    const float ksinv = 1./ k;
    if(nq >= inpaths.size()) {
        RUNTIME_ERROR(ks::sprintf("Wrong number of query/references. (ip size: %zu, nq: %zu\n", inpaths.size(), nq).data());
    }
    size_t nr = inpaths.size() - nq;
    float *arr = static_cast<float *>(std::malloc(nr * nq * sizeof(float)));
    if(!arr) throw std::bad_alloc();
#if TIMING
    auto start = std::chrono::high_resolution_clock::now();
#endif
    std::future<void> write_future, fmt_future;
    std::array<ks::string, 2> buffers;
    for(auto &b: buffers) b.resize(4 * nr);
    for(size_t qi = nr; qi < inpaths.size(); ++qi) {
        size_t qind =  qi - nr;
        switch(result_type) {


#define dist_sim(x, y) dist_index(similarity(x, y), ksinv)
#define fulldist_sim(x, y) full_dist_index(similarity(x, y), ksinv)
#define fullcont_sim(x, y) full_containment_dist(containment_index(x, y), ksinv)
#define cont_sim(x, y) containment_dist(containment_index(x, y), ksinv)
#define sym_cont_dist(x, y) containment_dist(symmetric_containment_func(x, y), ksinv)
#define DO_LOOP(func)\
                for(size_t j = 0; j < nr; ++j) {\
                    arr[qind * nr + j] = func(hlls[j], hlls[qi]);\
                }
            case MASH_DIST:
                #pragma omp parallel for schedule(dynamic)
                DO_LOOP(dist_sim);
                break;
            case FULL_MASH_DIST:
                #pragma omp parallel for schedule(dynamic)
                DO_LOOP(fulldist_sim);
                break;
            case JI:
                #pragma omp parallel for schedule(dynamic)
                DO_LOOP(similarity);
                break;
            case SIZES:
                #pragma omp parallel for schedule(dynamic)
                DO_LOOP(us::intersection_size);
                break;
            case CONTAINMENT_INDEX:
                #pragma omp parallel for schedule(dynamic)
                DO_LOOP(containment_index);
                break;
            case CONTAINMENT_DIST:
                #pragma omp parallel for schedule(dynamic)
                DO_LOOP(cont_sim);
                break;
            case FULL_CONTAINMENT_DIST:
                #pragma omp parallel for schedule(dynamic)
                DO_LOOP(fullcont_sim);
                break;
            case SYMMETRIC_CONTAINMENT_INDEX:
                #pragma omp parallel for schedule(dynamic)
                DO_LOOP(symmetric_containment_func)
                break;
            case SYMMETRIC_CONTAINMENT_DIST:
                #pragma omp parallel for schedule(dynamic)
                DO_LOOP(sym_cont_dist)
                break;
            default: throw std::runtime_error("Value not found");
#undef DO_LOOP
//#undef dist_sim
//#undef cont_sim
//#undef fulldist_sim
//#undef fullcont_sim
        }
        switch(emit_fmt) {
            case BINARY:
                if(write_future.valid()) write_future.get();
                write_future = std::async(std::launch::async, [ptr=arr + (qi - nr) * nq, nb=sizeof(float) * nr](const int fn) {
                    if(unlikely(::write(fn, ptr, nb) != ssize_t(nb))) RUNTIME_ERROR("Error writing to binary file");
                }, ::fileno(ofp));
                break;
            case UT_TSV: case UPPER_TRIANGULAR: default:
                // RUNTIME_ERROR(std::string("Illegal output format. numeric: ") + std::to_string(int(emit_fmt)));
            case FULL_TSV:
                if(fmt_future.valid()) fmt_future.get();
                fmt_future = std::async(std::launch::async, [nr,qi,ofp,ind=qi-nr,&inpaths,use_scientific,arr,&buffers,&write_future]() {
                    auto &buffer = buffers[qi & 1];
                    buffer += inpaths[qi];
                    const char *fmt = use_scientific ? "\t%e": "\t%f";
                    float *aptr = arr + ind * nr;
                    for(size_t i = 0; i < nr; ++i) {
                        buffer.sprintf(fmt, aptr[i]);
                    }
                    buffer.putc_('\n');
                    if(write_future.valid()) write_future.get();
                    write_future = std::async(std::launch::async, [ofp,&buffer]() {buffer.flush(::fileno(ofp));});
                });
                break;
        }
    }
    if(fmt_future.valid()) fmt_future.get();
    if(write_future.valid()) write_future.get();
#if TIMING
    auto end = std::chrono::high_resolution_clock::now();
#endif
    std::free(arr);
}

static const char *executable = nullptr;
static std::string get_executable() {
    return executable ? std::string(executable): "unspecified";
}
void main_usage(char **argv);
void dist_usage(const char *arg);
void sketch_usage(const char *arg);
void sketch_by_seq_usage(const char *arg);
void flatten_usage();
void union_usage [[noreturn]] (char *ex);
int sketch_main(int argc, char *argv[]);
int dist_main(int argc, char *argv[]);
int print_binary_main(int argc, char *argv[]);
int mkdist_main(int argc, char *argv[]);
int flatten_main(int argc, char *argv[]);
int setdist_main(int argc, char *argv[]);
int hll_main(int argc, char *argv[]);
int union_main(int argc, char *argv[]);
int view_main(int argc, char *argv[]);
int sketch_by_seq_main(int argc, char *argv[]);
int dist_by_seq_main(int argc, char *argv[]);
}

#endif /* DASHING_H__ */
