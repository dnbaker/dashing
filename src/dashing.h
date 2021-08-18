#ifndef DASHING_H__
#define DASHING_H__
#include <omp.h>
#include "sketch/bbmh.h"
#include "sketch/mh.h"
#include "sketch/hmh.h"
#include "sketch/mult.h"
#include "sketch/hk.h"
#include "sketch/bf.h"
#include "bonsai/encoder.h"
#include "khset/khset.h"
#include "distmat/distmat.h"
#include <sstream>
#include "getopt.h"
#include <sys/stat.h>
#include "substrs.h"
#include "khset64.h"
#include "enums.h"

#ifndef _OPENMP
#error("Need OpenMP")
#endif

#if __cplusplus >= 201703L && __cpp_lib_execution
#include <execution>
#endif
#ifndef BUFFER_FLUSH_SIZE
#define BUFFER_FLUSH_SIZE (1u << 18)
#endif

#define LO_ARG(LONG, SHORT) {LONG, required_argument, 0, SHORT},
#define LO_NO(LONG, SHORT) {LONG, no_argument, 0, SHORT},
#define LO_FLAG(LONG, SHORT, VAR, VAL) {LONG, no_argument, (int *)&VAR, VAL},

#define SHARED_OPTS \
    LO_FLAG("wj-exact", 145, gargs.exact_weighted, true)\
    LO_FLAG("use-wide-hll", 144, sketch_type, WIDE_HLL) \
    LO_FLAG("defer-hll", 146, gargs.defer_hll_creation, true)\
    /*LO_FLAG("use-hyperminhash", 147, sketch_type, HYPERMINHASH)*/\


using BKHash64 = sketch::minhash::BottomKHasher<sketch::WangHash, uint64_t>;

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
    /*LO_ARG("mkdist", 1337)*/\
    LO_FLAG("use-range-minhash", 128, sketch_type, RANGE_MINHASH)\
    LO_FLAG("use-full-khash-sets", 130, sketch_type, FULL_KHASH_SET)\
    LO_FLAG("use-full-hash-sets", 1000, sketch_type, FULL_KHASH_SET)\
    LO_FLAG("use-hash-sets", 1000, sketch_type, FULL_KHASH_SET)\
    LO_FLAG("hash-sets", 10001, sketch_type, FULL_KHASH_SET)\
    LO_FLAG("use-full-sets", 1000, sketch_type, FULL_KHASH_SET)\
    LO_FLAG("full-containment-dist", 133, result_type, FULL_CONTAINMENT_DIST) \
    LO_FLAG("use-bloom-filter", 134, sketch_type, BLOOM_FILTER)\
    LO_FLAG("use-nthash", 136, enct, NTHASH)\
    LO_FLAG("containment-index", 131, result_type, CONTAINMENT_INDEX) \
    LO_FLAG("containment-dist", 132, result_type, CONTAINMENT_DIST) \
    LO_FLAG("mash-dist", 'M', result_type, MASH_DIST)\
    LO_FLAG("symmetric-containment-index", 137, result_type, SYMMETRIC_CONTAINMENT_INDEX) \
    LO_FLAG("symmetric-containment-dist", 138, result_type, SYMMETRIC_CONTAINMENT_DIST) \
    LO_FLAG("use-cyclic-hash", 139, enct, CYCLIC)\
    LO_ARG("wj-cm-sketch-size", 140)\
    LO_ARG("wj-cm-nhashes", 141)\
    LO_FLAG("wj", 142, weighted_jaccard, true)\
    LO_ARG("nearest-neighbors", 143)\
    SHARED_OPTS \
    LO_ARG("nperbatch", 148)\
    {0,0,0,0}\
};



namespace bns {
struct khset64_t;
using sketch::mh::HyperLogLogHasher;
using sketch::HyperMinHash;
int flatten_all(const std::vector<std::string> &fpaths, const std::string outpath, std::vector<unsigned> &k_values);
namespace detail {void sort_paths_by_fsize(std::vector<std::string> &paths);}
size_t posix_fsizes(const std::string &path, const char sep=FNAME_SEP);
using namespace sketch;
using namespace hll;
using sketch::BBitMinHasher;
using option_struct = struct option;
using sketch::WangHash;
using CRMFinal = mh::FinalCRMinHash<uint64_t, uint32_t>;
using RMFinal = mh::FinalRMinHash<uint64_t, sketch::common::Allocator<uint64_t>>;
template<typename BaseHash>
struct SeededHash {
    BaseHash wh_;
    const uint64_t seed_;
    SeededHash(uint64_t seed): seed_(seed) {}
    uint64_t operator()(uint64_t x) const {return wh_(x ^ seed_);}
};


#if DASHING_USE_HK
#define DASHING_COUNTING_SKETCH ::sketch::hk::HeavyKeeper<6, 10, SeededHash<sketch::common::WangHash>>
#else
#define DASHING_COUNTING_SKETCH ::sketch::ccm_t
#endif
using CountingSketch = DASHING_COUNTING_SKETCH;

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
    return ji ? -std::log(2. * ji / (1. + ji)) * ksinv: 1.;
}

template<typename FType1, typename FType2,
         typename=typename std::enable_if<
            std::is_floating_point<FType1>::value && std::is_floating_point<FType2>::value
          >::type
         >
typename std::common_type<FType1, FType2>::type containment_dist(FType1 containment, FType2 ksinv) {
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
    UNRECOVERABLE_ERROR("This should only be called on the finalized sketches");
}

#define CONTAIN_OVERLOAD_FAIL(x)\
template<>\
inline double containment_index<x>(const x &b, const x &a) {\
    UNRECOVERABLE_ERROR(std::string("Containment index not implemented for ") + __PRETTY_FUNCTION__);\
}\
template<>\
inline std::array<double, 3> set_triple<x>(const x &b, const x &a) {\
    UNRECOVERABLE_ERROR(std::string("set_triple not implemented for ") + __PRETTY_FUNCTION__);\
}
CONTAIN_OVERLOAD_FAIL(RMFinal)
CONTAIN_OVERLOAD_FAIL(bf::bf_t)
CONTAIN_OVERLOAD_FAIL(wj::WeightedSketcher<RMFinal>)
CONTAIN_OVERLOAD_FAIL(wj::WeightedSketcher<bf::bf_t>)
CONTAIN_OVERLOAD_FAIL(CRMFinal)
using wjRMFinal = wj::WeightedSketcher<RMFinal, wj::ExactCountingAdapter>;
using wjBFFinal = wj::WeightedSketcher<bf::bf_t, wj::ExactCountingAdapter>;
CONTAIN_OVERLOAD_FAIL(wjRMFinal)
CONTAIN_OVERLOAD_FAIL(wjBFFinal)
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
    COUNTING_BB_MINHASH,
    WIDE_HLL,
    HYPERMINHASH,
    HMH = HYPERMINHASH
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
    "HMH/HyperMinHash"
};



struct GlobalArgs {
    uint32_t weighted_jaccard_cmsize = 22;
    uint32_t weighted_jaccard_nhashes = 10;
    uint32_t bbnbits = 16;
    uint32_t number_neighbors = 0; // set to 0 signifies that the option is not activated.
    size_t nperbatch = 16;
    bool exact_weighted = false;
    bool defer_hll_creation = false;
    void show() const {
        std::fprintf(stderr, "Global Arguments: %u wjcm, %u wjnh, %u bbits %u nn\n", weighted_jaccard_cmsize, weighted_jaccard_nhashes, bbnbits, number_neighbors);
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
#define FINAL_OVERLOAD(x) \
template<> struct FinalSketch<x> { \
    using final_type = typename x::final_type;}
#define FINAL_OVERLOAD2(x, y) \
template<> struct FinalSketch<x, y> { \
    using final_type = typename x, y::final_type;}
FINAL_OVERLOAD(mh::CountingRangeMinHash<uint64_t>);
FINAL_OVERLOAD(mh::RangeMinHash<uint64_t>);
FINAL_OVERLOAD(BKHash64);
FINAL_OVERLOAD(mh::BBitMinHasher<uint64_t>);
FINAL_OVERLOAD(WideHyperLogLogHasher<>);
FINAL_OVERLOAD(HyperLogLogHasher<>);
FINAL_OVERLOAD(SuperMinHashType);
FINAL_OVERLOAD(CBBMinHashType);
FINAL_OVERLOAD(bf::bf_t);
FINAL_OVERLOAD(hll::hll_t);
FINAL_OVERLOAD(khset64_t);
FINAL_OVERLOAD(HyperMinHash);
FINAL_OVERLOAD(wj::WeightedSketcher<bf::bf_t>);
FINAL_OVERLOAD(wj::WeightedSketcher<hll::hll_t>);
FINAL_OVERLOAD(wj::WeightedSketcher<khset64_t>);
FINAL_OVERLOAD(wj::WeightedSketcher<mh::CountingRangeMinHash<uint64_t>>);
FINAL_OVERLOAD(wj::WeightedSketcher<mh::RangeMinHash<uint64_t>>);
FINAL_OVERLOAD(wj::WeightedSketcher<BKHash64>);
FINAL_OVERLOAD(wj::WeightedSketcher<mh::BBitMinHasher<uint64_t>>);
FINAL_OVERLOAD(wj::WeightedSketcher<SuperMinHashType>);
FINAL_OVERLOAD(wj::WeightedSketcher<WideHyperLogLogHasher<>>);
FINAL_OVERLOAD(wj::WeightedSketcher<HyperLogLogHasher<>>);
FINAL_OVERLOAD(wj::WeightedSketcher<CBBMinHashType>);
FINAL_OVERLOAD(wj::WeightedSketcher<HyperMinHash>);
FINAL_OVERLOAD2(wj::WeightedSketcher<bf::bf_t, wj::ExactCountingAdapter>);
FINAL_OVERLOAD2(wj::WeightedSketcher<hll::hll_t, wj::ExactCountingAdapter>);
FINAL_OVERLOAD2(wj::WeightedSketcher<khset64_t, wj::ExactCountingAdapter>);
FINAL_OVERLOAD2(wj::WeightedSketcher<mh::CountingRangeMinHash<uint64_t>, wj::ExactCountingAdapter>);
FINAL_OVERLOAD2(wj::WeightedSketcher<mh::RangeMinHash<uint64_t>, wj::ExactCountingAdapter>);
FINAL_OVERLOAD2(wj::WeightedSketcher<BKHash64, wj::ExactCountingAdapter>);
FINAL_OVERLOAD2(wj::WeightedSketcher<mh::BBitMinHasher<uint64_t>, wj::ExactCountingAdapter>);
FINAL_OVERLOAD2(wj::WeightedSketcher<SuperMinHashType, wj::ExactCountingAdapter>);
FINAL_OVERLOAD2(wj::WeightedSketcher<WideHyperLogLogHasher<>, wj::ExactCountingAdapter>);
FINAL_OVERLOAD2(wj::WeightedSketcher<HyperLogLogHasher<>, wj::ExactCountingAdapter>);
FINAL_OVERLOAD2(wj::WeightedSketcher<CBBMinHashType, wj::ExactCountingAdapter>);
FINAL_OVERLOAD2(wj::WeightedSketcher<HyperMinHash, wj::ExactCountingAdapter>);
template<typename T>struct SketchFileSuffix {static constexpr const char *suffix = ".sketch";};
#define SSS(type, suf) \
    template<> struct SketchFileSuffix<type> {static constexpr const char *suffix = suf;};\
    template<> struct SketchFileSuffix<wj::WeightedSketcher<type>> {static constexpr const char *suffix = ".wj" suf;};\
    template<> struct SketchFileSuffix<wj::WeightedSketcher<type, wj::ExactCountingAdapter>> {static constexpr const char *suffix = ".wj.exact" suf;}

SSS(mh::CountingRangeMinHash<uint64_t>, ".crmh");
SSS(mh::RangeMinHash<uint64_t>, ".rmh");
SSS(BKHash64, ".rmh");
SSS(khset64_t, ".khs");
SSS(bf::bf_t, ".bf");
SSS(mh::BBitMinHasher<uint64_t>, ".bmh");
SSS(WideHyperLogLogHasher<>, ".whll");
SSS(SuperMinHashType, ".bbs");
SSS(CBBMinHashType, ".cbmh");
SSS(HyperMinHash, ".hmh");
SSS(hll::hll_t, ".hll");
SSS(HyperLogLogHasher<>, ".hll");
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


using HLLH = HyperLogLogHasher<>;
template<typename T> struct SketchEnum;
template<> struct SketchEnum<hll::hll_t> {static constexpr Sketch value = HLL;};
template<> struct SketchEnum<HLLH> {static constexpr Sketch value = HLL;};
template<> struct SketchEnum<bf::bf_t> {static constexpr Sketch value = BLOOM_FILTER;};
template<> struct SketchEnum<mh::RangeMinHash<uint64_t>> {static constexpr Sketch value = RANGE_MINHASH;};
template<> struct SketchEnum<BKHash64> {static constexpr Sketch value = RANGE_MINHASH;};
template<> struct SketchEnum<mh::CountingRangeMinHash<uint64_t>> {static constexpr Sketch value = COUNTING_RANGE_MINHASH;};
template<> struct SketchEnum<mh::BBitMinHasher<uint64_t>> {static constexpr Sketch value = BB_MINHASH;};
template<> struct SketchEnum<CBBMinHashType> {static constexpr Sketch value = COUNTING_BB_MINHASH;};
template<> struct SketchEnum<khset64_t> {static constexpr Sketch value = FULL_KHASH_SET;};
template<> struct SketchEnum<SuperMinHashType> {static constexpr Sketch value = BB_SUPERMINHASH;};
template<> struct SketchEnum<WideHyperLogLogHasher<>> {static constexpr Sketch value = WIDE_HLL;};
template<> struct SketchEnum<HyperMinHash> {static constexpr Sketch value = HYPERMINHASH;};


template<> struct SketchEnum<wj::WeightedSketcher<hll::hll_t>> {static constexpr Sketch value = HLL;};
template<> struct SketchEnum<wj::WeightedSketcher<HLLH>> {static constexpr Sketch value = HLL;};
template<> struct SketchEnum<wj::WeightedSketcher<bf::bf_t>> {static constexpr Sketch value = BLOOM_FILTER;};
template<> struct SketchEnum<wj::WeightedSketcher<mh::RangeMinHash<uint64_t>>> {static constexpr Sketch value = RANGE_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<BKHash64>> {static constexpr Sketch value = RANGE_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<mh::CountingRangeMinHash<uint64_t>>> {static constexpr Sketch value = COUNTING_RANGE_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<mh::BBitMinHasher<uint64_t>>> {static constexpr Sketch value = BB_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<CBBMinHashType>> {static constexpr Sketch value = COUNTING_BB_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<khset64_t>> {static constexpr Sketch value = FULL_KHASH_SET;};
template<> struct SketchEnum<wj::WeightedSketcher<SuperMinHashType>> {static constexpr Sketch value = BB_SUPERMINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<WideHyperLogLogHasher<>>> {static constexpr Sketch value = WIDE_HLL;};
template<> struct SketchEnum<wj::WeightedSketcher<HyperMinHash>> {static constexpr Sketch value = HYPERMINHASH;};

template<> struct SketchEnum<wj::WeightedSketcher<hll::hll_t, wj::ExactCountingAdapter>> {static constexpr Sketch value = HLL;};
template<> struct SketchEnum<wj::WeightedSketcher<HLLH, wj::ExactCountingAdapter>> {static constexpr Sketch value = HLL;};
template<> struct SketchEnum<wj::WeightedSketcher<bf::bf_t, wj::ExactCountingAdapter>> {static constexpr Sketch value = BLOOM_FILTER;};
template<> struct SketchEnum<wj::WeightedSketcher<mh::RangeMinHash<uint64_t>, wj::ExactCountingAdapter>> {static constexpr Sketch value = RANGE_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<BKHash64, wj::ExactCountingAdapter>> {static constexpr Sketch value = RANGE_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<mh::CountingRangeMinHash<uint64_t>, wj::ExactCountingAdapter>> {static constexpr Sketch value = COUNTING_RANGE_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<mh::BBitMinHasher<uint64_t>, wj::ExactCountingAdapter>> {static constexpr Sketch value = BB_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<CBBMinHashType, wj::ExactCountingAdapter>> {static constexpr Sketch value = COUNTING_BB_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<khset64_t, wj::ExactCountingAdapter>> {static constexpr Sketch value = FULL_KHASH_SET;};
template<> struct SketchEnum<wj::WeightedSketcher<SuperMinHashType, wj::ExactCountingAdapter>> {static constexpr Sketch value = BB_SUPERMINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<WideHyperLogLogHasher<>, wj::ExactCountingAdapter>> {static constexpr Sketch value = WIDE_HLL;};
template<> struct SketchEnum<wj::WeightedSketcher<HyperMinHash, wj::ExactCountingAdapter>> {static constexpr Sketch value = HYPERMINHASH;};


template<typename T>
inline void set_estim_and_jestim(T &x, hll::EstimationMethod estim, hll::JointEstimationMethod jestim) {}

template<typename Hashstruct>
inline void set_estim_and_jestim(hll::hllbase_t<Hashstruct> &h, hll::EstimationMethod estim, hll::JointEstimationMethod jestim) {
    h.set_estim(estim);
    h.set_jestim(jestim);
}
template<typename T> inline T construct(size_t ssarg);
template<typename T, bool is_weighted> struct Constructor;
template<typename T>
inline T construct(size_t ssarg) {
    Constructor<T, wj::is_weighted_sketch<T>::value> constructor;
    return constructor.create(ssarg);
}

template<typename T> struct Constructor<T, false> {
    static auto create(size_t ssarg) {
        return T(ssarg);
    }
};
template<> struct Constructor<BBitMinHasher<uint64_t>, false> {
    static auto create(size_t ssarg) {
        return BBitMinHasher<uint64_t>(ssarg, gargs.bbnbits);
    }
};
template<typename T> struct Constructor<T, true> {
    static auto create(size_t ssarg) {
        using base_type = typename T::base_type;
        using cm_type = typename T::cm_type;
        return T(cm_type(16, gargs.weighted_jaccard_cmsize, gargs.weighted_jaccard_nhashes), construct<base_type>(ssarg));
    }
};
template<template<typename> typename Adapter> struct Constructor<Adapter<BBitMinHasher<uint64_t>>, true> {
    using Type = Adapter<BBitMinHasher<uint64_t>>;
    static auto create(size_t ssarg) {
        using base_type = BBitMinHasher<uint64_t> /*typename Type::base_type*/;
        using cm_type = typename Adapter<BBitMinHasher<uint64_t>>::cm_type;
        return Type(cm_type(16, gargs.weighted_jaccard_cmsize, gargs.weighted_jaccard_nhashes), construct<base_type>(ssarg));
    }
};


template<typename T>
inline double cardinality_estimate(T &x) {
    return x.cardinality_estimate();
}
template<> inline double cardinality_estimate(hll::hll_t &x) {return x.report();}
template<> inline double cardinality_estimate(mh::FinalBBitMinHash &x) {return x.est_cardinality_;}
template<> inline double cardinality_estimate(mh::FinalDivBBitMinHash &x) {return x.est_cardinality_;}
template<> inline double cardinality_estimate(sketch::HyperMinHash &x) {return x.getcard();}

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
US_DEC(sketch::HyperMinHash)
US_DEC(CRMFinal)
US_DEC(khset64_t)
#undef US_DEC
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

template<typename ST>
float result_cmp(const ST &lhs, const ST &rhs, EmissionType result_type, double ksinv) {
    double ret = std::numeric_limits<double>::infinity();
    switch(result_type) {
        case FULL_MASH_DIST: case MASH_DIST: case JI: {
            ret = similarity<const ST>(lhs, rhs);
            if(result_type == MASH_DIST) ret = dist_index(ret, ksinv);
            else if(result_type == FULL_MASH_DIST) ret = full_dist_index(ret, ksinv);
        } break;
        case SYMMETRIC_CONTAINMENT_DIST: case SYMMETRIC_CONTAINMENT_INDEX: case SIZES: case FULL_CONTAINMENT_DIST: case CONTAINMENT_INDEX: case CONTAINMENT_DIST: {
            const auto triple = set_triple(lhs, rhs);
            ret = triple[2];
            if(result_type == SYMMETRIC_CONTAINMENT_INDEX || result_type == SYMMETRIC_CONTAINMENT_DIST) {
                ret /= (std::min(triple[0], triple[1]) + triple[2]);
                if(result_type == SYMMETRIC_CONTAINMENT_DIST) ret = containment_dist(ret, ksinv);
            } else if(result_type == FULL_CONTAINMENT_DIST || result_type == CONTAINMENT_DIST || result_type == CONTAINMENT_INDEX) {
                ret /= (triple[0] + triple[1] + triple[2]);
                if(result_type == CONTAINMENT_DIST) ret = containment_dist(ret, ksinv);
                else if(result_type == FULL_CONTAINMENT_DIST) ret = full_containment_dist(ret, ksinv);
            } // else, result_type is (SIZES), and we return ret
        } break;
        default: __builtin_unreachable();
    }
    return static_cast<float>(ret);
}

template<typename SketchType> inline hll::hll_t &get_hll(SketchType &s);
template<> inline hll::hll_t &get_hll<hll::hll_t>(hll::hll_t &s) {return s;}
template<typename SketchType> inline const hll::hll_t &get_hll(const SketchType &s);
template<> inline const hll::hll_t &get_hll<hll::hll_t>(const hll::hll_t &s) {return s;}

template<typename SketchType>
struct est_helper {
    const Spacer                      &sp_;
    const std::vector<std::string> &paths_;
    std::mutex                         &m_;
    const u64                          np_;
    const bool                      canon_;
    void                            *data_;
    std::vector<SketchType>         &hlls_;
    kseq_t                            *ks_;
};

template<typename SketchType, typename ScoreType=score::Lex>
void est_helper_fn(void *data_, long index, int tid) {
    est_helper<SketchType> &h(*(est_helper<SketchType> *)(data_));
    fill_lmers<ScoreType, SketchType>(h.hlls_[tid], h.paths_[index], h.sp_, h.canon_, h.data_, h.ks_ + tid);
}

template<typename SketchType, typename ScoreType=score::Lex>
void fill_sketch(SketchType &ret, const std::vector<std::string> &paths,
              unsigned k, uint16_t w, const spvec_t &spaces, bool canon=true,
              void *data=nullptr, int num_threads=1, u64 np=23, kseq_t *ks=nullptr) {
    // Default to using all available threads if num_threads is negative.
    if(num_threads < 0) {
        num_threads = std::thread::hardware_concurrency();
        LOG_INFO("Number of threads was negative and has been adjusted to all available threads (%i).\n", num_threads);
    }
    const Spacer space(k, w, spaces);
    if(num_threads <= 1) {
        LOG_DEBUG("Starting serial\n");
        for(u64 i(0); i < paths.size(); fill_lmers<ScoreType, SketchType>(ret, paths[i++], space, canon, data, ks));
    } else {
        LOG_DEBUG("Starting parallel\n");
        std::mutex m;
        KSeqBufferHolder kseqs(num_threads);
        std::vector<SketchType> sketches;
        while(sketches.size() < (unsigned)num_threads) sketches.emplace_back(ret.clone());
        est_helper<SketchType> helper{space, paths, m, np, canon, data, sketches, kseqs.data()};
        kt_for(num_threads, &est_helper_fn<SketchType, ScoreType>, &helper, paths.size());
        auto &rhll = get_hll(ret);
        for(auto &sketch: sketches) rhll += get_hll(sketch);
    }

}
template<typename ScoreType=score::Lex>
hll::hll_t make_hll(const std::vector<std::string> &paths,
                unsigned k, uint16_t w, spvec_t spaces, bool canon=true,
                void *data=nullptr, int num_threads=1, u64 np=23, kseq_t *ks=false, hll::EstimationMethod estim=hll::EstimationMethod::ERTL_MLE, uint16_t jestim=hll::JointEstimationMethod::ERTL_JOINT_MLE, bool clamp=true) {
    hll::hll_t master(np, estim, (hll::JointEstimationMethod)jestim, 1, clamp);
    fill_sketch<hll::hll_t, ScoreType>(master, paths, k, w, spaces, canon, data, num_threads, np, ks);
    return master;
}
template<typename ScoreType=score::Lex>
u64 estimate_cardinality(const std::vector<std::string> &paths,
                            unsigned k, uint16_t w, spvec_t spaces, bool canon,
                            void *data=nullptr, int num_threads=-1, u64 np=23, kseq_t *ks=nullptr, hll::EstimationMethod estim=hll::EstimationMethod::ERTL_MLE) {
    auto tmp(make_hll<ScoreType>(paths, k, w, spaces, canon, data, num_threads, np, ks, estim));
    return tmp.report();
}


template<typename SketchType>
void partdist_loop(std::FILE *ofp, SketchType *hlls, const std::vector<std::string> &inpaths, const bool use_scientific, const unsigned k, const EmissionType result_type, EmissionFormat emit_fmt, const size_t buffer_flush_size,
                   size_t nq)
{
    const float ksinv = 1./ k;
    if(nq >= inpaths.size()) {
        UNRECOVERABLE_ERROR(ks::sprintf("Wrong number of query/references. (ip size: %zu, nq: %zu\n", inpaths.size(), nq).data());
    }
    size_t nr = inpaths.size() - nq;
#if TIMING
    auto start = std::chrono::high_resolution_clock::now();
#endif
    std::future<void> write_future, fmt_future;
    std::array<ks::string, 2> buffers;
    for(auto &b: buffers) b.resize(4 * nr);
    for(size_t qi = nr; qi < inpaths.size(); ++qi) {
        auto &hq = hlls[qi];
        std::unique_ptr<float[]> arr(new float[nr]);
        OMP_PFOR_DYN
        for(size_t j = 0; j < nr; ++j) {
            arr[j] = result_cmp(hlls[j], hq, result_type, ksinv);
        }
        switch(emit_fmt) {
            case BINARY:
                if(write_future.valid()) write_future.get();
                write_future = std::async(std::launch::async, [arr=std::move(arr), nr,ofp]() {
                    if(unlikely(std::fwrite(arr.get(), sizeof(float), nr, ofp) != nr))
                    UNRECOVERABLE_ERROR("Error writing to binary file");
                });
                break;
            case UT_TSV: case UPPER_TRIANGULAR: default:
                // UNRECOVERABLE_ERROR(std::string("Illegal output format. numeric: ") + std::to_string(int(emit_fmt)));
            case FULL_TSV:
                if(fmt_future.valid()) fmt_future.get();
                fmt_future = std::async(std::launch::async, [nr,qi,ofp,ind=qi-nr,&inpaths,use_scientific,arr=std::move(arr),&buffers,&write_future]() {
                    auto &buffer = buffers[qi & 1];
                    buffer += inpaths[qi];
                    for(size_t i = 0; i < nr; ++i)
                        buffer.sprintf("\t%g", arr[i]);
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
    std::fprintf(stderr, "partdist (%zu by %zu) took %gms\n", nr, nq, std::chrono::duration<double, std::milli>(end - start).count());
#endif
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
int fold_main(int argc, char *argv[]);
int card_main(int argc, char *argv[]);
int panel_main(int argc, char *argv[]);
int dist_main(int argc, char *argv[]);
int print_binary_main(int argc, char *argv[]);
int mkdist_main(int argc, char *argv[]);
int flatten_main(int argc, char *argv[]);
int hll_main(int argc, char *argv[]);
int union_main(int argc, char *argv[]);
int view_main(int argc, char *argv[]);
int sketch_by_seq_main(int argc, char *argv[]);
int dist_by_seq_main(int argc, char *argv[]);
}

#endif /* DASHING_H__ */
