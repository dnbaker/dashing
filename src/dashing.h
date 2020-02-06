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

#if __cplusplus >= 201703L && __cpp_lib_execution
#include <execution>
#endif
#ifdef FNAME_SEP
#pragma message("Not: FNAME_SEP already defined. [not default \"' '\"]")
#else
#define FNAME_SEP ' '
#endif
#ifndef BUFFER_FLUSH_SIZE
#define BUFFER_FLUSH_SIZE (1u << 18)
#endif

#define LO_ARG(LONG, SHORT) {LONG, required_argument, 0, SHORT},
#define LO_NO(LONG, SHORT) {LONG, no_argument, 0, SHORT},
#define LO_FLAG(LONG, SHORT, VAR, VAL) {LONG, no_argument, (int *)&VAR, VAL},


namespace bns {
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

static int flatten_all(const std::vector<std::string> &fpaths, size_t nk, const std::string outpath) {
    if(fpaths.empty()) RUNTIME_ERROR("no fpaths, see usage.");
    // TODO: adapt this to pack an asymmetric comparison result into a flattened blaze distance matrix.
    std::vector<dm::DistanceMatrix<float>> dms;
    dms.reserve(nk);
    for(const auto &fp: fpaths)
        dms.emplace_back(fp.data());
    const uint64_t ne = dms.front().num_entries();
    assert(std::accumulate(dms.begin() + 1, dms.end(), true,
           [ne](bool val, const auto &x) {return val && x.num_entries() == ne;}));
    float *outp = static_cast<float *>(std::malloc(nk * ne * sizeof(float)));
    if(!outp) {
        std::fprintf(stderr, "Allocation of %zu bytes failed\n", size_t(nk * ne * sizeof(float))); return 1;
    }

    static constexpr uint64_t NB = 4096;
    #pragma omp parallel for
    for(size_t i = 0; i < ((NB - 1) + ne) / NB; ++i) {
        auto spos = i * NB, espos = std::min((i + 1) * NB, ne);
        auto destp = outp + spos * nk;//, endp = std::min(destp + nk, outp + ne);
        do {for(auto j = 0u; j < nk;*destp++ = dms[j++][spos]);} while(++spos < espos);
    }
    std::FILE *ofp = fopen(outpath.data(), "wb");
    if(!ofp) return 2;
    std::fwrite(&ne, sizeof(ne), 1, ofp);
    std::fwrite(outp, nk * ne, sizeof(float), ofp);
    std::fclose(ofp);
    std::free(outp);
    return 0;
}
// enums
//
enum EmissionFormat: unsigned {
    UT_TSV = 0,
    BINARY   = 1,
    UPPER_TRIANGULAR = 2,
    PHYLIP_UPPER_TRIANGULAR = 2,
    FULL_TSV = 3,
    JSON = 4,
    NEAREST_NEIGHBOR_TABLE = 8,
    NEAREST_NEIGHBOR_BINARY =  NEAREST_NEIGHBOR_TABLE | BINARY
};

enum Sketch: int {
    HLL,
    BLOOM_FILTER,
    RANGE_MINHASH,
    FULL_KHASH_SET,
    COUNTING_RANGE_MINHASH,
    BB_MINHASH,
    BB_SUPERMINHASH,
    COUNTING_BB_MINHASH, // TODO make this work.
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
};

enum sketching_method: int {
    EXACT = 0,
    CBF   = 1,
    BY_FNAME = 2
};


struct GlobalArgs {
    size_t weighted_jaccard_cmsize = 22;
    size_t weighted_jaccard_nhashes = 8;
    uint32_t bbnbits = 16;
    uint32_t number_neighbors = 0;
};

static GlobalArgs gargs;
enum EmissionType {
    MASH_DIST = 0,
    JI        = 1,
    SIZES     = 2,
    FULL_MASH_DIST = 3,
    FULL_CONTAINMENT_DIST = 4,
    CONTAINMENT_INDEX = 5,
    CONTAINMENT_DIST = 6,
    SYMMETRIC_CONTAINMENT_INDEX = 7,
    SYMMETRIC_CONTAINMENT_DIST = 8,
};

enum NNType {
    // Nearest Neighbor Type
    DECREASING = 0,
    INCREASING = 1,
    DIST_MEASURE = DECREASING,
    SIMILARITY_MEASURE = INCREASING
};

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
FINAL_OVERLOAD(SuperMinHashType);
FINAL_OVERLOAD(CBBMinHashType);
FINAL_OVERLOAD(wj::WeightedSketcher<mh::CountingRangeMinHash<uint64_t>>);
FINAL_OVERLOAD(wj::WeightedSketcher<mh::RangeMinHash<uint64_t>>);
FINAL_OVERLOAD(wj::WeightedSketcher<mh::BBitMinHasher<uint64_t>>);
FINAL_OVERLOAD(wj::WeightedSketcher<SuperMinHashType>);
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

enum EncodingType {
    BONSAI,
    NTHASH,
    RK,
    CYCLIC
};


template<typename T> struct SketchEnum;
template<> struct SketchEnum<hll::hll_t> {static constexpr Sketch value = HLL;};
template<> struct SketchEnum<bf::bf_t> {static constexpr Sketch value = BLOOM_FILTER;};
template<> struct SketchEnum<mh::RangeMinHash<uint64_t>> {static constexpr Sketch value = RANGE_MINHASH;};
template<> struct SketchEnum<mh::CountingRangeMinHash<uint64_t>> {static constexpr Sketch value = COUNTING_RANGE_MINHASH;};
template<> struct SketchEnum<mh::BBitMinHasher<uint64_t>> {static constexpr Sketch value = BB_MINHASH;};
template<> struct SketchEnum<CBBMinHashType> {static constexpr Sketch value = COUNTING_BB_MINHASH;};
template<> struct SketchEnum<khset64_t> {static constexpr Sketch value = FULL_KHASH_SET;};
template<> struct SketchEnum<SuperMinHashType> {static constexpr Sketch value = BB_SUPERMINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<hll::hll_t>> {static constexpr Sketch value = HLL;};
template<> struct SketchEnum<wj::WeightedSketcher<bf::bf_t>> {static constexpr Sketch value = BLOOM_FILTER;};
template<> struct SketchEnum<wj::WeightedSketcher<mh::RangeMinHash<uint64_t>>> {static constexpr Sketch value = RANGE_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<mh::CountingRangeMinHash<uint64_t>>> {static constexpr Sketch value = COUNTING_RANGE_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<mh::BBitMinHasher<uint64_t>>> {static constexpr Sketch value = BB_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<CBBMinHashType>> {static constexpr Sketch value = COUNTING_BB_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<khset64_t>> {static constexpr Sketch value = FULL_KHASH_SET;};
template<> struct SketchEnum<wj::WeightedSketcher<SuperMinHashType>> {static constexpr Sketch value = BB_SUPERMINHASH;};

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


#define US_DEC(type) \
template<> inline double union_size<type> (const type &a, const type &b) { \
    return a.union_size(b); \
}

US_DEC(RMFinal)
US_DEC(CRMFinal)
US_DEC(khset64_t)
US_DEC(hll::hllbase_t<>)
template<> inline double union_size<mh::FinalBBitMinHash> (const mh::FinalBBitMinHash &a, const mh::FinalBBitMinHash &b) {
    return (a.est_cardinality_ + b.est_cardinality_ ) / (1. + a.jaccard_index(b));
}
//template<> inline double union_size<mh::FinalBBitMinHash> (const mh::FinalBBitMinHash &a, const mh::FinalBBitMinHash &b) {
//    return (a.est_cardinality_ + b.est_cardinality_ ) / (1. + a.jaccard_index(b));
//}
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
                DO_LOOP(us::union_size);
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
