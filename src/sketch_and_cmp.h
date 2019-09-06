#pragma once
#include "dashing.h"

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
            case SYMMETRIC_CONTAINMENT_DIST:\
                perform_core_op(dists, nsketches, hlls, [ksinv](const auto &x, const auto &y) { \
                    const auto triple = set_triple(x, y);\
                    auto ret = triple[2] / (std::min(triple[0], triple[1]) + triple[2]);\
                    auto di = dist_index(ret, ksinv);\
                    assert(di <= dist_index(triple[2] / (triple[0] + triple[2]), ksinv));\
                    assert(di <= dist_index(triple[2] / (triple[1] + triple[2]), ksinv));\
                    return di;\
                }, i);\
                break;\
            case SYMMETRIC_CONTAINMENT_INDEX:\
                perform_core_op(dists, nsketches, hlls, [&](const auto &x, const auto &y) {\
                    const auto triple = set_triple(x, y);\
                    auto ret = triple[2] / (std::min(triple[0], triple[1]) + triple[2]);\
                    assert(ret >= triple[2] / (std::max(triple[0], triple[1]) + triple[2]) || triple[1] == 0. || triple[0] == 0.);\
                    return ret;\
                }, i);\
                break;\
            default: __builtin_unreachable();\
        } } while(0)



namespace bns {
using namespace sketch;
using namespace hll;

template<typename SketchType, typename T, typename Func>
INLINE void perform_core_op(T &dists, size_t nhlls, SketchType *hlls, const Func &func, size_t i) {
    auto &h1 = hlls[i];
    #pragma omp parallel for schedule(dynamic)
    for(size_t j = i + 1; j < nhlls; ++j)
        dists[j - i - 1] = func(hlls[j], h1);
    h1.free();
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
template<typename SketchType>
std::string make_fname(const char *path, size_t sketch_p, int wsz, int k, int csz, const std::string &spacing,
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
static size_t bytesl2_to_arg(int nblog2, Sketch sketch) {
    switch(sketch) {
        case HLL: return nblog2;
        case BLOOM_FILTER: return nblog2 + 3; // 8 bits per byte
        case RANGE_MINHASH: return size_t(1) << (nblog2 - 3); // 8 bytes per minimizer
        case COUNTING_RANGE_MINHASH: return (size_t(1) << (nblog2)) / (sizeof(uint64_t) + sizeof(uint32_t));
        case BB_MINHASH:
            return nblog2 - std::floor(std::log2(gargs.bbnbits / 8));
        case BB_SUPERMINHASH:
            return size_t(1) << (nblog2 - int(std::log2(gargs.bbnbits / 8)));
        case FULL_KHASH_SET: return 16; // Reserve hash set size a bit. Mostly meaningless, resizing as necessary.
        default: {
            char buf[128];
            std::sprintf(buf, "Sketch %s not yet supported.\n", (size_t(sketch) >= (sizeof(sketch_names) / sizeof(char *)) ? "Not such sketch": sketch_names[sketch]));
            RUNTIME_ERROR(buf);
            return -1337;
        }
    }

}

template<typename SketchType>
void partdist_loop(std::FILE *ofp, SketchType *hlls, const std::vector<std::string> &inpaths, const bool use_scientific,
                   const unsigned k, const EmissionType result_type, EmissionFormat emit_fmt, int nthreads, const size_t buffer_flush_size,
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
                for(size_t j = 0; j < nr; ++j) {
                    auto tmp = set_triple(hlls[j], hlls[qi]);
                    arr[qind * nr + j] = tmp[2] / (std::min(tmp[0], tmp[1]) + tmp[2]);
                }
                break;
            case SYMMETRIC_CONTAINMENT_DIST:
                #pragma omp parallel for schedule(dynamic)
                for(size_t j = 0; j < nr; ++j) {
                    auto tmp = set_triple(hlls[j], hlls[qi]);
                    arr[qind * nr + j] = dist_index(tmp[2] / (std::min(tmp[0], tmp[1]) + tmp[2]), ksinv);
                }
                break;
            default: throw std::runtime_error("Value not found");
#undef DO_LOOP
#undef dist_sim
#undef cont_sim
#undef fulldist_sim
#undef fullcont_sim
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
} // partdist_loop


template<typename SketchType>
void dist_loop(std::FILE *ofp, SketchType *hlls, const std::vector<std::string> &inpaths, const bool use_scientific, const unsigned k, const EmissionType result_type, EmissionFormat emit_fmt, int nthreads, const size_t buffer_flush_size=BUFFER_FLUSH_SIZE, size_t nq=0) {
    if(nq) {
        partdist_loop<SketchType>(ofp, hlls, inpaths, use_scientific, k, result_type, emit_fmt, nthreads, buffer_flush_size, nq);
        return;
    }
    if(!is_symmetric(result_type)) {
        char buf[1024];
        std::sprintf(buf, "Can't perform symmetric distance comparisons with a symmetric method (%s/%d). To perform an asymmetric distance comparison between a given set and itself, provide the same list of filenames to both -Q and -F.\n", emt2str(result_type), int(result_type));
        RUNTIME_ERROR(buf);
    }
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
            //LOG_DEBUG("Finished chunk %zu of %zu\n", i + 1, nsketches);
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
            std::fprintf(stderr, "Writing to file\n");
            dm.write(ofp);
        }
    }
} // dist_loop
template<typename SketchType>
void dist_sketch_and_cmp(const std::vector<std::string> &inpaths, std::vector<CountingSketch> &cms, KSeqBufferHolder &kseqs, std::FILE *ofp, std::FILE *pairofp,
                         Spacer sp,
                         unsigned ssarg, unsigned mincount, EstimationMethod estim, JointEstimationMethod jestim, bool cache_sketch, EmissionType result_type, EmissionFormat emit_fmt,
                         bool presketched_only, unsigned nthreads, bool use_scientific, std::string suffix, std::string prefix, bool canon, bool entropy_minimization, std::string spacing,
                         size_t nq=0, EncodingType enct=BONSAI)
{
    // nq -- number of queries
    //       for convenience, we will perform our comparisons (all-p) against (all-q) [remainder]
    //       and use the same guts for all portions of the process
    //       except for the final comparison and output.
    assert(nq <= inpaths.size());
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

    std::atomic<uint32_t> ncomplete;
    ncomplete.store(0);
    const unsigned k = sp.k_;
    const unsigned wsz = sp.w_;
    RollingHasher<uint64_t> rolling_hasher(k, canon);
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
            const std::string fpath(make_fname<SketchType>(path.data(), sketch_size, wsz, k, sp.c_, spacing, suffix, prefix, enct));
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
                if(entropy_minimization) {
                    FILL_SKETCH_MIN(score::Entropy);
                } else {
                    FILL_SKETCH_MIN(score::Lex);
                }
                CONST_IF(samesketch) {
                    if(cache_sketch && !isf) sketch.write(fpath);
                } else if(cache_sketch) final_sketches[i].write(fpath);
            }
        }
        ++ncomplete; // Atomic
    }
    _Pragma("omp parallel for")
    for(size_t i = 0; i < sketches.size(); ++i) {
        sketch_finalize(final_sketches[i]);
    }
    kseqs.free();
    ks::string str("#Path\tSize (est.)\n");
    assert(str == "#Path\tSize (est.)\n");
    str.resize(BUFFER_FLUSH_SIZE);
    {
        const int fn(fileno(ofp));
        for(size_t i(0); i < sketches.size(); ++i) {
            double card;
            CONST_IF(samesketch) card = cardinality_estimate(sketches[i]);
            else                 card = cardinality_estimate(final_sketches[i]);
            str.sprintf("%s\t%zu\n", inpaths[i].data(), size_t(card));
            if(str.size() >= BUFFER_FLUSH_SIZE) str.flush(fn);
        }
        str.flush(fn);
    }
    if(ofp != stdout) std::fclose(ofp);
    str.clear();
    if(emit_fmt == UT_TSV) {
        str.sprintf("##Names\t");
        for(size_t i = 0; i < inpaths.size() - nq; ++i)
            str.sprintf("%s\t", inpaths[i].data());
        str.back() = '\n';
        str.write(fileno(pairofp)); str.free();
    } else if(emit_fmt == UPPER_TRIANGULAR) { // emit_fmt == UPPER_TRIANGULAR
        std::fprintf(pairofp, "%zu\n", inpaths.size());
        std::fflush(pairofp);
    }
    dist_loop<final_type>(pairofp, final_sketches, inpaths, use_scientific, k, result_type, emit_fmt, nthreads, BUFFER_FLUSH_SIZE, nq);
    CONST_IF(!samesketch) {
#if __cplusplus >= 201703L
        std::destroy_n(
#  if __cpp_lib_execution
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
} // dist_sketch_and_cmp

#undef CORE_ITER

} // namespace bns
