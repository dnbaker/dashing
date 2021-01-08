#pragma once
#include "dashing.h"


using ::sketch::hll::EstimationMethod;
using ::sketch::hll::JointEstimationMethod;
using bns::EmissionFormat;
using bns::EmissionType;
using bns::Spacer;
using bns::KSeqBufferHolder;
using namespace sketch;

namespace bns {


template<typename FType, typename=typename std::enable_if<std::is_floating_point<FType>::value>::type>
size_t submit_emit_dists(int pairfi, const FType *ptr, u64 hs, size_t index, ks::string &str, const std::vector<std::string> &inpaths, EmissionFormat emit_fmt, bool use_scientific, const size_t buffer_flush_size=BUFFER_FLUSH_SIZE) {
    if(emit_fmt & BINARY) {
        const ssize_t nbytes = sizeof(FType) * (hs - index - 1);
        LOG_DEBUG("Writing %zd bytes for %zu items\n", nbytes, (hs - index - 1));
        ssize_t i = ::write(pairfi, ptr, nbytes);
        if(i != nbytes) {
            std::fprintf(stderr, "WARNING: written %zd bytes instead of expected %zd. Something is likely wrong.\n", i, nbytes);
        }
    } else {
        auto &strref = inpaths[index];
        str += strref;
        if(emit_fmt == UT_TSV) {
            constexpr const char *fmt = "\t%.6g";
            {
                u64 k;
                for(k = 0; k < index + 1;  ++k, kputsn_("\t-", 2, reinterpret_cast<kstring_t *>(&str)));
                for(k = 0; k < hs - index - 1; str.sprintf(fmt, ptr[k++]));
            }
        } else { // emit_fmt == UPPER_TRIANGULAR
            constexpr const char *fmt = "\t%.6g";
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
void dist_loop(std::FILE *&ofp, std::string ofpname, SketchType *sketches, const std::vector<std::string> &inpaths, const bool use_scientific, const unsigned k, const EmissionType result_type, EmissionFormat emit_fmt, int, const size_t buffer_flush_size, size_t nq);
using namespace sketch;
using namespace hll;
static size_t bytesl2_to_arg(int nblog2, Sketch sketch) {
    switch(sketch) {
        case HLL: case WIDE_HLL: return nblog2;
        case BLOOM_FILTER: return nblog2 + 3; // 8 bits per byte
        case RANGE_MINHASH: return size_t(1) << (nblog2 - 3); // 8 bytes per minimizer
        case COUNTING_RANGE_MINHASH: return (size_t(1) << (nblog2)) / double(sizeof(uint64_t) + sizeof(uint32_t));
        case BB_MINHASH:
            return nblog2 - std::floor(std::log2(gargs.bbnbits / 8));
        case BB_SUPERMINHASH:
            return size_t(1) << (nblog2 - int(std::log2(gargs.bbnbits / 8)));
        case FULL_KHASH_SET: return 16; // Reserve hash set size a bit. Mostly meaningless, resizing as necessary.
        case HYPERMINHASH: {
            switch(gargs.bbnbits) {
                case 8: return nblog2;
                case 16: return nblog2 - 1;
                case 32: return nblog2 - 2;
                case 64: return nblog2 - 3;
                default: {
                    if(gargs.bbnbits < 8) {gargs.bbnbits = 8; return nblog2;}
                    if(gargs.bbnbits < 16) {gargs.bbnbits = 16; return nblog2 - 1;}
                    if(gargs.bbnbits < 32) {gargs.bbnbits = 32; return nblog2 - 2;}
                    gargs.bbnbits = 64;
                    return nblog2 - 3;
                }
            }
        }
        default: {
            char buf[128];
            std::sprintf(buf, "Sketch %s not yet supported.\n", (size_t(sketch) >= (sizeof(sketch_names) / sizeof(char *)) ? "Not such sketch": sketch_names[sketch]));
            UNRECOVERABLE_ERROR(buf);
            return -1337;
        }
    }

}

template<typename SketchType>
void dist_by_seq(std::vector<std::string> &labels, std::string datapath,
                 std::FILE *pairofp, std::string outpath, int k,
                 EstimationMethod estim, JointEstimationMethod jestim, EmissionType result_type, EmissionFormat emit_fmt,
                 unsigned nthreads, std::string otherpath)
{
    gzFile sfp = gzopen(datapath.data(), "rb");
    if(!sfp) throw sketch::ZlibError(std::string("Failed to open file at ") + datapath);
    std::vector<SketchType> sketches;
    std::vector<std::string> qnames;
    sketches.reserve(labels.size());
    while(sketches.size() < labels.size()) {
        sketches.emplace_back(sfp);
        set_estim_and_jestim(sketches.back(), estim, jestim);
    }
    if(otherpath.size()) {
        qnames = get_paths((otherpath + ".names").data());
        if(qnames.empty()) UNRECOVERABLE_ERROR("Can't compare with empty qnames");
        sketches.reserve(qnames.size() + sketches.size());
        for(size_t i = 0; i < qnames.size(); ++i) {
            sketches.emplace_back(sfp);
            set_estim_and_jestim(sketches.back(), estim, jestim);
        }
        labels.insert(labels.end(), qnames.begin(), qnames.end());
    } else if(!is_symmetric(result_type)) {
        UNRECOVERABLE_ERROR("Can't perform asymmetric comparison without query paths");
    }
    const size_t nq = qnames.size();
    gzclose(sfp);
    ks::string str;
    if(emit_fmt == UT_TSV && !nq) {
        size_t nq = 0;
        str.sprintf("##Names\t");
        for(size_t i = 0; i < labels.size() - nq; ++i) {
            str.sprintf("%s\t", labels[i].data());
        }
        str.back() = '\n';
        str.flush(fileno(pairofp));
    } else if(emit_fmt == UPPER_TRIANGULAR) { // emit_fmt == UPPER_TRIANGULAR
        std::fprintf(pairofp, "%zu\n", labels.size());
        std::fflush(pairofp);
    }
    dist_loop<SketchType>(pairofp, outpath, sketches.data(), labels, /* use_scientific=*/ true, k, result_type, emit_fmt, nthreads, BUFFER_FLUSH_SIZE, nq);
}


template<typename SketchType>
void size_sketch_and_emit(std::vector<std::string> &inpaths, std::vector<CountingSketch> &cms, KSeqBufferHolder &kseqs, std::FILE *ofp,
                         Spacer sp,
                         unsigned ssarg, unsigned mincount, EncodingType enct, EstimationMethod estim, JointEstimationMethod jestim, bool cache_sketch, bool emit_binary,
                         bool use_scientific,
                         bool presketched_only, unsigned nthreads, std::string suffix, std::string prefix, bool canon, std::string spacing)
{
    // nq -- number of queries
    //       for convenience, we will perform our comparisons (all-p) against (all-q) [remainder]
    //       and use the same guts for all portions of the process
    //       except for the final comparison and output.
    assert(nq <= inpaths.size());
    using final_type = typename FinalSketch<SketchType>::final_type;
    static_assert(!sketch::wj::is_weighted_sketch<final_type>::value, "Can't have a weighted sketch be a final type");
    auto buf(std::make_unique<uint8_t[]>(inpaths.size() * sizeof(SketchType)));
    auto bufp = &buf[0];
    auto sketches = (SketchType *)bufp;
    const size_t npaths = inpaths.size();
    const uint32_t sketch_size = bytesl2_to_arg(ssarg, SketchEnum<SketchType>::value);
    OMP_PFOR
    for(size_t i = 0; i < npaths; ++i) {
        new(sketches + i) SketchType(construct<SketchType>(sketch_size));
        set_estim_and_jestim(sketches[i], estim, jestim);
    }
    static constexpr bool samesketch = std::is_same<SketchType, final_type>::value;
    final_type *final_sketches;
    std::unique_ptr<std::vector<final_type>> raii_final_sketches;

    std::atomic<uint32_t> ncomplete;
    ncomplete.store(0);
    const unsigned k = sp.k_;
    const unsigned wsz = sp.w_;
    RollingHasher<uint64_t> rolling_hasher(k, canon);
    if(npaths == 1 && presketched_only) {
        raii_final_sketches.reset(new std::vector<final_type>);
        gzFile ifp = gzopen(inpaths[0].data(), "rb");
        if(!ifp) UNRECOVERABLE_ERROR("Failed to open file.");
        for(;;) {
            try {
                //TD<decltype(raii_final_sketches)> td;
                raii_final_sketches->emplace_back(ifp);
                set_estim_and_jestim(raii_final_sketches->back(), estim, jestim);
            } catch(...) {
                break;
            }
        }
        final_sketches = raii_final_sketches->data();
        while(inpaths.size() < raii_final_sketches->size())
            inpaths.emplace_back(std::to_string(inpaths.size()));
        inpaths[0] = "0";
    } else {
        final_sketches =
            samesketch ? reinterpret_cast<final_type *>(sketches)
                       : static_cast<final_type *>(std::malloc(sizeof(*final_sketches) * npaths));
        OMP_PFOR_DYN
        for(size_t i = 0; i < npaths; ++i) {
            const std::string &path(inpaths[i]);
            auto &sketch = sketches[i];
            if(presketched_only)  {
                CONST_IF(samesketch) {
                    sketch.read(path);
                    set_estim_and_jestim(sketch, estim, jestim); // HLL is the only type that needs this, and it's the same
                } else {
                    new(final_sketches + i) final_type(path.data()); // Read from path
                }
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
                    Encoder<score::Lex> enc(nullptr, 0, sp, nullptr, canon);
                    if(cms.empty()) {
                        auto &h = sketch;
                        if(enct == BONSAI) for_each_substr([&](const char *s) {enc.for_each([&](u64 kmer){h.addh(kmer);}, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);
                        else if(enct == NTHASH) for_each_substr([&](const char *s) {enc.for_each_hash([&](u64 kmer){h.addh(kmer);}, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);
                        else for_each_substr([&](const char *s) {rolling_hasher.for_each_hash([&](u64 kmer){h.addh(kmer);}, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);
                    } else {
                        CountingSketch &cm = cms.at(tid);
                        const auto lfunc = [&](u64 kmer){if(cm.addh(kmer) >= mincount) sketch.addh(kmer);};
                        if(enct == BONSAI)      for_each_substr([&](const char *s) {enc.for_each(lfunc, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);
                        else if(enct == NTHASH) for_each_substr([&](const char *s) {enc.for_each_hash(lfunc, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);
                        else                    for_each_substr([&](const char *s) {rolling_hasher.for_each_hash(lfunc, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);
                        cm.clear();
                    }
                    CONST_IF(!samesketch) new(final_sketches + i) final_type(std::move(sketch)); 
                    CONST_IF(samesketch) {
                        if(cache_sketch && !isf) sketch.write(fpath);
                    } else if(cache_sketch) final_sketches[i].write(fpath);
                }
            }
            ++ncomplete; // Atomic
        }
    }
    OMP_PFOR
    for(size_t i = 0; i < npaths; ++i) {
        try {
            sketch_finalize(final_sketches[i]);
        } catch(const std::exception &ex) {
            std::cerr << "Failed to finalize sketch " << i << " for value " << inpaths[i] << ".\n"
            << "msg: " << ex.what() << '\n';
        }
    }
    kseqs.free();
    auto fbuf = std::make_unique<float[]>(npaths);
    OMP_PFOR
    for(size_t i = 0; i < npaths; ++i) {
        CONST_IF(samesketch) fbuf[i] = cardinality_estimate(sketches[i]);
        else                 fbuf[i] = cardinality_estimate(final_sketches[i]);
    }
    if(emit_binary) {
        if(std::fwrite(&fbuf[0], sizeof(float), npaths, ofp) != npaths) {
            throw std::system_error(std::ferror(ofp), std::system_category(), "Failed to write cardinality estimates to file");
        }
    } else {
        ks::string str("#Path\tSize (est.)\n");
        str.resize(BUFFER_FLUSH_SIZE);
        {
            const int fn(fileno(ofp));
            constexpr const char *scinotstr = "%s\t%0.12g\n";
            constexpr const char *stdnotstr = "%s\t%0.8f\n";
            const char *const ptr = use_scientific ? scinotstr: stdnotstr;
            for(size_t i(0); i < npaths; ++i) {
                str.sprintf(ptr, inpaths[i].data(), fbuf[i]);
                if(str.size() >= BUFFER_FLUSH_SIZE) str.flush(fn);
            }
            str.flush(fn);
        }
    }
    if(ofp != stdout) std::fclose(ofp);
    CONST_IF(!samesketch) {
        if(!raii_final_sketches) {
            OMP_PFOR
            for(size_t i = 0; i < npaths; ++i) {
                using T = typename std::decay<decltype(final_sketches[0])>::type;
                final_sketches[i].~T();
            }
            std::free(final_sketches);
        }
    }
    OMP_PFOR
    for(size_t i = 0; i < npaths; ++i) {
        sketches[i].~SketchType();
    }
} // size_sketch_and_emit


template<typename SketchType>
void dist_sketch_and_cmp(std::vector<std::string> &inpaths, std::vector<CountingSketch> &cms, KSeqBufferHolder &kseqs, std::FILE *ofp, std::FILE *&pairofp, std::string outpath,
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
    static_assert(!sketch::wj::is_weighted_sketch<final_type>::value, "Can't have a weighted sketch be a final type");
    std::vector<SketchType> sketches;
    sketches.reserve(inpaths.size());
    const uint32_t sketch_size = bytesl2_to_arg(ssarg, SketchEnum<SketchType>::value);
    while(sketches.size() < inpaths.size()) {
        sketches.emplace_back(construct<SketchType>(sketch_size));
        set_estim_and_jestim(sketches.back(), estim, jestim);
    }
    static constexpr bool samesketch = std::is_same<SketchType, final_type>::value;
    final_type *final_sketches;
    std::unique_ptr<std::vector<final_type>> raii_final_sketches;

    std::atomic<uint32_t> ncomplete;
    ncomplete.store(0);
    const unsigned k = sp.k_;
    const unsigned wsz = sp.w_;
    RollingHasher<uint64_t> rolling_hasher(k, canon);
    if(inpaths.size() == 1 && presketched_only) {
        raii_final_sketches.reset(new std::vector<final_type>);
        gzFile ifp = gzopen(inpaths[0].data(), "rb");
        if(!ifp) UNRECOVERABLE_ERROR("Failed to open file.");
        for(;;) {
            try {
                //TD<decltype(raii_final_sketches)> td;
                raii_final_sketches->emplace_back(ifp);
                set_estim_and_jestim(raii_final_sketches->back(), estim, jestim);
            } catch(...) {
                break;
            }
        }
        final_sketches = raii_final_sketches->data();
        while(inpaths.size() < raii_final_sketches->size())
            inpaths.emplace_back(std::to_string(inpaths.size()));
        inpaths[0] = "0";
    } else {
        final_sketches =
            samesketch ? reinterpret_cast<final_type *>(sketches.data())
                       : static_cast<final_type *>(std::malloc(sizeof(*final_sketches) * inpaths.size()));
        OMP_PFOR_DYN
        for(size_t i = 0; i < sketches.size(); ++i) {
            const std::string &path(inpaths[i]);
            auto &sketch = sketches[i];
            if(presketched_only)  {
                CONST_IF(samesketch) {
                    sketch.read(path);
                    set_estim_and_jestim(sketch, estim, jestim); // HLL is the only type that needs this, and it's the same
                } else {
                    //TD<final_type> td;
                    new(final_sketches + i) final_type(path.data()); // Read from path
                }
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
                    Encoder<score::Lex> enc(nullptr, 0, sp, nullptr, canon);
                    if(cms.empty()) {
                        auto &h = sketch;
                        if(enct == BONSAI) for_each_substr([&](const char *s) {enc.for_each([&](u64 kmer){h.addh(kmer);}, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);
                        else if(enct == NTHASH) for_each_substr([&](const char *s) {enc.for_each_hash([&](u64 kmer){h.addh(kmer);}, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);
                        else for_each_substr([&](const char *s) {rolling_hasher.for_each_hash([&](u64 kmer){h.addh(kmer);}, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);
                    } else {
                        CountingSketch &cm = cms.at(tid);
                        const auto lfunc = [&](u64 kmer){if(cm.addh(kmer) >= mincount) sketch.addh(kmer);};
                        if(enct == BONSAI)      for_each_substr([&](const char *s) {enc.for_each(lfunc, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);
                        else if(enct == NTHASH) for_each_substr([&](const char *s) {enc.for_each_hash(lfunc, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);
                        else                    for_each_substr([&](const char *s) {rolling_hasher.for_each_hash(lfunc, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);
                        cm.clear();
                    }
                    CONST_IF(!samesketch) new(final_sketches + i) final_type(std::move(sketch)); 
                    CONST_IF(samesketch) {
                        if(cache_sketch && !isf) sketch.write(fpath);
                    } else if(cache_sketch) final_sketches[i].write(fpath);
                }
            }
            ++ncomplete; // Atomic
        }
    }
    OMP_PFOR
    for(size_t i = 0; i < sketches.size(); ++i) {
        try {
            sketch_finalize(final_sketches[i]);
        } catch(const std::exception &ex) {
            std::cerr << "Failed to finalize sketch " << i << " for value " << inpaths[i] << ".\n"
            << "msg: " << ex.what() << '\n';
        }
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
    if(emit_fmt == UT_TSV && !nq) {
        str.sprintf("##Names\t");
        for(size_t i = 0; i < inpaths.size() - nq; ++i)
            str.sprintf("%s\t", inpaths[i].data());
        str.back() = '\n';
        str.write(fileno(pairofp)); str.free();
    } else if(emit_fmt == UPPER_TRIANGULAR) { // emit_fmt == UPPER_TRIANGULAR
        std::fprintf(pairofp, "%zu\n", inpaths.size());
        std::fflush(pairofp);
    }
    if(emit_fmt & NEAREST_NEIGHBOR_TABLE) {
        std::fprintf(stderr, "[%s] About to make nn table with result type = %s. Number inpaths: %zu. \n", __PRETTY_FUNCTION__, emt2str(result_type), inpaths.size());
        nndist_loop(pairofp, final_sketches, inpaths, k, result_type, emit_fmt, nq);
    } else {
        dist_loop<final_type>(pairofp, outpath, final_sketches, inpaths, use_scientific, k, result_type, emit_fmt, nthreads, BUFFER_FLUSH_SIZE, nq);
    }
    CONST_IF(!samesketch) {
        if(!raii_final_sketches) {
#if __cplusplus >= 201703L
            std::destroy_n(final_sketches, inpaths.size());
#else
            std::for_each(final_sketches, final_sketches + inpaths.size(), [](auto &sketch) {
                using T = typename std::decay<decltype(sketch)>::type;
                sketch.~T();
            });
#endif
            std::free(final_sketches);
        }
    }
} // dist_sketch_and_cmp
#define DECSKETCHCMP(DS) \
template void ::bns::dist_sketch_and_cmp<DS>(std::vector<std::string> &inpaths, std::vector<::bns::CountingSketch> &cms, KSeqBufferHolder &kseqs, std::FILE *ofp, std::FILE *&pairofp,\
                   std::string,\
                   Spacer sp,\
                   unsigned ssarg, unsigned mincount, EstimationMethod estim, JointEstimationMethod jestim, bool cache_sketch, EmissionType result_type, EmissionFormat emit_fmt,\
                   bool presketched_only, unsigned nthreads, bool use_scientific, std::string suffix, std::string prefix, bool canon, bool entropy_minimization, std::string spacing,\
                   size_t nq, EncodingType enct);\
template void ::bns::dist_sketch_and_cmp<sketch::wj::WeightedSketcher<DS>>(std::vector<std::string> &inpaths, std::vector<::bns::CountingSketch> &cms, KSeqBufferHolder &kseqs, std::FILE *ofp, std::FILE *&pairofp,\
                   std::string,\
                   Spacer sp,\
                   unsigned ssarg, unsigned mincount, EstimationMethod estim, JointEstimationMethod jestim, bool cache_sketch, EmissionType result_type, EmissionFormat emit_fmt,\
                   bool presketched_only, unsigned nthreads, bool use_scientific, std::string suffix, std::string prefix, bool canon, bool entropy_minimization, std::string spacing,\
                   size_t nq, EncodingType enct);\
template void ::bns::dist_sketch_and_cmp<sketch::wj::WeightedSketcher<DS, wj::ExactCountingAdapter>>(std::vector<std::string> &inpaths, std::vector<::bns::CountingSketch> &cms, KSeqBufferHolder &kseqs, std::FILE *ofp, std::FILE *&pairofp,\
                   std::string,\
                   Spacer sp,\
                   unsigned ssarg, unsigned mincount, EstimationMethod estim, JointEstimationMethod jestim, bool cache_sketch, EmissionType result_type, EmissionFormat emit_fmt,\
                   bool presketched_only, unsigned nthreads, bool use_scientific, std::string suffix, std::string prefix, bool canon, bool entropy_minimization, std::string spacing,\
                   size_t nq, EncodingType enct);


enum SketchFlags {
    SKIP_CACHED  = 1,
    CANONICALIZE = 2,
    ENTROPY_MIN  = 4
};

template<typename SketchType>
INLINE void sketch_core(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp,
                        const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix,
                        std::vector<CountingSketch> &cms, EstimationMethod estim, JointEstimationMethod jestim,
                        KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing,
                        int sketchflags, uint32_t mincount, EncodingType enct, std::string output_file)
{
    const auto canon = sketchflags & CANONICALIZE, skip_cached = sketchflags & SKIP_CACHED, entropy_minimization = sketchflags & ENTROPY_MIN;
    std::vector<SketchType> sketches;
    uint32_t sketch_size = bytesl2_to_arg(ssarg, SketchEnum<SketchType>::value);
    if(output_file.empty()) {
        while(sketches.size() < (u32)nthreads) sketches.push_back(construct<SketchType>(sketch_size)), set_estim_and_jestim(sketches.back(), estim, jestim);
    } else {
        while(sketches.size() < inpaths.size()) sketches.push_back(construct<SketchType>(sketch_size)), set_estim_and_jestim(sketches.back(), estim, jestim);
    }
    std::vector<std::string> fnames(nthreads);
    RollingHasher<uint64_t> rolling_hasher(k, canon);

    if(entropy_minimization)
        UNRECOVERABLE_ERROR("Unsupported option: entropy_minimization.");
    gzFile outputfp = nullptr;
    if(output_file.size()) {
        if((outputfp = gzopen((output_file + ".labels.gz").data(), "w")) == nullptr) UNRECOVERABLE_ERROR("Failed to write sequence labels to file");
        for(const auto &path: inpaths) {
            gzwrite(outputfp, path.data(), path.size());
            gzputc(outputfp, '\n');
        }
        gzclose(outputfp);
        LOG_DEBUG("Wrote labels to file\n");
    }
    {
        if(outputfp) {
            if((outputfp = gzopen(output_file.data(), "w")) == nullptr)
                UNRECOVERABLE_ERROR("Failed to write sketches to file");
#ifndef NDEBUG
            else std::fprintf(stderr, "Opened file at %s\n", output_file.data());
#endif
        }
        OMP_PFOR_DYN
        for(size_t i = 0; i < inpaths.size(); ++i) {
            const int tid = omp_get_thread_num();
            std::string &fname = fnames[tid];
            fname = make_fname<SketchType>(inpaths[i].data(), sketch_size, wsz, k, sp.c_, spacing, suffix, prefix, enct);
            LOG_DEBUG("fname: %s from %s with tid = %d\n", fname.data(), inpaths[i].data(), tid);
            unsigned sketchind = outputfp ? unsigned(i): tid;
            auto &h = sketches[sketchind];
            if(skip_cached && isfile(fname)) {
                if(outputfp) h.read(fname);
                else continue;
            }
            Encoder<bns::score::Lex> enc(nullptr, 0, sp, nullptr, canon);
            const auto &path = inpaths[i];
            if(use_filter.size() && use_filter[i]) {
                auto &cm = cms[tid];
                auto lfunc = [&](u64 kmer){if(cm.addh(kmer) >= mincount) h.add(kmer);};
                auto hlfunc = [&](u64 kmer){if(cm.addh(kmer) >= mincount) h.addh(kmer);};
#define FOR_EACH_FUNC(wrapperfunc) for_each_substr([&](const char *s) {wrapperfunc(lfunc, s, &kseqs[tid]);}, path, FNAME_SEP)
#define FOR_EACH_HASH_FUNC(wrapperfunc) for_each_substr([&](const char *s) {wrapperfunc(hlfunc, s, &kseqs[tid]);}, path, FNAME_SEP)
                if(enct == NTHASH)      FOR_EACH_FUNC(enc.for_each_hash);
                else if(enct == BONSAI) FOR_EACH_HASH_FUNC(enc.for_each);
                else                    FOR_EACH_FUNC(rolling_hasher.for_each_hash);
                cm.clear();
            } else {
                auto lfunc = [&](u64 kmer){h.add(kmer);};
                auto hlfunc = [&](u64 kmer){h.addh(kmer);};
                if(enct == NTHASH)      FOR_EACH_FUNC(enc.for_each_hash);
                else if(enct == BONSAI) FOR_EACH_HASH_FUNC(enc.for_each);
                else                    FOR_EACH_FUNC(rolling_hasher.for_each_hash);
            }
#undef FOR_EACH_FUNC
            sketch_finalize(h);
            if(!outputfp) h.write(fname.data()), h.clear();
        }
        if(outputfp) {
#ifndef NDEBUG
            size_t i = 0;
            auto it = inpaths.cbegin();
#endif
            for(const auto &sketch: sketches) {
#ifndef NDEBUG
                std::fprintf(stderr, "Writing sketch %zu/%zu from path %s to file\n", i++, inpaths.size(), (it++)->data());
#endif
                sketch.write(outputfp);
            }
            gzclose(outputfp);
        }
    }
}
template<typename SketchType>
INLINE void sketch_by_seq_core(uint32_t ssarg, uint32_t nthreads, const Spacer &sp,
                        const std::string &inpath, const std::string &outpath,
                        CountingSketch *cm, EstimationMethod estim, JointEstimationMethod jestim,
                        bool use_filter,
                        int sketchflags,
                        uint32_t mincount, EncodingType enct)
{
    const auto canon = sketchflags & CANONICALIZE;
    const uint32_t sketch_size = bytesl2_to_arg(ssarg, SketchEnum<SketchType>::value);
    const auto k = sp.k_;
    SketchType working_sketch(construct<SketchType>(sketch_size));
    set_estim_and_jestim(working_sketch, estim, jestim);
    RollingHasher<uint64_t> rolling_hasher(k, canon);
    Encoder<bns::score::Lex> enc(nullptr, 0, sp, nullptr, canon);
    gzFile fp = gzopen(inpath.data(), "rb");
    if(!fp) throw ZlibError(std::string("Failed to open file for reading at ") + inpath);
    gzFile ofp = gzopen(outpath.data(), "wb");
    if(!ofp) throw ZlibError(std::string("Failed to open file for writing at ") + outpath);
    std::string namepath = outpath == "/dev/stdout"
        ? std::string("stdout.names")
        : outpath + ".names";
    std::FILE *nameofp = fopen(namepath.data(), "w");
    if(!nameofp) UNRECOVERABLE_ERROR(std::string("Failed to open file for writing at ") + namepath);
    fprintf(nameofp, "#k=%d:Names for sequences sketched\n", k);
    kseq_t *ks = kseq_init(fp);
    auto &h = working_sketch; // just alias for less typing
    auto add = [&](u64 kmer) {h.addh(kmer);};
    auto cadd = [&](u64 kmer){if(cm->addh(kmer) >= mincount) h.addh(kmer);};
    std::vector<std::string> seqnames;
    while(kseq_read(ks) >= 0) {
        if(use_filter) {
            if(enct == NTHASH)
                enc.for_each_hash(cadd, ks->seq.s, ks->seq.l);
            else if(enct == BONSAI)
                enc.for_each(cadd, ks->seq.s, ks->seq.l);
            else
                rolling_hasher.for_each_hash(cadd, ks->seq.s, ks->seq.l);
            cm->clear();
        } else {
            if(enct == NTHASH)
                enc.for_each_hash(add, ks->seq.s, ks->seq.l);
            else if(enct == BONSAI)
                enc.for_each(add, ks->seq.s, ks->seq.l);
            else
                rolling_hasher.for_each_hash(add, ks->seq.s, ks->seq.l);
        }
        try {
            sketch_finalize(h);
        } catch(const std::exception &ex) {
            std::cerr << "Failed to finalize sketch  for sequence " << ks->seq.s << ".\n"
            << "msg: " << ex.what() << '\n';
        }
        if(std::fwrite(ks->name.s, 1, ks->name.l, nameofp) != ks->name.l) std::fprintf(stderr, "Warning: error in writing sequence name\n");
        std::fputc('\n', nameofp);
        h.write(ofp);
        h.clear();
    }
    kseq_destroy(ks);
    gzclose(fp);
    gzclose(ofp);
    std::fclose(nameofp);
}


using validx_t = std::pair<float, uint32_t>;

template<typename Cmp>
INLINE void lock_update(float val, validx_t *ptr,
                        unsigned nneighbors,
                        size_t j, std::mutex &mut, const Cmp &cmp)
{
    if(cmp(val, ptr->first)) {
        std::lock_guard<std::mutex> lg(mut);
        if(cmp(val, ptr->first)) { // after getting the lock, check again
            std::pop_heap (ptr, ptr + nneighbors, cmp);
            ptr[nneighbors - 1] = {val, j};
            std::push_heap(ptr, ptr + nneighbors, cmp);
        }
    }
}
template<typename Cmp>
INLINE void lockfree_update(float val, validx_t *ptr,
                            unsigned nneighbors,
                            size_t j, const Cmp &cmp)
{
    if(cmp(val, ptr->first)) {
        std::pop_heap (ptr, ptr + nneighbors, cmp);
        ptr[nneighbors - 1] = validx_t(val, j);
        std::push_heap(ptr, ptr + nneighbors, cmp);
    }
}


template<typename SketchType, typename Cmp, typename Func>
void perform_nns(validx_t *neighbors,
                 SketchType *sketches, const std::vector<std::string> &inpaths,
                 const unsigned k, const EmissionType result_type,
                 size_t nq,
                 unsigned nneighbors,
                 const Cmp &cmp, const Func &func) {
    float default_value = std::numeric_limits<float>::max();
    if(std::is_same<Cmp, std::greater<>>::value)
        default_value = -default_value;
    const size_t n = nq ? nq: inpaths.size();
    std::fprintf(stderr, "default value: %g\n", default_value);
    OMP_PFOR
    for(size_t i = 0; i < n; ++i)
        std::fill_n(&neighbors[i * nneighbors], nneighbors,
                  validx_t(default_value, uint32_t(-1)));
    if(nq == 0) {
        auto mutexes = std::make_unique<std::mutex[]>(n);
        OMP_PFOR_DYN
        for(size_t i = 0; i < n; ++i) {
            auto lhptr = &neighbors[i * nneighbors];
            const auto &h1 = sketches[i];
            for(size_t j = i + 1; j < n; ++j) {
                auto rhptr = &neighbors[j * nneighbors];
                const float distval = func(sketches[j], h1);
                lock_update(distval, lhptr, nneighbors, j, mutexes[i], cmp);
                lock_update(distval, rhptr, nneighbors, i, mutexes[j], cmp);
            }
        }
    } else {
        const size_t npaths = inpaths.size(), nr = npaths - nq;
        OMP_PFOR_DYN
        for(size_t qi = nr; qi < npaths; ++qi) {
            const size_t qind = qi - nr;
            const auto srcptr = &neighbors[qind * nneighbors];
            const auto &h1 = sketches[qi];
            for(size_t j = 0; j < nr; ++j) {
                lockfree_update(func(sketches[j], h1), srcptr, nneighbors, j, cmp);
            }
        }
    }
    LOG_DEBUG("Finished loop, now sorting\n");
    OMP_PFOR
    for(size_t i = 0; i < n; ++i) {
        auto start = neighbors + (i * nneighbors), end = start + nneighbors;
        std::sort(start, end, cmp);
#if VERBOSE_AF
        for(auto it = start; it < end; ++it) {
            std::fprintf(stderr, "%zu/%u/%g\n", i, it->second, it->first);
        }
#endif
    }
#ifndef NDEBUG
    if(std::find_if(neighbors, neighbors + nneighbors * n, [](auto x) {return x.second == uint32_t(-1);}) != neighbors + nneighbors * n)
        throw std::runtime_error("Empty\n");
#endif
}

template<bool MT=true, typename SketchType, typename T, typename Func>
inline void perform_core_op(T &dists, size_t nsketches, SketchType *sketches, const Func &func, size_t i) {
    auto &h1 = sketches[i];
#define compute_j(j) do {dists[j - i - 1] = func(sketches[j], h1);} while(0)
    if(MT) {
        OMP_PFOR_DYN
        for(size_t j = i + 1; j < nsketches; ++j) compute_j(j);
    } else {
        for(size_t j = i + 1; j < nsketches; ++j) compute_j(j);
    }
#undef compute_j
}

template<typename SketchType>
void nndist_loop(std::FILE *ofp, SketchType *sketches,
               const std::vector<std::string> &inpaths,
               const unsigned k, const EmissionType result_type, EmissionFormat emit_fmt,
               size_t nq) {
    gargs.show();
    unsigned nneighbors = gargs.number_neighbors;
    size_t npairs = (nq ? nq: inpaths.size());
    size_t possible_num_neighbors = nq ? inpaths.size() - nq: inpaths.size();
    if(nneighbors > possible_num_neighbors) {
        std::fprintf(stderr, "Only reporting %zu rather than %u neighbors due to their being only that many sets.\n", possible_num_neighbors, nneighbors);
        nneighbors = possible_num_neighbors;
    }
    const size_t qoffset = nq ? inpaths.size() - nq: size_t(0);
    const size_t ntups = npairs * nneighbors;
    auto neighbors = std::make_unique<validx_t[]>(ntups);
    std::fprintf(stderr, "made %zu pairs with %zu nq, %u neighbors and %zu tups\n", npairs, nq, nneighbors, ntups);
    const double ksinv = 1./ k;
    auto call_cmp = [result_type, ksinv](const auto &x, const auto &y) {return result_cmp(x, y, result_type, ksinv);};

    if(emt2nntype(result_type) == SIMILARITY_MEASURE) {
        std::fprintf(stderr, "Performing nn under similarity measure\n");
        perform_nns(neighbors.get(), sketches, inpaths, k, result_type, nq, nneighbors, std::greater<>(), call_cmp);
    } else {
        std::fprintf(stderr, "Performing nn under dissimilarity measure\n");
        perform_nns(neighbors.get(), sketches, inpaths, k, result_type, nq, nneighbors, std::less<>(), call_cmp);
    }
    if(emit_fmt & BINARY) {
        uint32_t n = inpaths.size();
        std::fwrite(&n, sizeof(n), 1, ofp);
        n = nneighbors;
        std::fwrite(&n, sizeof(n), 1, ofp);
        size_t nb = nneighbors * inpaths.size();
        if(unlikely(std::fwrite(neighbors.get(), sizeof(validx_t), nb, ofp) != nb))
            UNRECOVERABLE_ERROR("Failed to write neighbors to disk (binary)\n");
    } else {
        std::fprintf(ofp, "#File\tNeighbor ID:distance\t...\n");
        int nt = 1;
#ifdef _OPENMP
        _Pragma("omp parallel")
        {
            _Pragma("omp single")
            nt = omp_get_num_threads();
        }
#endif
        std::vector<ks::string> kstrs;
        while(kstrs.size() < unsigned(nt))
            kstrs.emplace_back(1ull<<10);
        std::fprintf(stderr, "Made buffer, ones per thread\n");
        const size_t npaths = inpaths.size();
        OMP_PFOR
        for(size_t i = 0; i < npaths; ++i) {
            auto tid = OMP_ELSE(omp_get_thread_num(), 0);
            auto &buf(kstrs[tid]);
            const validx_t *nptr = &neighbors[i * nneighbors];
            assert(i * nneighbors + nneighbors < npairs);
            size_t nameind = i + qoffset;
            assert(nameind < inpaths.size());
            buf += inpaths[nameind];
            for(unsigned j = 0; j < nneighbors; ++j) {
                assert(nptr + j < neighbors.get() + ntups || std::fprintf(stderr, "i = %zu, j = %zu, offset = %zu\n", i, j, i * nneighbors + j));
                const auto ndat = nptr[j];
                buf.putc_('\t');
                buf.putw_(ndat.second);
                buf.putc_(':');
                buf.sprintf("%g", ndat.first);
            }
            buf.putc_('\n');
            if(buf.size() > (1u << 14)) {
                OMP_CRITICAL
                {
                    buf.flush(ofp);
                }
            }
        }
        for(auto buf: kstrs) buf.flush(ofp);
    }
}

template<typename SketchType>
void dist_loop(std::FILE *&ofp, std::string ofp_name, SketchType *sketches, const std::vector<std::string> &inpaths, const bool use_scientific, const unsigned k, const EmissionType result_type, EmissionFormat emit_fmt, int, const size_t buffer_flush_size, size_t nq) {
    if(nq) {
        partdist_loop<SketchType>(ofp, sketches, inpaths, use_scientific, k, result_type, emit_fmt, buffer_flush_size, nq);
        return;
    }
    if(!is_symmetric(result_type)) {
        char buf[1024];
        std::sprintf(buf, "Can't perform symmetric distance comparisons with a symmetric method (%s/%d). To perform an asymmetric distance comparison between a given set and itself, provide the same list of filenames to both -Q and -F.\n", emt2str(result_type), int(result_type));
        UNRECOVERABLE_ERROR(buf);
    }
    const float ksinv = 1./ k;
    const int pairfi = fileno(ofp);
    const size_t nsketches = inpaths.size();
    std::future<size_t> submitter;
    auto cmp = [ksinv,result_type](const auto &x, const auto &y) {return result_cmp(x, y, result_type, ksinv);};
    if((emit_fmt & BINARY) == 0) {
        std::array<std::vector<float>, 2> dps;
        dps[0].resize(nsketches - 1);
        dps[1].resize(std::max(ssize_t(nsketches) - 2, ssize_t(1)));
        ks::string str;
        for(size_t i = 0; i < nsketches; ++i) {
            std::vector<float> &dists = dps[i & 1];
            perform_core_op<true>(dists, nsketches, sketches, cmp, i);
            //LOG_DEBUG("Finished chunk %zu of %zu\n", i + 1, nsketches);
            if(i) submitter.get();
            submitter = std::async(std::launch::async, submit_emit_dists<float>,
                                   pairfi, dists.data(), nsketches, i,
                                   std::ref(str), std::ref(inpaths), emit_fmt, use_scientific, buffer_flush_size);
        }
        submitter.get();
    } else {
        const float defv = static_cast<float>(emt2nntype(result_type) == SIMILARITY_MEASURE);
        // 1. for the diagonal for similarity, 0 for dissimilarity
        if(emit_fmt == FULL_TSV || ::isatty(::fileno(ofp))) {
            dm::DistanceMatrix<float> dm(nsketches, defv);
            OMP_PFOR
            for(size_t i = 0; i < nsketches - 1; ++i) {
                auto span = dm.row_span(i);
                auto &dists = span.first;
                perform_core_op<false>(dists, nsketches, sketches, cmp, i);
            }
            if(emit_fmt == FULL_TSV)
                dm.printf(ofp, use_scientific, &inpaths);
            else
                dm.write(ofp);
        } else {
            // Resize file
            ::ftruncate(::fileno(ofp), 1 + sizeof(uint64_t) + ((nsketches * (nsketches - 1)) >> 1));
            std::fclose(ofp);
            ofp = nullptr;
            std::fprintf(stderr, "Setting up distance matrix on disk, with %zu sketches\n", nsketches);
            // Modify in-place
            dm::DistanceMatrix<float> dm(ofp_name.data(), nsketches, defv);
            std::fprintf(stderr, "Created dm. %zu elem and %zu entrie\n", size_t(dm.nelem()), size_t(dm.num_entries()));
            dm::parallel_fill(dm, nsketches, [&cmp,sketches](size_t i, size_t j) {return cmp(sketches[i], sketches[j]);});
            float *dmp = dm.data();
            if(int rc = ::madvise(static_cast<void *>(dmp), dm.num_entries() * sizeof(float), MADV_SEQUENTIAL))
                std::fprintf(stderr, "Note: madvise call returned %d/%s, inessential\n", rc, std::strerrr(rc));
            for(size_t i = 0; i < nsketches - 1; ++i) {
                auto span = dm.row_span(i);
                std::fprintf(stderr, "Row %zu has ptr %p and %zu\n", i, (void *)span.first, span.second);
                std::unique_ptr<float[]> subrow(new float[span.second]);
                auto sp = subrow.get();
                auto &s1 = sketches[i];
                OMP_PFOR_DYN
                for(size_t j = i + 1; j < nsketches; ++j) {
                    sp[j - i - 1] = cmp(sketches[j], s1);
                }
                if(i) submitter.get();
                submitter = std::async(std::launch::async, [n=span.second,ofp,sp=std::move(subrow),&dmp]() -> size_t {
                    std::memcpy(dmp, sp.get(), n * sizeof(float));
                    dmp += n;
                    assert(dmp <= dm.data() + dm.num_entries());
                    return n * sizeof(float);
                });
            }
        }
    }
}
#define DECSKETCHCORE(DS) \
  template void sketch_core<DS>(uint32_t ssarg, uint32_t nthreads,\
                                uint32_t wsz, uint32_t k, const Spacer &sp,\
                                const std::vector<std::string> &inpaths,\
                                const std::string &suffix,\
                                const std::string &prefix, std::vector<CountingSketch> &counting_sketches,\
                                EstimationMethod estim, JointEstimationMethod jestim,\
                                KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing,\
                                int sketchflags, uint32_t mincount, EncodingType enct, std::string s);\
  template void sketch_core<wj::WeightedSketcher<DS>>(uint32_t ssarg, uint32_t nthreads,\
                                uint32_t wsz, uint32_t k, const Spacer &sp,\
                                const std::vector<std::string> &inpaths,\
                                const std::string &suffix,\
                                const std::string &prefix, std::vector<CountingSketch> &counting_sketches,\
                                EstimationMethod estim, JointEstimationMethod jestim,\
                                KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing,\
                                int sketchflags, uint32_t mincount, EncodingType enct, std::string s);\
  template void sketch_core<wj::WeightedSketcher<DS, wj::ExactCountingAdapter>>(uint32_t ssarg, uint32_t nthreads,\
                                uint32_t wsz, uint32_t k, const Spacer &sp,\
                                const std::vector<std::string> &inpaths,\
                                const std::string &suffix,\
                                const std::string &prefix, std::vector<CountingSketch> &counting_sketches,\
                                EstimationMethod estim, JointEstimationMethod jestim,\
                                KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing,\
                                int sketchflags, uint32_t mincount, EncodingType enct, std::string s);\


} // namespace bns
