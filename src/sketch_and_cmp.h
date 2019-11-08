#pragma once
#include "dashing.h"

#define FILL_SKETCH_MIN(MinType)  \
    {\
        Encoder<MinType> enc(nullptr, 0, sp, nullptr, canon);\
        if(cms.empty()) {\
            auto &h = sketch;\
            if(enct == BONSAI) for_each_substr([&](const char *s) {enc.for_each([&](u64 kmer){h.addh(kmer);}, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);\
            else if(enct == NTHASH) for_each_substr([&](const char *s) {enc.for_each_hash([&](u64 kmer){h.addh(kmer);}, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);\
            else for_each_substr([&](const char *s) {rolling_hasher.for_each_hash([&](u64 kmer){h.addh(kmer);}, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);\
        } else {\
            CountingSketch &cm = cms.at(tid);\
            const auto lfunc = [&](u64 kmer){if(cm.addh(kmer) >= mincount) sketch.addh(kmer);};\
            if(enct == BONSAI)      for_each_substr([&](const char *s) {enc.for_each(lfunc, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);\
            else if(enct == NTHASH) for_each_substr([&](const char *s) {enc.for_each_hash(lfunc, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);\
            else                    for_each_substr([&](const char *s) {rolling_hasher.for_each_hash(lfunc, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);\
            cm.clear();\
        }\
        CONST_IF(!samesketch) new(final_sketches + i) final_type(std::move(sketch)); \
    }

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
void dist_loop(std::FILE *ofp, SketchType *hlls, const std::vector<std::string> &inpaths, const bool use_scientific, const unsigned k, const EmissionType result_type, EmissionFormat emit_fmt, int nthreads, const size_t buffer_flush_size, size_t nq);
using namespace sketch;
using namespace hll;
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
void dist_by_seq(std::vector<std::string> &labels, std::string datapath,
                 std::FILE *pairofp, int k,
                 EstimationMethod estim, JointEstimationMethod jestim, EmissionType result_type, EmissionFormat emit_fmt,
                 unsigned nthreads, std::string otherpath) {
    gzFile sfp = gzopen(datapath.data(), "rb");
    if(!sfp) throw "up";
    std::vector<SketchType> sketches;
    std::vector<std::string> qnames;
    sketches.reserve(labels.size());
    while(sketches.size() < labels.size()) {
        sketches.emplace_back(sfp);
        set_estim_and_jestim(sketches.back(), estim, jestim);
    }
    if(otherpath.size()) {
        qnames = get_paths((otherpath + ".names").data());
        if(qnames.empty()) RUNTIME_ERROR("Can't compare with empty qnames");
        sketches.reserve(qnames.size() + sketches.size());
        for(size_t i = 0; i < qnames.size(); ++i) {
            sketches.emplace_back(sfp);
            set_estim_and_jestim(sketches.back(), estim, jestim);
        }
        labels.insert(labels.end(), qnames.begin(), qnames.end());
    } else if(!is_symmetric(result_type)) {
        RUNTIME_ERROR("Can't perform asymmetric comparison without query paths");
    }
    const size_t nq = qnames.size();
    gzclose(sfp);
    ks::string str;
    if(emit_fmt == UT_TSV) {
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
    dist_loop<SketchType>(pairofp, sketches.data(), labels, /* use_scientific=*/ true, k, result_type, emit_fmt, nthreads, BUFFER_FLUSH_SIZE, nq);
}


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
#define DECSKETCHCMP(DS) \
template void ::bns::dist_sketch_and_cmp<DS>(const std::vector<std::string> &inpaths, std::vector<::bns::CountingSketch> &cms, KSeqBufferHolder &kseqs, std::FILE *ofp, std::FILE *pairofp,\
                   Spacer sp,\
                   unsigned ssarg, unsigned mincount, EstimationMethod estim, JointEstimationMethod jestim, bool cache_sketch, EmissionType result_type, EmissionFormat emit_fmt,\
                   bool presketched_only, unsigned nthreads, bool use_scientific, std::string suffix, std::string prefix, bool canon, bool entropy_minimization, std::string spacing,\
                   size_t nq, EncodingType enct);\
template void ::bns::dist_sketch_and_cmp<sketch::wj::WeightedSketcher<DS>>(const std::vector<std::string> &inpaths, std::vector<::bns::CountingSketch> &cms, KSeqBufferHolder &kseqs, std::FILE *ofp, std::FILE *pairofp,\
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
                        int sketchflags, uint32_t mincount, EncodingType enct)
{
    const auto canon = sketchflags & CANONICALIZE, skip_cached = sketchflags & SKIP_CACHED, entropy_minimization = sketchflags & ENTROPY_MIN;
    std::vector<SketchType> sketches;
    uint32_t sketch_size = bytesl2_to_arg(ssarg, SketchEnum<SketchType>::value);
    while(sketches.size() < (u32)nthreads) sketches.push_back(construct<SketchType>(sketch_size)), set_estim_and_jestim(sketches.back(), estim, jestim);
    std::vector<std::string> fnames(nthreads);
    RollingHasher<uint64_t> rolling_hasher(k, canon);

    if(entropy_minimization)
        throw std::runtime_error("Removed.");
    #pragma omp parallel for schedule(dynamic)
    for(size_t i = 0; i < inpaths.size(); ++i) {
        const int tid = omp_get_thread_num();
        std::string &fname = fnames[tid];
        fname = make_fname<SketchType>(inpaths[i].data(), sketch_size, wsz, k, sp.c_, spacing, suffix, prefix, enct);
        LOG_DEBUG("fname: %s from %s\n", fname.data(), inpaths[i].data());
        if(skip_cached && isfile(fname)) continue;
        Encoder<bns::score::Lex> enc(nullptr, 0, sp, nullptr, canon);
        auto &h = sketches[tid];
        if(use_filter.size() && use_filter[i]) {
            auto &cm = cms[tid];
            if(enct == NTHASH) {
                for_each_substr([&](const char *s) {enc.for_each_hash([&](u64 kmer){if(cm.addh(kmer) >= mincount) h.add(kmer);}, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);
            } else if(enct == BONSAI) {
                for_each_substr([&](const char *s) {enc.for_each([&](u64 kmer){if(cm.addh(kmer) >= mincount) h.addh(kmer);}, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);
            } else {
                for_each_substr([&](const char *s) {rolling_hasher.for_each_hash([&](u64 kmer){if(cm.addh(kmer) >= mincount) h.addh(kmer);}, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);
            }
            cm.clear();  
        } else {
            if(enct == NTHASH) {
                for_each_substr([&](const char *s) {enc.for_each_hash([&](u64 kmer){h.add(kmer);}, inpaths[i].data(), &kseqs[tid]);}, inpaths[i], FNAME_SEP);
            } else if(enct == BONSAI) {
                for_each_substr([&](const char *s) {enc.for_each([&](u64 kmer){h.addh(kmer);}, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);
            } else {
                for_each_substr([&](const char *s) {rolling_hasher.for_each_hash([&](u64 kmer){h.addh(kmer);}, s, &kseqs[tid]);}, inpaths[i], FNAME_SEP);
            }
        }
        sketch_finalize(h);
        h.write(fname.data());
        h.clear();
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
    if(!nameofp) throw std::runtime_error(std::string("Failed to open file for writing at ") + namepath);
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
        sketch_finalize(h);
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

template<typename SketchType>
void dist_loop(std::FILE *ofp, SketchType *hlls, const std::vector<std::string> &inpaths, const bool use_scientific, const unsigned k, const EmissionType result_type, EmissionFormat emit_fmt, int nthreads, const size_t buffer_flush_size, size_t nq) {
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
        dps[1].resize(std::max(ssize_t(nsketches) - 2, ssize_t(1)));
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
            dm.write(ofp);
        }
    }
}
#define DECSKETCHCORE(DS) template void sketch_core<DS>(uint32_t ssarg, uint32_t nthreads,\
                                uint32_t wsz, uint32_t k, const Spacer &sp,\
                                const std::vector<std::string> &inpaths,\
                                const std::string &suffix,\
                                const std::string &prefix, std::vector<CountingSketch> &counting_sketches,\
                                EstimationMethod estim, JointEstimationMethod jestim,\
                                KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing,\
                                int sketchflags, uint32_t mincount, EncodingType enct);


} // namespace bns
