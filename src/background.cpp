#include "include/bonsai/encoder.h"
#include "src/substrs.h"
#include "background.h"
#ifdef _OPENMP
#include <omp.h>
#endif
namespace bns {

static inline auto background_match(unsigned lhid, unsigned rhid, const float *nucfreqs) {
    auto lhp = nucfreqs + (lhid * 4), rhp = nucfreqs + (rhid * 4);
#if __SSE2__
    auto mulval = _mm_mul_ps(_mm_loadu_ps(lhp), _mm_loadu_ps(rhp));
#if 0
    mulval = _mm_hadd_ps(mulval, mulval);
    mulval = _mm_hadd_ps(mulval, mulval);
    return _mm_cvtss_f32(mulval);
#else
    const float *p = reinterpret_cast<const float *>(&mulval);
    return p[0] + p[1] + p[2] + p[3];
#endif
#else
    std::fprintf(stderr, "background: %g %g %g %g %g %g %g %g\n",
        lhp[0], lhp[1], lhp[2], lhp[3],
        rhp[0], rhp[1], rhp[2], rhp[3]);
    return lhp[0] * rhp[0] + lhp[1] * rhp[1] +
           lhp[2] * rhp[2] + lhp[3] * rhp[3];
#endif
}


sumstats
nuc_freqs(const std::vector<std::string> &fpaths)
{
    sumstats rettuple(fpaths.size());
    auto &genome_sizes = rettuple.sizes();
    auto &ret = rettuple.freqs();
    auto &numseqs = rettuple.numseqs();
    const size_t np = fpaths.size();
    unsigned nt = 1;
    #pragma omp parallel
    {
        #pragma omp single
        nt = omp_get_num_threads();
    }
    KSeqBufferHolder kseqs(nt);
    #pragma omp parallel for
    for(size_t i = 0; i < np; ++i) {
        uint64_t arr [128]{0};
        size_t nkmer_windows = 0;
        unsigned nseq = 0;
        for_each_substr([&](const char *s) {
            std::fprintf(stderr, "string: %s\n", s);
            const int tid = omp_get_thread_num();
            gzFile fp = gzopen(s, "rb");
            if(!fp) {std::fprintf(stderr, "Missing file %s\n", s); std::exit(1);}
            kseq_t *ks = &kseqs[tid];
            while(kseq_read(ks) >= 0) {
                std::fprintf(stderr, "seqlen: %zu. name: %s\n", ks->seq.l, ks->name.s);
                nkmer_windows += ks->seq.l;
                ++nseq;
                for(size_t j = 0; j < ks->seq.l; ++arr[ks->seq.s[j++]]);
            }
            gzclose(fp);
        }, fpaths[i], FNAME_SEP);
        double total_valid = arr['A'] + arr['a'] + arr['t'] + arr['T'] +  arr['c'] + arr['C'] + arr['g'] + arr['G'];
        std::fprintf(stderr, "total valid: %g\n", total_valid);
        auto ptr = ret.data() + 4 * i;
        ptr[0] = (arr['A'] + arr['a']) / total_valid;
        ptr[1] = (arr['c'] + arr['C']) / total_valid;
        ptr[2] = (arr['g'] + arr['G']) / total_valid;
        ptr[3] = (arr['t'] + arr['T']) / total_valid;
        genome_sizes[i] = nkmer_windows;
        numseqs[i] = nseq;
    }
    return rettuple;
}

double jukes_cantor_p(const std::vector<unsigned> &ks, const float *values,
                      double background, double kmean, uint64_t lhsz, uint64_t rhsz,
                      unsigned lhns, unsigned rhns) {
    const size_t nk = ks.size();
    std::unique_ptr<double[]> lvals(new double[nk]);
    for(unsigned i = 0; i < nk; ++i) {
        const auto k = ks[i];
        auto corrected_size = [k](auto len, auto ns) {return len - uint64_t(k - 1) * ns;};
        auto lhc = corrected_size(lhsz, lhns), rhc = corrected_size(rhsz, rhns);
        lvals[i] = std::log(values[i] - std::pow(background, k) * 4 * lhc * rhc);
    }
    double lmean = std::accumulate(lvals.get(), &lvals[nk], 0.) / nk;
    double num = 0., denom = 0.;
    for(unsigned i = 0; i < nk; ++i) {
        auto kdiff = (ks[i] - kmean);
        denom += kdiff * kdiff;
        num += kdiff * (lvals[i] - lmean);
    }
    const double slope = num / denom;
    double p = std::exp(slope);
    return p;
}


dm::DistanceMatrix<float> mkmat2jcdistmat(std::string packedmatrix, EmissionType emtype,
                                          const uint64_t *genome_sizes,
                                          const float *corrected_background_rates,
                                          const unsigned *nseqs,
                                          bool reassign_to_distance)
{
    switch(emtype) {
        case MASH_DIST:
        case FULL_MASH_DIST:
        case SYMMETRIC_CONTAINMENT_INDEX:
        case SYMMETRIC_CONTAINMENT_DIST:
        case FULL_CONTAINMENT_DIST:
        case CONTAINMENT_INDEX:
        case JI:
        case CONTAINMENT_DIST:
            throw std::invalid_argument("EmissionType must be SIZES (intersection size)");
        case SIZES: std::fprintf(stderr, "Using intersection size\n"); break;
    }
    uint64_t number_sets, number_entries;
    gzFile ifp = gzopen(packedmatrix.data(), "rb");
    if(!ifp) UNRECOVERABLE_ERROR(std::string("Could not open file at ") + packedmatrix);
    uint32_t nk;
    gzread(ifp, &nk, sizeof(nk));
    std::vector<unsigned> k_values(nk);
    gzread(ifp, &number_entries, sizeof(number_entries));
    gzread(ifp, &number_sets, sizeof(number_sets));
    gzread(ifp, k_values.data(), k_values.size() * sizeof(unsigned));
    std::fprintf(stderr, "nk: %d\n", nk);
    for(const auto v: k_values) {
        std::fprintf(stderr, "k: %u\n", v);
    }
    double kmean = std::accumulate(k_values.begin(), k_values.end(), 0.) / k_values.size();
    assert(((number_sets * (number_sets - 1)) >> 1) == number_entries);
    dm::DistanceMatrix<float> ret(number_sets);
    std::fprintf(stderr, "num entries (at end): %zu\n", ret.num_entries());
    std::vector<float> value_set(nk);
    for(size_t i = 0; i < number_sets; ++i) {
        //auto rowptr = ret.row_ptr(i);
        const uint64_t genome_size1 = genome_sizes[i];
        const unsigned ns1 = nseqs ? nseqs[i]: 1u;
        // TODO: batch process and parallelize
        for(size_t j = i + 1; j < number_sets; ++j) {
            const unsigned ns2 = nseqs ? nseqs[j]: 1u;
            const uint64_t genome_size2 = genome_sizes[j];
            const double background = corrected_background_rates ? background_match(i, j, corrected_background_rates): .25;
            std::fprintf(stderr, "background rate: %g\n", background);
            int nb;
            if((nb = gzread(ifp, value_set.data(), nk * sizeof(float))) != int(nk * sizeof(float))) {
                auto msg = ks::sprintf("Failed to read from disk at %zu/%zu. Requested %d, got %d \n", i, j, int(nk * sizeof(float)), nb);
                throw std::runtime_error(msg);
            }
            auto jcp = jukes_cantor_p(k_values, value_set.data(), background, kmean, genome_size1, genome_size2, ns1, ns2);
            ret(i, j) = jcp; // TODO: turn this into distances?
        }
    }
    if(reassign_to_distance) {
        size_t i = 0;
#if 0
        const auto one = _mm_set1_ps(1), negft = _mm_set1_ps(-4. / 3);
        OMP_PFOR
        for(; i < ret.num_entries() / 4; ++i) {
            float *ptr = ret.data() + 4 * i;
            _mm_storeu_ps(ptr, _mm_add_ps(_mm_mul_ps(negft, _mm_sub_ps(one, _mm_loadu_ps(ptr))), one));
            for(float *pp = ptr + 4; ptr < pp; *ptr = std::log(*ptr), ++ptr);
        }
        i *= 4;
#endif
        for(;i < ret.num_entries(); ++i) {
            ret.data()[i] = std::log(1 - (4. / 3) * (1. - ret.data()[i]));
        }
    }
    gzclose(ifp);
    return ret;
}

} // namespace bns
