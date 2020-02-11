#include "include/bonsai/encoder.h"
#include "src/substrs.h"
namespace bns {

static inline auto background_match(unsigned lhid, unsigned rhid, const float *nucfreqs) {
    auto lhp = nucfreqs + (lhid * 4), rhp = nucfreqs + (rhid * 4);
    return lhp[0] * rhp[0] + lhp[1] * rhp[1] +
           lhp[2] * rhp[2] + lhp[3] * rhp[3];
}

std::pair<std::vector<float>, std::vector<uint64_t>> nuc_freqs(const std::vector<std::string> &fpaths, unsigned nt, unsigned k) {
    std::vector<uint64_t> genome_sizes(fpaths.size());
    std::vector<float> ret(4 * fpaths.size());
    const size_t np = fpaths.size();
    KseqBufferHolder kseqs(nt);
    OMP_PFOR
    for(size_t i = 0; i < np; ++i) {
        uint64_t arr [128]{0};
        size_t nkmer_windows = 0;
        for_each_substr([&](const char *s) {
            const int tid = OMP_ELSE(omp_get_thread_num(), 0);
            gzFile fp = gzopen(s, "rb");
            if(!fp) {std::fprintf(stderr, "Missing file %s\n", s); std::exit(1);}
            kseq_t *ks = &kseqs[i];
            kseq_assign(ks, fp);
            nkmer_windows += ks->seq.l;
            while(kseq_read(ks) >= 0) {
                for(size_t j = 0; j < ks->seq.l; ++arr[ks->seq.s[j++]]);
            }
            gzclose(fp);
            genome_sizes[i] = nkmer_windows;
        }, fpaths[i], FNAME_SEP);
        float cparr[4];
        double total_valid = arr['A'] + arr['a'] + arr['t'] + arr['T'] +  arr['c'] + arr['C'] + arr['g'] + arr['G'];
        cparr[0] = arr['A'] + arr['a'] / total_valid;
        cparr[1] = arr['c'] + arr['C'] / total_valid;
        cparr[2] = arr['g'] + arr['G'] / total_valid;
        cparr[3] = arr['t'] + arr['T'] / total_valid;
        std::memcpy(ret.data() + (4 * i), cparr, sizeof(cparr));
    }
    return std::make_pair(nuc_freqs, genome_sizes);
}

double jukes_cantor_p(const std::vector<unsigned> &ks, const std::vector<float> &values,
                         double background, bool is_intersection_size, uint64_t lhsz, uint64_t rhsz) {
    double kmean = std::accumulate(ks.begin(), ks.end(), 0.) / ks.size();
    std::unique_ptr<double> lvals(new double[ks.size()]);
    for(unsigned i = 0; i < ks.size(); ++i) {
        lvals[i] = std::log(values[i] - std::pow(background, k) * 4 * lhsz * rhsz);
    }
    double lmean = std::accumulate(lvals.get(), &lvals[ks.szez()], 0.) / ks.size();
    double num = 0., denom = 0.;
    for(unsigned i = 0; i < ks.size(); ++i) {
        auto kdiff = (ks[i] - kmean);
        denom += kdiff * kdiff;
        num += kdiff * (lvals[i] - lmean);
    }
    const double slope = num / denom;
    double p = std::exp(slope);
    return p;
}

INLINE double jcp2dist(double p) {
    return -0.75 * std::log1p(-4. / 3. * (1 - p));
}

dm::DistanceMatrix<float> mkmat2jcdistmat(std::string packedmatrix, EmissionType emtype, const float *corrected_background_rates=nullptr,
                                          const uint64_t *genome_sizes=nullptr)
{
    switch(emtype) {
        case MASH_DIST:
        case FULL_MASH_DIST:
        case FULL_CONTAINMENT_DIST:
        case CONTAINMENT_INDEX:
        case CONTAINMENT_DIST:
            throw std::invalid_argument("EmissionType must be JI or SIZES");
    }
    const bool is_intersection_size = (emtype == SIZES);
    uint64_t number_sets, number_entries;
    gzFile ifp = gzopen(packedmatrix.data(), "rb");
    if(!ifp) RUNTIME_ERROR(std::string("Could not open file at ") + packedmatrix);
    uint32_t nk;
    gzread(ifp, &nk, sizeof(nk));
    std::vector<unsigned> k_values(nk);
    gzread(ifp, &number_entries, sizeof(number_entries));
    gzread(ifp, &number_sets, sizeof(number_sets));
    gzread(ifp, k_values.data(), k_values.size() * sizeof(unsigned));
    assert(((number_sets * (number_sets - 1)) >> 1) == number_entries);
    dm::DistanceMatrix<float> ret(number_sets);
    std::vector<float> value_set(nk);
    for(size_t i = 0; i < number_sets; ++i) {
        auto rowptr = ret.row_ptr(i);
        const uint64_t genome_size1 = genome_size1[i];
        for(size_t j = i + 1; j < number_sets; ++j) {
            const uint64_t genome_size2 = genome_sizes[j];
            double background = corrected_background_rates ? background_match(i, j): .25;
            gzread(ifp, value_set.data(), nk * sizeof(float));
            float jcp = jukes_cantor_p(k_values, value_set, background, is_intersection_size, genome_size1, genome_size2);
            ret(i, j) = jcp2dist(jcp); // TODO: leave this as p?
        }
    }
    gzclose(ifp);
}

} // namespace bns
