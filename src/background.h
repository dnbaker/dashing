#ifndef DASHING_BACKGROUND_H__
#define DASHING_BACKGROUND_H__
#include <vector>
#include <string>
#include "distmat/distmat.h"
#include "enums.h"


namespace bns {
struct sumstats: public std::tuple<std::vector<float>, std::vector<uint64_t>, std::vector<unsigned>> {
    auto &sizes() {return std::get<1>(*this);}
    auto &freqs() {return std::get<0>(*this);}
    auto &numseqs() {return std::get<2>(*this);}
    sumstats(size_t n) {
        resize(n);
    }
    void resize(size_t n) {
        sizes().resize(n);
        freqs().resize(n * 4);
        numseqs().resize(n);
    }
};
sumstats nuc_freqs(const std::vector<std::string> &fpaths);

dm::DistanceMatrix<float> mkmat2jcdistmat(std::string packedmatrix, EmissionType emtype,
                                          const uint64_t *genome_sizes,
                                          const float *corrected_background_rates=nullptr,
                                          const unsigned *nseqs=nullptr,
                                          bool reassign_to_distance=false);

static INLINE double jcp2dist(double p) {
    return -0.75 * std::log1p(-4. / 3. * (1 - p));
}

} // namespace bns
#endif /* DASHING_BACKGROUND_H__ */
