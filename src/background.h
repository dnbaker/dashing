#ifndef DASHING_BACKGROUND_H__
#define DASHING_BACKGROUND_H__
#include <vector>
#include <string>
#include "distmat/distmat.h"
#include "enums.h"


namespace bns {
std::vector<float> nuc_freqs(const std::vector<std::string> &fpaths, unsigned nt);

static inline auto background_match(unsigned lhid, unsigned rhid, const std::vector<float> &nucfreqs) {
    auto lhp = nucfreqs.data() + (lhid * 4), rhp = nucfreqs.data() + (rhid * 4);
    return lhp[0] * rhp[0] + lhp[1] * rhp[1] +
           lhp[2] * rhp[2] + lhp[3] * rhp[3];
}

#if 0
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
#endif

dm::DistanceMatrix<float> mkmat2jcdistmat(std::string packedmatrix, EmissionType emtype);


} // namespace bns
#endif /* DASHING_BACKGROUND_H__ */
