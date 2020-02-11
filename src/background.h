#ifndef DASHING_BACKGROUND_H__
#define DASHING_BACKGROUND_H__
#include <vector>
#include <string>
#include "distmat/distmat.h"
#include "enums.h"


namespace bns {
std::pair<std::vector<float>, std::vector<uint64_t>> nuc_freqs(const std::vector<std::string> &fpaths, unsigned nt, unsigned k);


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
