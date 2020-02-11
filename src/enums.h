#ifndef EMISSION_FMT_H__
#define EMISSION_FMT_H__

#ifdef FNAME_SEP
#pragma message("Not: FNAME_SEP already defined. [not default \"' '\"]")
#else
#define FNAME_SEP ' '
#endif

namespace bns {


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

enum sketching_method: int {
    EXACT = 0,
    CBF   = 1,
    BY_FNAME = 2
};

enum EncodingType {
    BONSAI,
    NTHASH,
    RK,
    CYCLIC
};

enum NNType {
    // Nearest Neighbor Type
    DECREASING = 0,
    INCREASING = 1,
    DIST_MEASURE = DECREASING,
    SIMILARITY_MEASURE = INCREASING
};

}
#endif
