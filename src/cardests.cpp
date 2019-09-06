#include "dashing.h"

namespace bns {
template<> double cardinality_estimate(hll::hll_t &x) {return x.report();}
template<> double cardinality_estimate(mh::FinalBBitMinHash &x) {return x.est_cardinality_;}
template<> double cardinality_estimate(mh::FinalDivBBitMinHash &x) {return x.est_cardinality_;}
}
