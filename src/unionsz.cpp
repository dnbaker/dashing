#include "dashing.h"

namespace bns {namespace us {
template<> double union_size<khset64_t> (const khset64_t &a, const khset64_t &b) {return a.union_size(b);}
template<> double union_size<hll::hllbase_t<>> (const hll::hllbase_t<> &a, const hll::hllbase_t<> &b) {return a.union_size(b);}
template<> double union_size<RMFinal> (const RMFinal &a, const RMFinal &b) {return a.union_size(b);}
template<> double union_size<CRMFinal> (const CRMFinal &a, const CRMFinal &b) {return a.union_size(b);}
template<> double union_size<mh::FinalBBitMinHash> (const mh::FinalBBitMinHash &a, const mh::FinalBBitMinHash &b) {
    return (a.est_cardinality_ + b.est_cardinality_ ) / (1. + a.jaccard_index(b));
}
}} // bns::us
