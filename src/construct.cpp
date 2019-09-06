#include "dashing.h"

namespace bns {

template<> mh::BBitMinHasher<uint64_t> construct<mh::BBitMinHasher<uint64_t>>(size_t p) {return mh::BBitMinHasher<uint64_t>(p, gargs.bbnbits);}

}
