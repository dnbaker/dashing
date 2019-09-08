#include "sketch_and_cmp.h"
namespace bns {
using CBM = mh::CountingBBitMinHasher<uint64_t, uint16_t>;
DECSKETCHCORE(CBM)
}
