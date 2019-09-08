#include "sketch_and_cmp.h"
namespace bns {
template<> void sketch_core<hll::hll_t>(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp, const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix, std::vector<CountingSketch> &counting_sketches, EstimationMethod estim, JointEstimationMethod jestim, KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing, bool skip_cached, bool canon, uint32_t mincount, bool entropy_minimization, EncodingType enct);
}
