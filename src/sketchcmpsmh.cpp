#include "sketch_and_cmp.h"
template void bns::dist_sketch_and_cmp<bns::SuperMinHashType>(const std::vector<std::string> &inpaths, std::vector<CountingSketch> &cms, KSeqBufferHolder &kseqs, std::FILE *ofp, std::FILE *pairofp,
                     Spacer sp,
                     unsigned ssarg, unsigned mincount, EstimationMethod estim, JointEstimationMethod jestim, bool cache_sketch, EmissionType result_type, EmissionFormat emit_fmt,
                     bool presketched_only, unsigned nthreads, bool use_scientific, std::string suffix, std::string prefix, bool canon, bool entropy_minimization, std::string spacing,
                     size_t nq, EncodingType enct);
template void bns::dist_sketch_and_cmp<sketch::wj::WeightedSketcher<bns::SuperMinHashType>>(const std::vector<std::string> &inpaths, std::vector<CountingSketch> &cms, KSeqBufferHolder &kseqs, std::FILE *ofp, std::FILE *pairofp,
                     Spacer sp,
                     unsigned ssarg, unsigned mincount, EstimationMethod estim, JointEstimationMethod jestim, bool cache_sketch, EmissionType result_type, EmissionFormat emit_fmt,
                     bool presketched_only, unsigned nthreads, bool use_scientific, std::string suffix, std::string prefix, bool canon, bool entropy_minimization, std::string spacing,
                     size_t nq, EncodingType enct);
