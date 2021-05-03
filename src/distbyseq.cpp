#include "dashing.h"
#include "sketch_and_cmp.h"

namespace bns {
#define CALL_DIST(sketchtype) \
    dist_sketch_and_cmp<sketchtype>(inpaths, cms, kseqs, ofp, pairofp, pairofp_path, sp, sketch_size,\
                                    mincount, estim, jestim, cache_sketch, result_type,\
                                    emit_fmt, presketched_only, nthreads,\
                                    use_scientific, suffix, prefix, canon,\
                                    entropy_minimization, spacing, nq, enct)

#define CALL_DIST_WEIGHTED(sketchtype) CALL_DIST(sketch::wj::WeightedSketcher<sketchtype>)
#define CALL_DIST_WEIGHTED_EXACT(sketchtype) do {\
        using local_sketch = sketch::wj::WeightedSketcher<sketchtype, sketch::wj::ExactCountingAdapter>;\
        CALL_DIST(local_sketch); \
    } while(0)

#define CALL_DIST_BOTH(sketchtype) do {\
    if(gargs.exact_weighted || weighted_jaccard) {\
        if(gargs.exact_weighted) {\
            CALL_DIST_WEIGHTED_EXACT(sketchtype);\
        } else {\
            CALL_DIST_WEIGHTED(sketchtype);\
        }\
    } else CALL_DIST(sketchtype);\
} while(0)

void dist_by_seq_usage(const char *s=bns::executable) {
    std::fprintf(stderr, "Usage: %s <flags> -n [namefile] input_file\n"
                         "-p\t threads [1]\n"
                         "-o\t output path [/dev/stdout]\n"
                         "-b\t emit binary output\n"
                         "-U\t emit PHYLIP Upper Triangular output\n"
                         "-T\t emit full TSV format\n"
                         "data structures\n"
                         "-8\tb-bit minhash\n-B\tBloom Filter\n-C\tCounting Range MinHash\n-r\tRange MinHash\n\n"
                         "HLL options\n"
                         "Estimation methods - default MLE\n"
                         "-J\tJoint MLE\n-E\tOriginal Flajolet\n-I\tErtl Improved\n"
                         "Result options\n\n"
                         "Default: Jaccard Index\n"
                         "--containment-index\t Emit containment index\n"
                         "--containment-dist\t Emit containment distnace\n"
                         "--mash-dist, -M\t Emit mash distance"
                         "--symmetric-containment-index\tEmit symmetric containment index\n"
                         "--symmetric-containment-dist\tEmit symmetric containment distance\n"
                         "--sizes\tEmit intersection sizes\n"
                      , s);
    std::exit(EXIT_FAILURE);
}


int dist_by_seq_main(int argc, char *argv[]) {
    int c;
    std::string outpath = "/dev/stdout";
    hll::EstimationMethod estim = hll::EstimationMethod::ERTL_MLE;
    hll::JointEstimationMethod jestim = static_cast<hll::JointEstimationMethod>(hll::EstimationMethod::ERTL_MLE);
    std::string namefile, otherpath;
    EmissionFormat emit_fmt = UT_TSV;
    EmissionType result_type(JI);
    Sketch sketch_type = HLL;
    int k = -1;
    int nthreads = 1;
    static option_struct dbs_options [] {
        LO_FLAG("containment-index", 131, result_type, CONTAINMENT_INDEX) \
        LO_FLAG("containment-dist", 132, result_type, CONTAINMENT_DIST) \
        LO_FLAG("mash-dist", 'M', result_type, MASH_DIST)\
        LO_FLAG("symmetric-containment-index", 137, result_type, SYMMETRIC_CONTAINMENT_INDEX) \
        LO_FLAG("symmetric-containment-dist", 138, result_type, SYMMETRIC_CONTAINMENT_DIST) \
        LO_FLAG("sizes", 'Z', result_type, SIZES)\
        {0,0,0,0}
    };
    int option_index;
    while((c = getopt_long(argc, argv, "q:o:k:n:p:EIJMBbS8KCTrh?", dbs_options, &option_index)) >= 0) {
        switch(c) {
            case 'B': sketch_type = BLOOM_FILTER; break;
            case 'S':
            case '8': sketch_type = BB_MINHASH; break;
            case 'K': sketch_type = FULL_KHASH_SET; break;
            //case 'C': sketch_type = COUNTING_RANGE_MINHASH; break;
            case 'r': sketch_type = RANGE_MINHASH; break;
            case 'p': nthreads = std::atoi(optarg); break;
            case 'o': outpath = optarg; break;
            case 'E': jestim   = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ORIGINAL); break;
            case 'I': jestim   = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ERTL_IMPROVED); break;
            case 'J': jestim   = hll::JointEstimationMethod::ERTL_JOINT_MLE; break;
            case 'm': jestim   = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ERTL_MLE); LOG_WARNING("Note: ERTL_MLE is default. This flag is superfluous.\n"); break;
            case 'k': k = std::atoi(optarg); break;
            case 'n': namefile = optarg; break;
            case 'q': otherpath = optarg; break;
            case 'b': emit_fmt = BINARY; break;
            case 'T': emit_fmt = FULL_TSV; break;
            case 'U': emit_fmt = UPPER_TRIANGULAR; break;
            case 'h': dist_by_seq_usage(bns::executable);
        }
    }
    if(optind + 1 != argc || namefile.empty())
        dist_by_seq_usage(bns::executable);
    auto labels = get_paths(namefile.data());
    if(k <= 0) {
        std::ifstream reader(namefile);
        std::string line;
        std::getline(reader, line);
        auto tmp = std::atoi(line.data() + 3);
        if(tmp <= 0) tmp = 31; // Just guess
        k = tmp;
    }
    std::FILE *ofp = std::fopen(outpath.data(), "wb");
#define DBS(sketch)  \
    dist_by_seq<sketch>(labels, argv[optind], ofp, outpath, \
                        k, estim, jestim, result_type, emit_fmt, nthreads, otherpath); \
    break

    switch(sketch_type) {
        case HLL: DBS(hll::hll_t);
        case WIDE_HLL: DBS(whll::wh119_t);
        case FULL_KHASH_SET:
             DBS(khset64_t);
        case BLOOM_FILTER: DBS(bf_t);
        case BB_MINHASH: case BB_SUPERMINHASH:
            DBS(FinalBBitMinHash);
        case RANGE_MINHASH:
            DBS(RMFinal);
        default:
            UNRECOVERABLE_ERROR("Unexpected value " + std::to_string(int(sketch_type)));
    }
#undef DBS
    std::fclose(ofp);
    return EXIT_SUCCESS;
}

#undef DIST_LONG_OPTS

} // namespace bns
