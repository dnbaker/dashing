#include "dashing.h"
#include "sketch_and_cmp.h"

namespace bns {
#define DISTEXT(sketchtype) \
     extern template void  \
                     dist_sketch_and_cmp<sketchtype>(std::vector<std::string> &inpaths, std::vector<CountingSketch> &cms, KSeqBufferHolder &kseqs, std::FILE *ofp, std::FILE *&pairofp, std::string outpath,\
                         Spacer sp,\
                         unsigned ssarg, unsigned mincount, EstimationMethod estim, JointEstimationMethod jestim, bool cache_sketch, EmissionType result_type, EmissionFormat emit_fmt,\
                         bool presketched_only, unsigned nthreads, bool use_scientific, std::string suffix, std::string prefix, bool canon, bool entropy_minimization, std::string spacing,\
                         size_t nq, EncodingType enct);\
     extern template void  \
                     dist_sketch_and_cmp<sketch::wj::WeightedSketcher<sketchtype>>(std::vector<std::string> &inpaths, std::vector<CountingSketch> &cms, KSeqBufferHolder &kseqs, std::FILE *ofp, std::FILE *&pairofp, std::string outpath,\
                         Spacer sp,\
                         unsigned ssarg, unsigned mincount, EstimationMethod estim, JointEstimationMethod jestim, bool cache_sketch, EmissionType result_type, EmissionFormat emit_fmt,\
                         bool presketched_only, unsigned nthreads, bool use_scientific, std::string suffix, std::string prefix, bool canon, bool entropy_minimization, std::string spacing,\
                         size_t nq, EncodingType enct);

DISTEXT(HyperLogLogHasher<>)
DISTEXT(bf::bf_t)
//DISTEXT(mh::RangeMinHash<uint64_t>)
DISTEXT(BKHash64)
DISTEXT(mh::CountingRangeMinHash<uint64_t>)
DISTEXT(khset64_t)
//DISTEXT(SuperMinHashType)
DISTEXT(mh::BBitMinHasher<uint64_t>)
#undef DISTEXT
int dist_main(int argc, char *argv[]) {
    int wsz(0), k(31), sketch_size(10), use_scientific(false), co, cache_sketch(false),
        nthreads(1), mincount(5), nhashes(1), cmsketchsize(-1);
    int canon(true), presketched_only(false), entropy_minimization(false),
         avoid_fsorting(false), weighted_jaccard(false);
    Sketch sketch_type = HLL;
    EmissionFormat emit_fmt = UT_TSV;
    EncodingType enct = BONSAI;
    EmissionType result_type(JI);
    hll::EstimationMethod estim = hll::EstimationMethod::ERTL_MLE;
    hll::JointEstimationMethod jestim = static_cast<hll::JointEstimationMethod>(hll::EstimationMethod::ERTL_MLE);
    std::string spacing, paths_file, suffix, prefix, pairofp_labels, pairofp_path;
    FILE *ofp(stdout), *pairofp(stdout);
    sketching_method sm = EXACT;
    std::vector<std::string> querypaths;
    uint64_t seedseedseed = 1337u;
    if(argc == 1) dist_usage(bns::executable);
    int option_index;
    DIST_LONG_OPTS
    while((co = getopt_long(argc, argv, "n:Q:P:x:F:c:p:o:s:w:O:S:k:=:t:R:D:8TgazlICbMEeHJhZBNyUmqW?", dist_long_options, &option_index)) >= 0) {
        switch(co) {
            case 1337: case 'D': /* does nothing */ break;
            case '8': sketch_type = BB_MINHASH; break;
            case 'B': gargs.bbnbits = std::atoi(optarg);   break;
            case 'F': paths_file = optarg;              break;
            case 'P': prefix     = optarg;              break;
            case 'U': emit_fmt   =  UPPER_TRIANGULAR;   break;
            case 'l': result_type = FULL_MASH_DIST;     break;
            case 'T': emit_fmt = FULL_TSV;              break;
            case 'Q': querypaths = get_paths(optarg); break;
            case 'R': seedseedseed = std::strtoull(optarg, nullptr, 10); break;
            case 'E': jestim   = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ORIGINAL); break;
            case 'I': jestim   = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ERTL_IMPROVED); break;
            case 'J': jestim   = hll::JointEstimationMethod::ERTL_JOINT_MLE; break;
            case 'm': jestim   = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ERTL_MLE); LOG_WARNING("Note: ERTL_MLE is default. This flag is superfluous.\n"); break;
            case 'S': sketch_size = std::atoi(optarg);  break;
            case 'e': use_scientific = true; break;
            case 'C': canon = false; break;
            case 'b': emit_fmt = BINARY;                break;
            case 'c': mincount = std::atoi(optarg);     break;
            case 'g': entropy_minimization = true; LOG_WARNING("Entropy-based minimization is probably theoretically ill-founded, but it might be of practical value.\n"); break;
            case 'k': k        = std::atoi(optarg);           break;
            case 'M': result_type = MASH_DIST; break;
            case 'o': if((ofp = fopen(optarg, "w")) == nullptr) LOG_EXIT("Could not open file at %s for writing.\n", optarg); break;
            case 'p': nthreads = std::atoi(optarg);     break;
            case 'q': nhashes  = std::atoi(optarg);     break;
            case 't': cmsketchsize = std::atoi(optarg); break;
            case 's': spacing  = optarg;                break;
            case 'w': wsz      = std::atoi(optarg);         break;
            case 'x': suffix   = optarg;                 break;
            case 'O': if((pairofp = fopen(optarg, "wb")) == nullptr)
                          LOG_EXIT("Could not open file at %s for writing.\n", optarg);
                      pairofp_labels = std::string(optarg) + ".labels";
                      pairofp_path = optarg;
                      break;
            case 140:
                gargs.weighted_jaccard_cmsize  = std::atoi(optarg); weighted_jaccard = true; break;
            case 141:
                gargs.weighted_jaccard_nhashes = std::atoi(optarg); weighted_jaccard = true; break;
            case 148:
                gargs.nperbatch = std::max(std::strtoull(optarg, nullptr, 10), static_cast<unsigned long long>(1)); break;
            case 143:
                gargs.number_neighbors = std::atoi(optarg);
                BNS_REQUIRE(gargs.number_neighbors > 0);
                std::fprintf(stderr, "Calculating %u nearest neighbors\n", gargs.number_neighbors);
                break;
            case 145:
                gargs.exact_weighted = true; break;
            case 'W': cache_sketch = true; break;
            // Should also be set by getopt, but users are reporting that it does not.
            case 'h': case '?': dist_usage(bns::executable);
        }
    }
    if(k > 32 && enct == BONSAI)
        UNRECOVERABLE_ERROR("k must be <= 32 for non-rolling hashes.");
    if(k > 32 && spacing.size())
        UNRECOVERABLE_ERROR("kmers must be unspaced for k > 32");
    if(nthreads < 0) nthreads = 1;
    if(gargs.number_neighbors > 0) {
        gargs.show();
        emit_fmt = EmissionFormat(unsigned(emit_fmt) | NEAREST_NEIGHBOR_TABLE);
    }
#if !NDEBUG
    std::fprintf(stderr, "paths file: %s/%zu\n", paths_file.data(), paths_file.size());
#endif
    std::vector<std::string> inpaths(paths_file.size() ? get_paths(paths_file.data())
                                                       : std::vector<std::string>(argv + optind, argv + argc));
    if(inpaths.empty())
        std::fprintf(stderr, "No paths. See usage.\n"), dist_usage(bns::executable);
    omp_set_num_threads(nthreads);
    Spacer sp(k, wsz, parse_spacing(spacing.data(), k));
    size_t nq = querypaths.size();
    if(nq == 0 && !is_symmetric(result_type)) {
        querypaths = inpaths;
        nq = querypaths.size();
        LOG_WARNING("Note: No query files provided, but an asymmetric distance was requested. Switching to a query/reference format with all references as queries.\n"
                    "In the future, this will throw an error.\nYou must provide query and reference paths (-Q/-F) to calculate asymmetric distances.\n");
    }
    if(!presketched_only && !avoid_fsorting) {
        detail::sort_paths_by_fsize(inpaths);
        detail::sort_paths_by_fsize(querypaths);
    }
    inpaths.reserve(inpaths.size() + querypaths.size());
    for(auto &p: querypaths)
        inpaths.push_back(std::move(p));
    {
        decltype(querypaths) tmp;
        std::swap(tmp, querypaths);
    }
    std::vector<CountingSketch> cms;
    KSeqBufferHolder kseqs(nthreads);
    switch(sm) {
        case CBF: case BY_FNAME: {
            if(cmsketchsize < 0) {
                cmsketchsize = 20;
                LOG_WARNING("CM Sketch size not set. Defaulting to 20, 1048576 entries per table\n");
            }
            cms.reserve(nthreads);
            while(cms.size() < static_cast<unsigned>(nthreads))
                cms.emplace_back(cmsketchsize, nhashes, 1.08, (cms.size() ^ seedseedseed) * 1337uL);
            break;
        }
        case EXACT: default: break;
    }
    if(enct == NTHASH)
        std::fprintf(stderr, "Use nthash's rolling hash for kmers. This comes at the expense of reversibility\n");
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

    switch(sketch_type) {
        case BB_MINHASH:      CALL_DIST_BOTH(mh::BBitMinHasher<uint64_t>); break;
        case HLL:      if(gargs.defer_hll_creation) CALL_DIST_BOTH(HyperLogLogHasher<>);
                       else                         CALL_DIST_BOTH(hll::hll_t);
        break;
        case WIDE_HLL:        CALL_DIST_BOTH(sketch::WideHyperLogLogHasher<>); break;
        case RANGE_MINHASH:   CALL_DIST_BOTH(BKHash64); break;
        case BLOOM_FILTER:    CALL_DIST_BOTH(bf::bf_t); break;
        case FULL_KHASH_SET:  CALL_DIST_BOTH(khset64_t); break;
        default: {
                char buf[128];
                std::sprintf(buf, "Sketch %s not yet supported.\n", (size_t(sketch_type) >= (sizeof(sketch_names) / sizeof(char *)) ? "Not such sketch": sketch_names[sketch_type]));
                UNRECOVERABLE_ERROR(buf);
        }
    }

    std::future<void> label_future;
    if(emit_fmt == BINARY) {
        if(pairofp_labels.empty()) pairofp_labels = "unspecified";
        label_future = std::async(std::launch::async, [&inpaths](const std::string &labels) {
            std::FILE *fp = std::fopen(labels.data(), "wb");
            if(fp == nullptr) UNRECOVERABLE_ERROR(std::string("Could not open file at '") + labels + "' for writing");
            for(const auto &path: inpaths) std::fwrite(path.data(), path.size(), 1, fp), std::fputc('\n', fp);
            std::fclose(fp);
        }, pairofp_labels);
    }
    if(pairofp  && pairofp != stdout) std::fclose(pairofp);
    if(label_future.valid()) label_future.get();
    return EXIT_SUCCESS;
} // dist_main

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
    std::fprintf(stderr, "Writing to %s\n", outpath.data());
    std::FILE *ofp = std::fopen(outpath.data(), "wb");
    std::fprintf(stderr, "Prepared everything. Now calling dist_by_seq\n");
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
