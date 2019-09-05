#include "dashing.h"
#include "sketch_and_cmp.h"

namespace bns {
#define DIST_LONG_OPTS \
static option_struct dist_long_options[] = {\
    LO_FLAG("full-tsv", 'T', emit_fmt, FULL_TSV)\
    LO_FLAG("emit-binary", 'b', emit_fmt, BINARY)\
    LO_FLAG("phylip", 'U', emit_fmt, UPPER_TRIANGULAR)\
    LO_FLAG("no-canon", 'C', canon, false)\
    LO_FLAG("by-entropy", 'g', entropy_minimization, true) \
    LO_FLAG("use-bb-minhash", '8', sketch_type, BB_MINHASH)\
    LO_FLAG("full-mash-dist", 'l', result_type, FULL_MASH_DIST)\
    LO_FLAG("mash-dist", 'M', result_type, MASH_DIST)\
    LO_FLAG("countmin", 'y', sm, CBF)\
    LO_FLAG("sketch-by-fname", 'N', sm, BY_FNAME)\
    LO_FLAG("sizes", 'Z', result_type, SIZES)\
    LO_FLAG("use-scientific", 'e', use_scientific, true)\
    LO_FLAG("cache-sketches", 'W', cache_sketch, true)\
    LO_FLAG("presketched", 'H', presketched_only, true)\
    LO_FLAG("avoid-sorting", 'n', avoid_fsorting, true)\
    LO_ARG("out-sizes", 'o') \
    LO_ARG("query-paths", 'Q') \
    LO_ARG("out-dists", 'O') \
    LO_ARG("bbits", 'B')\
    LO_ARG("original", 'E')\
    LO_ARG("improved", 'I')\
    LO_ARG("ertl-joint-mle", 'J')\
    LO_ARG("ertl-mle", 'm')\
    LO_ARG("paths", 'F')\
    LO_ARG("prefix", 'P')\
    LO_ARG("nhashes", 'q')\
    LO_ARG("seed", 'R')\
    LO_ARG("sketch-size", 'S')\
    LO_ARG("bbits", 'B')\
    LO_ARG("original", 'E')\
    LO_ARG("improved", 'I')\
    LO_ARG("ertl-joint-mle", 'J')\
    LO_ARG("ertl-mle", 'm')\
    LO_ARG("paths", 'F')\
    LO_ARG("prefix", 'P')\
    LO_ARG("nhashes", 'q')\
    LO_ARG("seed", 'R')\
    LO_ARG("sketch-size", 'S')\
    LO_ARG("kmer-length", 'k')\
    LO_ARG("min-count", 'c')\
    LO_ARG("nthreads", 'p')\
    LO_ARG("cm-sketch-size", 't')\
    LO_ARG("spacing", 's')\
    LO_ARG("window-size", 'w')\
    LO_ARG("suffix", 'x')\
    LO_ARG("help", 'h')\
    LO_FLAG("use-range-minhash", 128, sketch_type, RANGE_MINHASH)\
    LO_FLAG("use-counting-range-minhash", 129, sketch_type, COUNTING_RANGE_MINHASH)\
    LO_FLAG("use-full-khash-sets", 130, sketch_type, FULL_KHASH_SET)\
    LO_FLAG("containment-index", 131, result_type, CONTAINMENT_INDEX) \
    LO_FLAG("containment-dist", 132, result_type, CONTAINMENT_DIST) \
    LO_FLAG("full-containment-dist", 133, result_type, FULL_CONTAINMENT_DIST) \
    LO_FLAG("use-bloom-filter", 134, sketch_type, BLOOM_FILTER)\
    LO_FLAG("use-super-minhash", 135, sketch_type, BB_SUPERMINHASH)\
    LO_FLAG("use-nthash", 136, enct, NTHASH)\
    LO_FLAG("symmetric-containment-index", 137, result_type, SYMMETRIC_CONTAINMENT_INDEX) \
    LO_FLAG("symmetric-containment-dist", 138, result_type, SYMMETRIC_CONTAINMENT_DIST) \
    LO_FLAG("use-cyclic-hash", 139, enct, NTHASH)\
    LO_ARG("wj-cm-sketch-size", 140)\
    LO_ARG("wj-cm-nhashes", 141)\
    LO_FLAG("wj", 142, weighted_jaccard, true)\
    {0,0,0,0}\
};

int dist_main(int argc, char *argv[]) {
    int wsz(0), k(31), sketch_size(10), use_scientific(false), co, cache_sketch(false),
        nthreads(1), mincount(5), nhashes(4), cmsketchsize(-1);
    int canon(true), presketched_only(false), entropy_minimization(false),
         avoid_fsorting(false), weighted_jaccard(false);
    Sketch sketch_type = HLL;
         // bool sketch_query_by_seq(true);
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
    if(argc == 1) dist_usage(*argv);
    int option_index;
    DIST_LONG_OPTS
    while((co = getopt_long(argc, argv, "n:Q:P:x:F:c:p:o:s:w:O:S:k:=:t:R:8TgDazlICbMEeHJhZBNyUmqW?", dist_long_options, &option_index)) >= 0) {
        switch(co) {
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
            case 'm': jestim   = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ERTL_MLE); LOG_WARNING("Note: ERTL_MLE is default. This flag is redundant.\n"); break;
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
            case 'W': cache_sketch = true; break;
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
            case 'h': case '?': dist_usage(*argv);
        }
    }
    if(k > 32 && enct == BONSAI)
        RUNTIME_ERROR("k must be <= 32 for non-rolling hashes.");
    if(k > 32 && spacing.size())
        RUNTIME_ERROR("kmers must be unspaced for k > 32");
    if(nthreads < 0) nthreads = 1;
    std::vector<std::string> inpaths(paths_file.size() ? get_paths(paths_file.data())
                                                       : std::vector<std::string>(argv + optind, argv + argc));
    if(inpaths.empty())
        std::fprintf(stderr, "No paths. See usage.\n"), dist_usage(*argv);
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
    std::vector<sketch::cm::ccm_t> cms;
    KSeqBufferHolder kseqs(nthreads);
    switch(sm) {
        case CBF: case BY_FNAME: {
            if(cmsketchsize < 0) {
                cmsketchsize = 20;
                LOG_WARNING("CM Sketch size not set. Defaulting to 20, 1048576 entries per table\n");
            }
            unsigned nbits = std::log2(mincount) + 1;
            cms.reserve(nthreads);
            while(cms.size() < static_cast<unsigned>(nthreads))
                cms.emplace_back(nbits, cmsketchsize, nhashes, (cms.size() ^ seedseedseed) * 1337uL);
            break;
        }
        case EXACT: default: break;
    }
    if(enct == NTHASH)
        std::fprintf(stderr, "Use nthash's rolling hash for kmers. This comes at the expense of reversibility\n");
#define CALL_DIST(sketchtype) \
    dist_sketch_and_cmp<sketchtype>(inpaths, cms, kseqs, ofp, pairofp, sp, sketch_size,\
                                    mincount, estim, jestim, cache_sketch, result_type,\
                                    emit_fmt, presketched_only, nthreads,\
                                    use_scientific, suffix, prefix, canon,\
                                    entropy_minimization, spacing, nq, enct)

#define CALL_DIST_WEIGHTED(sketchtype) CALL_DIST(sketch::wj::WeightedSketcher<sketchtype>)

#define CALL_DIST_BOTH(sketchtype) do {if(weighted_jaccard) CALL_DIST_WEIGHTED(sketchtype); else CALL_DIST(sketchtype);} while(0)

    switch(sketch_type) {
        case BB_MINHASH:      CALL_DIST_BOTH(mh::BBitMinHasher<uint64_t>); break;
        case BB_SUPERMINHASH: CALL_DIST_BOTH(SuperMinHashType); break;
        case HLL:             CALL_DIST_BOTH(hll::hll_t); break;
        case RANGE_MINHASH:   CALL_DIST_BOTH(mh::RangeMinHash<uint64_t>); break;
        case BLOOM_FILTER:    CALL_DIST_BOTH(bf::bf_t); break;
        case FULL_KHASH_SET:  CALL_DIST_BOTH(khset64_t); break;
        case COUNTING_RANGE_MINHASH:
                              CALL_DIST_BOTH(mh::CountingRangeMinHash<uint64_t>); break;
        default: {
                char buf[128];
                std::sprintf(buf, "Sketch %s not yet supported.\n", (size_t(sketch_type) >= (sizeof(sketch_names) / sizeof(char *)) ? "Not such sketch": sketch_names[sketch_type]));
                RUNTIME_ERROR(buf);
        }
    }

    std::future<void> label_future;
    if(emit_fmt == BINARY) {
        if(pairofp_labels.empty()) pairofp_labels = "unspecified";
        label_future = std::async(std::launch::async, [&inpaths](const std::string &labels) {
            std::FILE *fp = std::fopen(labels.data(), "wb");
            if(fp == nullptr) RUNTIME_ERROR(std::string("Could not open file at ") + labels);
            for(const auto &path: inpaths) std::fwrite(path.data(), path.size(), 1, fp), std::fputc('\n', fp);
            std::fclose(fp);
        }, pairofp_labels);
    }
    if(pairofp != stdout) std::fclose(pairofp);
    if(label_future.valid()) label_future.get();
    return EXIT_SUCCESS;
} // dist_main

#undef DIST_LONG_OPTS

} // namespace bns
