#include "dashing.h"
#include "sketch_and_cmp.h"

#define CALL_SIZE_MAIN(sketchtype) \
 size_sketch_and_emit<sketchtype>(inpaths, cms, kseqs, ofp, sp, sketch_size, mincount, enct, estim, jestim, cache_sketch, emit_binary, use_scientific, \
                                  presketched_only, nthreads, prefix, suffix, canon, spacing)

namespace bns {
int card_main(int argc, char *argv[]) {
    int wsz(0), k(31), sketch_size(10), use_scientific(false), co, cache_sketch(false),
        nthreads(1), mincount(5), nhashes(1), cmsketchsize(-1);
    int canon(true), presketched_only(false), entropy_minimization(false),
         avoid_fsorting(false), weighted_jaccard(false);
    Sketch sketch_type = HLL;
    EmissionType result_type = JI;
    EmissionFormat emit_fmt = (EmissionFormat)0;
    hll::EstimationMethod estim = hll::EstimationMethod::ERTL_MLE;
    hll::JointEstimationMethod jestim = static_cast<hll::JointEstimationMethod>(hll::EstimationMethod::ERTL_MLE);
    std::string spacing, paths_file, suffix, prefix, pairofp_labels;
    FILE *ofp(stdout);
    sketching_method sm = EXACT;
    EncodingType enct = BONSAI;
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
            case 'l': result_type = FULL_MASH_DIST;     break;
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
            case 140:
                gargs.weighted_jaccard_cmsize  = std::atoi(optarg); weighted_jaccard = true; break;
            case 141:
                gargs.weighted_jaccard_nhashes = std::atoi(optarg); weighted_jaccard = true; break;
            case 143:
                gargs.number_neighbors = std::atoi(optarg);
                BNS_REQUIRE(gargs.number_neighbors > 0);
                std::fprintf(stderr, "gargs.number_neighbors: %u\n", gargs.number_neighbors);
                break;
            case 145:
                gargs.exact_weighted = true; break;
            case 'W': cache_sketch = true; break;
            // Should also be set by getopt, but users are reporting that it does not.
            case 'h': case '?': dist_usage(bns::executable);
        }
    }
    const bool emit_binary = emit_fmt & BINARY;
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

    switch(sketch_type) {
        //case BB_MINHASH:      CALL_SIZE_MAIN(mh::BBitMinHasher<uint64_t>); break;
        //case BB_SUPERMINHASH: CALL_SIZE_MAIN(SuperMinHashType); break;
        case HLL:      if(gargs.defer_hll_creation) CALL_SIZE_MAIN(HyperLogLogHasher<>);
                       else                         CALL_SIZE_MAIN(hll::hll_t);
        break;
        case RANGE_MINHASH:   CALL_SIZE_MAIN(mh::RangeMinHash<uint64_t>); break;
        case BLOOM_FILTER:    CALL_SIZE_MAIN(bf::bf_t); break;
        case FULL_KHASH_SET:  CALL_SIZE_MAIN(khset64_t); break;
        default: {
                char buf[128];
                std::sprintf(buf, "Sketch %s not yet supported.\n", (size_t(sketch_type) >= (sizeof(sketch_names) / sizeof(char *)) ? "Not such sketch": sketch_names[sketch_type]));
                UNRECOVERABLE_ERROR(buf);
        }
    }

    std::future<void> label_future;
    if((int)emit_fmt & BINARY) {
        if(pairofp_labels.empty()) pairofp_labels = "unspecified";
        label_future = std::async(std::launch::async, [&inpaths](const std::string &labels) {
            std::FILE *fp = std::fopen(labels.data(), "wb");
            if(fp == nullptr) UNRECOVERABLE_ERROR(std::string("Could not open file at '") + labels + "' for writing");
            for(const auto &path: inpaths) std::fwrite(path.data(), path.size(), 1, fp), std::fputc('\n', fp);
            std::fclose(fp);
        }, pairofp_labels);
    }
    if(label_future.valid()) label_future.get();
    return EXIT_SUCCESS;
} // size_main

} // namespace bns
