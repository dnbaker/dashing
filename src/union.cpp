#include "dashing.h"
namespace bns {
using sketch::hll_t;

void show(const std::vector<std::string> &p) {
    for(const auto &x: p) std::fprintf(stderr, "%s ", x.data());
    std::fputc('\n', stderr);
}

template<typename T>
void union_core(std::vector<std::string> &paths, gzFile ofp, size_t nthreads) {
    // Read from disk
    if(paths.size() < 1) {
        std::fprintf(stderr, "require >= 1 paths. See usage.\n");
        std::exit(1);
    }
    T u(paths.front().data());
    for(size_t i = 1; i < paths.size(); ++i)
        u += T(paths[i]);
    u.write(ofp);
}

int union_main(int argc, char *argv[]) {
    if(std::find_if(argv, argc + argv,
                    [](const char *s) {return std::strcmp(s, "--help") == 0 || std::strcmp(s, "-h") == 0;})
       != argc + argv)
        union_usage(*argv);
    bool compress = false;
    int compression_level = 6, nthreads = 1;
    const char *opath = "/dev/stdout";
    std::vector<std::string> paths;
    Sketch sketch_type = HLL;
    for(int c;(c = getopt(argc, argv, "p:b:o:F:zZ:h?")) >= 0;) {
        switch(c) {
            case 'h': union_usage(*argv);
            case 'Z': compression_level = std::atoi(optarg); [[fallthrough]];
            case 'z': compress = true; break;
            case 'o': opath = optarg; break;
            case 'F': paths = get_paths(optarg); break;
            case 'r': sketch_type = RANGE_MINHASH; break;
            case 'H': sketch_type = FULL_KHASH_SET; break;
            case 'b': sketch_type = BLOOM_FILTER; break;
            case 'p': nthreads = std::atoi(optarg); break;
        }
    }
    if(argc == optind && paths.empty()) union_usage(*argv);
    std::for_each(argv + optind, argv + argc, [&](const char *s){paths.emplace_back(s);});
    char mode[6];
    if(compress && compression_level)
        std::sprintf(mode, "wb%d", compression_level % 23);
    else
        std::sprintf(mode, "wT");
    gzFile ofp = gzopen(opath, mode);
    if(!ofp) throw std::runtime_error(std::string("Could not open file at ") + opath);
    using sketch::whll::wh119_t;
    switch(sketch_type) {
        case HLL: union_core<hll::hll_t>(paths, ofp, nthreads); break;
        case WIDE_HLL: union_core<wh119_t>(paths, ofp, nthreads); break;
        case BLOOM_FILTER: union_core<bf::bf_t>(paths, ofp, nthreads); break;
        case FULL_KHASH_SET: union_core<khset64_t>(paths, ofp, nthreads); break;
        case RANGE_MINHASH: union_core<mh::FinalRMinHash<uint64_t>>(paths, ofp, nthreads); break;
        default: throw NotImplementedError(ks::sprintf("Union not implemented for %s\n", sketch_names[sketch_type]).data());
    }
    gzclose(ofp);
    return 0;
}

}
