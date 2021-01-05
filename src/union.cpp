#include "dashing.h"
namespace bns {
using sketch::hll_t;

void show(const std::vector<std::string> &p) {
    for(const auto &x: p) std::fprintf(stderr, "%s ", x.data());
    std::fputc('\n', stderr);
}

template<typename T>
void union_core(std::vector<std::string> &paths, gzFile ofp) {
    T ret(paths.back().data());
    paths.pop_back();
    for(const auto &path: paths) {
        T tmp(path.data());
        ret += tmp;
    }
    ret.write(ofp);
}

template<>
void union_core<sketch::hll_t>(std::vector<std::string> &paths, gzFile ofp) {
    sketch::hll_t fs(paths[0]);
    fs.sum();
    char *p = new char[sizeof(sketch::hll_t) * (paths.size() - 1)];
    hll_t *oh = reinterpret_cast<hll_t *>(p);
    OMP_PFOR
    for(size_t i = 1; i < paths.size(); ++i) {
        new(oh + i - 1) hll_t(paths[i].data());
    }
    for(size_t i = 1; i < paths.size(); ++i) {
        fs += oh[i - 1];
    }
    for(size_t i = 1; i < paths.size(); ++i) oh[i].~hll_t();
    fs.write(ofp);
}

int union_main(int argc, char *argv[]) {
    if(std::find_if(argv, argc + argv,
                    [](const char *s) {return std::strcmp(s, "--help") == 0 || std::strcmp(s, "-h") == 0;})
       != argc + argv)
        union_usage(*argv);
    bool compress = false;
    int compression_level = 6;
    const char *opath = "/dev/stdout";
    std::vector<std::string> paths;
    Sketch sketch_type = HLL;
    for(int c;(c = getopt(argc, argv, "b:o:F:zZ:h?")) >= 0;) {
        switch(c) {
            case 'h': union_usage(*argv);
            case 'Z': compression_level = std::atoi(optarg); [[fallthrough]];
            case 'z': compress = true; break;
            case 'o': opath = optarg; break;
            case 'F': paths = get_paths(optarg); break;
            case 'r': sketch_type = RANGE_MINHASH; break;
            case 'H': sketch_type = FULL_KHASH_SET; break;
            case 'b': sketch_type = BLOOM_FILTER; break;
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
        case HLL: union_core<hll::hll_t>(paths, ofp); break;
        case WIDE_HLL: union_core<wh119_t>(paths, ofp); break;
        case BLOOM_FILTER: union_core<bf::bf_t>(paths, ofp); break;
        case FULL_KHASH_SET: union_core<khset64_t>(paths, ofp); break;
        case RANGE_MINHASH: union_core<mh::FinalRMinHash<uint64_t>>(paths, ofp); break;
        default: throw NotImplementedError(ks::sprintf("Union not implemented for %s\n", sketch_names[sketch_type]).data());
    }
    gzclose(ofp);
    return 0;
}

}
