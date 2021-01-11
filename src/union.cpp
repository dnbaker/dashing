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
    size_t np = paths.size();
    char *p = new char[sizeof(T) * np] + (sizeof(T) - 1);
    char *ap = p;
    while(reinterpret_cast<uint64_t>(ap) % sizeof(T)) ++ap;
    T *oh = reinterpret_cast<T *>(ap);
    new(oh) T(paths.front().data());
    paths.pop_back();
    if(paths.empty()) return;

    OMP_PFOR
    for(size_t i = 0; i < std::min(nthreads, np); ++i) {
        new(oh) T(paths[i].data());
    }
    if(np > nthreads) {
        OMP_PFOR
        for(size_t i = nthreads; i < np; ++i) {
            oh[omp_get_thread_num()] += T(paths[i].data());
        }
    }
    if(nthreads > 1) {
        const size_t nloops = std::ceil(std::log2(nthreads));
        for(size_t i = 0; i < nloops; ++i) {
            size_t step = 1ull << i, block_size = step << 1;
            size_t nblocks = (np + block_size - 1) / block_size;
            OMP_PFOR
            for(size_t j = 0; j < nblocks; ++j) {
                auto base = step * j, oidx = base + step;
                if(oidx < nthreads) oh[base] += oh[oidx];
            }
        }
    }
    oh[0].write(ofp);
    OMP_PFOR
    for(size_t i = 0; i < nthreads; ++i) {
        oh[i].~T();
    }
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
