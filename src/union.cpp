#include "dashing.h"
namespace bns {
using sketch::hll_t;

void show(const std::vector<std::string> &p) {
    for(const auto &x: p) std::fprintf(stderr, "%s ", x.data());
    std::fputc('\n', stderr);
}

template<typename T>
void reduce(T *x, size_t n) {
    for(size_t i = 1; i < n; ++i)
        *x += x[i];
}

template<typename T>
void par_reduce(T *x, size_t n) {
    const unsigned int ln = static_cast<int>(std::ceil(std::log2(n)));
    for(size_t i = 0; i < ln; ++i) {
        const size_t step_size = 1 << i;
        const size_t sweep_size = (i + 1) << 1;
        const size_t nsweeps = (n + (sweep_size - 1)) / sweep_size;
        OMP_PFOR
        for(size_t j = 0; j < nsweeps; ++j) {
            const auto lh = j * sweep_size, rh = lh + step_size;
            if(rh < n)
                x[lh] += x[rh];
        }
    }
}

template<typename T>
void union_core(std::vector<std::string> &paths, gzFile ofp, size_t nthreads) {
    if(paths.size() < 1) {
        std::fprintf(stderr, "require >= 1 paths. See usage.\n");
        std::exit(1);
    }
    T *items = nullptr;
    const size_t cap = std::min(paths.size(), nthreads);
    if(posix_memalign((void **)&items, 64, sizeof(T) * cap))
        throw std::bad_alloc();
    OMP_PFOR
    for(size_t i = 0; i < cap; ++i)
        new(&items[i]) T(paths[i].data());
    if(cap < paths.size()) {
        OMP_PFOR
        for(size_t i = cap; i < paths.size(); ++i) {
            items[omp_get_thread_num()] += T(paths[i].data());
        }
    }
    par_reduce(items, cap);
    items[0].write(ofp);
    OMP_PFOR
    for(size_t i = 0; i < cap; ++i) items[i].~T();
    std::free(items);
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
