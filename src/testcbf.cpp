#include <cassert>
#include "bonsai/hll/cbf.h"
#include "bonsai/kspp/ks.h"
#include "bonsai/hll/aesctr.h"
#include "omp.h"
#include "bonsai/tinythreadpp/source/fast_mutex.h"
#include <set>

void usage() {
    std::fprintf(stderr, "TODO: Write this usage menu\n");
    std::exit(EXIT_FAILURE);
}

using namespace bf;

static const std::vector<unsigned> DEFAULT_BFS {10, 12, 14, 16, 18, 20, 22};
static const std::vector<unsigned> DEFAULT_NHASHES {1, 2, 4, 8, 16};
static const std::vector<unsigned> DEFAULT_NBFS {8, 16, 32};
static const std::vector<uint64_t> DEFAULT_SIZES {
    10ull,
    1024ull,
    2048ull,
    4096ull,
    8192ull,
    16384ull,
    32768ull,
    65536ull,
    131072ull,
    262144ull,
    524288ull,
    1048576ull,
    2097152ull,
    4194304ull,
    8388608ull,
    16777216ull,
    33554432ull,
    67108864ull,
    134217728ull // Powers of 2 from 10 to 28
};

struct LockSmith {
    // Simple lock-holder to avoid writing to the same file twice.
    tthread::fast_mutex &m_;
    LockSmith(tthread::fast_mutex &m): m_(m) {m_.lock();}
    ~LockSmith() {m_.unlock();}
};

template<typename T, typename=std::enable_if_t<std::is_integral_v<T>>>
unsigned fastl2(T val) {
    return sizeof(val) * CHAR_BIT - __builtin_clz(val) - 1;
}

int main(int argc, char *argv[]) {
    if(argc == 1) usage();
    // Approximate counting bloom filter
    uint64_t seedseedseed = 1337; // Wheeeee!
    std::vector<uint64_t> sizes;
    std::vector<unsigned> bfsizes;
    std::vector<unsigned> nhashes;
    std::vector<unsigned> nbfs;
    int c, nthreads = 1;
    std::FILE *ofp = stdout;
    while((c = getopt(argc, argv, "f:b:z:n:o:p:s:S:H:dh")) >= 0) {
        switch(c) {
            case 'S': seedseedseed = std::strtoull(optarg, nullptr, 10); break;
            case 's': sizes.emplace_back(std::strtoull(optarg, nullptr, 10)); break;
            case 'd': sizes = DEFAULT_SIZES, bfsizes = DEFAULT_BFS, nhashes = DEFAULT_NHASHES, nbfs = DEFAULT_NBFS; break;
            case 'p': nthreads = std::atoi(optarg); break;
            case 'n': nbfs.push_back(std::atoi(optarg)); break;
            case 'H': nhashes.push_back(std::atoi(optarg)); break;
            case 'z': bfsizes.push_back(std::atoi(optarg)); break;
            case 'o':
                if(ofp != stdout) throw std::runtime_error("Can't assign ofp twice! This would memory leak.");
                if((ofp = std::fopen(optarg, "w")) == nullptr)
                    throw std::runtime_error("Could not open file at "s + optarg);
                break;
            case 'h': usage(); break;
        }
    }
    tthread::fast_mutex lock;
    nthreads = std::max(nthreads, 1);
    omp_set_num_threads(nthreads);
    const int fn = fileno(ofp);
    if(sizes.size() * nbfs.size() * bfsizes.size() * nhashes.size() == 0)
        sizes = DEFAULT_SIZES, bfsizes = DEFAULT_BFS, nhashes = DEFAULT_NHASHES, nbfs = DEFAULT_NBFS;
    size_t total_number = sizes.size() * nbfs.size() * bfsizes.size() * nhashes.size();
    std::vector<std::tuple<unsigned, unsigned, unsigned, uint64_t>> combs;
    combs.reserve(total_number);
    for(const auto nh: nhashes)
        for(const auto bfs: bfsizes)
            for(const auto nbf: nbfs)
                for(const auto size: sizes)
                    combs.emplace_back(nh, bfs, nbf, size);
    std::fprintf(ofp, "Size\tNumber of hashes\tBloom Filter size (log2)\tNumber of bloom filters\t[...]: Occupancy rates, then popcounts of each filter.\n");
    std::fflush(ofp);
    #pragma omp parallel for
    for(size_t i = 0; i < total_number; ++i) {
        auto [nh, bfsize, nbfs, size] = combs[i];
        const uint64_t seed = WangHash()(i);
        aes::AesCtr<uint64_t, 8> gen(seed);
        std::vector<uint64_t> buf(size);
        while(buf.size() < size) buf.emplace_back(gen());
        bf::cbf_t cbf(nbfs, bfsize, nh, seedseedseed + 666 * 777 * (i + 1)); // Convolution of bad and good luck.
        for(const auto val: buf) cbf.addh(val);
        gen.seed(seed);
        std::vector<uint64_t> counts(nbfs + 1);
        for(const auto val: buf) {
            unsigned count = cbf.est_count(val);
            if(count == 0) {
                std::fputs(ks::sprintf("Value which should be at least 1 is missing. (%" PRIu64 ")", val).data(), stderr);
                ++counts[0];
                continue;
            }
            if(fastl2(count) > counts.size()) throw std::runtime_error(ks::sprintf("Count of %u found, which is too big for cbf of %u elements", count, nbfs).data());
            ++counts[fastl2(cbf.est_count(val))];
        }
        std::fflush(stderr);
        assert(std::accumulate(counts.begin(), counts.end(), 0) == size);
        ks::string outstr = ks::sprintf("%" PRIu64 "\t%u\t%u\t%u", size, nh, bfsize, nbfs);
        for(size_t i(0); i < counts.size(); ++i) outstr.sprintf("\t%u|%" PRIu64 "", 1<<i, counts[i]);
        outstr.sprintf("\t<popcount/m>\t");
        for(const auto &bf: cbf) {
            outstr.sprintf("%u/%zu\t", bf.popcnt(), bf.m());
        }
        outstr.putc_('\n');
        {
            LockSmith he_who_holds_the_keys(lock);
            outstr.write(fn);
        }
        buf.clear();
        //for(auto &bf: cbf) std::fprintf(stderr, "vals: %s\n", bf.print_vals().data());
        //for(auto &bf: cbf) std::fprintf(stderr, "seeds: %s\n", bf.seedstring().data());
    }
    return EXIT_SUCCESS;
}
