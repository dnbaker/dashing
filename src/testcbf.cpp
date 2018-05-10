#include <cassert>
#include "bonsai/hll/cbf.h"
#include "bonsai/kspp/ks.h"
#include "bonsai/hll/aesctr.h"
#include "omp.h"
#include "bonsai/tinythreadpp/source/fast_mutex.h"
#include <set>

void usage() {
    std::fprintf(stderr, "Purpose: Tests error rates in cbf as a function of filter size.\n"
                         "-S\tSet seed for all seeds for seeds.\n"
                         "-s\tAdd a size to the numbers of unique items to add.\n"
                         "-d\tUse default sets of values for all permutations.\n"
                         "-p\tSet number of threads to use.\n"
                         "-n\tAdd a number to the set of bloom filters to use for each cbf.\n"
                         "-H\tAdd a number to the set of the number of hash functions to use per filter.\n"
                         "-z\tAdd a number to the set of the sizees of bloom filters to use.\n"
                         "-o\tSet output file [stdout]\n"
                         "-h\tThis help menu\n");
    std::exit(EXIT_FAILURE);
}

using namespace bf;

static const std::vector<unsigned> DEFAULT_BFS {10, 12, 14, 16, 18, 20, 22};
static const std::vector<unsigned> DEFAULT_NHASHES {4, 8, 16};
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
    tthread::fast_mutex mutex;
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
    std::vector<ks::string> strings;
    std::vector<std::vector<uint64_t>> countvecs;
    std::vector<std::vector<uint64_t>> bufs;
    while(strings.size() < (unsigned)nthreads) strings.emplace_back(65536u);
    while(countvecs.size() < (unsigned)nthreads) countvecs.emplace_back(16);
    while(bufs.size() < (unsigned)nthreads) bufs.emplace_back(10000);
    #pragma omp parallel for
    for(size_t i = 0; i < total_number; ++i) {
        const auto [nh, bfsize, nbfs, size] = combs[i];
        const uint64_t seed = WangHash()(i);
        aes::AesCtr<uint64_t, 8> gen(seed);
        bf::cbf_t cbf(nbfs, bfsize, nh, seedseedseed + 666 * 777 * (i + 1)); // Convolution of bad and good luck.
        const auto tid = omp_get_thread_num();
        std::vector<uint64_t> &buf(bufs[tid]);
        buf.resize(size);
        for(size_t i(0); i < buf.size(); buf[i] = gen(), cbf.addh(buf[i]), ++i);
        std::vector<uint64_t> &counts(countvecs[tid]);
        counts.resize(nbfs);
        std::memset(counts.data(), 0, sizeof(counts[0]) * counts.size());
        for(const auto val: buf) {
            unsigned result = cbf.est_count(val);
            if(result == 0) {
                std::fputs(ks::sprintf("Value which should be at least 1 is missing. (%" PRIu64 ")", val).data(), stderr);
                ++counts[0];
                continue;
            }
            result = fastl2(result);
            if(__builtin_expect(result > counts.size(), 0))
                throw std::runtime_error(ks::sprintf("log2(count) %u found, which is too big for cbf of %u elements", result, nbfs).data());
            ++counts[result];
        }
        assert(std::accumulate(counts.begin(), counts.end(), 0ull) == size);
        ks::string &outstr = strings[tid];
        outstr.sprintf("%" PRIu64 "\t%u\t%u\t%u", size, nh, bfsize, nbfs);
        for(size_t i(0); i < counts.size(); ++i) outstr.sprintf("\t%u|%" PRIu64 "", 1<<i, counts[i]);
        outstr.sprintf("\t<popcount/m>\t");
        for(const auto &bf: cbf) outstr.sprintf("%u/%zu\t", bf.popcnt(), bf.m());
        outstr.putc_('\n');
        if(outstr.size() > 1u << 16) {
            LockSmith he_who_holds_the_keys(mutex);
            outstr.write(fn);
            outstr.clear();
        }
        //std::fill(std::begin(counts), std::end(counts), 0ull);
        //for(auto &bf: cbf) std::fprintf(stderr, "vals: %s\n", bf.print_vals().data());
        //for(auto &bf: cbf) std::fprintf(stderr, "seeds: %s\n", bf.seedstring().data());
    }
    for(auto &ks: strings) ks.write(fn);
    return EXIT_SUCCESS;
}
