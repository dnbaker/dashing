#include "hll/hll.h"
#include "bonsai/include/aesctr.h"
#include <random>
#include "cppitertools/product.hpp"
#include "klib/kthread.h"
#include <getopt.h>

void usage() {
    std::fprintf(stderr, "Usage:\ntestsat <opts>\n"
                         "-r\tAdd a random number size\n"
                         "-s\tAdd a sketch size\n"
                         "-p\tSet number of threads [1]\n"
                         "-n\tSet number of iterations. [Default: 50]\n"
                         "Purpose: test the saturation of HyperLogLogs, including how often they are within expected bounds at various cardinalities and sketch sizes.\n");
    std::exit(EXIT_FAILURE);
}

struct kth_t {
    const std::vector<size_t>               &rns_;
    const std::vector<unsigned>             &ss_;
    std::vector<std::vector<std::uint64_t>> &bufs_;
};

void func(void *data_, long index, int tid) {
    aes::AesCtr<uint64_t, 8> gen(index);
    kth_t & data(*(kth_t *)data_);
    auto &buf = data.bufs_[tid];
    auto &ss = data.ss_;
    const size_t rn = data.rns_[index];
    buf.clear();
}

int main(int argc, char *argv[]) {
    if(argc == 1) usage();
    size_t default_buf_size = 1 << 20;
    std::vector<size_t> rnum_sizes;
    std::vector<unsigned> sketch_sizes;
    int c, nthreads(1);
    while((c = getopt(argc, argv, "r:s:p:b:h?")) >= 0) {
        switch(c) {
            case 'r':
                rnum_sizes.emplace_back(static_cast<size_t>(std::strtoull(optarg, nullptr, 10))); break;
            case 's':
                sketch_sizes.emplace_back(std::strtoul(optarg, nullptr, 10)); break;
            case 'p': nthreads = std::atoi(optarg); break;
            case 'b': default_buf_size = static_cast<size_t>(std::strtoull(optarg, nullptr, 10)); break;
            case 'h': case '?': usage();
        }
    }
    if(rnum_sizes.empty() || sketch_sizes.empty()) {
        std::fprintf(stderr, "Error: at least one each of sketch sizes and rnum sizes should be provided.\n");
        usage();
    }
    std::vector<std::vector<uint64_t>> bufs;
    std::fill_n(std::back_inserter(bufs), nthreads, std::vector<uint64_t>(default_buf_size)); 
    kth_t data{rnum_sizes, sketch_sizes, bufs};
    kt_for(nthreads, &func, (void *)&data, rnum_sizes.size());
    return EXIT_SUCCESS;
}
