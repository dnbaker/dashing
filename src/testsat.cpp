#include "hll/hll.h"
#include "bonsai/include/aesctr.h"
#include "bonsai/kspp/ks.h"
#include "bonsai/include/util.h"
#include <random>
#include <mutex>
#include "cppitertools/product.hpp"
#include "klib/kthread.h"
#include <getopt.h>

void usage() {
    std::fprintf(stderr, "Usage:\ntestsat <opts>\n"
                         "-r\tAdd a random number size\n"
                         "-s\tAdd a sketch size\n"
                         "-p\tSet number of threads [1]\n"
                         "-n\tSet number of iterations. [Default: 250]\n"
                         "Purpose: test the saturation of HyperLogLogs, including how often they are within expected bounds at various cardinalities and sketch sizes.\n");
    std::exit(EXIT_FAILURE);
}

template<typename MutexType>
class LockSmith {
    std::lock_guard<MutexType> lock_;
public:
    LockSmith(MutexType &m): lock_(m) {}
};

struct acc {
    size_t             exact_;
    size_t     nwithinbounds_;
    size_t            nabove_;
    size_t            nbelow_;
    std::vector<double> ests_;
    acc(size_t val): exact_(val), nwithinbounds_(0), nabove_(0), nbelow_(0) {}
    acc(const acc &other) = default;
    acc(acc &&other)      = default;
    double mean_bias() const {
        return std::accumulate(std::begin(ests_), std::end(ests_), 0., [&](auto a, auto b) {return a + b - exact_;}) / ests_.size();
    }
    double mean_error() const {
        return std::accumulate(std::begin(ests_), std::end(ests_), 0., [&](auto a, auto b) {return a + std::abs(b - exact_);}) / ests_.size();
    }
    double mse() const {
        return std::accumulate(std::begin(ests_), std::end(ests_), 0., [&](auto a, auto b) {return a + (b - exact_) * (b - exact_);}) / ests_.size();
    }
    double sde() const {
        const double mb = mean_error();
        return std::sqrt(std::accumulate(std::begin(ests_), std::end(ests_), 0., [&](auto a, auto b) {return a + ((b - exact_) - mb) * ((b - exact_) - mb);}) / (ests_.size() - 1));
    }
    double sdb() const {
        const double mb = mean_bias();
        return std::sqrt(std::accumulate(std::begin(ests_), std::end(ests_), 0., [&](auto a, auto b) {return a + ((b - exact_) - mb) * ((b - exact_) - mb);}) / (ests_.size() - 1));
    }
    std::pair<double, double> conf95() {
        std::sort(std::begin(ests_), std::end(ests_));
        return std::make_pair(ests_[ests_.size() * 0.025], ests_[ests_.size() * 0.975]);
    };
    void add(double est, unsigned ss) {
        ests_.push_back(est);
        const auto experr = (1.03896 / std::sqrt(1ull << ss) * exact_);
        const auto err = static_cast<double>(exact_) - est;
        nwithinbounds_ += std::abs(err) < experr;
        nabove_ += (err < 0);
        nbelow_ += (err > 0);
    }
    void clear() {ests_.clear(); nwithinbounds_ = nabove_ = nbelow_ = 0}
};

struct kth_t {
    const std::vector<size_t>                &rns_;
    const std::vector<unsigned>               &ss_;
    std::vector<std::vector<std::uint64_t>> &bufs_;
    std::vector<std::vector<hll::hll_t>>    &hlls_;
    std::vector<ks::string>              &strings_;
    std::mutex                                 &m_;
    const size_t                            niter_;
    std::FILE                                 *fp_;
};

void func(void *data_, long index, int tid) {
    aes::AesCtr<uint64_t, 8> gen(index); // Seed with job index
    kth_t &data(*(kth_t *)data_);
    auto &buf = data.bufs_[tid];
    const auto &ss = data.ss_;
    auto &hlls(data.hlls_[tid]);
    const size_t rn = data.rns_[index];
    std::vector<acc> accumulators;
    accumulators.reserve(hlls.size());
    std::generate_n(std::back_emplacer(accumulators), hlls.size(), [rn](){return acc(rn);});
    // Consider just using the vanilla method, which involves more labor-intensive copying but is simpler.
    if(rn > buf.size()) buf.resize(rn);
    for(size_t inum(0); inum < data.niter_; ++inum) {
        size_t i(0);
        while(i < rn * 8 / gen.BUFSIZE) {
            gen.generate_new_values();
            std::memcpy(buf.data() + (i++ * gen.BUFSIZE / 8), gen.buf(), gen.BUFSIZE);
        }
        for(i *= gen.BUFSIZE / 8; i<rn; buf[i++] = gen());
        for(size_t i(0); i < rn; ++i) for(auto &h: hlls) h.addh(buf[i]);
        for(size_t i(0); i < hlls.size(); ++i) accumulators[i].add(hlls[i].report(), ss[i]);
        //for(auto &h: hlls) std::fprintf(stderr, "Estimated %lf with exact %zu\n", h.report(), rn);
        for(auto &h: hlls) h.clear();
    }
    auto &ks(data.strings_[tid]);
    for(size_t i(0); i < hlls.size(); ++i) {
        auto c95 = accumulators[i].conf95();
        ks.sprintf("%u\t%zu\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf|%lf\n",
                   ss[i], rn, accumulators[i].mean_bias(), accumulators[i].mean_error(),
                   accumulators[i].mse(), (1.03896 / std::sqrt(1ull << ss[i]) * rn), static_cast<double>(accumulators[i].nwithinbounds_) / data.niter_,
                   static_cast<double>(accumulators[i].nabove_) / data.niter_, static_cast<double>(accumulators[i].nbelow_) / data.niter_,
                   accumulators[i].sde(), accumulators[i].sdb(), c95.first, c95.second);
    }
    for(auto &a: accumulators) a.clear();
    {
        LockSmith lock(data.m_);
        ks.write(fileno(data.fp_));
    }
    ks.clear();
}

std::vector<size_t> DEFAULT_RNUMS {
    1ull << 10, 1ull << 12, 1ull << 14, 1ull << 16, 1ull << 18, 1ull << 20, 1ull << 22, 1ull << 24, 1ull << 26, 1ull << 28, 1ull << 30, 1ull << 31, 1ull << 32, 1ull << 33
};

std::vector<unsigned> DEFAULT_SIZES {
    10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20
};

int main(int argc, char *argv[]) {
    if(argc == 1) usage();
    size_t default_buf_size = 1 << 20;
    size_t niter = 250;
    std::vector<size_t> rnum_sizes;
    std::vector<unsigned> sketch_sizes;
    int c, nthreads(1);
    std::FILE *fp = stdout;
    while((c = getopt(argc, argv, "o:n:r:s:p:b:dh?")) >= 0) {
        switch(c) {
            case 'n': niter = std::strtoull(optarg, nullptr, 10); break;
            case 'r':
                rnum_sizes.emplace_back(static_cast<size_t>(std::strtoull(optarg, nullptr, 10))); break;
            case 's':
                sketch_sizes.emplace_back(std::strtoul(optarg, nullptr, 10)); break;
            case 'p': nthreads = std::atoi(optarg); break;
            case 'b': default_buf_size = static_cast<size_t>(std::strtoull(optarg, nullptr, 10)); break;
            case 'd': rnum_sizes = DEFAULT_RNUMS; sketch_sizes = DEFAULT_SIZES; break;
            case 'o': fp = std::fopen(optarg, "w"); break;
            case 'h': case '?': usage();
        }
    }
    if(rnum_sizes.empty() || sketch_sizes.empty()) {
        std::fprintf(stderr, "Error: at least one each of sketch sizes and rnum sizes should be provided.\n");
        usage();
    }
    std::vector<ks::string> kstrings;
    kstrings.reserve(nthreads);
    std::generate_n(std::back_emplacer(kstrings), nthreads, []{return ks::string(1024);});
    std::vector<std::vector<uint64_t>> bufs;
    std::vector<std::vector<hll::hll_t>> hlls;
    std::generate_n(std::back_emplacer(hlls), nthreads, [&] {
        std::vector<hll::hll_t> ret;
        std::generate_n(std::back_emplacer(ret), sketch_sizes.size(), [&]{return hll::hll_t(sketch_sizes[ret.size()]);});
        return ret;
    });
    std::fill_n(std::back_emplacer(bufs), nthreads, std::vector<uint64_t>(default_buf_size));
    std::mutex m;
    kth_t data{rnum_sizes, sketch_sizes, bufs, hlls, kstrings, m, niter, fp};
    std::fprintf(fp, "#Sketch size (log2)\tExact size\tMean bias\tMean error\tMean squared error\tTheoretical Mean Error\tFraction within bounds\tFraction Overestimated\tFraction Underestimated\tError Std Deviation\tBias Std Deviation\t95%% confidence interval\n");
    std::fflush(fp);
    kt_for(nthreads, &func, (void *)&data, rnum_sizes.size());
    if(fp != stdout) std::fclose(fp);
    return EXIT_SUCCESS;
}
