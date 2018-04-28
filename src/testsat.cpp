#include "hll/hll.h"
#include "bonsai/include/aesctr.h"
#include "bonsai/kspp/ks.h"
#include "bonsai/include/util.h"
#include <random>
#include <mutex>
#include <limits>
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
    double max() const {
        return ests_.size() ? *std::max_element(ests_.begin(), ests_.end(), [&](auto a, auto b){return std::abs(a - exact_) < std::abs(b - exact_);}): std::numeric_limits<double>::quiet_NaN();
    }
    std::pair<double, double> conf(unsigned pct, bool sort=true) {
        if(sort) std::sort(std::begin(ests_), std::end(ests_));
        const double start = (1. - (pct * 0.01)) / 2., end = 1. - start;
        return std::make_pair(ests_[ests_.size() * start], ests_[ests_.size() * end]);
    }
#define CONF_N(pct) std::pair<double, double> conf##pct(bool sort=true) {return conf(pct, sort);}
    CONF_N(90)
    CONF_N(95)
    CONF_N(80)
    CONF_N(10)
    CONF_N(25)
    CONF_N(50)
#undef CONF_N
    void add(double est, unsigned ss) {
        ests_.push_back(est);
        const auto experr = (1.03896 / std::sqrt(1ull << ss) * exact_);
        const auto err = static_cast<double>(exact_) - est;
        nwithinbounds_ += std::abs(err) < experr;
        nabove_ += (err < 0);
        nbelow_ += (err > 0);
    }
    void clear() {ests_.clear(); nwithinbounds_ = nabove_ = nbelow_ = 0;}
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
    const bool         emit_full_percentile_range_;
};

struct kth_jacc_t: public kth_t {
    template<typename... Args>
    kth_jacc_t(const std::vector<double> &jaccs, Args &&... args):
        kth_t(std::forward<Args>(args)...), ohlls_(this->hlls_), obufs_(this->bufs_), jaccs_(jaccs) {}
    std::vector<std::vector<hll::hll_t>>     ohlls_; // Actually copies
    std::vector<std::vector<std::uint64_t>>  obufs_; // Actually copies
    const std::vector<double>            &jaccs_;
};

template<typename T>
void random_buff(T *data, size_t nelem, aes::AesCtr<T, 8> &gen) {
    using AesType = typename aes::AesCtr<T, 8>;
    using ResultType = typename AesType::result_type;
    size_t i(0);
    while(i < nelem * sizeof(ResultType) / gen.BUFSIZE) {
        gen.generate_new_values();
        std::memcpy(data + (i++ * gen.BUFSIZE / sizeof(ResultType)), gen.buf(), gen.BUFSIZE);
    }
    for(i *= gen.BUFSIZE / sizeof(ResultType);i < nelem; data[i++] = gen());
}

void jacc_func(void *data_, long index, int tid) {
    aes::AesCtr<uint64_t, 8> gen(index); // Seed with job index
    kth_jacc_t &data(*(kth_jacc_t *)data_);
    auto &buf = data.bufs_[tid];
    auto &obuf = data.obufs_[tid];
    const auto &ss = data.ss_;
    auto &hlls(data.hlls_[tid]);
    auto &ohlls(data.ohlls_[tid]);
    const size_t rn = data.rns_[index % data.rns_.size()];
    const double ji = data.jaccs_[index / data.rns_.size()];
    if(rn > buf.size()) buf.resize(rn);
    uint64_t hv;
    hll::WangHash hf;
    const uint64_t us = rn;
    const uint64_t is = rn * ji;
    const uint64_t unique_size = (us - is) / 2;
    const double exact_ji = static_cast<double>(is) / static_cast<double>(us);
    std::vector<std::uint64_t> shared_buf(is); // The only heap allocation in the core of the program. These could be cached.
    std::vector<acc> accumulators;
    accumulators.reserve(hlls.size());
    std::generate_n(std::back_emplacer(accumulators), hlls.size(), [rn](){return acc(rn);});
    // Core loop
    for(size_t inum(0); inum < data.niter_; ++inum) {
        random_buff(shared_buf.data(), shared_buf.size(), gen);
        for(const auto &v: shared_buf) {
            hv = hf(v);
            for(auto &h: hlls) h.add(hv);
            for(auto &oh: ohlls) oh.add(hv);
        }
        random_buff(buf.data(), unique_size, gen);
        for(size_t i(0); i < unique_size;) {
            hv = hf(buf[i++]);
            for(auto &h: hlls) h.add(hv);
        }
        random_buff(obuf.data(), unique_size, gen);
        for(size_t i(0); i < unique_size;) {
            hv = hf(obuf[i++]);
            for(auto &oh: ohlls) oh.add(hv);
        }
        for(size_t i(0); i < hlls.size(); ++i) {
            throw std::runtime_error("NotImplementedError: TODO @dnbaker -- write statistical accumulation code.");
        }
        for(auto &h: hlls) h.clear();
        for(auto &oh: hlls) oh.clear();
    }
    auto &ks(data.strings_[tid]);
#if 0
    for(size_t i(0); i < hlls.size(); ++i) {
        auto c95 = accumulators[i].conf95();
        auto c90 = accumulators[i].conf90(false);
        auto c80 = accumulators[i].conf80(false);
        auto c50 = accumulators[i].conf50(false);
        auto c10 = accumulators[i].conf10(false);
        ks.sprintf("%u\t%zu\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf|%lf\t",
                   ss[i], rn, accumulators[i].mean_bias(), accumulators[i].mean_error(),
                   accumulators[i].mse(), (1.03896 / std::sqrt(1ull << ss[i]) * rn), static_cast<double>(accumulators[i].nwithinbounds_) / data.niter_,
                   static_cast<double>(accumulators[i].nabove_) / data.niter_, static_cast<double>(accumulators[i].nbelow_) / data.niter_,
                   accumulators[i].sde(), accumulators[i].sdb(), c95.first, c95.second);
        if(data.emit_full_percentile_range_)
            ks.sprintf("%lf|%lf\t%lf|%lf\t%lf|%lf\t%lf|%lf\t", c90.first, c90.second, c80.first, c80.second, c50.first, c50.second, c10.first, c10.second);
        ks.sprintf("%lf\n", accumulators[i].max());
    }
    for(auto &a: accumulators) a.clear();
#endif
    {
        LockSmith lock(data.m_);
        ks.write(fileno(data.fp_));
    }
    ks.clear();
}
void card_func(void *data_, long index, int tid) {
    aes::AesCtr<uint64_t, 8> gen(index); // Seed with job index
    kth_t &data(*(kth_t *)data_);
    auto &buf = data.bufs_[tid];
    const auto &ss = data.ss_;
    auto &hlls(data.hlls_[tid]);
    const size_t rn = data.rns_[index];
    std::vector<acc> accumulators;
    accumulators.reserve(hlls.size());
    std::generate_n(std::back_emplacer(accumulators), hlls.size(), [rn](){return acc(rn);});
    if(rn > buf.size()) buf.resize(rn);
    uint64_t hv;
    hll::WangHash hf;
    for(size_t inum(0); inum < data.niter_; ++inum) {
        random_buff(buf.data(), rn, gen);
        for(size_t i(0); i < rn;) {
            hv = hf(buf[i++]);
            for(auto &h: hlls) h.add(hv);
        }
        for(size_t i(0); i < hlls.size(); ++i) accumulators[i].add(hlls[i].report(), ss[i]);
        //for(auto &h: hlls) std::fprintf(stderr, "Estimated %lf with exact %zu\n", h.report(), rn);
        for(auto &h: hlls) h.clear();
    }
    auto &ks(data.strings_[tid]);
    for(size_t i(0); i < hlls.size(); ++i) {
        auto c95 = accumulators[i].conf95();
        auto c90 = accumulators[i].conf90(false);
        auto c80 = accumulators[i].conf80(false);
        auto c50 = accumulators[i].conf50(false);
        auto c10 = accumulators[i].conf10(false);
        ks.sprintf("%u\t%zu\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf|%lf\t",
                   ss[i], rn, accumulators[i].mean_bias(), accumulators[i].mean_error(),
                   accumulators[i].mse(), (1.03896 / std::sqrt(1ull << ss[i]) * rn), static_cast<double>(accumulators[i].nwithinbounds_) / data.niter_,
                   static_cast<double>(accumulators[i].nabove_) / data.niter_, static_cast<double>(accumulators[i].nbelow_) / data.niter_,
                   accumulators[i].sde(), accumulators[i].sdb(), c95.first, c95.second);
        if(data.emit_full_percentile_range_)
            ks.sprintf("%lf|%lf\t%lf|%lf\t%lf|%lf\t%lf|%lf\t", c90.first, c90.second, c80.first, c80.second, c50.first, c50.second, c10.first, c10.second);
        ks.sprintf("%lf\n", accumulators[i].max());
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
    10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22
};

std::vector<double> DEFAULT_JACCS {
    0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0
};

int main(int argc, char *argv[]) {
    if(argc == 1) usage();
    size_t default_buf_size = 1 << 20;
    size_t niter = 50;
    std::vector<size_t> rnum_sizes;
    std::vector<unsigned> sketch_sizes;
    bool emit_full_percentile_range(false), jaccard_exploration(false);
    int c, nthreads(1);
    std::vector<double> jaccs;
    std::FILE *fp = stdout;
    while((c = getopt(argc, argv, "o:n:r:s:p:b:J:jPdh?")) >= 0) {
        switch(c) {
            case 'n': niter = std::strtoull(optarg, nullptr, 10); break;
            case 'r':
                rnum_sizes.emplace_back(static_cast<size_t>(std::strtoull(optarg, nullptr, 10))); break;
            case 's':
                sketch_sizes.emplace_back(std::strtoul(optarg, nullptr, 10)); break;
            case 'p': nthreads = std::atoi(optarg); break;
            case 'J': jaccs.emplace_back(std::atof(optarg)); break;
            case 'j': jaccard_exploration = true; break;
            case 'b': default_buf_size = static_cast<size_t>(std::strtoull(optarg, nullptr, 10)); break;
            case 'd': rnum_sizes = DEFAULT_RNUMS; sketch_sizes = DEFAULT_SIZES; break;
            case 'o': fp = std::fopen(optarg, "w"); break;
            case 'P': emit_full_percentile_range = true; break;
            case 'h': case '?': usage();
        }
    }
    if(rnum_sizes.empty() || sketch_sizes.empty()) {
        std::fprintf(stderr, "Error: at least one each of sketch sizes and rnum sizes should be provided.\n");
        usage();
    }
    {
        aes::AesCtr<uint64_t, 8> gen(1337);
        std::vector<uint64_t> rvals(10000);
        random_buff(rvals.data(), rvals.size(), gen);
        if(std::find(rvals.begin(), rvals.end(), 0) != rvals.end()) throw std::runtime_error("ZOMG");
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
    kth_t data{rnum_sizes, sketch_sizes, bufs, hlls, kstrings, m, niter, fp, emit_full_percentile_range};
    if(!jaccard_exploration) {
        std::fprintf(fp, "#Sketch size (log2)\tExact size\tMean bias\tMean error\tMean squared error\tTheoretical Mean Error\tFraction within bounds\tFraction Overestimated\tFraction Underestimated\tError Std Deviation\tBias Std Deviation\t95%% Interval");
        if(emit_full_percentile_range)
            std::fprintf(fp, "90%% Interval\t80%% Interval\t50%% Interval\t10%% Interval\t");
        std::fprintf(fp, "Max Error\n");
        std::fflush(fp);
        kt_for(nthreads, &card_func, (void *)&data, rnum_sizes.size());
    } else {
        if(jaccs.empty()) jaccs = DEFAULT_JACCS;
        kth_jacc_t jacc_data(jaccs, data);
        std::fprintf(fp, "#Sketch size (log2)\tExact Union Size\tExact Size1\tExact Size2\tMean Est Union Size\tMean Est Intersection Size\tMean Est Size1\tMean Est Size2\n"
                         "Exact JI\tMean Est JI\tMean Error US\tMean Error IS\tMean Error JI\tMean Bias JI\tMean Bias IS\t90%% range\n");
        std::fflush(fp);
        kt_for(nthreads, &jacc_func, (void *)&jacc_data, rnum_sizes.size() * jaccs.size());
        throw std::runtime_error("NotImplementedError.");
    }
    if(fp != stdout) std::fclose(fp);
    return EXIT_SUCCESS;
}
