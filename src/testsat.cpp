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
                         "Purpose: test the saturation of HyperLogLogs, including how often they are within expected bounds at various cardinalities and sketch sizes.\n"
                         "If -j is specified, instead experiments over ranges of jaccard indices.\n"
                         "-J\tAdd a jaccard index to use\n"
                         "-l\t[int] Cover the full range over [0,1] in incremenets of 1./<value>\n"
                );
    std::exit(EXIT_FAILURE);
}

template<typename MutexType>
class LockSmith {
    std::lock_guard<MutexType> lock_;
public:
    LockSmith(MutexType &m): lock_(m) {}
};

template<typename T>
double arr_mean(const std::vector<T> &vec) {
    return std::accumulate(std::cbegin(vec), std::cend(vec), 0., [](auto a, auto b){return a + b;}) / static_cast<double>(std::size(vec));
}

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
        kth_t(std::forward<Args>(args)...), ohlls_(this->hlls_), jaccs_(jaccs) {}
    std::vector<std::vector<hll::hll_t>>     ohlls_; // Actually copies
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

struct jacc_acc {
    const uint64_t us_; // union size
    const uint64_t is_; // insersection size
    const double   ji_; // jaccard index
    std::vector<double> jis_;  // Jaccard indices
    std::vector<double> unions_;
    std::vector<double> isns_; // Insersections
    std::vector<double> s1s_; // Size 1s
    std::vector<double> s2s_; // Size 2s
    jacc_acc(uint64_t us, uint64_t is, double ji, size_t niter): us_(us), is_(is), ji_(ji)
    {
        jis_.reserve(niter);
        unions_.reserve(niter);
        isns_.reserve(niter);
        s1s_.reserve(niter);
        s2s_.reserve(niter);
    }
    void add(double ji, double us, double ins, double s1, double s2) {
        jis_.push_back(ji);
        unions_.push_back(us);
        isns_.push_back(ins);
        s1s_.push_back(s1);
        s2s_.push_back(s2);
    }
};


#define GEN_ADD(sketchvec, num) \
    do {\
        isleft = num;\
        while(isleft > NPERBUF) { \
            gen.generate_new_values(); \
            for(const auto &val: gen.template view<uint64_t>()) { \
                hv = hf(val); \
                for(auto &h: sketchvec) h.add(hv); \
            } \
            isleft -= NPERBUF; \
        }\
        while(isleft--) {\
            hv = hf(gen());\
            for(auto &h: sketchvec) h.add(hv);\
        }\
    } while(0)

void jacc_func(void *data_, long index, int tid) {
    aes::AesCtr<uint64_t, 8> gen(index); // Seed with job index
    kth_jacc_t &data(*(kth_jacc_t *)data_);
    const auto &ss = data.ss_;
    auto &hlls(data.hlls_[tid]);
    auto &ohlls(data.ohlls_[tid]);
    const size_t rn = data.rns_[index % data.rns_.size()];
    const double ji = data.jaccs_[index / data.rns_.size()];
    uint64_t hv;
    hll::WangHash hf;
    const uint64_t us = rn;
    const uint64_t is = rn * ji;
    const uint64_t unique_size = (us - is) / 2;
    const double exact_ji = static_cast<double>(is) / static_cast<double>(us);
    // This is the only heap allocation in the core loop. Caching these could help threadscaling.
    std::vector<jacc_acc> accumulators;
    accumulators.reserve(hlls.size());
    std::generate_n(std::back_emplacer(accumulators), hlls.size(), [=](){return jacc_acc(us, is, exact_ji, data.niter_);});
    assert(std::accumulate(std::begin(accumulators), std::end(accumulators), true, [](auto a, auto b) {return a && b.unions_.size() == 0;}));
    // Core loop
    for(size_t inum(0); inum < data.niter_; ++inum) {
        size_t isleft = is;
        static constexpr size_t NPERBUF = gen.BUFSIZE / sizeof(uint64_t);
        while(isleft > NPERBUF) {
            gen.generate_new_values();
            for(const auto &val: gen.template view<uint64_t>()) {
                hv = hf(val);
                for(auto &h: hlls) h.add(hv);
                for(auto &oh: ohlls) oh.add(hv);
            }
            isleft -= NPERBUF;
        }
        while(isleft--) {
            hv = hf(gen());
            for(auto &h: hlls) h.add(hv);
            for(auto &oh: ohlls) oh.add(hv);
        }
        GEN_ADD(hlls, unique_size);
        GEN_ADD(ohlls, unique_size);

        for(size_t i(0); i < hlls.size(); ++i) {
            hlls[i].sum(); ohlls[i].sum();
            double sz1(hlls[i].report());
            double sz2(ohlls[i].report());
            double est_us(hll::union_size(hlls[i], ohlls[i]));
            double isn(sz1 + sz2 - est_us);
            accumulators[i].add(isn / est_us, est_us, isn, sz1, sz2);
        }
        for(auto &h: hlls) h.clear();
        for(auto &oh: hlls) oh.clear();
    }
    auto &ks(data.strings_[tid]);
    double mv, iv, jv;
    for(size_t i(0); i < hlls.size(); ++i) {
        auto &a(accumulators[i]);
        ks.sprintf("%u\t%lf\t%" PRIu64 "\t", ss[i], exact_ji, us); // sketch size
        mv = arr_mean(a.unions_);
        // Mean US, Mean Error, Mean Bias
        ks.sprintf("%lf\t%lf\t%lf\t", mv, std::accumulate(a.unions_.begin(), a.unions_.end(), 0., [&](auto a, auto b){return a + std::abs(b - us);}) / data.niter_, mv - us);
        iv = arr_mean(a.isns_);
        // Mean IS, IS Error, IS Bias
        ks.sprintf("%lf\t%lf\t%lf\t", iv, std::accumulate(a.isns_.begin(), a.isns_.end(), 0., [&](auto a, auto b){return a + std::abs(b - is);}) / data.niter_, iv - is);
        // Mean JI, JI Error, JI Bias
        jv = arr_mean(a.jis_);
        ks.sprintf("%lf\t%lf\t%lf\n", jv, std::accumulate(a.jis_.begin(), a.jis_.end(), 0., [&](auto a, auto b){return a + std::abs(b - is);}) / data.niter_, jv - exact_ji);
    }
    {
        LockSmith lock(data.m_);
        ks.write(fileno(data.fp_));
    }
    ks.clear();
}
void card_func(void *data_, long index, int tid) {
    aes::AesCtr<uint64_t, 8> gen(index); // Seed with job index
    kth_t &data(*(kth_t *)data_);
    const auto &ss = data.ss_;
    auto &hlls(data.hlls_[tid]);
    const size_t rn = data.rns_[index];
    std::vector<acc> accumulators;
    accumulators.reserve(hlls.size());
    std::generate_n(std::back_emplacer(accumulators), hlls.size(), [rn](){return acc(rn);});
    uint64_t hv;
    hll::WangHash hf;
    size_t isleft;
    for(size_t inum(0); inum < data.niter_; ++inum) {
        static constexpr size_t NPERBUF = gen.BUFSIZE / sizeof(uint64_t);
        GEN_ADD(hlls, rn);
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

#undef GEN_ADD

std::vector<size_t> DEFAULT_RNUMS {
    1ull << 10, 1ull << 12, 1ull << 14, 1ull << 16, 1ull << 18, 1ull << 20, 1ull << 22, 1ull << 24, 1ull << 26, 1ull << 28, 1ull << 30, 1ull << 31, 1ull << 32, 1ull << 33
};

std::vector<unsigned> DEFAULT_SIZES {
    10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22
};

std::vector<double> DEFAULT_JACCS {
/*
 * import numpy as np
 * vals = np.random.random((25,))
 * vals.sort()
 * print(f"    {', '.join(map(str, vals))}")
 */
    0.006722968595635481, 0.017796518494720748, 0.07826379255402727, 0.10464927948391956, 0.12273725109390998, 0.12361422045187453, 0.12545888740551148, 0.12570606740917067, 0.12575697804761177, 0.16168503173635285, 0.1859235765405335, 0.24223370475397865, 0.28467302368862824, 0.3019017718325645, 0.3171949327372131, 0.346532512061814, 0.3895198723519192, 0.4254367781126247, 0.484925922588482, 0.5265848548958305, 0.5367582354724071, 0.5698976125206885, 0.7336933923750625, 0.8532835057417715, 0.9413034301399099
};

template<typename FloatType=float, typename=std::enable_if_t<std::is_arithmetic_v<float>>>
std::vector<FloatType> linspace(size_t ndiv) {
    std::vector<FloatType> ret; ret.reserve(ndiv);
    std::generate_n(std::back_inserter(ret), ndiv, [&]{return static_cast<FloatType>(ret.size()) / ndiv;});
    return ret;
}


int main(int argc, char *argv[]) {
    using bns::ForPool;
    if(argc == 1) usage();
    size_t niter = 50;
    std::vector<size_t> rnum_sizes;
    std::vector<unsigned> sketch_sizes;
    bool emit_full_percentile_range(false), jaccard_exploration(false);
    int c, nthreads(1);
    std::vector<double> jaccs;
    std::FILE *fp = stdout;
    while((c = getopt(argc, argv, "l:o:n:r:s:p:b:J:jPdh?")) >= 0) {
        switch(c) {
            case 'n': niter = std::strtoull(optarg, nullptr, 10); break;
            case 'r':
                rnum_sizes.emplace_back(static_cast<size_t>(std::strtoull(optarg, nullptr, 10))); break;
            case 's':
                sketch_sizes.emplace_back(std::strtoul(optarg, nullptr, 10)); break;
            case 'p': nthreads = std::atoi(optarg); break;
            case 'J': jaccs.emplace_back(std::atof(optarg)); break;
            case 'j': jaccard_exploration = true; break;
            case 'l': jaccs = std::move(linspace<double>(std::atoi(optarg))); break;
            case 'd': rnum_sizes = DEFAULT_RNUMS; sketch_sizes = DEFAULT_SIZES; jaccs = std::move(linspace<double>(100)); break;
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
    std::vector<std::vector<hll::hll_t>> hlls;
    std::generate_n(std::back_emplacer(hlls), nthreads, [&] {
        std::vector<hll::hll_t> ret;
        std::generate_n(std::back_emplacer(ret), sketch_sizes.size(), [&]{return hll::hll_t(sketch_sizes[ret.size()]);});
        return ret;
    });
    std::mutex m;
    kth_t data{rnum_sizes, sketch_sizes, hlls, kstrings, m, niter, fp, emit_full_percentile_range};
    ForPool pool(nthreads);
    if(!jaccard_exploration) {
        std::fprintf(fp, "#Sketch size (log2)\tExact size\tMean bias\tMean error\tMean squared error\tTheoretical Mean Error\tFraction within bounds\tFraction Overestimated\tFraction Underestimated\tError Std Deviation\tBias Std Deviation\t95%% Interval");
        if(emit_full_percentile_range)
            std::fprintf(fp, "90%% Interval\t80%% Interval\t50%% Interval\t10%% Interval\t");
        std::fprintf(fp, "Max Error\n");
        std::fflush(fp);
        pool.forpool(&card_func, (void *)&data, rnum_sizes.size());
    } else {
        if(jaccs.empty()) {
            std::fprintf(stderr, "Error: Some jaccard indices must be specified, either by successive -J{float} calls, -l{ndiv} (e.g., linspace(0, 1, ndiv + 1)[:-1]), or -d {default_jaccs}\n");
            usage();
        }
        kth_jacc_t jacc_data(jaccs, data);
        std::fprintf(fp, "#Sketch size\tExact JI\tExact Union Size\tMean Est US\tMean US error\tMean US Bias\tMean IS\tMean error IS\tMean IS bias\tMean JI\tMean JI error\tMean JI bias\n");
        std::fflush(fp);
        pool.forpool(&jacc_func, (void *)&jacc_data, rnum_sizes.size() * jaccs.size());
    }
    if(fp != stdout) std::fclose(fp);
    return EXIT_SUCCESS;
}
