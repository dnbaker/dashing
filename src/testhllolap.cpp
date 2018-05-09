#include <array>
#include "hll/hll.h"
#include "aesctr.h"
#include "omp.h"
#include "util.h"
#include <thread>
using namespace bns;
using namespace hll;
using namespace hll::detail;

template<typename HllType1, typename HllType2>
ks::string olap_gen(size_t set1, size_t set2, double olap_frac, size_t niter, HllType1 &h1, HllType2 &h2) {
    const auto olap = size_t(olap_frac * std::min(set1, set2));
    const double ji = double(olap) / (set1 + set2 - olap);
    assert(olap <= set1 && olap <= set2);
    aes::AesCtr<std::uint64_t, 8> gen((set1) * (olap + 3) * (set2 + 7) * (olap + 63) + h1.m());
    ks::string ret;
    for(size_t iternum(0); iternum < niter; ++iternum) {
        size_t i;
        for(i = 0; i < olap; ++i) {
            auto val = gen();
            h1.addh(val); h2.addh(val);
        }
        for(i = olap; i < set1; ++i) h1.addh(gen());
        for(i = olap; i < set2; ++i) h2.addh(gen());
    }
    ret.sprintf("%lf\t", ji); // True jaccard
    return ret;
}

struct comb_t {
    size_t size1; size_t size2; double olap_frac; size_t ss;
};

int main(int argc, char *argv[]) {
    unsigned nthreads = argc > 1 ? std::strtoul(argv[1], nullptr, 10):  8;
    size_t niter      = argc > 2 ? std::strtoull(argv[2], nullptr, 10): 25;
    std::vector<size_t> size1s       {12, 14, 16, 18, 20, 24, 26};
    std::vector<size_t> size2s       {12, 14, 16, 18, 20, 24, 26};
    std::vector<size_t> sketch_sizes {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    std::vector<double> olap_fracs;
    {
        unsigned div      = argc > 3 ? std::strtoul(argv[3], nullptr, 10): 20;
        double divinv = 1./div;
        for(double val = 0.; val <= 1.; olap_fracs.push_back(val), val += divinv);
    }
    std::fprintf(stdout, "#Mean error hlf\tMean error hll\tMean error hlf median\t""Mean error strength borrowing\t"
                          "Mean diffs hlf\tMean diffs hll\tMean diffs hlf med\tMean diffs strength borrowing\t"
                          "Mean fraction off (hll)\tmean frac off (hlf borrow)\tMean Ertl ML diff\tMean Ertl ML error\t"
                          "sketch size l2\tNumber of subfilters\tnelem\tError using mean of both methods\n");
    std::fflush(stdout);
    omp_set_num_threads(nthreads);
    std::vector<comb_t> combs;
    for(const auto s1: size1s)
        for(const auto s2: size2s)
            for(const auto frac: olap_fracs)
                for(const auto ss: sketch_sizes)
                    combs.emplace_back(comb_t{s1, s2, frac, ss});
    #pragma omp parallel for
    for(unsigned i = 0; i < combs.size(); ++i) {
        hll_t h1(combs[i].ss), h2(combs[i].ss);
        auto str = olap_gen(combs[i].size1, combs[i].size2, combs[i].olap_frac, niter, h1, h2);
        {
            #pragma omp critical
            ::write(STDOUT_FILENO, str.data(), str.size());
        }
    }
}
