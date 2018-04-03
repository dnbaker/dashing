#define ENABLE_HLL_DEVELOP 1
#include <array>
#include "hll/hll.h"
#include "aesctr.h"
#include "omp.h"
#include "util.h"
#include <thread>
using namespace bns;
using namespace hll;
using namespace hll::detail;

template<typename SketchType>
double fprate(const std::unordered_set<u64> &oset, const SketchType &sketch) {
    size_t nfp = 0;
    for(const auto &el: oset) nfp += sketch.may_contain(el);
    std::fprintf(stderr, "%zu of %zu might be contained.\n", nfp, oset.size());
    return static_cast<double>(nfp) / oset.size();
}

template<typename Hasher>
std::string calculate_errors(size_t ss, size_t nfiltl2, size_t niter, size_t nelem) {
    aes::AesCtr<std::uint64_t, 8> gen((ss + 1) * (nfiltl2 + 3) * (niter + 7) * (nelem + 63));
    std::array<double,8> mdiffs {0.,0.,0.,0.,0.,0.};
    if((int)ss - (int)nfiltl2 < 6) {
        std::fprintf(stderr, "Can't do this many.\n");
        return std::string{};
    }
    hlfbase_t<hllbase_t<Hasher>> hlf(1 << nfiltl2, 137, ss - nfiltl2);
    hllbase_t<Hasher> hll(ss);
    chlf_t<Hasher>   chlf(nfiltl2, ERTL_MLE, ERTL_JOINT_MLE, ss, ss * nelem + niter);
    double frac = 0., fracborrow = 0.;
    double ediff = 0., eabsdiff = 0., oabsdiff = 0;
    double cdiff = 0., cabsdiff = 0.;
    double fphll = 0., fphlf = 0., fpchlf = 0.;
    std::unordered_set<u64> set, oset;
    for(size_t i(0); i < niter; ++i) {
        while(set.size() < nelem) set.emplace(gen());
        while(oset.size() < nelem) {
            u64 val;
            do {val = gen();} while(set.find(val) != set.end());
            oset.insert(val);
        }
        for(const auto &val: set) {
            hll.addh(val);
            hlf.addh(val);
            chlf.addh(val);
        }
        mdiffs[0] += std::abs(nelem - hlf.report());
        mdiffs[1] += std::abs(nelem - hll.report());
        mdiffs[2] += std::abs(nelem - hlf.med_report());
        mdiffs[3] += std::abs(nelem - hlf.chunk_report());
        mdiffs[4] += nelem - hlf.report();
        mdiffs[5] += nelem - hll.report();
        mdiffs[6] += nelem - hlf.med_report();
        mdiffs[7] += nelem - hlf.chunk_report();
        frac += nelem / hll.report();
        fracborrow += nelem / hlf.chunk_report();
        ediff += nelem - ertl_ml_estimate(hll);
        eabsdiff += std::abs(nelem - ertl_ml_estimate(hll));
        auto ctmp = nelem - chlf.chunk_report();
        cdiff += ctmp;
        cabsdiff += std::abs(ctmp);
        double omidd = (ertl_ml_estimate(hll) + hll.report()) * 0.5;
        oabsdiff += std::abs(nelem - omidd);
        fphll += fprate(oset, hll);
        fphlf += fprate(oset, hlf);
        fpchlf += fprate(oset, chlf);
        hlf.clear();
        hll.clear();
        chlf.clear();
        set.clear();
        oset.clear();
    }
    frac /= niter;
    fracborrow /= niter;
    ediff /= niter;
    eabsdiff /= niter;
    oabsdiff /= niter;
    cabsdiff /= niter;
    cdiff /= niter;
    for(auto &md: mdiffs) md *= 1./niter;
    fphll /= niter;
    fphlf /= niter;
    fpchlf /= niter;
    char buf[512];
    std::sprintf(buf, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%zu\t%zu\t%zu\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", mdiffs[0], mdiffs[1], mdiffs[2], mdiffs[3],
                 mdiffs[4], mdiffs[5], mdiffs[6], mdiffs[7], frac, fracborrow, ediff, eabsdiff, ss, size_t(1ull << nfiltl2), nelem, oabsdiff, cdiff, cabsdiff, fphll, fphlf, fpchlf);
    return std::string(buf);
}

struct comb_t {
    size_t size; size_t nelem; size_t nfiltl2;
};

int main(int argc, char *argv[]) {
    unsigned nthreads = argc > 1 ? std::strtoul(argv[1], nullptr, 10): 8;
    size_t niter      = argc > 2 ? std::strtoull(argv[2], nullptr, 10): 50;
    std::string prefix = argc > 3 ? argv[3]: "hlfout";
    std::vector<size_t> sizes    {10, 11, 12, 14, 16, 18, 20};
    std::vector<size_t> nelems   {1 << 18, 1 << 20, 1 << 14, 1 << 12, 1 << 24};
    std::vector<size_t> nfiltl2s {1, 2, 3, 4, 5, 6};
    omp_set_num_threads(nthreads);
    std::vector<comb_t> combs;
    for(auto size: sizes)
        for(auto nelem: nelems)
            for(auto nfiltl2: nfiltl2s)
                combs.emplace_back(comb_t{size, nelem, nfiltl2});
    std::FILE *ofp = std::fopen((prefix + ".wang.out").data(), "w");
    if(ofp == nullptr) throw 1;
    std::fprintf(ofp, "#Mean error hlf\tMean error hll\tMean error hlf median\t""Mean error strength borrowing\t"
                       "Mean diffs hlf\tMean diffs hll\tMean diffs hlf med\tMean diffs strength borrowing\t"
                       "Mean fraction off (hll)\tmean frac off (hlf borrow)\tMean Ertl ML diff\tMean Ertl ML error\t"
                       "sketch size l2\tNumber of subfilters\tnelem\tError using mean of both methods\tMean chlf bias\tMean chlf error\thll fp containment\thlf fp containment\tchlf fp containment\n");
    std::fflush(ofp);
    LOG_INFO("Okay about to do first loop\n");
    #pragma omp parallel for
    for(unsigned i = 0; i < combs.size(); ++i) {
       auto str = calculate_errors<WangHash>(combs[i].size, combs[i].nfiltl2, niter, combs[i].nelem);
       {
            #pragma omp critical
            ::write(fileno(ofp), str.data(), str.size());
       }
    }
    std::fclose(ofp);
    ofp = std::fopen((prefix + ".mur.out").data(), "w");
    if(ofp == nullptr) throw 1;
    std::fprintf(ofp, "#Mean error hlf\tMean error hll\tMean error hlf median\t""Mean error strength borrowing\t"
                       "Mean diffs hlf\tMean diffs hll\tMean diffs hlf med\tMean diffs strength borrowing\t"
                       "Mean fraction off (hll)\tmean frac off (hlf borrow)\tMean Ertl ML diff\tMean Ertl ML error\t"
                       "sketch size l2\tNumber of subfilters\tnelem\tError using mean of both methods\tMean chlf bias\tMean chlf error\thll fp containment\thlf fp containment\tchlf fp containment\n");
    std::fflush(ofp);
    #pragma omp parallel for
    for(unsigned i = 0; i < combs.size(); ++i) {
       auto str = calculate_errors<MurFinHash>(combs[i].size, combs[i].nfiltl2, niter, combs[i].nelem);
       {
            #pragma omp critical
            ::write(fileno(ofp), str.data(), str.size());
       }
    }
    std::fclose(ofp);
    ofp = std::fopen((prefix + ".clhash.out").data(), "w");
    if(ofp == nullptr) throw 1;
    std::fprintf(ofp, "#Mean error hlf\tMean error hll\tMean error hlf median\t""Mean error strength borrowing\t"
                       "Mean diffs hlf\tMean diffs hll\tMean diffs hlf med\tMean diffs strength borrowing\t"
                       "Mean fraction off (hll)\tmean frac off (hlf borrow)\tMean Ertl ML diff\tMean Ertl ML error\t"
                       "sketch size l2\tNumber of subfilters\tnelem\tError using mean of both methods\tMean chlf bias\tMean chlf error\thll fp containment\thlf fp containment\tchlf fp containment\n");
    std::fflush(ofp);
    #pragma omp parallel for
    for(unsigned i = 0; i < combs.size(); ++i) {
       auto str = calculate_errors<clhasher>(combs[i].size, combs[i].nfiltl2, niter, combs[i].nelem);
       {
            #pragma omp critical
            ::write(fileno(ofp), str.data(), str.size());
       }
    }
    std::fclose(ofp);
}
