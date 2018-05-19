#define ENABLE_HLL_DEVELOP 1
#include <array>
#include "hll/hll.h"
#include "hll/aesctr/aesctr.h"
#include "omp.h"
#include "util.h"
#include <thread>
#include <unordered_set>
#include "bonsai/bonsai/include/logutil.h"
#include "bonsai/clhash/include/clhash.h"
using u64 = std::uint64_t;
using namespace sketch;
using namespace hll;
using namespace hll::detail;

template<typename SketchType>
double fprate(const std::unordered_set<u64> &oset, const SketchType &sketch) {
    size_t nfp = 0;
    for(const auto &el: oset) nfp += sketch.may_contain(el);
#if !NDEBUG
    std::fprintf(stderr, "%zu of %zu might be contained.\n", nfp, oset.size());
#endif
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
        // Generate and insert values.
        while(__builtin_expect(nelem - set.size() > gen.BUFSIZE, 1)) {
            for(const auto &val: gen.template view<std::uint64_t>()) if(set.find(val) != set.end()) set.insert(val);
            gen.generate_new_values();
        }
        while(__builtin_expect(nelem - oset.size() > gen.BUFSIZE, 1)) {
            for(const auto &val: gen.template view<std::uint64_t>()) oset.insert(val);
            gen.generate_new_values();
        }
        while(oset.size() < nelem) oset.insert(gen());
        while(set.size() < nelem)  set.insert(gen());
        for(const auto &val: set) hll.addh(val), hlf.addh(val), chlf.addh(val);


        // Woof.
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
        // Clear structures.
        hlf.clear(), hll.clear(), chlf.clear(), set.clear(), oset.clear();
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

void usage() {
    std::fprintf(stderr, "Usage:\ntesthlf <nthreads=8>, <niter=50>, <prefix='hlfout'>\n");
    std::exit(EXIT_FAILURE);
}

int main(int argc, char *argv[]) {
    if(auto it = std::find_if(argv, argv + argc, [](const auto x){return std::strcmp(x, "-h") == 0;}); it != argv + argc) usage();
    unsigned nthreads = argc > 1 ? std::strtoul(argv[1], nullptr, 10): 8;
    size_t niter      = argc > 2 ? std::strtoull(argv[2], nullptr, 10): 50;
    std::string prefix = argc > 3 ? argv[3]: "hlfout";
    std::vector<size_t> sizes    {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    std::vector<size_t> nelems   {1 << 18, 1 << 20, 1 << 14, 1 << 12, 1 << 24, 1 << 26};
    std::vector<size_t> nfiltl2s {1, 2, 3, 4, 5, 6};
    omp_set_num_threads(nthreads);
    std::vector<comb_t> combs;
    for(auto size: sizes)
        for(auto nelem: nelems)
            for(auto nfiltl2: nfiltl2s)
                combs.emplace_back(comb_t{size, nelem, nfiltl2});
    std::FILE *ofp = std::fopen((prefix + ".wang.out").data(), "w");
    if(ofp == nullptr) throw std::runtime_error("Could not open file at "s + prefix + ".wang.out" + " for reading");
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
    if((ofp = std::fopen((prefix + ".mur.out").data(), "w")) == nullptr) throw std::runtime_error("Could not open file at "s + prefix + ".mur.out" + " for reading");
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
    if((ofp = std::fopen((prefix + ".clhash.out").data(), "w")) == nullptr)
        throw std::runtime_error("Could not open file at "s + prefix + ".clhash.out" + " for reading");
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
