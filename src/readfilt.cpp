#include <omp.h>
#include "bonsai/util.h"
#include "bonsai/database.h"
#include "bonsai/bitmap.h"
#include "sketch/sparse.h"
#include <map>
#include "getopt.h"

using namespace sketch;
using namespace hll;

int usage() {
    std::fprintf(stderr, "readfilt <flags> in.fq [in2.fq]\n-f\tFraction cutoff (0.5)\n-s\tPath to HLL\n-o\toutput [stdout]\n-k\tSet kmer [21]\n");
    return 1;
}

inline int emit(kseq_t *ks1, kseq_t *ks2, std::FILE *ofp, double ci, const std::array<double, 3> &res) {
    if(ks1->qual.s) {
        const char *comment = ks1->comment.s ? ks1->comment.s: "";
        if(unlikely(std::fprintf(ofp, "@%s %s|CI:%lf|%lf|%lf|%lf|\n%s\n+\n%s\n", ks1->name.s, comment, ci, res[0], res[1], res[2], ks1->seq.s, ks1->qual.s) < 0)) return -1;
        if(ks2) {
            comment = ks2->comment.s ? ks2->comment.s: "";
            if(unlikely(std::fprintf(ofp, "@%s %s|%lf\n%s\n+\n%s\n", ks2->name.s, comment, ci, ks2->seq.s, ks2->qual.s) < 0)) return -1;
        }
    } else {
        const char *comment = ks1->comment.s ? ks1->comment.s: "";
        std::fprintf(ofp, ">%s %s|%lf\n%s\n", ks1->name.s, comment, ci, ks1->seq.s);
        if(ks2) {
            comment = ks2->comment.s ? ks2->comment.s: "";
            std::fprintf(ofp, ">%s %s|%lf\n%s\n", ks2->name.s, comment, ci, ks2->seq.s);
        }
    }
    return 0;
}

int main(int argc, char *argv[]) {
    if(argc == 1) return usage();
    std::string hllpath;
    int c, k = 21;
    bool canon = true;
    const char *opath = nullptr;
    double frac_cutoff = 0.5;
    while((c = getopt(argc, argv, "Ch?k:s:f:")) >= 0) {
        switch(c) {
            case 's': hllpath = optarg; break;
            case 'f': frac_cutoff = std::atof(optarg); break;
            case 'k': k = std::atoi(optarg); break;
            case 'o': opath = optarg; break;
            case 'C': canon = false; break;
            case 'h': case '?': return usage();
        }
    }
    std::FILE *ofp = opath ? std::fopen(opath, "w"): stdout;
    std::vector<std::string> inputs(argv + optind, argv + argc);
    if(inputs.empty()) throw "a party";
    if(inputs.size() > 2) throw std::string("WOOO");
    std::fprintf(stderr, "Processing %zu files (%s, %s) with sketch from %s\n", inputs.size(), inputs[0].data(), inputs.size() > 1 ? inputs[1].data(): "single-end", hllpath.data());
    hll_t hll(hllpath);
    hll.csum();
    gzFile ifp1 = gzopen(inputs[0].data(), "rb"), ifp2 = inputs.size() > 1 ? gzopen(inputs[1].data(), "rb"): nullptr;
    if(ifp1 == nullptr) throw 1;
    kseq_t *ks = kseq_init(ifp1), *ks2 = ifp2 ? kseq_init(ifp2): nullptr;
    int rc;
    int p = hll.p();
    std::fprintf(stderr, "Querying with sketch size = %d\n", p);
    bns::Encoder<> enc(k, canon);
    size_t nread = 0;
#if USE_SPARSE
    if(p > 26) {
        throw std::runtime_error("Sparse representation only supported for p <= 26");
    }
    const auto hllhist = hll::detail::sum_counts(hll.core());
    sparse::SparseHLL<> qhll(hll.p());
    std::map<uint32_t, uint8_t> rmap;
#else
    hll_t qhll = hll.clone();
#endif
    while(likely((rc = kseq_read(ks)) >= 0)) {
        if((++nread & 0xFFFu) == 0)
            std::fprintf(stderr, "processed %zu reads\n", nread);
#if USE_SPARSE
        qhll.clear();
#endif
        auto func = [&](uint64_t kmer) {
#if USE_SPARSE
            kmer = hll.hash(kmer);
            auto pos = kmer >> (64 - p);
            uint8_t v = clz(((kmer << 1)|1) << (p - 1)) + 1;
            rmap[pos] = std::max(rmap[pos], v);
#else
            qhll.addh(kmer);
#endif
        };
        enc.for_each(func, ks->seq.s, ks->seq.l);
        if(ks2) {
            if(unlikely((rc = kseq_read(ks2)) < 0)) {
                std::fprintf(stderr, "Warning: mismatched numbers of reads between paired-end files. Error code: %d\n", rc);
                break;
            }
            enc.for_each(func, ks2->seq.s, ks2->seq.l);
        }
#if USE_SPARSE
        auto vals = sparse::pair_query(rmap, hll, &hllhist);
#if VERBOSE_AF
        std::fprintf(stderr, "vals: %lf/%lf/%lf\n", vals[0], vals[1], vals[2]);
#endif
#if !NDEBUG
        qhll.fill_from_pairs(rmap.begin(), rmap.end());
        rmap.clear();
        auto vals = qhll.query(hll, &hllhist);
        double ci = vals[2] / (vals[0] + vals[2]);
#else
        std::fprintf(stderr, "vals: only in query: %lf. only in reference: %lf. shared: %lf\n", vals[0], vals[1], vals[2]);
        double ci = vals[2] / (vals[0] + vals[2]);
#endif
#else
        qhll.csum();
        auto vals = ertl_joint(qhll, hll);
        double ci = vals[2] / (vals[0] + vals[2]);
#endif
        if(ci >= frac_cutoff)
            emit(ks, ks2, ofp, ci, vals);
#if 0
        else {
            std::fprintf(stderr, "ci (with vals %lf, %lf, %lf) of %lf failed\n", vals[0], vals[1], vals[2], ci);
        }
#endif
#if USE_SPARSE
        qhll.clear();
#else
        qhll.reset();
#endif
    }
    gzclose(ifp1);
    if(ifp2) gzclose(ifp2);
    kseq_destroy(ks);
    if(ks2) kseq_destroy(ks2);
    if(ofp != stdout) std::fclose(ofp);
}
