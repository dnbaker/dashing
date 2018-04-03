#include "tinythreadpp/source/fast_mutex.h"
#include <fstream>
#include <omp.h>
#include "bonsai/bonsai/include/util.h"
#include "bonsai/bonsai/include/database.h"
#include "bonsai/bonsai/include/bitmap.h"
#include "bonsai/bonsai/include/setcmp.h"
//#include "klib/kthread.h"
#include <sstream>

namespace bns {

void main_usage(char **argv) {
    std::fprintf(stderr, "Usage: %s <subcommand> [options...]. Use %s <subcommand> for more options. [Subcommands: sketch, dist, setdist, hll.]\n",
                 *argv, *argv);
    std::exit(EXIT_FAILURE);
}

void dist_usage(const char *arg) {
    std::fprintf(stderr, "Usage: %s <opts> [genomes if not provided from a file with -F]\n"
                         "Flags:\n"
                         "-h/-?: Usage\n"
                         "-k\tSet kmer size [31]\n"
                         "-p\tSet number of threads [1]\n"
                         "-s\tadd a spacer of the format <int>x<int>,<int>x<int>,"
                         "..., where the first integer corresponds to the space "
                         "between bases repeated the second integer number of times\n"
                         "-w\tSet window size [max(size of spaced kmer, [parameter])]\n"
                         "-S\tSet sketch size [16, for 2**16 bytes each]\n"
                         "-c:\tCache sketches/use cached sketches\n"
                         "-H:\tTreat provided paths as pre-made sketches.\n"
                         "-C:\tDo not canonicalize. [Default: canonicalize]\n"
                         "-P\tSet prefix for sketch file locations [empty]\n"
                         "-x\tSet suffix in sketch file names [empty]\n"
                         "-o\tOutput for genome size estimates [stdout]\n"
                         "-O\tOutput for genome distance matrix [stdout]\n"
                         "-e\tEmit in scientific notation\n"
                         "-f\tReport results as float. (Only important for binary format.) This halves the memory footprint at the cost of precision loss.\n"
                         "-F\tGet paths to genomes from file rather than positional arguments\n"
                , arg);
    std::exit(EXIT_FAILURE);
}

// Usage, utilities
void sketch_usage(const char *arg) {
    std::fprintf(stderr, "Usage: %s <opts> [genomes if not provided from a file with -F]\n"
                         "Flags:\n"
                         "-h/-?:\tUsage\n"
                         "-k\tSet kmer size [31]\n"
                         "-p\tSet number of threads [1]\n"
                         "-s\tadd a spacer of the format <int>x<int>,<int>x<int>,"
                         "..., where the first integer corresponds to the space "
                         "between bases repeated the second integer number of times\n"
                         "-w\tSet window size [max(size of spaced kmer, [parameter])]\n"
                         "-S\tSet sketch size [16, for 2**16 bytes each]\n"
                         "-F\tGet paths to genomes from file rather than positional arguments\n"
                         "-b:\tBatch size [16 genomes]\n"
                         "-c:\tCache sketches/use cached sketches\n"
                         "-C:\tDo not canonicalize. [Default: canonicalize]\n"
                         "-P\tSet prefix for sketch file locations [empty]\n"
                         "-x\tSet suffix in sketch file names [empty]\n"
                         "-E\tUse Flajolet, not Ertl, quantitation method for hll. [Default: Ertl]\n"
                , arg);
    std::exit(EXIT_FAILURE);
}

std::string hll_fname(const char *path, size_t sketch_p, int wsz, int k, int csz, const std::string &spacing, const std::string &suffix="", const std::string &prefix="") {
    std::string ret(prefix);
    {
        const char *p;
        if(ret.size() && (p = strrchr(get_cstr(path), '/')))
            ret += std::string(p);
        else
            ret += get_cstr(path);
    }
    ret += ".w";
    ret + std::to_string(std::max(csz, wsz));
    ret += ".";
    ret += std::to_string(k);
    ret += ".spacing";
    ret += spacing;
    ret += '.';
    if(suffix.size()) {
        ret += "suf";
        ret += suffix;
        ret += '.';
    }
    ret += std::to_string(sketch_p);
    ret += ".hll";
    return ret;
}


namespace detail {

struct kt_sketch_helper {
    std::vector<hll::hll_t>  &hlls_; // HyperLogLog scratch space
    std::vector<kseq_t>     &kseqs_;
    const int bs_, sketch_size_, kmer_size_, window_size_, csz_;
    const spvec_t sv_;
    std::vector<std::vector<std::string>> &ssvec_; // scratch string vector
    const std::string &suffix_;
    const std::string &prefix_;
    const std::string &spacing_;
    const bool skip_cached_;
    const bool canon_;
    const hll::EstimationMethod estim_;
    const bool write_to_dev_null_;
    const bool write_gz_;
};

void kt_for_helper(void  *data_, long index, int tid) {
    kt_sketch_helper &helper(*(kt_sketch_helper *)data_);
    hll::hll_t &hll(helper.hlls_[tid]);
    std::string fname;
    for(size_t i(helper.bs_ * index); i < std::min((uint64_t)(helper.bs_ * (index + 1)), (uint64_t)(helper.ssvec_.size())); ++i) {
        std::vector<std::string> &scratch_stringvec(helper.ssvec_[i]);
        fname = hll_fname(scratch_stringvec[0].data(), helper.sketch_size_, helper.window_size_, helper.kmer_size_, helper.csz_, helper.spacing_, helper.suffix_, helper.prefix_);
        if(helper.write_gz_) fname += ".gz";
        if(helper.skip_cached_ && isfile(fname)) continue;
        fill_hll(hll, scratch_stringvec, helper.kmer_size_, helper.window_size_, helper.sv_, helper.canon_, nullptr, 1, helper.sketch_size_, &helper.kseqs_[tid]); // Avoid allocation fights.
        hll.set_estim(helper.estim_);
        hll.write(helper.write_to_dev_null_ ? "/dev/null": fname.data(), helper.write_gz_);
        hll.clear();
    }
}

}

// Main functions
int sketch_main(int argc, char *argv[]) {
    int wsz(-1), k(31), sketch_size(16), skip_cached(false), co, nthreads(1), bs(16);
    bool canon(true), write_to_dev_null(false), write_gz(false), clamp(true);
    hll::EstimationMethod estim = hll::EstimationMethod::ERTL_MLE;
    hll::JointEstimationMethod jestim = hll::JointEstimationMethod::ERTL_JOINT_MLE;
    std::string spacing, paths_file, suffix, prefix;
    while((co = getopt(argc, argv, "P:F:c:p:x:s:S:k:w:jLzEDIcCeh?")) >= 0) {
        switch(co) {
            case 'L': clamp = false; break;
            case 'b': bs = std::atoi(optarg); break;
            case 'k': k = std::atoi(optarg); break;
            case 'x': suffix = optarg; break;
            case 'p': nthreads = std::atoi(optarg); break;
            case 'P': prefix = optarg; break;
            case 's': spacing = optarg; break;
            case 'E': jestim = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ORIGINAL); break;
            case 'I': jestim = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ERTL_IMPROVED); break;
            case 'J': jestim = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ERTL_MLE); break;
            case 'S': sketch_size = std::atoi(optarg); break;
            case 'w': wsz = std::atoi(optarg); break;
            case 'c': skip_cached = true; break;
            case 'F': paths_file = optarg; break;
            case 'C': canon = false; break;
            case 'D': write_to_dev_null = true; break;
            case 'z': write_gz = true; break;
            case 'h': case '?': sketch_usage(*argv);
        }
    }
    LOG_DEBUG("Sketch size: %i\n", sketch_size);
    LOG_DEBUG("Using %zu threads\n", nthreads);
    omp_set_num_threads(nthreads);
    spvec_t sv(parse_spacing(spacing.data(), k));
    Spacer sp(k, wsz, sv);
    std::vector<std::vector<std::string>> ivecs;
    {
        std::vector<std::string> inpaths(paths_file.size() ? get_paths(paths_file.data())
                                                           : std::vector<std::string>(argv + optind, argv + argc));
        for(const auto &el: inpaths) ivecs.emplace_back(std::vector<std::string>{el});
    }
    std::vector<hll::hll_t> hlls;
    std::vector<kseq_t>    kseqs;
    while(hlls.size() < (unsigned)nthreads) hlls.emplace_back(sketch_size, estim, jestim, 1, clamp), kseqs.emplace_back(kseq_init_stack());
    for(auto &hll: hlls) hll.set_clamp(clamp);
    assert(hlls[0].size() == ((1ull << sketch_size)));
    if(wsz < sp.c_) wsz = sp.c_;
    if(ivecs.size() == 0) {
        std::fprintf(stderr, "No paths. See usage.\n");
        sketch_usage(*argv);
    }
    if(ivecs.size() / (unsigned)(nthreads) > (unsigned)bs) bs = (ivecs.size() / (nthreads) / 2);
    detail::kt_sketch_helper helper {hlls, kseqs, bs, sketch_size, k, wsz, (int)sp.c_, sv, ivecs, suffix, prefix, spacing, skip_cached, canon, estim, write_to_dev_null, write_gz};
    kt_for(nthreads, detail::kt_for_helper, &helper, ivecs.size() / bs + (ivecs.size() % bs != 0));
    for(auto &kseq: kseqs) {
        kseq_destroy_stack(kseq);
    }
    LOG_DEBUG("Finished sketching\n");
    return EXIT_SUCCESS;
}

struct LockSmith {
    // Simple lock-holder to avoid writing to the same file twice.
    tthread::fast_mutex &m_;
    LockSmith(tthread::fast_mutex &m): m_(m) {m_.lock();}
    ~LockSmith() {m_.unlock();}
};

template<typename FType, typename=std::enable_if_t<std::is_floating_point_v<FType>>>
size_t submit_emit_dists(const int pairfi, const FType *ptr, u64 hs, size_t index, ks::string &str, const std::vector<std::string> &inpaths, bool write_binary, bool use_scientific, const size_t buffer_flush_size=1ull<<18) {
    if(write_binary) {
       ::write(pairfi, ptr, sizeof(FType) * hs);
    } else {
        const char *const fmt(use_scientific ? "%e\t": "%f\t");
        str += inpaths[index];
        str.putc_('\t');
        {
            u64 k;
            for(k = 0; k < index + 1;  ++k, str.putsn_("-\t", 2));
            for(k = 0; k < hs - index - 1; str.sprintf(fmt, ptr[k++]));
        }
        str.back() = '\n';
        if(str.size() >= 1 << 18) str.write(pairfi), str.clear();
    }
    return index;
}

template<typename FType1, typename FType2,
         typename=std::enable_if_t<std::is_floating_point_v<FType1> && std::is_floating_point_v<FType2>>>
FType2 dist_index(FType1 ji, FType2 ksinv) {
    // Adapter from Mash https://github.com/Marbl/Mash
    return -std::log(2. * ji / (1. + ji)) * ksinv;
}

template<typename FType, typename=std::enable_if_t<std::is_floating_point_v<FType>>>
void dist_loop(const int pairfi, std::vector<hll::hll_t> &hlls, const std::vector<std::string> &inpaths, const bool use_scientific, const unsigned k, const bool emit_jaccard, bool write_binary, const size_t buffer_flush_size=1ull<<18) {
    std::array<std::vector<FType>, 2> dps;
    dps[0].resize(hlls.size() - 1);
    dps[1].resize(hlls.size() - 2);
    ks::string str;
    const FType ksinv = 1./ k;
    std::future<size_t> submitter;
    for(size_t i = 0; i < hlls.size(); ++i) {
        hll::hll_t &h1(hlls[i]); // TODO: consider working backwards and pop_back'ing.
        std::vector<FType> &dists = dps[i & 1];
        if(emit_jaccard) {
            #pragma omp parallel for schedule(dynamic)
            for(size_t j = i + 1; j < hlls.size(); ++j) {
                dists[j - i - 1] = jaccard_index(hlls[j], h1);
            }
        } else {
            #pragma omp parallel for schedule(dynamic)
            for(size_t j = i + 1; j < hlls.size(); ++j) {
                dists[j - i - 1] = dist_index(jaccard_index(hlls[j], h1), ksinv);
            }
        }
        h1.free();
        LOG_DEBUG("Finished chunk %zu of %zu\n", i + 1, hlls.size());
#if !NDEBUG
        if(i) LOG_DEBUG("Finished writing row %zu\n", submitter.get());
#else
        if(i) submitter.get();
#endif
        submitter = std::async(std::launch::async, submit_emit_dists<FType>, pairfi, dists.data(), hlls.size(), i, std::ref(str), std::ref(inpaths), write_binary, use_scientific, buffer_flush_size);
    }
    submitter.get();
    if(!write_binary) str.write(pairfi), str.clear();
}

enum CompReading: unsigned {
    UNCOMPRESSED,
    GZ,
    AUTODETECT
};

int dist_main(int argc, char *argv[]) {
    int wsz(-1), k(31), sketch_size(16), use_scientific(false), co, cache_sketch(false), nthreads(1);
    bool canon(true), presketched_only(false), write_binary(false), emit_jaccard(true), emit_float(false), clamp(true);
    hll::EstimationMethod estim = hll::EstimationMethod::ERTL_MLE;
    hll::JointEstimationMethod jestim = hll::JointEstimationMethod::ERTL_JOINT_MLE;
    std::string spacing, paths_file, suffix, prefix;
    CompReading reading_type = UNCOMPRESSED;
    FILE *ofp(stdout), *pairofp(stdout);
    omp_set_num_threads(1);
    while((co = getopt(argc, argv, "P:x:F:c:p:o:s:w:O:S:k:azfJICbMEeHh?")) >= 0) {
        switch(co) {
            case 'z': reading_type = GZ; break;
            case 'a': reading_type = AUTODETECT; break;
            case 'C': canon = false; break;
            case 'c': clamp = false; break;
            case 'E': jestim = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ORIGINAL); break;
            case 'F': paths_file = optarg; break;
            case 'H': presketched_only = true; break;
            case 'I': jestim = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ERTL_IMPROVED); break;
            case 'm': jestim = (hll::JointEstimationMethod)(estim = hll::EstimationMethod::ERTL_MLE); break;
            case 'J': emit_jaccard = false; break;
            case 'O': pairofp = fopen(optarg, "wb"); if(pairofp == nullptr) LOG_EXIT("Could not open file at %s for writing.\n", optarg); break;
            case 'P': prefix = optarg; break;
            case 'S': sketch_size = std::atoi(optarg); break;
            case 'W': cache_sketch = true;  break;
            case 'b': write_binary = true; break;
            case 'e': use_scientific = true; break;
            case 'f': emit_float = true; break;
            case 'k': k = std::atoi(optarg); break;
            case 'o': ofp = fopen(optarg, "w"); if(ofp == nullptr) LOG_EXIT("Could not open file at %s for writing.\n", optarg); break;
            case 'p': nthreads = std::atoi(optarg); break;
            case 's': spacing = optarg; break;
            case 'w': wsz = std::atoi(optarg); break;
            case 'x': suffix = optarg; break;
            case 'h': case '?': dist_usage(*argv);
        }
    }
    spvec_t sv(parse_spacing(spacing.data(), k));
    Spacer sp(k, wsz, sv);
    std::vector<std::string> inpaths(paths_file.size() ? get_paths(paths_file.data())
                                                       : std::vector<std::string>(argv + optind, argv + argc));
    if(inpaths.size() == 0) {
        std::fprintf(stderr, "No paths. See usage.\n");
        dist_usage(*argv);
    }
    omp_set_num_threads(nthreads);
    std::vector<hll::hll_t> hlls(inpaths.size(), presketched_only ? hll::hll_t(): hll::hll_t(sketch_size, estim, jestim, 1, clamp));
    {
        // Scope to force deallocation of scratch_vv.
        std::vector<std::vector<std::string>> scratch_vv(nthreads, std::vector<std::string>{"empty"});
        if(wsz < sp.c_) wsz = sp.c_;
        #pragma omp parallel for
        for(size_t i = 0; i < hlls.size(); ++i) {
            const std::string &path(inpaths[i]);
            static const std::string suf = ".gz";
            if(presketched_only) {
                hlls[i].read(path);
            } else {
                const std::string fpath(hll_fname(path.data(), sketch_size, wsz, k, sp.c_, spacing, suffix, prefix));
                if(cache_sketch && isfile(fpath)) {
                    LOG_DEBUG("Sketch found at %s with size %zu, %u\n", fpath.data(), size_t(1ull << sketch_size), sketch_size);
                    hlls[i].read(fpath);
                } else {
                    // By reserving 256 character, we make it probably that no allocation is necessary in this section.
                    std::vector<std::string> &scratch_stringvec(scratch_vv[omp_get_thread_num()]);
                    scratch_stringvec[0] = inpaths[i];
                    fill_hll(hlls[i], scratch_stringvec, k, wsz, sv, canon, nullptr, 1, sketch_size); // Avoid allocation fights.
                    if(cache_sketch) hlls[i].write(fpath, (reading_type == GZ ? 1: reading_type == AUTODETECT ? std::equal(suf.rbegin(), suf.rend(), fpath.rbegin()): false));
                }
            }
        }
    }
    ks::string str("#Path\tSize (est.)\n");
    assert(str == "#Path\tSize (est.)\n");
    str.resize(1 << 18);
    {
        const int fn(fileno(ofp));
        for(size_t i(0); i < hlls.size(); ++i) {
            str.sprintf("%s\t%lf\n", inpaths[i].data(), hlls[i].report());
            if(str.size() >= 1 << 18) str.write(fn), str.clear();
        }
        str.write(fn), str.clear();
    }
    // TODO: Emit overlaps and symmetric differences.
    if(ofp != stdout) std::fclose(ofp);
    str.clear();
    if(write_binary) {
        const size_t hs(hlls.size());
        std::fwrite(&hs, sizeof(hs), 1, pairofp);
    } else {
        str.sprintf("##Names \t");
        for(const auto &path: inpaths) str.sprintf("%s\t", path.data());
        str.back() = '\n';
        str.write(fileno(pairofp)); str.free();
    }
    if(emit_float)
        dist_loop<float>(fileno(pairofp), hlls, inpaths, k, use_scientific, emit_jaccard, write_binary);
    else
        dist_loop<double>(fileno(pairofp), hlls, inpaths, k, use_scientific, emit_jaccard, write_binary);
    if(pairofp != stdout) std::fclose(pairofp);
    return EXIT_SUCCESS;
}

int setdist_main(int argc, char *argv[]) {
    int wsz(-1), k(31), use_scientific(false), co;
    bool canon(true);
    unsigned bufsize(1 << 18);
    std::string spacing, paths_file;
    FILE *ofp(stdout), *pairofp(stdout);
    omp_set_num_threads(1);
    while((co = getopt(argc, argv, "F:c:p:o:O:S:B:k:CMeh?")) >= 0) {
        switch(co) {
            case 'B': std::stringstream(optarg) << bufsize; break;
            case 'k': k = std::atoi(optarg); break;
            case 'p': omp_set_num_threads(std::atoi(optarg)); break;
            case 's': spacing = optarg; break;
            case 'C': canon = false; break;
            case 'w': wsz = std::atoi(optarg); break;
            case 'F': paths_file = optarg; break;
            case 'o': ofp = fopen(optarg, "w"); break;
            case 'O': pairofp = fopen(optarg, "w"); break;
            case 'e': use_scientific = true; break;
            case 'h': case '?': dist_usage(*argv);
        }
    }
    std::vector<char> rdbuf(bufsize);
    spvec_t sv(parse_spacing(spacing.data(), k));
    Spacer sp(k, wsz, sv);
    std::vector<std::string> inpaths(paths_file.size() ? get_paths(paths_file.data())
                                                       : std::vector<std::string>(argv + optind, argv + argc));
    std::vector<khash_t(all) *> hashes;
    while(hashes.size() < inpaths.size()) hashes.emplace_back((khash_t(all) *)calloc(sizeof(khash_t(all)), 1));
    const size_t nhashes(hashes.size());
    if(wsz < sp.c_) wsz = sp.c_;
    if(inpaths.size() == 0) {
        std::fprintf(stderr, "No paths. See usage.\n");
        dist_usage(*argv);
    }
    #pragma omp parallel for
    for(size_t i = 0; i < hashes.size(); ++i) {
        const char *path(inpaths[i].data());
        khash_t(all) *hash(hashes[i]);
        fill_set_genome<score::Lex>(path, sp, hash, i, nullptr, canon);
    }
    LOG_DEBUG("Filled genomes. Now analyzing data.\n");
    ks::string str;
    str.sprintf("#Path\tSize (est.)\n");
    {
        const int fn(fileno(ofp));
        for(size_t i(0); i < hashes.size(); ++i) {
            str.sprintf("%s\t%zu\n", inpaths[i].data(), kh_size(hashes[i]));
            if(str.size() > 1 << 17) str.write(fn), str.clear();
        }
        str.write(fn), str.clear();
    }
    // TODO: Emit overlaps and symmetric differences.
    if(ofp != stdout) std::fclose(ofp);
    std::vector<double> dists(nhashes - 1);
    str.clear();
    str.sprintf("##Names \t");
    for(auto &path: inpaths) str.sprintf("%s\t", path.data());
    str.back() = '\n';
    str.write(fileno(pairofp)); str.free();
    setvbuf(pairofp, rdbuf.data(), _IOLBF, rdbuf.size());
    const char *const fmt(use_scientific ? "\t%e": "\t%f");
    for(size_t i = 0; i < hashes.size(); ++i) {
        auto &h1(hashes[i]);
        size_t j;
        #pragma omp parallel for
        for(j = i + 1; j < hashes.size(); ++j)
            dists[j - i - 1] = jaccard_index(hashes[j], h1);
        for(j = 0; j < i + 1; ++j) fputc('\t', pairofp), fputc('-', pairofp);
        for(j = 0; j < hashes.size() - i - 1; ++j)
            fprintf(pairofp, fmt, dists[j]);
        fputc('\n', pairofp);
        khash_destroy(h1), h1 = nullptr;
        // Delete data as soon as we don't need it.
    }
    return EXIT_SUCCESS;
}

int hll_main(int argc, char *argv[]) {
    int c, wsz(-1), k(31), num_threads(-1), sketch_size(24);
    bool canon(true);
    std::string spacing, paths_file;
    std::ios_base::sync_with_stdio(false);
    if(argc < 2) {
        usage: LOG_EXIT("Usage: %s <opts> <paths>\nFlags:\n"
                        "-k:\tkmer length (Default: 31. Max: 31)\n"
                        "-w:\twindow size (Default: -1)  Must be -1 (ignored) or >= kmer length.\n"
                        "-s:\tspacing (default: none). format: <value>x<times>,<value>x<times>,...\n"
                        "   \tOmitting x<times> indicates 1 occurrence of spacing <value>\n"
                        "-S:\tsketch size (default: 24). (Allocates 2 << [param] bytes of memory per HyperLogLog.\n"
                        "-p:\tnumber of threads.\n"
                        "-F:\tPath to file which contains one path per line\n"
                        , argv[0]);
    }
    while((c = getopt(argc, argv, "Cw:s:S:p:k:tfh?")) >= 0) {
        switch(c) {
            case 'C': canon = false; break;
            case 'h': case '?': goto usage;
            case 'k': k = std::atoi(optarg); break;
            case 'p': num_threads = std::atoi(optarg); break;
            case 's': spacing = optarg; break;
            case 'S': sketch_size = std::atoi(optarg); break;
            case 'w': wsz = std::atoi(optarg); break;
            case 'F': paths_file = optarg; break;
        }
    }
    if(wsz < k) wsz = k;
    std::vector<std::string> inpaths(paths_file.size() ? get_paths(paths_file.data())
                                                       : std::vector<std::string>(argv + optind, argv + argc));
    spvec_t sv(parse_spacing(spacing.data(), k));
    LOG_INFO("Processing %zu paths with %i threads\n", inpaths.size(), num_threads);
    const double est(estimate_cardinality<score::Lex>(inpaths, k, wsz, sv, canon, nullptr, num_threads, sketch_size));
    std::fprintf(stdout, "Estimated number of unique exact matches: %lf\n", est);
    return EXIT_SUCCESS;
}

} // namespace bns

using namespace bns;

int main(int argc, char *argv[]) {
    if(argc == 1) main_usage(argv);
    if(std::strcmp(argv[1], "sketch") == 0) return sketch_main(argc - 1, argv + 1);
    else if(std::strcmp(argv[1], "dist") == 0) return dist_main(argc - 1, argv + 1);
    else if(std::strcmp(argv[1], "setdist") == 0) return setdist_main(argc - 1, argv + 1);
    else if(std::strcmp(argv[1], "hll") == 0) return hll_main(argc - 1, argv + 1);
    else throw std::runtime_error(std::string("Invalid subcommand ") + argv[1] + " provided.");
}
