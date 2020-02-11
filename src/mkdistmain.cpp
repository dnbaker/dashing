#include "dashing.h"
#include "kspp/ks.h"
#include "background.h"

namespace bns {

void mkdist_usage() {
    std::fprintf(stderr, "mkdist: --multik <outpref>,<start>,<end> [optional: ,<step>], plus all dist usage options\n");
    dist_usage(bns::executable);
}

int mkdist_main(int argc, char *argv[]) {
    std::string path_to_dashing="dashing";
    for(const auto s: {"-h", "--help"}) {
        if(std::find_if(argv, argv + argc, [s](auto x) {return std::strcmp(x, s) == 0;}) != argv + argc)
            mkdist_usage();
    }
    auto it = std::find_if(argv, argv + argc, [](auto x) {
        return std::strcmp("--multik", x) == 0;
    });
    std::string _outpref;
    if(it == argv + argc || argv + argc == ++it) {
        std::fprintf(stderr, "required: --multik <outpref>,<start>,<end> [optional: ,<step>]\nIf dashing is not in $PATH, specify it with -D<path_to_dashing>\n");
        mkdist_usage();
    }
    int e, step, s;
    auto pos = std::strchr(*it, ',');
    if(!pos) {
        std::fprintf(stderr, "Ill-formatted");
        mkdist_usage();
    }
    _outpref = std::string(*it, pos - *it);
    s = std::atoi(++pos);
    const char *outpref = _outpref.data();
    pos = std::strchr(pos, ',');
    e = std::atoi(++pos);
    if((pos = std::strchr(pos, ','))) {
        step = std::atoi(++pos);
        assert(step > 0 ? e > s: e < s);
    } else step = e > s ? 1: -1;
    std::fprintf(stderr, "step: %d\n", step);
    std::pair<std::string, std::string> ea(std::string(_outpref.size() + 32, 0), std::string(_outpref.size() + 32, 0));
    //std::vector<char *> args(argv, argv + argc);
    std::vector<std::string> fpaths;
    size_t itind = it - 1 - argv;
    argv[0] = const_cast<char *>("dist");
    argv[itind] = &ea.first[0];
    argv[itind + 1] = &ea.second[0];
    assert(!std::strcmp("--multik", argv[itind]));
    std::vector<unsigned> k_values;
    std::vector<std::string> fnames;
    bool binary = false;
    std::vector<char *> cargs(argv, argv + argc);
    for(int c;(c = getopt(cargs.size(), cargs.data(), "n:Q:P:x:F:c:p:o:s:w:O:S:k:=:t:R:8TgazlICbMEeHJhZBNyUmqW?D:-")) >= 0;) {
        // Ignore all other options, pass along to dist_main, which ignores the -D option
        if(c == 'D')
            path_to_dashing = const_cast<const char *>(optarg);
        else if(c == 'F') {
            fnames = get_paths(optarg); break;
        } else if(c == 'l' || c == 'M' || c == 133 || c == 131 || c == 132 || c == 137 || c == 138) {
            throw std::runtime_error("multi-k needs to use intersection size currently");
        } else if(c == 'b') binary = true;
    }
    if(optind == argc && fnames.empty()) {
        RUNTIME_ERROR("fnames must be specified by -F or as positional arguments\n");
    }
    if(fnames.empty()) {
        fnames.assign(argv + optind, argv + argc);
    }
    for(int ind = s; (e > s ? ind < e: ind > e); ind += step) {
        k_values.push_back(ind);
        std::string sizes_name = std::string("-o_") + outpref + "_" + std::to_string(ind);
        std::sprintf(argv[itind], "-bO_%s_%d", outpref, ind);
        std::sprintf(argv[itind + 1], "-k%d", ind);
        fpaths.push_back(std::string("_") + std::string(outpref) + '_' + std::to_string(ind));
        assert(itind != size_t(argc));
        std::vector<char *> largs(argv, argv + argc);
        largs.push_back(&sizes_name[0]);
        ks::string cmd = path_to_dashing;
        for(const auto arg: largs)
            cmd.sprintf(" %s", arg);
#ifndef NDEBUG
        std::fprintf(stderr, "About to call: '%s'\n", cmd.data());
#endif
        int rc = std::system(cmd.data());
        POST_REQ(rc >= 0, std::strerror(rc));
        POST_REQ(WIFEXITED(rc), "non-zero exit status");
    }
    std::fprintf(stderr, "Finished distance matrix calculations. Now flattening\n");
    std::string outpath = outpref;
    outpath += ".bin";
    auto ret = flatten_all(fpaths, outpath, k_values);
    sumstats fstatistics = nuc_freqs(fpaths);
    auto kcdistmat = mkmat2jcdistmat(outpath, SIZES, fstatistics.sizes().data(),
                                     fstatistics.freqs().data(), fstatistics.numseqs().data(),
                                     true);
    if(binary) {
        kcdistmat.write(outpath + ".jcd");
    } else {
        std::FILE *ofp = std::fopen((outpath + ".jcd").data(), "wb");
        kcdistmat.write(ofp);
        fclose(ofp);
    }
    for(const auto &f: fpaths)
        std::system(ks::sprintf("rm %s %s.labels", f.data(), f.data()).data());
    return ret;
} // mkdist_main
} // namespace bns
