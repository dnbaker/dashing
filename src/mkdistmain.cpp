#include "dashing.h"

namespace bns {

template<typename T>
void print_vec(const T &x) {
    std::fprintf(stderr, "About to call: '");
    for(const auto v: x)
         std::fprintf(stderr, "%s ", v);
    std::fputs("'\n", stderr);
}
template<typename IT>
void print_vec(IT i, IT i2) {
    for(ssize_t in = 0; in < i2 - i; ++in) {
        std::fprintf(stderr, "%s/%d\n", *(i + in), int(in));
    }
    std::fprintf(stderr, "About to call: '");
    while(i != i2) 
         std::fprintf(stderr, "%s ", *i++);
    std::fputs("'\n", stderr);
}

int mkdist_main(int argc, char *argv[]) {
    auto it = std::find_if(argv, argv + argc, [](auto x) {
        return std::strcmp("--multik", x) == 0;}
    );
    std::string _outpref;
    if(it == argv + argc || argv + argc == ++it) {
        RUNTIME_ERROR("required: --multik <outpref>,<start>,<end> [optional: ,<step>]");
    }
    int e, step, s;
    auto pos = std::strchr(*it, ',');
    if(!pos) RUNTIME_ERROR("Ill-formatted");
    _outpref = std::string(*it, pos - *it);
    s = std::atoi(++pos);
    const char *outpref = _outpref.data();
    std::fprintf(stderr, "outpref: %s\n", _outpref.data());
    std::fprintf(stderr, "pos: %s. start: %d. end: %d\n", pos, s, std::atoi(pos + 1));
    pos = std::strchr(pos, ',');
    e = std::atoi(++pos);
    if((pos = std::strchr(pos, ','))) {
        step = std::atoi(++pos);
        std::fprintf(stderr, "step: %d\n", step);
        assert(step > 0 ? e > s: e < s);
    } else step = e > s ? 1: -1;
    std::fprintf(stderr, "step: %d\n", step);
    std::pair<std::string, std::string> ea(std::string(_outpref.size() + 32, 0), std::string(_outpref.size() + 32, 0));
    //std::vector<char *> args(argv, argv + argc);
    size_t nk = 0;
    std::vector<std::string> fpaths;
    size_t itind = it - 1 - argv;
    argv[0] = const_cast<char *>("dist");
    argv[itind] = &ea.first[0];
    argv[itind + 1] = &ea.second[0];
    assert(!std::strcmp("--multik", argv[itind]));
    for(int ind = s; (e > s ? ind < e: ind > e); ind += step) {
        //args[0] = const_cast<char *>("dist");
        std::sprintf(argv[itind], "-bO_%s_%d", outpref, ind);
        std::sprintf(argv[itind + 1], "-k%d", ind);
        fpaths.push_back(std::string(" ") + std::string(outpref) + ' ' + std::to_string(ind));
        assert(itind != size_t(argc));
        print_vec(argv, argv + argc);
        std::vector<char *> largs(argv, argv + argc);
        /*POST_REQ(*/dist_main(largs.size(), largs.data())/*, "non-zero exit status")*/;
        std::fprintf(stderr, "successfully called\n");
        ++nk;
    }
    std::string outpath = outpref;
    outpath += ".bin";
    return flatten_all(fpaths, nk, outpath);
}
} // namespace bns
