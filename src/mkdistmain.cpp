#include "dashing.h"

namespace bns {

int mkdist_main(int argc, char *argv[]) {
    auto it = std::find_if(argv, argv + argc, [](auto x) {return std::strcmp("--multik", x) == 0;});
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
    e = std::atoi(++pos);
    if((pos = std::strchr(pos, ','))) {
        step = std::atoi(++pos);
        assert(step > 0 ? e > s: e < s);
    } else step = e > s ? 1: -1;
    std::vector<std::string> ea(2);
    for(auto &i: ea) i.resize(_outpref.size() + 32);
    std::vector<char *> args(argv, argv + argc);
    size_t nk = 0;
    std::vector<std::string> fpaths;
    for(int ind = s; (e > s ? ind < e: ind > e); ind += step) {
        char buf[256]{0};
        ea[0] = std::string(buf, std::sprintf(buf, "-bO_%s_%d", outpref, ind));
        ea[1] = std::string(buf, std::sprintf(buf, "-k%d", ind));
        fpaths.push_back(std::string(buf, std::sprintf(buf, "_%s_%d", outpref, ind)));
        size_t j = std::find_if(argv, argv + argc, [](auto x) {return std::strcmp("--multik", x) == 0;}) - argv;
        assert(j != size_t(argc));
        args[j] = &ea[0][0];
        args[1] = &ea[1][0];
        dist_main(args.size(), args.data());
        ++nk;
    }
    std::string outpath = outpref;
    outpath += ".bin";
    return flatten_all(fpaths, nk, outpath);
}
} // namespace bns
