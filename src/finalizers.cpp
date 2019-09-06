#include "dashing.h"

namespace bns {
template<> void sketch_finalize<khset64_t>(khset64_t &x) {x.cvt2shs();}
#define CONTAIN_OVERLOAD_FAIL(x)\
template<>\
double containment_index<x>(const x &b, const x &a) {\
    RUNTIME_ERROR(std::string("Containment index not implemented for ") + __PRETTY_FUNCTION__);\
}\
template<>\
std::array<double, 3> set_triple<x>(const x &b, const x &a) {\
    RUNTIME_ERROR(std::string("set_triple not implemented for ") + __PRETTY_FUNCTION__);\
}
CONTAIN_OVERLOAD_FAIL(RMFinal)
CONTAIN_OVERLOAD_FAIL(bf::bf_t)
CONTAIN_OVERLOAD_FAIL(wj::WeightedSketcher<RMFinal>)
CONTAIN_OVERLOAD_FAIL(wj::WeightedSketcher<bf::bf_t>)
CONTAIN_OVERLOAD_FAIL(CRMFinal)
#undef CONTAIN_OVERLOAD_FAIL
namespace detail {
void sort_paths_by_fsize(std::vector<std::string> &paths) {
    if(paths.size() < 2) return;
    uint32_t *fsizes = static_cast<uint32_t *>(std::malloc(paths.size() * sizeof(uint32_t)));
    if(!fsizes) throw std::bad_alloc();
    #pragma omp parallel for
    for(size_t i = 0; i < paths.size(); ++i)
        fsizes[i] = posix_fsizes(paths[i].data());
    std::vector<path_size> ps(paths.size());
    #pragma omp parallel for
    for(size_t i = 0; i < paths.size(); ++i)
        ps[i] = path_size(paths[i], fsizes[i]);
    std::free(fsizes);
    std::sort(ps.begin(), ps.end(), [](const auto &x, const auto &y) {return x.size > y.size;});
    paths.clear();
    for(const auto &p: ps) paths.emplace_back(std::move(p.path));
}
} // detail
size_t posix_fsizes(const std::string &path, const char sep) {
    size_t ret = 0;
    for_each_substr([&ret](const char *s) {struct stat st; ::stat(s, &st); ret += st.st_size;}, path, sep);
    return ret;
}
} // bns
