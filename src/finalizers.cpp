#include "dashing.h"

namespace bns {
//template<> void sketch_finalize<khset64_t>(khset64_t &x) {x.cvt2shs();}
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
