#include <vector>
#include <string>

std::vector<float> nuc_freqs(const std::vector<std::string> &fpaths, unsigned nt);

static inline auto background_match(unsigned lhid, rhid, const std::vector<float> &nucfreqs) {
    auto lhp = nucfreqs.data() + (lhid * 4), rhp = nucfreqs.data() + (rhid * 4);
    return lhp[0] * rhp[0] + lhp[1] * rhp[1] +
           lhp[2] * rhp[2] + lhp[3] * rhp[3];
}
