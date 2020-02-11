#include "dashing.h"


namespace bns {
int flatten_all(const std::vector<std::string> &fpaths, size_t nk, const std::string outpath) {
    if(fpaths.empty()) RUNTIME_ERROR("no fpaths, see usage.");
    // TODO: adapt this to pack an asymmetric comparison result into a flattened blaze distance matrix.
    std::vector<dm::DistanceMatrix<float>> dms;
    dms.reserve(nk);
    for(const auto &fp: fpaths)
        dms.emplace_back(fp.data());
    const uint64_t ne = dms.front().num_entries();
    assert(std::accumulate(dms.begin() + 1, dms.end(), true,
           [ne](bool val, const auto &x) {return val && x.num_entries() == ne;}));
    float *outp = static_cast<float *>(std::malloc(nk * ne * sizeof(float)));
    if(!outp) {
        std::fprintf(stderr, "Allocation of %zu bytes failed\n", size_t(nk * ne * sizeof(float))); return 1;
    }

    static constexpr uint64_t NB = 4096;
    #pragma omp parallel for
    for(size_t i = 0; i < ((NB - 1) + ne) / NB; ++i) {
        auto spos = i * NB, espos = std::min((i + 1) * NB, ne);
        auto destp = outp + spos * nk;//, endp = std::min(destp + nk, outp + ne);
        do {for(auto j = 0u; j < nk;*destp++ = dms[j++][spos]);} while(++spos < espos);
    }
    std::FILE *ofp = fopen(outpath.data(), "wb");
    if(!ofp) return 2;
    uint64_t number_sets = fpaths.size();
    std::fwrite(&ne, sizeof(ne), 1, ofp);
    std::fwrite(&number_sets, sizeof(number_sets), 1, ofp);
    std::fwrite(outp, nk * ne, sizeof(float), ofp);
    std::fclose(ofp);
    std::free(outp);
    return 0;
}

}
