#include "dashing.h"


namespace bns {
int flatten_all(const std::vector<std::string> &fpaths, const std::string outpath, std::vector<unsigned> &k_values) {
    if(fpaths.empty()) RUNTIME_ERROR("no fpaths, see usage.");
    const size_t nk = k_values.size();
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
    uint32_t numk = k_values.size();
#if 0
    gzread(ifp, &nk, sizeof(nk));
    std::vector<unsigned> k_values(nk);
    gzread(ifp, &number_entries, sizeof(number_entries));
    gzread(ifp, &number_sets, sizeof(number_sets));
    gzread(ifp, k_values.data(), k_values.size() * sizeof(unsigned));
#endif
    if(std::fwrite(&numk, sizeof(numk), 1, ofp) != 1) throw std::runtime_error("Failed to write nk");
    if(std::fwrite(&ne, sizeof(ne), 1, ofp) != 1) throw std::runtime_error("Failed to write num entries");
    std::fwrite(&number_sets, sizeof(number_sets), 1, ofp);
    std::fprintf(stderr, "Wrote %u nk and %zu num e and %zu num sets\n", numk, ne, size_t(number_sets));
    if(std::fwrite(k_values.data(), sizeof(unsigned), k_values.size(), ofp) != k_values.size()) throw std::runtime_error("Wrong count of written values");
    size_t nwritten;
    if((nwritten = std::fwrite(outp, sizeof(float), nk * ne, ofp)) != nk * ne) {
        throw std::runtime_error("Failed to write nk * ne");
    }
    std::fprintf(stderr, "Wrote %zu items (%zu * %zu) [%zu bytes]\n", nk * ne, nk, ne, nk * ne * 4);
    std::fclose(ofp);
    std::free(outp);
    return 0;
}

}
