#include "include/bonsai/encoder.h"
#include "src/substrs.h"
namespace bns {

std::vector<float> nuc_freqs(const std::vector<std::string> &fpaths, unsigned nt) {
    std::vector<float> ret(4 * fpaths.size());
    const size_t np = fpaths.size();
    KseqBufferHolder kseqs(nt);
    OMP_PFOR
    for(size_t i = 0; i < np; ++i) {
        uint64_t arr [128]{0};
        size_t nbp = 0;
        kseq_t *ks =
        for_each_substr([arr,&nbp](const char *s) {
            const int tid = OMP_ELSE(omp_get_thread_num(), 0);
            gzFile fp = gzopen(s, "rb");
            if(!fp) {std::fprintf(stderr, "Missing file %s\n", s); std::exit(1);}
            kseq_t *ks = &kseqs[i];
            kseq_assign(ks, fp);
            while(kseq_read(ks) >= 0) {
                for(size_t j = 0; j < ks->seq.l; ++arr[ks-seq.s[j++]]);
            }
            gzclose(fp);
        }, fpaths[i], FNAME_SEP);
        float cparr[4];
        double total_valid = arr['A'] + arr['a'] + arr['t'] + arr['T'] +  arr['c'] + arr['C'] + arr['g'] + arr['G'];
        cparr[0] = arr['A'] + arr['a'] / total_valid;
        cparr[1] = arr['c'] + arr['C'] / total_valid;
        cparr[2] = arr['g'] + arr['G'] / total_valid;
        cparr[3] = arr['t'] + arr['T'] / total_valid;
        std::memcpy(ret.data() + (4 * i), cparr, sizeof(cparr));
    }
    return nuc_freqs;
}

dm::DistanceMatrix<float> mkmat2jcdistmat(std::string packedmatrix, EmissionType emtype) {
    switch(emtype) {
        case MASH_DIST:
        case FULL_MASH_DIST:
        case FULL_CONTAINMENT_DIST:
        case CONTAINMENT_INDEX:
        case CONTAINMENT_DIST:
            throw std::invalid_argument("EmissionType must be JI or SIZES");
    }
    const bool is_intersection_size = (emtype == SIZES);
    uint64_t number_sets, number_entries;
    gzFile ifp = gzopen(packedmatrix.data(), "rb");
    if(!ifp) RUNTIME_ERROR(std::string("Could not open file at ") + packedmatrix);
    gzread(ifp, &number_sets, sizeof(number_sets));
    gzread(ifp, &number_entries, sizeof(number_entries));
    assert(((number_sets * (number_sets - 1)) >> 1) == number_entries);
    dm::DistanceMatrix<float> ret(number_sets);
    for(size_t i = 0; i < number_sets; ++i) {
        auto rowptr = ret.row_ptr(i);
    }
    gzclose(ifp);
}

} // namespace bns
