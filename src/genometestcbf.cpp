#include <zlib.h>
#include "bonsai/hll/cbf.h"
#include "bonsai/bonsai/include/encoder.h"
using namespace bns;
void usage() {
    std::exit(1);
}

int main(int argc, char *argv[]) {
    if(argc == 1) {
        std::fprintf(stderr, "TODO: this\n");
        usage();
    }
    khash_t(64) *counts = kh_init(64);
    gzFile ofp = gzopen(argc > 2 ? argv[2]: "/dev/stdout", "wb7");
    Encoder enc(31);
    unsigned nh = 1;
    unsigned bfsize = 22;
    unsigned nbfs = 10;
    bf::cbf_t cbf(nbfs, bfsize, nh, 137);
    enc.for_each([&](const u64 val){
        khiter_t ki;
        ki = kh_get(64, counts, val);
        if(ki == kh_end(counts)) {
            int khr;
            ki = kh_put(64, counts, val, &khr);
            kh_val(counts, ki) = 1;
        } else {
            assert(kh_key(counts, ki) == val);
            ++kh_val(counts, ki);
        }
        cbf.addh(val);
    }, argv[1]);
    ks::string str(256);
    for(khiter_t ki(0); ki < kh_end(counts); ++ki) {
        if(!kh_exist(counts, ki)) continue;
        str.sprintf("%s\t%" PRIu64 "\t%u\n", enc.sp_.to_string(kh_key(counts, ki)).data(), kh_val(counts, ki), cbf.est_count(kh_key(counts, ki)));
        if(str.size() > (1<<16)) str.write(ofp), str.clear();
    }
    str.write(ofp), str.clear();
    kh_destroy(64, counts);
    gzclose(ofp);
}
