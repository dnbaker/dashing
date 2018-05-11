#include "bonsai/hll/cbf.h"
#include "bonsai/hll/aesctr.h"

int main() {
    bf::cbf_t cbf(16, 10, 4, 1337);
    aes::AesCtr<uint64_t, 8> gen;
}
