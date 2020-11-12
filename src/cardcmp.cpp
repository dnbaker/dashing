#include "dashing.h"
#include "sketch_and_cmp.h"

using namespace bns;
int main(int argc, char **argv) {
    bns::executable = "cardcmp";
    return card_main(argc, argv);
}
