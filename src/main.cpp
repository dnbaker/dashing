#include "dashing.h"
#include "sketch_and_cmp.h"
using namespace bns;
#if HAS_AVX_512
#  pragma message("Building with AVX512 support")
#elif __AVX2__
#  pragma message("Building with AVX2 support")
#elif __SSE4_1__
#  pragma message("Building with SSE4.1 support")
#else
#  pragma message("Building with no vectorization support [this will likely fail]")
#endif


void version_info(char *argv[]) {
    std::fprintf(stderr, "Dashing version: %s\n", DASHING_VERSION);
    std::exit(1);
}

int main(int argc, char *argv[]) {
    bns::executable = argv[0];
    const std::unordered_map<std::string, int (*)(int, char**)> submap
{
    {"sketch", sketch_main},
    {"union", union_main},
    {"setdist", setdist_main},
    {"dist", dist_main},
    {"cmp", dist_main},
    {"hll", hll_main},
    {"view", view_main},
    {"panel", panel_main},
<<<<<<< HEAD
    {"mkdist", mkdist_main},
=======
    //{"mkdist", mkdist_main},
>>>>>>> 522d52a7b0bf187d2329a9ae8b953be81adb23dd
    {"card", card_main},
    {"printmat", print_binary_main},
    {"dist_by_seq", dist_by_seq_main},
    {"sketch_by_seq", sketch_by_seq_main},
};
    
    std::fprintf(stderr, "Dashing version: %s\n", DASHING_VERSION);
    if(argc == 1) main_usage(argv);
    auto it = submap.find(argv[1]);
    if(it != submap.end()) return it->second(argc - 1, argv + 1);
    else {
        for(const char *const *p(argv + 1); *p; ++p) {
            std::string v(*p);
            std::transform(v.begin(), v.end(), v.begin(), [](auto c) {return std::tolower(c);});
            if(v == "-h" || v == "--help") main_usage(argv);
            if(v == "-v" || v == "--version") version_info(argv);
        }
        std::fprintf(stderr, "Usage: %s <subcommand> [options...]. Use %s <subcommand> for more options.\n"
                             "Subcommands:\tsketch\tdist\thll\tunion\tprintmat\tview\tflatten\t"
                             "sketch_by_seq\tdist_by_seq\tpanel"
                             "\ncmp is also now a synonym for dist, which will be deprecated in the future.\n"
                             //"mkdist, which is simply a wrapper around dist which runs the computation for multiple 'k' and organizes them into one result.\n"
                             "\npanel is experimental; please don't use it.\n"
                     , *argv, *argv);
        UNRECOVERABLE_ERROR(std::string("Invalid subcommand ") + argv[1] + " provided.");
    }
}
