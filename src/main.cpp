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
    std::fprintf(stderr, "Dashing version: %s\n", DASHING_VERSION);
    if(argc == 1) main_usage(argv);
    if(std::strcmp(argv[1], "sketch") == 0) return sketch_main(argc - 1, argv + 1);
    else if(std::strcmp(argv[1], "dist") == 0 || std::strcmp(argv[1], "cmp") == 0)
        return dist_main(argc - 1, argv + 1);
    else if(std::strcmp(argv[1], "union") == 0) return union_main(argc - 1, argv + 1);
    else if(std::strcmp(argv[1], "setdist") == 0) return setdist_main(argc - 1, argv + 1);
    else if(std::strcmp(argv[1], "hll") == 0) return hll_main(argc - 1, argv + 1);
    else if(std::strcmp(argv[1], "view") == 0) return view_main(argc - 1, argv + 1);
    else if(std::strcmp(argv[1], "mkdist") == 0) return mkdist_main(argc - 1, argv + 1);
#if 0
    else if(std::strcmp(argv[1], "flatten") == 0) return flatten_main(argc - 1, argv + 1);
#endif
    else if(std::strcmp(argv[1], "printmat") == 0) return print_binary_main(argc - 1, argv + 1);
    else if(std::strcmp(argv[1], "sketch_by_seq") == 0) return sketch_by_seq_main(argc - 1, argv + 1);
    else if(std::strcmp(argv[1], "dist_by_seq") == 0) return dist_by_seq_main(argc - 1, argv + 1);
    else {
        for(const char *const *p(argv + 1); *p; ++p) {
            std::string v(*p);
            std::transform(v.begin(), v.end(), v.begin(), [](auto c) {return std::tolower(c);});
            if(v == "-h" || v == "--help") main_usage(argv);
            if(v == "-v" || v == "--version") version_info(argv);
        }
        std::fprintf(stderr, "Usage: %s <subcommand> [options...]. Use %s <subcommand> for more options.\n"
                             "Subcommands:\nsketch\ndist\nhll\nunion\nprintmat\nview\nflatten\n"
                             "sketch_by_seq\ndist_by_seq\n"
                             "\ncmp is also now a synonym for dist, which will be deprecated in the future.\n"
                             "mkdist, which is simply a wrapper around dist which runs the computation for multiple 'k' and organizes them into one result.\n"
                     , *argv, *argv);
        UNRECOVERABLE_ERROR(std::string("Invalid subcommand ") + argv[1] + " provided.");
    }
}
