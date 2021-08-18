#include "dashing.h"
namespace bns {

template<typename SketchType>
void panel_query(const std::string &path, const std::vector<typename FinalSketch<SketchType>::final_type>, const std::vector<std::string> &labels, gzFile fp) {

}

#undef UNRECOVERABLE_ERROR
#define UNRECOVERABLE_ERROR(x) do {std::cout << x << '\n' << std::flush; return 1;} while(0)

template<typename SketchType>
int panel_core(std::string dbpath, std::string outpath, std::string inpath, bool filesinfile=false) {
    using final_type = typename FinalSketch<SketchType>::final_type;
    std::vector<std::string> labels;
    gzFile fp = gzopen((dbpath + ".labels").data(), "rb");
    if(!fp) UNRECOVERABLE_ERROR("Can't open labels file");
    {
        std::unique_ptr<char[]> buf(new char[10000]);
        for(char *s;(s = gzgets(fp, buf.get(), 10000));labels.push_back(s));
    }
    gzclose(fp);
    const size_t nref = labels.size();
    std::fprintf(stderr, "nref: %zu\n", nref);
    if(nref == 0) UNRECOVERABLE_ERROR("No references found");
    std::vector<final_type> references;
    references.reserve(nref);
    if((fp = gzopen(dbpath.data(), "rb")) == nullptr) UNRECOVERABLE_ERROR("Failed to open DB");
    // Maybe stop after refid = nref?
    size_t refid = 0;
    for(;;) {
        try {
            references.emplace_back(fp);
            std::fprintf(stderr, "reading refid %zu\n", ++refid);
        } catch(...) {
            break;
        }
    }
    gzclose(fp);
    if(references.size() != labels.size()) {
        std::string msg = "Failed to read the correct number of errors. Found: " + std::to_string(references.size()) + ". Expected: " + std::to_string(nref);
        UNRECOVERABLE_ERROR(msg);
    }
    if((fp = gzopen(outpath.data(), "wb")) == nullptr) UNRECOVERABLE_ERROR("Failed to open output file");
    if(filesinfile)
        for(const auto &path: bns::get_paths(inpath.data()))
            panel_query<SketchType>(path, references, labels, fp);
    else panel_query<SketchType>(inpath, references, labels, fp);
    gzclose(fp);
    return EXIT_SUCCESS;
}

int panel_main(int argc, char **argv) {
    bool filesinfile = false;
    std::string outpath = "/dev/stdout", databasepath;
    for(int c;(c = getopt(argc, argv, "p:o:Fh?") >= 0);) {
        switch(c) {
            case 'p': omp_set_num_threads(std::atoi(optarg)); break;
            case 'F': filesinfile = true; break;
            case 'o': outpath = optarg; break;
        }
    }
    if(argc != optind + 2)
        UNRECOVERABLE_ERROR("Failure");
    //databasepath = argv[optind];
    return panel_core<hll::hll_t>(argv[optind], outpath, argv[optind + 1], filesinfile);
}

}
