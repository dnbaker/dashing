#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <getopt.h>
#include "distmat/distmat.h"
#include "bonsai/kspp/ks.h"

#define FOREVER for(;;)

using namespace std::literals;

dm::DistanceMatrix<float, 1> parse_flashdans_matrix(const std::string &path) {
    std::ifstream ifs(path);
    std::string line;
    if(!std::getline(ifs, line)) throw std::runtime_error("Could not read a first line from the file at "s + path);
    std::vector<size_t> offsets;
    ks::split(line.data(), '\t', offsets);
    dm::DistanceMatrix<float, 1> ret(offsets.size() - 1);
#if !NDEBUG
    std::fprintf(stderr, "Size: %zu\n", offsets.size());
#endif
    size_t linenum = 0;
    while(std::getline(ifs, line)) {
        offsets.clear();
        ks::split(line.data(), '\t', offsets);
        for(size_t i(1); i < offsets.size(); ++i)
            if(line[offsets[i]] != '-')
                ret(linenum, i - 1) = std::atof(line.data() + offsets[i]);
        ++linenum;
        //if(linenum % 100 == 0) std::fprintf(stderr, "line number: %zu. Lines left to go: %zu\n", linenum, offsets.size() - 1);
    }
    return ret;
}

std::vector<float> linspace(size_t ndiv) {
    std::vector<float> ret; ret.reserve(ndiv + 1);
    for(size_t i(0); i < ndiv; ret.emplace_back(static_cast<float>(i++) / ndiv));
    return ret;
}

void usage(const char *ex) {
    std::fprintf(stderr, "%s -n <ndiv [default: 100]> -S <subset size [default: 40]> distmat.txt\n", ex);
    std::exit(EXIT_FAILURE);
}

int main(int argc, char *argv[]) {
    int c;
    unsigned ndiv = 100, npairs = 40;
    if(argc == 1) goto usage;
    while((c = getopt(argc, argv, "n:s:h?")) >= 0) {
        switch(c) {
            case 'n': ndiv = std::atoi(optarg); break;
            case 's': npairs = std::atoi(optarg); break;
            case 'h': case '?': usage: usage(*argv);
        }
    }
    ndiv = std::min(ndiv, static_cast<unsigned>(mat.size()));
    auto lins = linspace(ndiv);
    std::vector<size_t> pair_counts(ndiv, static_cast<size_t>(0));
    size_t i(0);
    std::fprintf(stderr, "Starting loop %zu\n", size_t(ndiv));
    while(std::accumulate(std::cbegin(pair_counts), std::cend(pair_counts), static_cast<size_t>(0), [](auto a, auto b) {return a + (b != 0);}) < ndiv && i < mat.size()) {
        for(size_t j(0); j < i; ++pair_counts[static_cast<size_t>(ndiv * mat(i, j++))]);
        ++i;
    }
    std::fprintf(stderr, "[W:%s:%s:%d] Currently not using subset_size (value: %u). i: %zu\n", __PRETTY_FUNCTION__, __FILE__, __LINE__, subset_size, i);
    auto tinymat = dm::DistanceMatrix<float, 1>(i);
    for(size_t j(0); j < i; ++j) {
        for(size_t k(0); k < i; ++k) {
            tinymat(k, j) = mat(k, j);
        }
    }
    tinymat.printf(stdout);
}
