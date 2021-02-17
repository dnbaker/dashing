#include "distmat.h"
#include <iostream>
#include <random>

template<typename T>
void test_other(const size_t n, size_t bs) {
    {
        dm::DistanceMatrix<T> mat(n);
        std::gamma_distribution<double> gamrock(137);
        const T val = 0.42;
        dm::parallel_fill(mat, n, [&](uint64_t x, uint64_t y) {return val;}, bs);
        //mat.printf(stdout);
        for(size_t i = 0; i < mat.rows(); ++i) {
            for(size_t j = i + 1; j < mat.rows(); ++j) {
                double v  = mat(i, j);
                //if(v != val) std::fprintf(stderr, "Val %g found at %zu/%zu, expected %g\n", double(v), i, j, double(val));
            }
        }
        for(size_t i = 0; i < mat.num_entries(); ++i) {
            if(mat.data()[i] != 0.42) std::fprintf(stderr, "index %zu/%zu is wrong\n", i, mat.num_entries());
        }
        auto v = std::all_of(mat.data(), mat.data() + mat.num_entries(), [&](auto x) {return x == val;});
        if(!v) std::fprintf(stderr, "Failed with %zu/%zu\n", n, bs);
    }
}

int main() {
    for(const size_t n: {280u, 800u, 10000u}) {
        test_other<double>(n, 4);
        test_other<double>(n, 40);
        test_other<double>(n, 140);
    }
}
