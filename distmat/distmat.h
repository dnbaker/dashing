#pragma once
#include <array>
#include <cassert>
#include <cinttypes>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <system_error>
#include <type_traits>
#include <vector>
#if ZWRAP_USE_ZSTD
#  include "zstd_zlibwrapper.h"
#else
#  include <zlib.h>
#endif
#include "unistd.h"
#include "./mio.hpp"

#ifndef INLINE
#  if defined(__GNUC__) || defined(__clang__)
#    define INLINE __attribute__((always_inline)) inline
#  elif defined(__CUDACC__)
#    define INLINE __forceinline__ inline
#  else
#    define INLINE inline
#  endif
#endif

namespace std {
#if __cplusplus < 201703L
template <class C>
constexpr auto size(const C& c) -> decltype(c.size())
{
    return c.size();
}

template <class T, std::size_t N>
constexpr std::size_t size(const T (&array)[N]) noexcept
{
    return N;
}
#endif
} // std


namespace dm {
using std::fputc;

template<typename T> inline std::string to_string(T x) {return std::to_string(x);}
template<> inline std::string to_string<__uint128_t>(__uint128_t num)
{
    std::string str;
    do {
        int digit = num % 10;
        num /= 10;
        str = std::to_string(digit) + str;
    } while(num);
    return str;
}
template<> inline std::string to_string<__int128_t>(__int128_t n) {
    std::string str;
    bool signbit;
    if(n < 0) signbit = 1, n = -n;
    else signbit = 0;
    do {
        int digit = n % 10;
        n /= 10;
        str = std::to_string(digit) + str;
    } while(n);
    if(signbit) str = std::string("-") + str;
    return str;
}
template<typename T> struct numeric_limits: public std::numeric_limits<T> {};

template<> struct numeric_limits<__uint128_t> {
    static constexpr __uint128_t max() {return __uint128_t(-1);}
    static constexpr __uint128_t min() {return __uint128_t(0);}
};
template<> struct numeric_limits<__int128_t> {
    static constexpr __int128_t max() {return  (__int128) (((unsigned __int128) 1 << ((__SIZEOF_INT128__ * __CHAR_BIT__) - 1)) - 1);}
    static constexpr __int128_t min() {return (-max() - 1);}
};

namespace more_magic {
template<typename ArithType>
struct MAGIC_NUMBER {
    const char *name() {
        throw std::runtime_error("NotImplemented");
        return "NOTIMPLEMENTED";
    }
};
enum MagicNumber: uint8_t {
    FLOAT,
    DOUBLE,
    UINT8_T,
    UINT16_T,
    UINT32_T,
    UINT64_T,
    UINT128_T,
    INT8_T,
    INT16_T,
    INT32_T,
    INT64_T,
    INT128_T
};
static constexpr const char *arr[] {
    "float",
    "double",
    "uint8_t",
    "uint16_t",
    "uint32_t",
    "uint64_t",
    "uint128_t",
    "int8_t",
    "int16_t",
    "int32_t",
    "int64_t",
    "int128_t",
};

#define DEC_MAGIC(type, STR, num) \
    template<> struct MAGIC_NUMBER<type> {\
        static constexpr MagicNumber magic_number = num;\
        static constexpr const char *name() {return "DM::" STR;}\
        static const std::string magic_name() {return std::string("DM::") + arr[magic_number];}\
    }

DEC_MAGIC(float, "float", FLOAT);
DEC_MAGIC(double,"double", DOUBLE);
DEC_MAGIC(uint8_t,"uint8_t", UINT8_T);
DEC_MAGIC(uint16_t,"uint16_t", UINT16_T);
DEC_MAGIC(uint32_t,"uint32_t", UINT32_T);
DEC_MAGIC(uint64_t,"uint64_t", UINT64_T);
DEC_MAGIC(__uint128_t,"uint128_t", UINT128_T);
DEC_MAGIC(int8_t,"int8_t", INT8_T);
DEC_MAGIC(int16_t,"int16_t", INT16_T);
DEC_MAGIC(int32_t,"int32_t", INT32_T);
DEC_MAGIC(int64_t,"int64_t", INT64_T);
DEC_MAGIC(__int128_t,"int128_t", INT128_T);

} // namespace more_magic

#undef DEC_MAGIC


/* *
 * DistanceMatrix holds an upper-triangular matrix.
 * You can access rows with row_span()
 * or individual entries with (i, j) notation (like Eigen/Blaze/&c.)
 * You can set the default value with set_default_value.
 *
*/
template<typename ArithType=float,
         size_t DefaultValue=0,
         bool force=false,
         typename=typename std::enable_if<std::is_arithmetic<ArithType>::value || std::is_same<ArithType,__uint128_t>::value || std::is_same<ArithType,__int128_t>::value || force>::type
         >
class DistanceMatrix {
    ArithType *data_;
    std::unique_ptr<ArithType> dup_;
    uint64_t  nelem_, num_entries_;
    ArithType default_value_;
    std::unique_ptr<mio::mmap_sink> mfbp_;

public:
    static constexpr const char *magic_string() {return more_magic::MAGIC_NUMBER<ArithType>::name();}
    static constexpr more_magic::MagicNumber magic_number() {return more_magic::MAGIC_NUMBER<ArithType>::magic_number;}
    using value_type = ArithType;
    using pointer_type = ArithType *;
    using const_pointer_type = const ArithType *;
    static constexpr ArithType DEFAULT_VALUE = static_cast<ArithType>(DefaultValue);
    void set_default_value(ArithType val) {default_value_ = val;}
    DistanceMatrix(size_t n, ArithType default_value=DEFAULT_VALUE, ArithType *prevdat=static_cast<ArithType *>(nullptr)):
        nelem_(n), num_entries_((nelem_ * (nelem_ - 1)) >> 1), default_value_(default_value)
    {
        if(prevdat) data_ = prevdat;
        else {
            data_ = new ArithType[num_entries_];
            dup_.reset(data_);
        }
    }
    DistanceMatrix(): DistanceMatrix(size_t(0), DEFAULT_VALUE) {}
    pointer_type       data()       {return data_;}
    const_pointer_type data() const {return data_;}
    DistanceMatrix(DistanceMatrix &&other) = default;
    DistanceMatrix(const char *path, size_t nelem=0, ArithType default_value=DEFAULT_VALUE, ArithType *prevdat=nullptr, bool forcestream=false):
        nelem_(nelem), default_value_(default_value)
    {
        std::fprintf(stderr, "path: %s. nelem: %zu. defv: %g\n", path, nelem_, double(default_value_));
        if(!forcestream && ::access(path, F_OK) == -1 && nelem > 0) {
            num_entries_ = (nelem_ * (nelem_ - 1)) >> 1;
            std::fprintf(stderr, "File does not exist, so create it. nelem = %zu, num entries = %zu\n", nelem_, num_entries_);
            // If file does not exist,
            // open a new file on disk and resize it.
            std::FILE *ofp = std::fopen(path, "wb");
            const size_t nb = 1 + sizeof(num_entries_) + sizeof(ArithType) * num_entries_;
            std::fputc(magic_number(), ofp);
            if(std::fwrite(&nelem_, sizeof(nelem_), 1, ofp) != 1) throw std::runtime_error("Failed to write nelem to disk");
            // Resize
            //std::fprintf(stderr, "Resized to %zu bytes for %zu entries\n", nb, num_entries_);
            ::ftruncate(::fileno(ofp), nb);
            std::fclose(ofp);
            mfbp_.reset(new mio::mmap_sink(path));
            data_ = reinterpret_cast<ArithType *>((*mfbp_).data() + 1 + sizeof(nelem_));
            std::fprintf(stderr, "Created new file of %zu bytes, %zu nelem, and %zu entries\n", nb, nelem_, num_entries_);
        } else read(path, prevdat, forcestream);
    }
    auto nelem() const {return nelem_;}
    DistanceMatrix(const DistanceMatrix &other, ArithType *prevdat=static_cast<ArithType *>(nullptr)):
            nelem_(other.nelem_),
            num_entries_(other.num_entries_),
            default_value_(other.default_value_)
    {
        if(prevdat) data_ = prevdat;
        else        data_ = new ArithType[num_entries_], dup_.reset(data_);
        std::memcpy(data_, other.data_, num_entries_ * sizeof(value_type));
    }
    size_t num_entries() const {return num_entries_;}
#define ARRAY_ACCESS(row, column) (((row) * (nelem_ * 2 - row - 1)) / 2 + column - (row + 1))
    INLINE size_t index(size_t row, size_t column) const {
        return row < column ? ARRAY_ACCESS(row, column): ARRAY_ACCESS(column, row);
    }
#undef ARRAY_ACCESS
    INLINE value_type &operator()(size_t row, size_t column) {
        if(__builtin_expect(row == column, 0)) return default_value_;
        return data_[index(row, column)];
    }
    INLINE const value_type &operator()(size_t row, size_t column) const {
        if(__builtin_expect(row == column, 0)) return default_value_;
        return data_[index(row, column)];
    }
    pointer_type row_ptr(size_t row) {
        auto ret = data_ + nelem_ * row - (row * (row + 1) / 2);
        return ret;
    }
    const_pointer_type row_ptr(size_t row) const {
        auto ret = data_ + nelem_ * row - (row * (row + 1) / 2);
        return ret;
    }
    std::pair<pointer_type, size_t> row_span(size_t i) {
        return std::make_pair(row_ptr(i), nelem_ - i - 1);
    }
    std::pair<const_pointer_type, size_t> row_span(size_t i) const {
        return std::make_pair(row_ptr(i), nelem_ - i - 1);
    }
    value_type &operator[](size_t index) {
       return data_[index];
    }
    const value_type &operator[](size_t index) const {
       return data_[index];
    }
    value_type &operator[](std::pair<size_t, size_t> &idx) {
       return data_[index(idx.first, idx.second)];
    }
    const value_type &operator[](std::pair<size_t, size_t> &idx) const {
       return data_[index(idx.first, idx.second)];
    }
    void resize(size_t new_size) {
        if(new_size == nelem_) return; // Already done! Aren't we fast?
        if(new_size < nelem_) throw std::runtime_error("NotImplemented: shrinking.");
        const auto nsz = (new_size * (new_size - 1)) >> 1;
        if(!dup_) throw std::runtime_error("Can't resize external data"); // Can't resize externl
        nelem_ = new_size;
        num_entries_ = nsz;
        dup_.reset(new ArithType[num_entries_]);
        std::fill_n(data_, num_entries_, static_cast<value_type>(-1));
        // Invalid -- to ensure data is calculated before use
    }
    auto begin() {return data_;}
    auto end()   {return data_ + num_entries_;}
    auto begin() const {return data_;}
    auto end()   const {return data_ + num_entries_;}
    void write(const std::string &path) const {
        this->write(path.data());
    }
    std::string to_string(bool use_scientific=false, const std::vector<std::string> *labels=nullptr) const {
        std::string ret;
        ret.reserve(size() * size() * 6 + (labels ? size_t(10 * labels->size()): size_t(0)));
        if(labels) {
            ret = "##Labels";
            for(const auto &l: *labels)
                ret += '\t', ret += l;
            ret += '\n';
        }
        for(size_t i = 0; i < size(); ++i) {
            if(labels)
                ret += labels->operator[](i), ret += '\t';
            for(size_t j = 0; j < size(); ++j) {
                ret += to_string(this->operator()(i, j)), ret += '\t';
            }
            ret.back() = '\n'; // Extra tab is now a newline
        }
        return ret;
    }
    void printf(gzFile fp, bool use_scientific=false, const std::vector<std::string> *labels=nullptr) const {
        if(labels) {
            gzprintf(fp, "#Names");
            for(const auto &s: *labels) {
                gzputc(fp, '\t');
                gzwrite(fp, s.data(), s.size());
            }
            gzputc(fp, '\n');
        }
        std::array<std::array<char, 5>, 2> fmts{{{'%','l','f','\n','\0'}, {'%','l','f','\t','\0'}}};
        if(use_scientific) fmts[0][2] = fmts[1][2] = 'e';
        for(unsigned i(0); i < nelem_; ++i) {
            if(labels) {
                const auto &label = (*labels)[i];
                gzwrite(fp, label.data(), label.size());
                gzputc(fp, '\t');
            }
            for(unsigned j(0); j < nelem_; ++j) {
                gzprintf(fp, fmts[j != nelem_ - 1].data(), static_cast<double>(this->operator()(i, j)));
            }
        }
    }
    void printf(std::FILE *fp, bool use_scientific=false, const std::vector<std::string> *labels=nullptr) const {
        if(labels) {
            fprintf(fp, "#Names");
            for(const auto &s: *labels) {
                fputc('\t', fp);
                fwrite(s.data(), s.size(), 1, fp);
            }
            fputc('\n', fp);
        }
        std::array<std::array<char, 5>, 2> fmts;
        fmts[0] = {'%','l','f','\n','\0'};
        fmts[1] = {'%','l','f','\t','\0'};
        if(use_scientific) fmts[0][2] = fmts[1][2] = 'e';
        for(unsigned i(0); i < nelem_; ++i) {
            if(labels) {
                const auto &label = labels->operator[](i);
                fwrite(label.data(), label.size(), 1, fp);
                fputc('\t', fp);
            }
            for(unsigned j(0); j < nelem_; ++j) {
                fprintf(fp, fmts[j != nelem_ - 1].data(), static_cast<double>(this->operator()(i, j)));
            }
        }
    }
    size_t write(const char *path, int compression_level=0) const {
        std::string fmt = compression_level ? std::string("wT"): (std::string("wb") + std::to_string(compression_level % 10));
        gzFile fp = gzopen(std::strcmp(path, "-") ? path: "/dev/stdout", fmt.data());
        if(!fp) throw std::runtime_error(std::string("Could not open file at ") + path);
        size_t ret = write(fp);
        gzclose(fp);
        return ret;
    }
    size_t write(gzFile fp) const {
        size_t ret = gzputc(fp, magic_number());
        ret += gzwrite(fp, &nelem_, sizeof(nelem_));
        ret += gzwrite(fp, data_, sizeof(ArithType) * num_entries_);
        return ret;
    }
    size_t write(std::FILE *fp) const {
        if(__builtin_expect(fputc(magic_number(), fp) != magic_number(), 0))
            throw std::system_error(std::ferror(fp), std::system_category(), "Failed to write magic number to file");
        std::fflush(fp);
        int fn = fileno(fp);
        size_t ret = 1;
        if(::write(fn, &nelem_, sizeof(nelem_)) != sizeof(nelem_)) {
            throw std::system_error(errno, std::system_category(), ::strerror(errno));
        }
        ret += sizeof(nelem_);
        const ssize_t nb =  sizeof(ArithType) * num_entries_;
        if(::write(fn, data_, nb) != nb) {
            throw std::system_error(errno, std::system_category(), ::strerror(errno));
        }
        ret += sizeof(ArithType) * num_entries_;
        return ret;
    }
    void read(const char *path, ArithType *prevdat=static_cast<ArithType *>(nullptr), bool forcestream=false) {
        // Else, open from file on disk
        //std::fprintf(stderr, "Starting read with path = %s, fs = %d, and nelem = %zu\n", path, forcestream, nelem_);
        using more_magic::MagicNumber;
        path = std::strcmp(path, "-") ? path: "/dev/stdin";
        std::FILE *fp = std::fopen(path, "r");
        if(fp == nullptr) throw std::runtime_error(std::string("Could not open file at ") + path);
        int fd = ::fileno(fp);
        const bool isstream = ::isatty(fd) || forcestream;
        const auto fc = std::fgetc(fp);
        std::ungetc(fc, fp);
        std::fclose(fp);
        fp = std::fopen(path, "r");
        fd = ::fileno(fp);
        gzFile gzfp = gzdopen(fd, "r");
        const MagicNumber magic = MagicNumber(gzgetc(gzfp));
        if(magic != magic_number()) {
            char buf[256];
            std::sprintf(buf, "Wrong magic number read from file (%d/%s), expected (%d/%s)\n", magic, more_magic::arr[magic], magic_number(), magic_string());
            throw std::runtime_error(buf);
        }
        gzread(gzfp, &nelem_, sizeof(nelem_));
        num_entries_ = ((nelem_ - 1) * nelem_) >> 1;
        // If this is streaming, or the file is compressed,
        // copy the memory out
        if(prevdat) data_ = prevdat;
        else {
            dup_.reset(new ArithType[num_entries_]);
            data_ = dup_.get();
        }
        gzread(gzfp, data_, sizeof(ArithType) * num_entries_);
        gzclose(gzfp);
        std::fclose(fp);
    }
    size_t size() const {return nelem_;}
    size_t rows() const {return nelem_;}
    size_t columns() const {return nelem_;}
#if 0
    DistanceMatrix &operator=(DistanceMatrix &&o) = default;
    DistanceMatrix &operator=(const DistanceMatrix &o) {
        if(!data_) throw std::runtime_error("data is null");
        std::memcpy(data_, o.data_, num_entries_ * sizeof(ArithType));
    }
#endif
    bool operator==(const DistanceMatrix &o) const {
        //std::fprintf(stderr, "data_: %p. o.data: %p. nelems: %zu, %zu\n", data_, o.data_, nelem_, o.nelem_);
        return nelem_ == o.nelem_ &&
            (data_ && o.data_ ? (std::memcmp(data_, o.data_, num_entries_ * sizeof(ArithType)) == 0)
                              : data_ == o.data_);
    }
};

template<typename T>
struct is_distance_matrix: public std::false_type {};
template<typename ArithType,
         size_t DefaultValue,
         bool force>
struct is_distance_matrix<DistanceMatrix<ArithType, DefaultValue, force>>:
    public std::true_type {};

#if __cplusplus >= 201703L
template<typename T>
constexpr bool is_distance_matrix_v = is_distance_matrix<T>::value;
#endif


template<typename ArithType=float,
         size_t DefaultValue=0,
         bool force=false>
inline std::ostream &operator<<(std::ostream &os, const DistanceMatrix<ArithType, DefaultValue, force> &m) {
    const size_t nr = m.size();
    for(size_t i = 0; i < nr; ++i) {
        auto rowspan = m.row_span(i);
        os << "Row " << i << "/" << nr << ":( ";
        for(size_t j = 0; j < rowspan.second; ++j) {
            os << rowspan.first[j] << '\t';
        }
        os << ")\n";
    }
    return os;
}



} // namespace dm
