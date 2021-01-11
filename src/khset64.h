#pragma once
#include "dashing.h"
#include "kxsort.h"

namespace bns {
struct khset64_t: public kh::khset64_t {
    // TODO: change to sorted hash sets for faster comparisons, implement parallel merge sort.
    using base = kh::khset64_t;
    using final_type = khset64_t;
    void addh(uint64_t v) {this->insert(v);}
    void add(uint64_t v) {this->insert(v);}
    double cardinality_estimate() const {
        return this->size();
    }
    khset64_t(): kh::khset64_t() {}
    //khset64_t(khset64_t &&o) = default;
    khset64_t(const khset64_t &o) = default;
    khset64_t(size_t reservesz): kh::khset64_t(reservesz) {}
    khset64_t(gzFile fp): kh::khset64_t(fp) {}
    khset64_t(std::string s) {
        this->n_occupied = this->n_buckets = 0;
        this->keys = 0;
        this->flags = 0;
        this->vals = 0;
        this->read(s);
    }
    void reset() {
        if(flags) {
            std::memset(flags, 0xaa, __ac_fsize(n_buckets) * sizeof(uint32_t));
            size_ref() = n_occupied = 0;
        }
    }
    khset64_t(khset64_t &&o) {
        this->n_occupied = o.n_occupied;
        this->n_buckets = o.n_buckets;
        o.n_occupied = o.n_buckets = o.size_ref() = 0;
        this->upper_bound =    o.upper_bound;
        this->keys = o.keys;   o.keys = nullptr;
        this->vals = nullptr;   o.vals = nullptr;
        this->flags = o.flags; o.flags = nullptr;
    }
    void cvt2shs() {
        if(this->flags == nullptr) return;
        uint64_t *newp = reinterpret_cast<uint64_t *>(this->keys);
        static_assert(sizeof(*newp) == sizeof(*this->keys), "same");
        size_t i = 0;
        for(size_t ki = 0; ki != this->n_buckets; ++ki) {
            if(kh_exist(this, ki))
                newp[i++] = kh_key(this, ki);
        }
        assert(i == this->n_occupied);
        std::free(this->flags);
        this->flags = nullptr;
        kx::radix_sort(newp, newp + i);
    }
    void read(const std::string &s) {read(s.data());}
    void read(const char *s) {
        gzFile fp = gzopen(s, "rb");
        this->read(fp);
        gzclose(fp);
    }
    void read(gzFile fp) {
        uint64_t nelem;
        if(gzread(fp, &nelem, sizeof(nelem)) != sizeof(nelem))
            UNRECOVERABLE_ERROR("Failure to read");
        static_assert(sizeof(*this->keys) == sizeof(uint64_t), "must be same");
        if((this->keys = static_cast<khint64_t *>(std::realloc(this->keys, nelem * sizeof(uint64_t)))) == nullptr)
            std::bad_alloc();
        if(gzread(fp, this->keys, nelem * sizeof(uint64_t)) != ssize_t(nelem * sizeof(uint64_t)))
            UNRECOVERABLE_ERROR("Failure to read");
    }
    void free() {
        auto ptr = reinterpret_cast<kh::khash_t(set64) *>(this);
        std::free(ptr->keys);
        std::free(ptr->flags);
        std::memset(ptr, 0, sizeof(*this));
    }
    void write(const std::string &s) const {write(s.data());}
    void write(const char *s) const {
        gzFile fp = gzopen(s, "wb");
        this->write(fp);
        gzclose(fp);
    }
    void write(gzFile fp) const {
        uint64_t nelem = this->n_occupied;
        std::vector<uint64_t> tmp(nelem);
        auto it = tmp.begin();
        for(khiter_t ki = 0; ki != this->n_buckets; ++ki)
            if(kh_exist(this, ki))
                *it++ = kh_key(this, ki);
        kx::radix_sort(tmp.begin(), tmp.end());
        if(gzwrite(fp, &nelem, sizeof(nelem)) != sizeof(nelem)) UNRECOVERABLE_ERROR("Failed to write khash set to disk.");
        if(gzwrite(fp, this->keys, sizeof(*this->keys) * nelem) != ssize_t(sizeof(*this->keys) * nelem))
            UNRECOVERABLE_ERROR("Failed to write khash set to disk.");
    }
    struct Counter
    {
      struct value_type { template<typename T> value_type(const T&) { } };
      void push_back(const value_type &x) { ++count; }
      Counter(): count(0) {}
      size_t count;
    };
    std::array<double, 3> full_set_comparison(const khset64_t &other) const {
        if(flags) UNRECOVERABLE_ERROR("flags must be null by now\n");
        if(other.flags) UNRECOVERABLE_ERROR("flags must be null by now (other)\n");
        Counter c;
        assert(std::is_sorted(this->keys, this->keys + this->n_occupied));
        assert(std::is_sorted(other.keys, other.keys + other.n_occupied));
        assert(c.count == 0);
        std::set_intersection(this->keys, this->keys + this->n_occupied,
                              other.keys, other.keys + other.n_occupied,
                              std::back_inserter(c));
        double is = c.count;
        return std::array<double, 3>{{this->n_occupied - is, other.n_occupied - is, is}};
    }
    size_t intersection_size(const khset64_t &other) const {
        return full_set_comparison(other)[2];
    }
    size_t isz(const khset64_t &o) const {return intersection_size(o);}
    double jaccard_index(const khset64_t &other) const {
        auto cmps = full_set_comparison(other);
        return double(cmps[2]) / (cmps[0] + cmps[1] + cmps[2]);
    }
    double containment_index(const khset64_t &other) const {
        auto cmps = full_set_comparison(other);
        return double(cmps[2]) / (cmps[0] + cmps[2]);
    }
    double symmetric_containment_index(const khset64_t &other) const {
        auto cmps = full_set_comparison(other);
        return double(cmps[2]) / (std::min(cmps[0], cmps[1]) + 1e-20 + cmps[2]);
    }
    khset64_t &operator=(const khset64_t &o) {
        static_cast<kh::khset64_t *>(this)->operator=(static_cast<kh::khset64_t>(o));
        return *this;
    }
    khset64_t &operator+=(const khset64_t &o) {
        throw NotImplementedError("This shouldn't be called.");
    }
    uint64_t union_size(const khset64_t &other) const {
        auto cmps = full_set_comparison(other);
        return cmps[0] + cmps[1] + cmps[2];
    }
};


} // bns
