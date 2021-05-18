#ifndef FOR_EACH_SUBSTR_H__
#define FOR_EACH_SUBSTR_H__
#include <string>
#include <cstring>


template<typename F>
void for_each_substr(const F &func, const std::string &s, const int sep=' ') {
    const char *p;
    if((p = std::strchr(s.data(), sep)) == nullptr) {
        func(s.data());
        return;
    }
    const char *p2 = s.data();
    std::string tmp(p2, p);
    for(;;) {
        func(tmp.data());
        std::swap(p2, ++p);
        if((p = std::strchr(p2, sep)) == nullptr) {
            func(p2);
            break;
        }
        tmp = std::string(p2, p);
        if(std::all_of(tmp.begin(), tmp.end(), [](auto x) {return std::isspace(x);})) break;
    }
}

#endif
