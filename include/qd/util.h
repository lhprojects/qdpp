#ifndef QD_UTILS_H
#define QD_UTILS_H
#include <string>

#include <iostream>
#include <stdint.h>
#include <stddef.h>
#include <limits.h>

namespace qd {

    inline void error(const char* msg);
    /* This routine is called whenever a fatal error occurs. */
    inline void error(const char* msg)
    {
        if (msg) { std::cerr << "ERROR " << msg << std::endl; }
    }

}
void append_expn(std::string &str, int expn);


// report error by return -1
// return 0 if everything is OK
constexpr inline int read_int(char const* s, int& e)
{
    int v = 0;
    int sign = 0;
    if (*s == '-') {
        sign = -1;
        s += 1;
    } else if (*s == '+') {
        sign = '+';
        s += 1;
    }

    for (; *s; ++s) {
        char ch = *s;
        if (ch >= '0' && ch <= '9') {
            if (v < INT_MIN / 10) {
                // overflow
                return -1;
            }
            v *= 10;
            int d = *s - '0';
            v += d;
            if (v < 0) {
                // overflow
                break;
            }
        }
    }

    if (sign == -1) {
        if (v == INT_MIN) {
            e = INT_MIN;
        } else if (v < 0) {
            // true overflow
            return -1;
        } else {
            e = -v;
        }
    } else {
        e = v;
    }
    return 0;

}

// report error by return INT_MIN
constexpr inline int to_int(char const* s)
{
    int a;
    if (read_int(s, a) < 0) {
        return INT_MIN;
    }
    return a;
}
#include "util.inl.h"
#endif
