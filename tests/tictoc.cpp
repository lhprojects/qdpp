/*
 * tests/tictoc.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2006
 *
 * Contains function used for timing.
 */

#include "tictoc.h"

void tic(tictoc *tv) {
    *tv = std::chrono::high_resolution_clock::now();
}

double toc(tictoc *tv) {
    typedef std::chrono::duration<double> seconds;

    tictoc tv2 = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::duration dur = tv2 - *tv;
    auto sds = std::chrono::duration_cast<seconds>(dur);
    return sds.count();
}

