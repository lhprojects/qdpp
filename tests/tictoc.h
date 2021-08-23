/*
 * tests/timer.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Contains function used for timing.
 */

#ifndef TICTOC_H__
#define TICTOC_H__

#include <chrono>

typedef std::chrono::high_resolution_clock::time_point tictoc;

void   tic(tictoc *tv);   /* start timing. */
// return duration in seconds
double toc(tictoc *tv);   /* stop  timing. */

#endif
