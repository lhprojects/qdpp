/*
 * src/fpu.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Contains functions to set and restore the round-to-double flag in the
 * control word of a x86 FPU.
 */

#include <qd/config.h>
#include <qd/fpu.h>

extern "C" {

void c_fpu_fix_start(unsigned int *old_cw) {
    fpu_fix_end(old_cw);
}

void c_fpu_fix_end(unsigned int *old_cw) {
    fpu_fix_end(old_cw);
}

#ifdef HAVE_FORTRAN

#define f_fpu_fix_start FC_FUNC_(f_fpu_fix_start, F_FPU_FIX_START)
#define f_fpu_fix_end   FC_FUNC_(f_fpu_fix_end,   F_FPU_FIX_END)

void f_fpu_fix_start(unsigned int* old_cw)
{
    fpu_fix_start(old_cw);
}

void f_fpu_fix_end(unsigned int* old_cw)
{
    fpu_fix_end(old_cw);
}

#endif

}
