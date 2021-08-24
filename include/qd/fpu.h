/*
 * include/fpu.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2001
 *
 * Contains functions to set and restore the round-to-double flag in the
 * control word of a x86 FPU.  The algorithms in the double-double and
 * quad-double package does not function with the extended mode found in
 * these FPU.
 */
#ifndef _QD_FPU_H
#define _QD_FPU_H

 /*
  * Set the round-to-double flag, and save the old control word in old_cw.
  * If old_cw is NULL, the old control word is not saved.
  */
void fpu_fix_start(unsigned int* old_cw);

/*
 * Restore the control word.
 */
void fpu_fix_end(unsigned int* old_cw);


#include "qd_config.h"

#ifdef X86
#ifdef  _WIN32
#include <float.h>
#else

#ifdef HAVE_FPU_CONTROL_H
#include <fpu_control.h>
#endif

#ifndef _FPU_GETCW
#define _FPU_GETCW(x) asm volatile ("fnstcw %0":"=m" (x));
#endif

#ifndef _FPU_SETCW
#define _FPU_SETCW(x) asm volatile ("fldcw %0": :"m" (x));
#endif

#ifndef _FPU_EXTENDED
#define _FPU_EXTENDED 0x0300
#endif

#ifndef _FPU_DOUBLE
#define _FPU_DOUBLE 0x0200
#endif

#endif
#endif /* X86 */


inline void fpu_fix_start(unsigned int* old_cw)
{
#ifdef X86
#ifdef _WIN32
#ifndef _M_X64 

    /* Win 32 MSVC */
    unsigned int cw = _control87(0, 0);
    _control87(0x00010000, 0x00030000);
    if (old_cw) {
        *old_cw = cw;
    }

#endif
#else
    /* Linux */
    volatile unsigned short cw, new_cw;
    _FPU_GETCW(cw);

    new_cw = (cw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
    _FPU_SETCW(new_cw);

    if (old_cw) {
        *old_cw = cw;
    }
#endif
#endif
}

inline void fpu_fix_end(unsigned int* old_cw)
{
#ifdef X86
#ifdef _WIN32
#ifndef _M_X64 

    /* Win 32 MSVC */
    if (old_cw) {
        _control87(*old_cw, 0xFFFFFFFF);
    }

#endif // #ifndef _M_X64 

#else
    /* Linux */
    if (old_cw) {
        int cw;
        cw = *old_cw;
        _FPU_SETCW(cw);
    }
#endif
#endif
}

#endif  /* _QD_FPU_H */
