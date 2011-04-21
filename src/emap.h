#ifndef EMAP_H
#define EMAP_H
#include <assert.h>

#define EMAP_COMMENT_CHARS "#;"
#define EMAP_FLTCMP_DELTA 1e-6

#include <math.h>

#if defined(EMAP_POINT_DOUBLE)
typedef double emap_float;

static inline emap_float emap_float_abs(emap_float n)
{
        return fabs(n);
}
#else
typedef float  emap_float;

static inline emap_float emap_float_abs(emap_float n)
{
        return fabsf(n);
}
#endif

#define EMAP_SUCCESS  0
#define EMAP_EFAULT   1
#define EMAP_EINVAL   2
#define EMAP_ESYSTEM  3
#define EMAP_ENODATA  4
#define EMAP_ENOBUF   5
#define EMAP_ENOMEM   6
#define EMAP_EUNKNOWN 255

#define EMAP_LINEBUFFER_SIZE 65536

#endif /* EMAP_H */
