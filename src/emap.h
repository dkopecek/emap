#ifndef EMAP_H
#define EMAP_H
#include <assert.h>

#define EMAP_COMMENT_CHARS "#;"
#define EMAP_FLTCMP_DELTA 1e-6
#include <math.h>

#if defined(EMAP_POINT_DOUBLE)
typedef double emap_float;
#define EMAP_FLTFMT "f"

static inline emap_float emap_float_abs(emap_float n)
{
        return fabs(n);
}

static inline emap_float emap_strtoflt(const char *str, char **save)
{
        return strtod(str, save);
}
#else
typedef float  emap_float;
#define EMAP_FLTFMT "f"

static inline emap_float emap_float_abs(emap_float n)
{
        return fabsf(n);
}

static inline emap_float emap_strtoflt(const char *str, char **save)
{
        return strtof(str, save);
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

#include "pointdb.h"

#ifdef _OPENMP
# include <pthread.h>
#endif

typedef struct {
#ifdef _OPENMP
        pthread_mutex_t mutex;
#endif
        unsigned int    cntinc; /**< count increment for reallocation */
        emap_point_t  **points; /**< all collected points */
        size_t          pcount;
        size_t          palloc;

        emap_point_t  **minima; /**< points representing local minima */
        size_t          mcount;

        emap_point_t  **transp; /**< points representing transition states */
        size_t          tcount;
} POIdb_t;

#define POI_FLAG_MINIMUM    0x00000001
#define POI_FLAG_TRANSITION 0x00000002

#endif /* EMAP_H */
