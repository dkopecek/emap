#ifndef EMAP_H
#define EMAP_H
#include <config.h>

#ifndef EMAP_VERSION
# define EMAP_VERSION PACKAGE_VERSION
#endif

enum emap_transfn {
  EMAP_TRANSFORM_NONE,
  EMAP_TRANSFORM_LOG,
  EMAP_TRANSFORM_LOG2,
  EMAP_TRANSFORM_LOG10
};

#include <assert.h>

#define EMAP_COMMENT_CHARS "#;"
#define EMAP_FLTCMP_DELTA 1e-6
#include <math.h>

#if defined(EMAP_POINT_DOUBLE)
typedef long double emap_float;
#define EMAP_FLTFMT ".10Lf"

#include <stdlib.h>

static inline emap_float emap_float_abs(emap_float n)
{
        return fabsl(n);
}

static inline emap_float emap_strtoflt(const char *str, char **save)
{
        return strtold(str, save);
}

static inline emap_float emap_log(emap_float n)
{
	return logl(n);
}

static inline emap_float emap_log2(emap_float n)
{
	return log2l(n);
}

static inline emap_float emap_log10(emap_float n)
{
	return log10l(n);
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

static inline emap_float emap_log(emap_float n)
{
	return logf(n);
}

static inline emap_float emap_log2(emap_float n)
{
	return log2f(n);
}

static inline emap_float emap_log10(emap_float n)
{
	return log10f(n);
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
