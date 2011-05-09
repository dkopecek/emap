#ifndef HELPERS_H
#define HELPERS_H

/**
 * Memory allocation helpers
 */
#include <stdlib.h>

/**
 * Allocate storage for a variable of type `T'.
 */
#define alloc_type(T)                           \
        ((T *) malloc(sizeof(T)))

/**
 * Allocate an array of variables of type `T'.
 */
#define alloc_array(T, count)                   \
        ((T *) malloc(sizeof(T) * (count)))

#define realloc_array(ptr, T, count)            \
        ((T *) realloc(ptr, sizeof(T) * (count)))

#define __lsym(s) s ## __LINE__

/**
 * Various useful functions, macros, ...
 */
#include <errno.h>

/**
 * Save errno, execute the block, restore errno.
 */
#define errno_protect \
        for(int __lsym(errno) = errno, __lsym(exec) = 1; __lsym(exec) == 1; __lsym(exec) = 0, errno = __lsym(errno))

#ifdef __GNUC__
# define __predict(expr, v) __builtin_expect (expr, v)
#else
# define __predict(expr, v) expr
#endif /* __GNUC__ */

/**
 * Output to stderr in debug mode, otherwise do nothing.
 */
#ifndef NDEBUG
# define dP(fmt, ...) fprintf(stderr, fmt, __VA_ARGS__)
#else
# define dP(fmt, ...) while(0)
#endif

#endif /* HELPERS_H */
