#ifndef POINTDB_H
#define POINTDB_H

#include <stdint.h>
#include "emap.h"

typedef struct {
        emap_float y;   /**< dependent variable */
        emap_float x[]; /**< inline array of independent variables */
} emap_point_t;

#define EMAP_POINT_SIZE(arity) (sizeof(emap_point_t) + (sizeof(emap_float) * (arity)))

typedef struct {
        uint16_t       flags; /**< flags */
        uint16_t       arity; /**< number of independent variables */
        size_t         count; /**< number of points in the array */
        emap_point_t  *point; /**< array of points */
        emap_point_t **psort; /**< array of sorted points by each independent variable */
} emap_pointdb_t;

static inline emap_point_t *emap_pointp(emap_pointdb_t *pdb, size_t pi)
{
        return (emap_point_t *)(((uint8_t *)(pdb->point)) + (EMAP_POINT_SIZE(pdb->arity) * pi));
}

static inline emap_point_t *emap_psortp(emap_pointdb_t *pdb, size_t xn, size_t xi)
{
        return (NULL);
}

#define EMAP_PDBF_INIT 0x0001 /**< set if properly initialized */
#define EMAP_PDBF_FREE 0x0002 /**< set if memory for the struct was dynamically allocated */
#define EMAP_PDBF_LOAD 0x0004 /**< set if some data points were loaded */
#define EMAP_PDBF_SORT 0x0008 /**< set if the data points are sorted */

#define EMAP_PDB_INITIALIZED(p)                         \
        do {                                            \
                if (!((p)->flags & EMAP_PDBF_INIT))     \
                        abort();                        \
        } while(0)

#define EMAP_PDB_LOADED(p)                              \
        do {                                            \
                if (!((p)->flags & EMAP_PDBF_LOAD))     \
                        abort();                        \
        } while(0)

#define EMAP_PDB_SORTED(p)                              \
        do {                                            \
                if (!((p)->flags & EMAP_PDBF_SORT))     \
                        abort();                        \
        } while(0)

emap_pointdb_t *emap_pointdb_init(emap_pointdb_t *pdb);
void emap_pointdb_free(emap_pointdb_t *pdb);
int emap_pointdb_load(emap_pointdb_t *pdb, const char *path, uint32_t y_n, const char *c_chars);
int emap_pointdb_sort(emap_pointdb_t *pdb);

#endif /* POINTDB_H */
