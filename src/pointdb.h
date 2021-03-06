#ifndef POINTDB_H
#define POINTDB_H

#include <stdbool.h>
#include <stdint.h>
#include "emap.h"
#include "helpers.h"
#include "rbt_common.h"

typedef struct {
        uint32_t   flags; /**< flags */
        uint32_t  *ptkey; /**< point key */
        uint32_t   line;  /**< line number in the original file */
        emap_float y;     /**< dependent variable */
        emap_float x[];   /**< inline array of independent vari`ables */
} emap_point_t;

#define EMAP_POINT_VISITED   0x00000001
#define EMAP_POINT_MINIMUM_L 0x00000010
#define EMAP_POINT_MINIMUM_G 0x00000020

/**
 * Atomically set flag `f' on point `p'.
 * @return true if the flag was set by this operation, false if the flag was already set
 */
static inline bool emap_point_setflag(emap_point_t *p, uint32_t f)
{
        register uint32_t o;
#ifndef NDEBUG
        /* Check that only one flag is being set */
        register uint32_t i, c;
        for (i = f, c = 0; i > 0; i>>=1)
                if (f & 1) ++c;
#endif
        o = __sync_fetch_and_or(&p->flags, f);

        return (o & f ? false : true);
}

/**
 * Atomically read flag `f' from point `p'.
 * @return true if the flag is set, false otherwise
 */
static inline bool emap_point_getflag(emap_point_t *p, uint32_t f)
{
#ifndef NDEBUG
        /* Check that only one flag is being read */
        register uint32_t i, c;
        for (i = f, c = 0; i > 0; i>>=1)
                if (f & 1) ++c;
#endif
        return (bool)(__sync_fetch_and_add(&p->flags, 0) & f);
}

#define EMAP_POINT_SIZE(arity) (sizeof(emap_point_t) + (sizeof(emap_float) * (arity)))

typedef struct {
        uint16_t       flags; /**< flags */
        uint16_t       arity; /**< number of independent variables */
        size_t         count; /**< number of points in the array */
        emap_point_t  *point; /**< array of points */
        emap_point_t **psort; /**< array of sorted points by each independent variable */
        rbt_t         *pindex; /**< point index for nearest neighbor search */
        uint32_t      *keymax; /**< maximum values of each key component */
	enum emap_transfn y_trans; /**< y transformation function */
} emap_pointdb_t;

static inline emap_point_t *emap_pointp(emap_pointdb_t *pdb, size_t pi)
{
        return (emap_point_t *)(((uint8_t *)(pdb->point)) + (EMAP_POINT_SIZE(pdb->arity) * pi));
}

static inline emap_point_t *emap_psortp(emap_pointdb_t *pdb, size_t xn, size_t xi)
{
        return ((emap_point_t **)(pdb->psort + xn * pdb->count))[xi];
}

#define EMAP_PDBF_INIT 0x0001 /**< set if properly initialized */
#define EMAP_PDBF_FREE 0x0002 /**< set if memory for the struct was dynamically allocated */
#define EMAP_PDBF_LOAD 0x0004 /**< set if some data points were loaded */
#define EMAP_PDBF_SORT 0x0008 /**< set if the data points are sorted */

#define EMAP_PDB_INITIALIZED(p)                                         \
        do {                                                            \
                if (__predict(!((p)->flags & EMAP_PDBF_INIT), 0))       \
                        abort();                                        \
        } while(0)

#define EMAP_PDB_LOADED(p)                                              \
        do {                                                            \
                if (__predict(!((p)->flags & EMAP_PDBF_LOAD), 0))       \
                        abort();                                        \
        } while(0)

#define EMAP_PDB_SORTED(p)                                              \
        do {                                                            \
                if (__predict(!((p)->flags & EMAP_PDBF_SORT), 0))       \
                        abort();                                        \
        } while(0)

emap_pointdb_t *emap_pointdb_init(emap_pointdb_t *pdb);
void emap_pointdb_free(emap_pointdb_t *pdb);
int emap_pointdb_load(emap_pointdb_t *pdb, const char *path, uint32_t y_n, uint32_t skip_x_n[], size_t nmemb, const char *c_chars);
int emap_pointdb_index(emap_pointdb_t *pdb);

void emap_pointdb_apply(emap_pointdb_t *pdb, void (*fn)(emap_pointdb_t *, emap_point_t *));

void emap_pointdb_apply_r(emap_pointdb_t *pdb, void (*fn)(emap_pointdb_t *, emap_point_t *, void *), void *fnarg);

uint32_t emap_pointnb_applyp(emap_pointdb_t *pdb, emap_point_t *p, uint32_t (*fn)(emap_pointdb_t *, emap_point_t *, emap_point_t *));

int  emap_point_keydistance(const emap_point_t *a, const emap_point_t *b, emap_pointdb_t *pdb);
bool emap_point_keyneighbor(const emap_point_t *a, const emap_point_t *b, emap_pointdb_t *pdb);

#endif /* POINTDB_H */
