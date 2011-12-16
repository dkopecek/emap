#ifndef PES_H
#define PES_H

#include <stddef.h>
#include "pointdb.h"

typedef struct emap_spoint emap_spoint_t;

typedef struct {
	size_t         dist;
	emap_spoint_t *pptr;
} emap_spdist_t;

struct emap_spoint {
        uint32_t       flags;
        emap_point_t  *cmaximum;  /**< pointer to a component with the highest energy */
        emap_point_t  *cminimum;  /**< pointer to a component with the lowest energy */
        emap_point_t **component; /**< surface point components */
        size_t         compcount; /**< number of components */
	emap_spdist_t *cPOI;      /**< POIs sorted by distance from this point */
	size_t         cPOIcount;
};

typedef struct {
        emap_spoint_t **tpoint;
        size_t          tcount;
        emap_spoint_t **mpoint;
        size_t          mcount;
} emap_surface_t;

int emap_spoint_ymaxcmp(const emap_spoint_t *a, const emap_spoint_t *b);
int emap_spoint_ymincmp(const emap_spoint_t *a, const emap_spoint_t *b);

size_t emap_spoint_mindistance(const emap_spoint_t *a, const emap_spoint_t *b, emap_pointdb_t *pdb);
size_t emap_spoint_mindistance_fast(const emap_spoint_t *a, const emap_spoint_t *b, emap_pointdb_t *pdb);

emap_spoint_t *emap_mindist_tpoint(const emap_spoint_t *a, const emap_spoint_t *b, emap_pointdb_t *pdb, emap_float El, emap_float Eh);

int emap_spoint_merge(emap_spoint_t *dst, const emap_spoint_t *src);

int emap_surface_init(emap_surface_t *es);
emap_surface_t *emap_surface_new(void);

#endif /* PES_H */
