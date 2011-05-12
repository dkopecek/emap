#ifndef PES_H
#define PES_H

#include <stddef.h>
#include "pointdb.h"

typedef struct {
        uint32_t       flags;
        emap_point_t  *cmaximum;  /**< pointer to a component with the highest energy */
        emap_point_t  *cminimum;  /**< pointer to a component with the lowest energy */
        emap_point_t **component; /**< surface point components */
        size_t         compcount; /**< number of components */
} emap_spoint_t;

typedef struct {
        emap_spoint_t **tpoint;
        size_t          tcount;
        emap_spoint_t **mpoint;
        size_t          mcount;
} emap_surface_t;

int emap_spoint_ymaxcmp(const emap_spoint_t *a, const emap_spoint_t *b);
int emap_spoint_ymincmp(const emap_spoint_t *a, const emap_spoint_t *b);
int emap_spoint_mindistance(const emap_spoint_t *a, const emap_spoint_t *b, emap_pointdb_t *pdb);

int emap_spoint_merge(emap_spoint_t *dst, const emap_spoint_t *src);

int emap_surface_init(emap_surface_t *es);
emap_surface_t *emap_surface_new(void);

#endif /* PES_H */
