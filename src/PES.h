#ifndef PES_H
#define PES_H

#include <stddef.h>
#include "pointdb.h"

typedef struct {
        emap_point_t  *cmaximum;  /**< pointer to a component with the highest energy */
        emap_point_t **component; /**< surface point components */
        size_t         compcount; /**< number of components */
} emap_spoint_t;

typedef struct {
        emap_spoint_t *tpoint;
        size_t         tcount;
        emap_spoint_t *mpoint;
        size_t         mcount;
} emap_surface_t;


int emap_surface_init(emap_surface_t *es);
emap_surface_t *emap_surface_new(void);

#endif /* PES_H */
