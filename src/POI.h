#ifndef POI_H
#define POI_H

#include "pointdb.h"
#include "PES.h"

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

uint32_t is_POI(emap_pointdb_t *pdb, emap_point_t *p, emap_point_t *n);
void     collect_POI(emap_pointdb_t *pdb, emap_point_t *p, void *arg);
emap_surface_t *POI_postprocess(emap_pointdb_t *pdb, POIdb_t *POIdb, bool progress);

#endif /* POI_H */
