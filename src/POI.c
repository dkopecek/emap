#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include "emap.h"
#include "helpers.h"
#include "pointdb.h"
#include "POI.h"
#include "PES.h"

uint32_t is_POI(emap_pointdb_t *pdb, emap_point_t *p, emap_point_t *n)
{
        (void)pdb;

        if (p->y < n->y)
                return POI_FLAG_MINIMUM;
        if (p->y > n->y)
                return POI_FLAG_TRANSITION;
        else
                return POI_FLAG_MINIMUM|POI_FLAG_TRANSITION;
}

void collect_POI(emap_pointdb_t *pdb, emap_point_t *p, void *arg)
{
        register POIdb_t *POIdb = (POIdb_t *)arg;
        register uint32_t POI;

        POI = emap_pointnb_applyp(pdb, p, is_POI);

        if (POI != 0) {
                if (POI & POI_FLAG_MINIMUM)
                        emap_point_setflag(p, POI_FLAG_MINIMUM);
                if (POI & POI_FLAG_TRANSITION)
                        emap_point_setflag(p, POI_FLAG_TRANSITION);
#ifdef _OPENMP
                if (__predict(pthread_mutex_lock(&POIdb->mutex) != 0, 0))
                        abort();
#endif
                /* realloc the POI array if needed */
                if (__predict(POIdb->palloc == POIdb->pcount, 0)) {
                        POIdb->palloc += POIdb->cntinc;
                        POIdb->points  = realloc_array(POIdb->points, emap_point_t *, POIdb->palloc);

                        if (POIdb->points == NULL)
                                abort();
                }

                POIdb->points[POIdb->pcount++] = p;

                if (POI & POI_FLAG_MINIMUM)
                        ++POIdb->mcount;
                if (POI & POI_FLAG_TRANSITION)
                        ++POIdb->tcount;
#ifdef _OPENMP
                if (__predict(pthread_mutex_unlock(&POIdb->mutex) != 0, 0))
                        abort();
#endif
        }
}

static int POI_cmp(const emap_point_t **a, const emap_point_t **b, emap_pointdb_t *pdb)
{
        if (emap_point_keydistance(*a, *b, pdb) <= pdb->arity)
                return (0);
        else
                return (1);
}

emap_surface_t *POI_postprocess(emap_pointdb_t *pdb, POIdb_t *POIdb)
{
        emap_surface_t *es;

        es = emap_surface_new();

        while (POIdb->pcount > 0) {
                emap_spoint_t *sp = alloc_type(emap_spoint_t);
                register size_t i, j;

                if (sp == NULL) {
#ifndef NDEBUG
                        fprintf(stderr, "DEBUG: can't allocate memory for a new surface point\n");
#endif
                        goto __fail;
                }

                sp->component = alloc_array(emap_point_t *, 1);
                sp->compcount = 1;

                sp->component[0] = sp->cmaximum = POIdb->points[POIdb->pcount - 1];
                array_remove(POIdb->points, &POIdb->pcount, POIdb->pcount - 1);

                for (i = 0; i < POIdb->pcount; ++i) {
                        for(j = 0; j < sp->compcount; ++j) {
                                /*
                                 * Merge neighbors
                                 */
                                if (emap_point_keyneighbor(POIdb->points[i],
                                                           sp->component[j], pdb))
                                {
#ifndef NDEBUG
                                        fprintf(stderr, "DEBUG: adding new point %p to component %p\n",
                                                POIdb->points[i], sp->component);
#endif
                                        /* update component points */
                                        sp->component = realloc_array(sp->component, emap_point_t *, ++sp->compcount);
                                        sp->component[sp->compcount - 1] = POIdb->points[i];

                                        /* update component maximum */
                                        if (sp->cmaximum->y < sp->component[sp->compcount - 1]->y)
                                                sp->cmaximum = sp->component[sp->compcount - 1];

                                        /* update the array of local minina */
                                        array_remove(POIdb->points, &POIdb->pcount, i);

                                        /* start over */
                                        i = 0;
                                        break;
                                }
                        }
                }

#ifndef NDEBUG
                fprintf(stderr, "DEBUG: new surface point %p with %zu components\n", sp, sp->compcount);
                for (i = 0; i < sp->compcount; ++i)
                        fprintf(stderr, "DEBUG:  y[%zu](%p) = %f\n", i, sp->component[i], sp->component[i]->y);
#endif
        }

        return (es);
__fail:
        if (es != NULL) {
                /* free */
        }

        return (NULL);
}
