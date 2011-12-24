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

#if 0
static int POI_cmp(const emap_point_t **a, const emap_point_t **b, emap_pointdb_t *pdb)
{
        if (emap_point_keydistance(*a, *b, pdb) <= pdb->arity)
                return (0);
        else
                return (1);
}
#endif

static int spoint_ymaxcmp(const emap_spoint_t **a, const emap_spoint_t **b)
{
        return emap_spoint_ymaxcmp(*a, *b);
}

static int spoint_ymincmp(const emap_spoint_t **a, const emap_spoint_t **b)
{
        return emap_spoint_ymincmp(*a, *b);
}

static int _spdist_cmp(const emap_spdist_t *a, const emap_spdist_t *b)
{
	if (a->pptr < b->pptr)
		return -1;
	if (a->pptr > b->pptr)
		return  1;
	return emap_spoint_ymincmp(a->pptr, b->pptr);
}

emap_surface_t *POI_postprocess(emap_pointdb_t *pdb, POIdb_t *POIdb, bool progress)
{
        emap_surface_t *es;
	register size_t i, j;
	emap_spdist_t *dist;
	emap_spoint_t *cmin;

        es = emap_surface_new();

        while (POIdb->pcount > 0) {
                emap_spoint_t *sp = alloc_type(emap_spoint_t);

                if (sp == NULL) {
#ifndef NDEBUG
                        fprintf(stderr, "DEBUG: can't allocate memory for a new surface point\n");
#endif
                        goto __fail;
                }

                sp->component = alloc_array(emap_point_t *, 1);
                sp->compcount = 1;

		sp->cPOIcount = 0;
		sp->cPOI = NULL;

                sp->component[0] = POIdb->points[POIdb->pcount - 1];
                sp->flags    = sp->component[0]->flags;
                sp->cmaximum = sp->cminimum = sp->component[0];
                array_remove(POIdb->points, &POIdb->pcount, POIdb->pcount - 1);

                for (i = 0; i < POIdb->pcount; ++i) {
			if (progress) {
				printf("\r[i] %3.u%% m=%zu t=%zu p=%zu",
				       (unsigned int)(((float)(i + 1)/(float)POIdb->pcount) * 100.0),
				       es->mcount, es->tcount, POIdb->pcount);
			}

                        for(j = 0; j < sp->compcount; ++j) {

                                /*
                                 * Merge neighbors
                                 */
                                if (emap_point_keyneighbor(POIdb->points[i],
                                                           sp->component[j], pdb) && (sp->flags & POIdb->points[i]->flags))
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
                                        /* update component minimum */
                                        else if (sp->cminimum->y > sp->component[sp->compcount - 1]->y)
                                                sp->cminimum = sp->component[sp->compcount - 1];

                                        /* update the array of local minina */
                                        array_remove(POIdb->points, &POIdb->pcount, i);

                                        /* start over */
                                        i = 0;
                                        break;
                                }
                        }
                }

		if (progress)
			printf("\n");

#ifndef NDEBUG
                fprintf(stderr, "DEBUG: new surface point %p with %zu components\n", sp, sp->compcount);
#endif
                sp->flags = UINT32_MAX; /**< set all bits */

                for (i = 0; i < sp->compcount; ++i) {
#ifndef NDEBUG
                        fprintf(stderr, "DEBUG:  y[%zu](%p) = %"EMAP_FLTFMT"\n", i, sp->component[i], sp->component[i]->y);
#endif
                        sp->flags &= sp->component[i]->flags;
                }

                assert((sp->flags & (POI_FLAG_MINIMUM|POI_FLAG_TRANSITION)) != (POI_FLAG_MINIMUM|POI_FLAG_TRANSITION));
#ifndef NDEBUG
                fprintf(stderr, "DEBUG:  the surface point is a %s\n",
                        sp->flags & POI_FLAG_MINIMUM ? "local minimum" : "transition point");
#endif
                if (sp->flags & POI_FLAG_MINIMUM) {
                        es->mpoint = realloc_array(es->mpoint, emap_spoint_t *, ++es->mcount);
                        es->mpoint[es->mcount - 1] = sp;
                } else {
                        es->tpoint = realloc_array(es->tpoint, emap_spoint_t *, ++es->tcount);
                        es->tpoint[es->tcount - 1] = sp;
                }
        }

        /**
         * Sort the minima so we can walk thru the array lineary
         * as the energy band moves up to higher energies.
         */
        qsort(es->mpoint, es->mcount, sizeof(emap_spoint_t *),
              (int(*)(const void *, const void *))spoint_ymaxcmp);

        /**
         * Sort the transition points by the lowest energy of a component
         * since we are interested in such points between two minima.
         */
        qsort(es->tpoint, es->tcount, sizeof(emap_spoint_t *),
              (int(*)(const void *, const void *))spoint_ymincmp);

#ifndef NDEBUG
        fprintf(stderr, "DEBUG: initializing the cPOI field for each minumum (mcount=%zu, tcount=%zu)\n", es->mcount, es->tcount);
#endif
	/**
	 * For each minumum, calculate the distance to each transition point
	 */
#pragma omp parallel for private(cmin, dist)
	for (i = 0; i < es->mcount; ++i) {
		dist = alloc_array(emap_spdist_t, es->tcount);
		cmin = es->mpoint[i];

		for (j = 0; j < es->tcount; ++j) {
			dist[j].pptr = es->tpoint[j];
			dist[j].dist = emap_spoint_mindistance(cmin, es->tpoint[j], pdb);
		}

		qsort(dist, es->tcount, sizeof(emap_spdist_t),
		      (int(*)(const void *, const void *))_spdist_cmp);

		cmin->cPOI = dist;
		cmin->cPOIcount = es->tcount;
	}

#ifndef NDEBUG
        fprintf(stderr, "DEBUG: mcount=%zu, tcount=%zu\n", es->mcount, es->tcount);
#endif
        return (es);
__fail:
        if (es != NULL) {
                /* free */
        }

        return (NULL);
}
