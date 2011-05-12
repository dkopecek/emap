#include <stdio.h>
#include <string.h>
#include "POI.h"
#include "PES.h"
#include "DG.h"

typedef struct {
        struct DG_node *dgnode;
        emap_spoint_t  *spoint;
        bool            spcopy;
        bool            firstm; /**< first merge of this basin? */
} SB_t;

DG_t *DG_create(emap_pointdb_t *pdb, emap_surface_t *es, emap_float dE)
{
        SB_t  *sb;
        register size_t i, j, k;
        size_t sbcount;
        emap_float El, Eh; /**< energy band boundaries - low, high */
        size_t Ti; /**< transition point array index */

        sbcount = es->mcount;
        sb      = alloc_array(SB_t, sbcount);

        /**
         * Initialize the lowest level of basins, they are equivalent to
         * all local minima of the abstract surface.
         */
#pragma omp parallel for
        for (i = 0; i < sbcount; ++i) {
                sb[i].spoint = es->mpoint[i];
                sb[i].spcopy = true;
                sb[i].firstm = true;
                sb[i].dgnode = alloc_type(struct DG_node);
                sb[i].dgnode->child = NULL;
                sb[i].dgnode->count = i;
                sb[i].dgnode->point = es->mpoint[i];
        }

        El   = sb[0].spoint->cmaximum->y;
        Eh   = El + dE;
        Ti   = 0;

        while (sbcount > 1) {
                /*
                 * Find the index of the lowest transition point with
                 * energy larger than or equal to `El'.
                 * XXX: implement as binary search since the array is sorted
                 */
                while(es->tpoint[Ti]->cminimum->y < El)
                        ++Ti;
#ifndef NDEBUG
                fprintf(stderr, "DEBUG: Ti = %zu\n", Ti);
#endif
                if (es->tpoint[Ti]->cminimum->y < Eh) {
                        /**
                         * Perform superbasin analysis, merge basins if they are connected
                         * by a pathway in the current energy band. Merging creates new DG
                         * nodes with edges to the basins that are being merged.
                         */
                        for (i = 0; i < sbcount; ++i) {
                        __restart:
                                for (j = i + 1; j < sbcount; ++j) {
                                        /*
                                         * find the minimal distance from `i' to `j' over
                                         * transition points with energy larger than or
                                         * equal to`El'.
                                         */
                                        emap_spoint_t *Tmin = NULL; /**< transition point of the minimal pathway distance */
                                        size_t         Dmin = SIZE_MAX; /**< minimal distance */
                                        size_t         Dcur; /**< for intermediate results */

                                        for (k = Ti; k < es->tcount; ++k) {
                                                Dcur = emap_spoint_mindistance(sb[i].spoint, es->tpoint[k], pdb) + \
                                                        emap_spoint_mindistance(sb[j].spoint, es->tpoint[k], pdb);

                                                if (Dcur < Dmin) {
                                                        Tmin    = es->tpoint[k];
                                                        Dmin    = Dcur;
                                                }
                                        }

                                        assert(Tmin != NULL);

                                        /**
                                         * If `Tmin' is in the current energy band we can merge the basins.
                                         */
                                        if (Tmin->cminimum->y < Eh) {
#ifndef NDEBUG
                                                fprintf(stderr, "DEBUG: found a transition point between sb[%zu] <-> sb[%zu] in the current energy band\n", i, j);
#endif
                                                /**
                                                 * merge sb[j] into sb[i]
                                                 */
                                                if (sb[i].firstm) {
                                                        struct DG_node *node = alloc_type(struct DG_node);

                                                        node->child = alloc_array(struct DG_node *, 2);
                                                        node->point = NULL;
                                                        node->count = 2;
                                                        node->child[0] = sb[i].dgnode;
                                                        node->child[1] = sb[j].dgnode;

                                                        sb[i].dgnode = node;
                                                        sb[i].firstm = false;
                                                } else {
                                                        sb[i].dgnode->child = realloc_array(sb[i].dgnode->child,
                                                                                            struct DG_node *, ++(sb[i].dgnode->count));
                                                        sb[i].dgnode->child[sb[i].dgnode->count - 1] = sb[j].dgnode;
                                                }

                                                if (sb[i].spcopy) {
                                                        emap_spoint_t *spcopy = alloc_type(emap_spoint_t);

                                                        spcopy->flags     = sb[i].spoint->flags;
                                                        spcopy->cmaximum  = sb[i].spoint->cmaximum;
                                                        spcopy->cminimum  = sb[i].spoint->cminimum;

                                                        spcopy->component = alloc_array(emap_point_t *, sb[i].spoint->compcount);
                                                        spcopy->compcount = sb[i].spoint->compcount;

                                                        memcpy(spcopy->component, sb[i].spoint->component,
                                                               sizeof(emap_point_t *) * spcopy->compcount);

                                                        sb[i].spcopy = false;
                                                        sb[i].spoint = spcopy;
                                                }

                                                emap_spoint_merge(sb[i].spoint, sb[j].spoint);
                                                array_remove(sb, &sbcount, j);
                                                goto __restart;
                                        }
                                }

                                /* reset the first merge flag for the next level */
                                sb[i].firstm = true;
                        }
                }

                El = Eh;
                Eh = El + dE;
#ifndef NDEBUG
                fprintf(stderr, "DEBUG: moving to energy band <%f, %f>\n", El, Eh);
#endif
        }

        return NULL;
}

int DG_write(DG_t *dg, const char *path)
{
        return (-1);
}
