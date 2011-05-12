#include "POI.h"
#include "PES.h"
#include "DG.h"

typedef struct {
        struct DG_node *dgnode;
        emap_spoint_t  *spoint;
} SB_t;

DG_t *DG_create(emap_surface_t *es, emap_float dE)
{
        SB_t  *sb;
        register size_t sbcount, i;

        emap_float El, Eh; /**< energy band boundaries - low, high */

        sbcount = es->mcount;
        sb      = alloc_array(SB_t, sbcount);

#pragma omp parallel for
        for (i = 0; i < sbcount; ++i) {
                sb[i].spoint = es->mpoint[i];
                sb[i].dgnode = alloc_type(struct DG_node);
                sb[i].dgnode->child = NULL;
                sb[i].dgnode->count = i;
                sb[i].dgnode->point = es->mpoint[i];
        }

        El = sb[0].spoint->cmaximum->y;
        Eh = El + dE;

        while (sbcount > 1) {
                /*
                 * 
                 */
        }

        return NULL;
}

int DG_write(DG_t *dg, const char *path)
{
        return (-1);
}
