#include <assert.h>
#include "PES.h"
#include "helpers.h"

int emap_spoint_ymaxcmp(const emap_spoint_t *a, const emap_spoint_t *b)
{
        assert(a != b);
        assert(a != NULL);
        assert(b != NULL);

        if (a->cmaximum->y < b->cmaximum->y)
                return (-1);
        if (a->cmaximum->y > b->cmaximum->y)
                return (1);
        else
                return (0);
}

int emap_spoint_ymincmp(const emap_spoint_t *a, const emap_spoint_t *b)
{
        assert(a != b);
        assert(a != NULL);
        assert(b != NULL);

        if (a->cminimum->y < b->cminimum->y)
                return (-1);
        if (a->cminimum->y > b->cminimum->y)
                return (1);
        else
                return (0);
}

int emap_spoint_mindistance(const emap_spoint_t *a, const emap_spoint_t *b, emap_pointdb_t *pdb)
{
        assert(a != b);
        assert(a != NULL);
        assert(b != NULL);

        size_t i, j, cur, min = SIZE_MAX;

#pragma omp parallel for private(j, cur)
        for (i = 0; i < a->compcount; ++i) {
                for (j = 0; j < b->compcount; ++j) {
                        cur = emap_point_keydistance(a->component[i], b->component[j], pdb);
#pragma omp critical
                        {
                                if (cur < min)
                                        min = cur;
                        }
                }
        }

        return (min);
}

int emap_surface_init(emap_surface_t *es)
{
        es->tpoint = NULL;
        es->tcount = 0;
        es->mpoint = NULL;
        es->mcount = 0;

        return (-1);
}

emap_surface_t *emap_surface_new(void)
{
        emap_surface_t *es = alloc_type(emap_surface_t);

        if (es == NULL)
                return (NULL);

        emap_surface_init(es);

        return (es);
}
