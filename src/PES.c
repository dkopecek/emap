#include "PES.h"
#include "helpers.h"

int emap_spoint_ycmp(const emap_spoint_t *a, const emap_spoint_t *b)
{
        if (a->cmaximum->y < b->cmaximum->y)
                return (-1);
        if (a->cmaximum->y > b->cmaximum->y)
                return (1);
        else
                return (0);
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
