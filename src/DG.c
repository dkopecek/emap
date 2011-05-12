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
        size_t sbcount;

        

        return NULL;
}

int DG_write(DG_t *dg, const char *path)
{
        return (-1);
}
