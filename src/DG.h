#ifndef DG_H
#define DG_H

#include "emap.h"
#include "POI.h"
#include "PES.h"

struct DG_node {
        struct DG_node **child; /**< pointer to the lower level superbasins or the leaf points/minima */
        size_t           count; /**< if child != NULL, the number of child nodes, otherwise the level of the minimum (0 - means that this is a global minimum */
        emap_surface_t  *point; /**< The abstract surface point representing a minimum */
};

typedef struct {
        struct DG_node *root;
} DG_t;

DG_t *DG_create(emap_surface_t *es, emap_float dE);
int DG_write(DG_t *dg, const char *path);

#endif /* DG_H */
