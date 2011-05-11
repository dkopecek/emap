#ifndef DG_H
#define DG_H

#include "POI.h"
#include "PES.h"

typedef struct {
        int foo;
} DG_t;

DG_t *DG_create(emap_surface_t *es);
int DG_write(DG_t *dg, const char *path);

#endif /* DG_H */
