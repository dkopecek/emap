#ifndef DG_H
#define DG_H

#include "emap.h"
#include "POI.h"
#include "PES.h"
#include "pointdb.h"

struct DG_node {
        struct DG_node **child; /**< pointer to the lower level superbasins or the leaf points/minima */
        size_t           count; /**< if child != NULL, the number of child nodes, otherwise the level of the minimum (0 - means that this is a global minimum */
        emap_spoint_t   *point; /**< The abstract surface point representing a minimum */
        uint32_t         id;
	emap_float       El;    /**< if child != NULL, the energy level of the basin, otherwise the energy of the minimum */
};

enum emap_etrans {
	EMAP_ETRANS_NONE,
	EMAP_ETRANS_LOG,
	EMAP_ETRANS_LOG2,
	EMAP_ETRANS_LOG10
};

typedef struct {
        struct DG_node *root;
	size_t mcount; /**< number of minima */
	emap_float dE;
	enum emap_etrans Etrans; /**< transformation method of the energy values */
	emap_float Emax;
	emap_float Emin;
} DG_t;

DG_t *DG_create(emap_pointdb_t *pdb, emap_surface_t *es, emap_float dE, bool progress);
int DG_write(DG_t *dg, const char *path);
int DG_write_emap(DG_t *dg, const char *path);

#endif /* DG_H */
