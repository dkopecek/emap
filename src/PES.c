#include <stdio.h>
#include <string.h>
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

static int spdist_bcmp(const emap_spoint_t *k, const emap_spdist_t *sp)
{
	if (k < sp->pptr)
		return -1;
	if (k > sp->pptr)
		return  1;
	return 0;
}

size_t emap_spoint_mindistance_fast(const emap_spoint_t *a, const emap_spoint_t *b, emap_pointdb_t *pdb)
{
        assert(a != b);
        assert(a != NULL);
        assert(b != NULL);
	assert(a->cPOI != NULL);

       	emap_spdist_t *res;

	(void)pdb;

	res = bsearch(b, a->cPOI, a->cPOIcount, sizeof(emap_spdist_t),
		      (int(*)(const void *, const void *))spdist_bcmp);

	assert(res != NULL);

	return (res->dist);
}

size_t emap_spoint_mindistance(const emap_spoint_t *a, const emap_spoint_t *b, emap_pointdb_t *pdb)
{
        assert(a != b);
        assert(a != NULL);
        assert(b != NULL);

        size_t i, j, cur, min = SIZE_MAX;

#pragma omp parallel for if(a->compcount > 8) shared(min) private(j, cur) schedule(dynamic,32)
        for (i = 0; i < a->compcount; ++i) {
                for (j = 0; j < b->compcount; ++j) {
                        cur = emap_point_keydistance(a->component[i], b->component[j], pdb);
#pragma omp critical
                        if (cur < min) min = cur;
                }
        }

        return (min);
}

emap_spoint_t *emap_mindist_tpoint(const emap_spoint_t *a, const emap_spoint_t *b, emap_pointdb_t *pdb, emap_float El, emap_float Eh)
{
	emap_spoint_t *Tmin = NULL;
	register size_t d = SIZE_MAX, cur_d;
	register size_t i, j;

	(void)pdb;
	(void)Eh;
	(void)El;

	assert(a != b);
	assert(a != NULL);
	assert(b != NULL);
	assert(a->cPOI != NULL);
	assert(b->cPOI != NULL);

	//cut_d = a->cPOI[0].dist + b->cPOI[0].dist;

	for (i = 0; i < a->cPOIcount; ++i) {
		if (d < a->cPOI[i].dist)
			break;

		if (a->cPOI[i].pptr->cmaximum->y < a->cminimum->y)
			continue;

		for (j = 0; j < b->cPOIcount; ++j) {
			if (d < b->cPOI[j].dist)
				break;

			if (a->cPOI[i].pptr == b->cPOI[j].pptr &&
			    b->cPOI[j].pptr->cmaximum->y >= Eh)
			{
				cur_d = a->cPOI[i].dist + b->cPOI[j].dist;

				if (cur_d >= d) {
					break;
				}

				if (cur_d < d) {
					//fprintf(stdout, "FOO!\n");
					d    = cur_d;
					Tmin = a->cPOI[i].pptr;
#ifndef NDEBUG
					fprintf(stderr, "DEBUG: HIT! i=%zu, j=%zu, dist=%zu, Tmin=%p, E=%"EMAP_FLTFMT"\n", i, j, d, Tmin, Tmin->cminimum->y);
#endif
				}
			}
		}
	}

	return (Tmin);
}

int emap_spoint_merge(emap_spoint_t *dst, const emap_spoint_t *src)
{
	register size_t i;

        dst->flags     &= src->flags;
        dst->component  = realloc_array(dst->component, emap_point_t *,
                                        dst->compcount + src->compcount);

        memcpy(dst->component + dst->compcount,
               src->component, sizeof(emap_point_t *) * src->compcount);

        dst->compcount += src->compcount;

        if (dst->cmaximum->y < src->cmaximum->y)
                dst->cmaximum = src->cmaximum;

        if (dst->cminimum->y > src->cminimum->y)
                dst->cminimum = src->cminimum;

	assert(src->cPOIcount == dst->cPOIcount);

	for (i = 0; i < dst->cPOIcount; ++i) {
		if (dst->cPOI[i].dist > src->cPOI[i].dist)
			dst->cPOI[i].dist = src->cPOI[i].dist;
	}

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
