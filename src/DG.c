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

DG_t *DG_create(emap_pointdb_t *pdb, emap_surface_t *es, emap_float dE, bool progress)
{
        DG_t *dg;
        SB_t *sb;
        register size_t i, j, k;
        size_t sbcount;
        emap_float El, Eh; /**< energy band boundaries - low, high */
        size_t Ti; /**< transition point array index */
	size_t Ti_min;
	size_t ncount;
        uint32_t next_id;
	int progress_out = 0;
	int merge_event = 0;

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
		sb[i].dgnode->El    = es->mpoint[i]->cminimum->y;
                sb[i].dgnode->id    = i;
		sb[i].dgnode->msum  = 1;
        }

	ncount  = sbcount;
        next_id = sbcount;

        El   = sb[0].spoint->cminimum->y;
        Eh   = El + dE;
        Ti   = 0;
	Ti_min = SIZE_MAX;

        while (sbcount > 1) {
                /*
                 * Find the index of the lowest transition point with
                 * energy larger than or equal to `El'.
                 * XXX: implement as binary search since the array is sorted
                 */
                while(es->tpoint[Ti]->cmaximum->y < El)
                        ++Ti;
#ifndef NDEBUG
                fprintf(stderr, "DEBUG: Ti = %zu\n", Ti);
#endif
                if (es->tpoint[Ti]->cmaximum->y < Eh) {
                        /**
                         * Perform superbasin analysis, merge basins if they are connected
                         * by a pathway in the current energy band. Merging creates new DG
                         * nodes with edges to the basins that are being merged.
                         */

                        for (i = 0; i < sbcount; ++i) {
				if (sb[i].spoint->cminimum->y > Eh)
					continue;
                        __restart:
                                for (j = i + 1; j < sbcount; ++j) {
					if (sb[j].spoint->cminimum->y > Eh)
						continue;
                                        /*
                                         * find the minimal distance from `i' to `j' over
                                         * transition points with energy larger than or
                                         * equal to`El'.
                                         */
                                        emap_spoint_t *Tmin = NULL; /**< transition point of the minimal pathway distance */
                                        size_t         Dmin = SIZE_MAX; /**< minimal distance */
                                        size_t         Dcur; /**< for intermediate results */
                                        size_t         T1max;
					//time_t t1, t2;

                                        T1max = Ti;

                                        while (T1max < es->tcount) {
                                                if (es->tpoint[T1max]->cmaximum->y < Eh)
                                                        ++T1max;
                                                else
                                                        break;
                                        }

                                        --T1max;

#define b_mindist(b1, b2, tp) (emap_spoint_mindistance_fast(b1, tp, pdb) + emap_spoint_mindistance(b2, tp, pdb))

#pragma omp parallel for if((T1max - Ti) >= 8) shared(Tmin, Dmin) private(k, Dcur) schedule(static)
                                        for (k = Ti; k <= T1max; ++k) {
						Dcur = b_mindist(sb[i].spoint, sb[j].spoint, es->tpoint[k]);
#pragma omp critical
                                                {
                                                        if (Dcur < Dmin) {
                                                                Tmin = es->tpoint[k];
                                                                Dmin = Dcur;
								Ti_min = k;
                                                        }
                                                }
                                        }

                                        for (k = T1max + 1; k < es->tcount; ++k) {
						Dcur = b_mindist(sb[i].spoint, sb[j].spoint, es->tpoint[k]);

                                                if (Dcur < Dmin) {
                                                        Tmin = es->tpoint[k];
                                                        Dmin = Dcur;
							Ti_min = k;

                                                        if (Tmin->cminimum->y >= Eh)
                                                                break;
                                                }
                                        }

					assert(Tmin != NULL);

                                        /**
                                         * If `Tmin' is in the current energy band we can merge the basins.
                                         */
                                        if (Tmin->cmaximum->y < Eh &&
					    Tmin->cmaximum->y > sb[i].spoint->cminimum->y &&
					    Tmin->cmaximum->y > sb[j].spoint->cminimum->y)
					{
#ifndef NDEBUG
                                                fprintf(stderr, "DEBUG: found a transition point between sb[%zu] <-> sb[%zu] in the current energy band (sbcount=%zu)\n", i, j, sbcount);
#endif
						merge_event = 1;

                                                /**
                                                 * merge sb[j] into sb[i]
                                                 */
                                                if (sb[i].firstm) {
                                                        struct DG_node *node = alloc_type(struct DG_node);

                                                        node->id    = next_id++;
                                                        node->child = alloc_array(struct DG_node *, 2);
                                                        node->point = NULL;
							node->El    = Eh;
                                                        node->count = 2;
                                                        node->child[0] = sb[i].dgnode;
                                                        node->child[1] = sb[j].dgnode;
							node->msum     = sb[i].dgnode->msum + sb[j].dgnode->msum;

                                                        sb[i].dgnode = node;
                                                        sb[i].firstm = false;

							++ncount;
                                                } else {
                                                        sb[i].dgnode->child = realloc_array(sb[i].dgnode->child,
                                                                                            struct DG_node *, ++(sb[i].dgnode->count));
                                                        sb[i].dgnode->child[sb[i].dgnode->count - 1] = sb[j].dgnode;
                                                        sb[i].dgnode->msum += sb[j].dgnode->msum;
                                                }

                                                if (sb[i].spcopy) {
                                                        emap_spoint_t *spcopy = alloc_type(emap_spoint_t);

                                                        spcopy->flags     = sb[i].spoint->flags;
                                                        spcopy->cmaximum  = sb[i].spoint->cmaximum;
                                                        spcopy->cminimum  = sb[i].spoint->cminimum;

                                                        spcopy->component = alloc_array(emap_point_t *, sb[i].spoint->compcount);
                                                        spcopy->compcount = sb[i].spoint->compcount;
                                                        spcopy->cPOI      = sb[i].spoint->cPOI;
                                                        spcopy->cPOIcount = sb[i].spoint->cPOIcount;

                                                        memcpy(spcopy->component, sb[i].spoint->component,
                                                               sizeof(emap_point_t *) * spcopy->compcount);

                                                        sb[i].spcopy = false;
                                                        sb[i].spoint = spcopy;
                                                }

                                                emap_spoint_merge(sb[i].spoint, sb[j].spoint);
                                                array_remove(sb, &sbcount, j);

                                                if (progress) {
							progress_out = 1;
                                                        printf("\r[i] %3.u%% E=<%"EMAP_FLTFMT", %"EMAP_FLTFMT"> TP=%zu SB=%zu           ",
                                                               (unsigned int)((1.0 - ((double)sbcount/(double)es->mcount)) * 100.0),
                                                               El, Eh, es->tcount - Ti - 1, sbcount);
						}
                                                goto __restart;
                                        }
                                }
				
                                /* reset the first merge flag for the next level */
                                sb[i].firstm = true;
				//printf("A!\n");
                        }
                } else if (es->tpoint[Ti]->cminimum->y > Eh) {
			Ti_min = Ti;
		}

		if (merge_event || Ti_min == SIZE_MAX) {
			El = Eh;
			Eh = El + dE;

			/*while (Eh == El) {
				dE = dE * 10;
				Eh = El + dE;
				}*/

			merge_event = 0;
			Ti_min = SIZE_MAX;
		} else {
			emap_spoint_t *Tmin;
			emap_float dE2;

			Tmin = es->tpoint[Ti_min];
		
			assert(Tmin->cmaximum->y > Eh);
			assert(((Tmin->cmaximum->y - Eh)/ dE) > 0.0);

			//dE2 = (roundl((Tmin->cmaximum->y - Eh)/ dE) - 1) * dE;

			El = Tmin->cmaximum->y - dE/2;
			Eh = El + dE;

			//El  = Eh + dE2;
			//Eh  = El + dE;

			/*
			  while (Eh == El) {
				dE = dE * 10;
				Eh = El + dE;
				}*/

			assert(El != Eh);

			Ti_min = SIZE_MAX;

#ifndef NDEBUG
			fprintf(stderr, "DEBUG: energy band skip: dE2=%"EMAP_FLTFMT"\n", dE2);
#endif
		}
#ifndef NDEBUG
                fprintf(stderr, "DEBUG: moving to energy band <%"EMAP_FLTFMT", %"EMAP_FLTFMT">\n", El, Eh);
#endif
		if (progress && progress_out) {
			progress_out = 0;
			//putchar('\n');
		}
        }

        dg = alloc_type(DG_t);
        dg->root = sb[0].dgnode;
	dg->Etrans = EMAP_ETRANS_NONE;
	dg->mcount = es->mcount;
	dg->ncount = ncount;
	dg->dE = dE;
	dg->Emax = Eh;
	dg->Emin = es->mpoint[0]->cminimum->y;

        return (dg);
}

static void DG_write_node(FILE *fp, struct DG_node *n)
{
        register size_t i;

        if (n == NULL)
                return;

        if (n->child != NULL) {
                for (i = 0; i < n->count; ++i)
                        fprintf(fp, "%u -- %u\n", n->id, n->child[i]->id);
                for (i = 0; i < n->count; ++i)
                        DG_write_node(fp, n->child[i]);
        } else {
                // what?
        }
}

int DG_write(DG_t *dg, const char *path)
{
        FILE *fp;

	if (path == NULL)
		fp = stdout;
	else {
		fp = fopen(path, "w");

		if (fp == NULL) {
#ifndef NDEBUG
			fprintf(stderr, "DEBUG: can't open \"%s\" for writing\n", path);
#endif
			return (-1);
		}
	}

        fprintf(fp, "graph DG {\n");
        fprintf(fp, "\t node  [decorate=false,shape=point,style=filled];\n");
        fprintf(fp, "\t edge  [decorate=false,dir=none,arrowhead=none,arrowtail=none];\n");

        DG_write_node(fp, dg->root);
        fprintf(fp, "}\n");

	if (fp != stdout)
		fclose(fp);

        return (0);
}

static void DG_write_node_emap(FILE *fp, struct DG_node *n, struct DG_node *p)
{
        register size_t i;

        if (n == NULL)
                return;

        if (n->child != NULL) {
		fprintf(fp, "B %u %"EMAP_FLTFMT" %u %zu %zu ",
			n->id, n->El, p != NULL ? p->id : n->id, n->msum, n->count);

                for (i = 0; i < n->count; ++i)
                        fprintf(fp, "%u%c", n->child[i]->id, i+1 == n->count ? '\n' : ' ');

                for (i = 0; i < n->count; ++i)
                        DG_write_node_emap(fp, n->child[i], n);
        } else {
		fprintf(fp, "M %u %"EMAP_FLTFMT" %u\n", n->id, n->El, p != NULL ? p->id : n->id);
        }
}

int DG_write_emap(DG_t *dg, const char *path)
{
	FILE *fp;
	char *Etrans;

	assert (dg != NULL);

	if (path == NULL)
		fp = stdout;
	else {
		fp = fopen(path, "w");
		if (fp == NULL) {
#ifndef NDEBUG
			fprintf(stderr, "DEBUG: can't open \"%s\" for writing\n", path);
#endif
			return (-1);
		}
	}

	switch(dg->Etrans) {
	case EMAP_ETRANS_LOG:
		Etrans = "log";
		break;
	case EMAP_ETRANS_LOG2:
		Etrans = "log2";
		break;
	case EMAP_ETRANS_LOG10:
		Etrans = "log10";
		break;
	case EMAP_ETRANS_NONE:
	default:
		Etrans = "none";
		break;
	}

	fprintf(fp, "#emap-"EMAP_VERSION"\n");
	fprintf(fp, "#fltfmt=%s\n", EMAP_FLTFMT);
	fprintf(fp, "#ncount=%zu\n", dg->ncount);
	fprintf(fp, "#mcount=%zu\n", dg->mcount);
	fprintf(fp, "#dE=%"EMAP_FLTFMT"\n", dg->dE);
	fprintf(fp, "#Etrans=%s\n", Etrans);
	fprintf(fp, "#Emax=%"EMAP_FLTFMT"\n", dg->Emax);
	fprintf(fp, "#EmaxID=%u\n", dg->root ? dg->root->id : 0);
	fprintf(fp, "#Emin=%"EMAP_FLTFMT"\n", dg->Emin);
	fprintf(fp, "#EminID=%u\n", 0); /**< we always assign 0 to the lowest/global minimum */
	fprintf(fp, "#!\n");
	fprintf(fp, "#! B <id> <E> <pid> <msum> <m> <cid> [<cid> ...]\n");
	fprintf(fp, "#! M <id> <E> <pid>\n");
	fprintf(fp, "#!\n");
	fprintf(fp, "#!   id - node identification\n");
	fprintf(fp, "#!    E - node identification\n");
	fprintf(fp, "#!  pid - parent id\n");
	fprintf(fp, "#! msum - sum of multiplicities in the subtree\n");
	fprintf(fp, "#!    m - multiplicity\n");
	fprintf(fp, "#!  cid - child id\n");
	fprintf(fp, "#!\n");
	fprintf(fp, "\n");

	DG_write_node_emap(fp, dg->root, NULL);

	if (fp != stdout)
		fclose(fp);

	return (0);
}
