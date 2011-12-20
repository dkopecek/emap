#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <errno.h>
#include "emap.h"
#include "helpers.h"
#include "pointdb.h"
#include "rbt_u32mem.h"

emap_pointdb_t *emap_pointdb_init(emap_pointdb_t *pdb)
{
        if (pdb == NULL)
                return (NULL);

        pdb->flags = EMAP_PDBF_INIT;
        pdb->arity = 0;
        pdb->count = 0;
        pdb->point = NULL;
        pdb->psort = NULL;
        pdb->pindex = NULL;
        pdb->keymax = NULL;

        return (pdb);
}

void emap_pointdb_free(emap_pointdb_t *pdb)
{
        if (pdb->point != NULL)
                free(pdb->point);
        if (pdb->psort != NULL)
                free(pdb->psort);
        if (pdb->flags & EMAP_PDBF_FREE)
                free(pdb);
}

static int uint32_cmp(const uint32_t *a, const uint32_t *b)
{
        return *a - *b;
}

int emap_pointdb_load(emap_pointdb_t *pdb, const char *path, uint32_t y_n, uint32_t skip_x_n[], size_t skip_n, const char *c_chars)
{
        fpos_t fl_pos;
        FILE  *fp;
        char   line_buffer[EMAP_LINEBUFFER_SIZE], *row, **tok, *toksave;
        register unsigned int lines, xi, pi;
        register unsigned int ts, i, l;
        unsigned int gi;
        emap_float res;

        if (pdb == NULL || path == NULL)
                return (EMAP_EFAULT);
        if (skip_n > 0 && skip_x_n == NULL)
                return (EMAP_EINVAL);

        EMAP_PDB_INITIALIZED(pdb);

        /*
         * Sort the skip_x_n array so that after we get the column count
         * we can sanitize the skip_n value.
         */
        qsort(skip_x_n, skip_n, sizeof(uint32_t), (int(*)(const void *, const void *))uint32_cmp);

        fp = fopen(path, "r");

        if (fp == NULL)
                return (EMAP_ESYSTEM);

        /*
         * Read in the first data row and determine parameters for subsequent reads
         */
        l = 0; /* track current line number in this variable */

        for(;;) {
                if (fgets(line_buffer, sizeof line_buffer, fp) == NULL) {
                        int r;

                        r = feof(fp) ? EMAP_ENODATA : EMAP_ESYSTEM;
                        fclose(fp);

                        return (r);
                }

                /* skip whitespace characters at the beginning */
                i = 0;
                while(isspace(line_buffer[i])) ++i;

                row = line_buffer + i;

                /* determine whether we've read a complete data row */
                if (strchr(row, '\n') == NULL) {
                        fclose(fp);
                        return (EMAP_ENOBUF);
                }

                ++l;

                /* skip comments */
                if (strchr(c_chars, *row) == NULL)
                        break;
#ifndef NDEBUG
                fprintf(stderr, "DEBUG: skipping comment line\n");
#endif
        }

#ifndef NDEBUG
        fprintf(stderr, "DEBUG: found a data row\n");
#endif

        ts  = 8;
        tok = alloc_array(char *, ts);
        toksave = NULL;

        for (i = 0;; ++i) {
                tok[i] = strtok_r(i == 0 ? row : NULL, "\t ", &toksave);

                if (tok[i] == NULL)
                        break;

                if (i+1 == ts) {
                        ts += 8;
                        tok = realloc(tok, sizeof(char *) * ts);
                }
        }

        /* there should be at least 2 columns */
        if (i < 2) {
                fclose(fp);
                free(tok);

                return (EMAP_ENODATA);
        }

        ts  = i;
        tok = realloc(tok, sizeof(char *) * ts);

#ifndef NDEBUG
        fprintf(stderr, "DEBUG: found %d columns\n", i);
#endif
        if (fgetpos(fp, &fl_pos) != 0) {
                fclose(fp);
                free(tok);

                return (EMAP_ESYSTEM);
        }

        lines = 2; /* count in an extra line, the last row may or may not end in '\n' */

        while(!feof(fp)) {
                if (fgetc(fp) == '\n')
                        ++lines;
        }
#ifndef NDEBUG
        fprintf(stderr, "DEBUG: estimated # of rows is ~%d\n", lines);
#endif
        if (fsetpos(fp, &fl_pos) != 0) {
                fclose(fp);
                free(tok);

                return (EMAP_ESYSTEM);
        }

        if (y_n == 0)
                y_n = ts;

        --y_n;

        if (bsearch(&y_n, skip_x_n, skip_n, sizeof(uint32_t),
                    (int(*)(const void *, const void *))uint32_cmp) != NULL)
        {
                /* one of the specified column matches the y column */
                --skip_n;
        }

        /*
         * Update the skip_x_n array size, ignore anything beyond
         * the number of tokens and ensure that at least one column
         * will be read as an independend variable.
         */
        if (skip_n > ts)
                skip_n = ts - 2;

        pdb->arity = ts - 1 - skip_n;
        pdb->count = lines;

#ifndef NDEBUG
        fprintf(stderr, "DEBUG: allocating %zu MiB of memory for data points (point = %zu bytes)\n",
                (EMAP_POINT_SIZE(pdb->arity) * lines)/(1024*1024), EMAP_POINT_SIZE(pdb->arity));
#endif
        pdb->point = malloc(EMAP_POINT_SIZE(pdb->arity) * lines);

        if (pdb->point == NULL) {
                fclose(fp);
                free(tok);

                return (EMAP_ENOMEM);
        }

#ifndef NDEBUG
        fprintf(stderr, "DEBUG: y_n = %u\n", y_n);
#endif
        pi = 0;

        do {
                /*
                 * convert
                 */
                for (gi = 0, xi = 0; gi < ts; ++gi) {
                        toksave = NULL;
                        errno   = 0;

                        /* convert string token to float/double */
                        if (gi == y_n) {
			  emap_pointp(pdb, pi)->y     = res = log(emap_strtoflt(tok[gi], &toksave));
                                emap_pointp(pdb, pi)->flags = 0;
                                emap_pointp(pdb, pi)->ptkey = alloc_array(uint32_t, pdb->arity);
                        } else {
                                if (bsearch(&gi, skip_x_n, skip_n, sizeof(uint32_t),
                                            (int(*)(const void *, const void *))uint32_cmp) == NULL)
                                        {
                                                emap_pointp(pdb, pi)->x[xi] = res = emap_strtoflt(tok[gi], &toksave);
                                                ++xi;
                                        }
                                else
                                        continue;
                        }

#if !defined(NDEBUG) && defined(EMAP_PRINT_POINTS)
                        fprintf(stderr, "\"%s\" => %"EMAP_FLTFMT"| ", tok[gi], res);
#endif
                        /* error check */
                        if ((res == 0 && toksave == tok[gi]) || errno == ERANGE) {
#ifndef NDEBUG
                                fprintf(stderr, "DEBUG: str2flt failed: xi=%u, gi=%u, tok=%s\n", xi, gi, tok[gi]);
#endif
                                /* TODO: free resources */
                                return (EMAP_EINVAL);
                        }
                }
#if !defined(NDEBUG) && defined(EMAP_PRINT_POINTS)
                fprintf(stderr, "\n");
#endif
                emap_pointp(pdb, pi)->line = l;

                /*
                 * checks
                 */
                if (++pi >= pdb->count) {
#ifndef NDEBUG
                        fprintf(stderr, "DEBUG: more data points than estimated!\n");
#endif
                        return (EMAP_EUNKNOWN);
                }

                /*
                 * prepare next token array
                 */
                for(;;) {
                        if (fgets(line_buffer, sizeof line_buffer, fp) != NULL) {
                                /* skip whitespace characters at the beginning */
                                i = 0;
                                while(isspace(line_buffer[i])) ++i;

                                row = line_buffer + i;

                                /* skip comments */
                                if (strchr(c_chars, *row) == NULL)
                                        break;
#ifndef NDEBUG
                                fprintf(stderr, "DEBUG: skipping comment/empty line\n");
#endif
                        } else {
                                if (feof(fp))
                                        goto _done;
                                else if (ferror(fp)) {
#ifndef NDEBUG
                                        fprintf(stderr, "DEBUG: error while reading line\n");
#endif
                                        return (EMAP_ESYSTEM);
                                } else {
#ifndef NDEBUG
                                        fprintf(stderr, "DEBUG: unknown error while reading line\n");
#endif
                                        return (EMAP_EUNKNOWN);
                                }
                        }
                }

                for (i = 0; i < ts; ++i) {
                        tok[i] = strtok_r(i == 0 ? row : NULL, "\t ", &toksave);

                        if (tok[i] == NULL)
                                break;
                }

                if (i != ts) {
#ifndef NDEBUG
                        fprintf(stderr, "DEBUG: number of columns changed! i=%u\n", i);
#endif
                        free(tok);
                        return (EMAP_EINVAL);
                }

                ++l;
        } while(!feof(fp));
_done:
        fclose(fp);
        free(tok);

        if (pdb->count > pi)
                pdb->point = realloc(pdb->point, EMAP_POINT_SIZE(pdb->arity) * pi);

        pdb->count  = pi;
        pdb->flags |= EMAP_PDBF_LOAD;

#ifndef NDEBUG
        fprintf(stderr, "DEBUG: loaded %zu data points from %s\n", pdb->count, path);
#endif
        return (EMAP_SUCCESS);
}

static int _pointcmp_stage1(const emap_point_t *a, const emap_point_t *b)
{
#if !defined(NDEBUG) && defined(EMAP_PRINT_POINTS)
        fprintf(stderr, "%"EMAP_FLTFMT" - %"EMAP_FLTFMT"\n", a->y, b->y);
#endif
        if (a->y > b->y)
                return 1;
        if (a->y < b->y)
                return -1;
        else
                return 0;
}

static int _pointcmp_stage2(const emap_point_t **a, const emap_point_t **b, void *arg)
{
        uint8_t xi = (uint8_t)(arg);

        if ((*a)->x[xi] > (*b)->x[xi])
                return 1;
        if ((*a)->x[xi] < (*b)->x[xi])
                return -1;
        else
                return 0;
}

int emap_pointdb_index(emap_pointdb_t *pdb)
{
        register unsigned int i, x, k;
        register uint8_t i8;

        if (pdb == NULL)
                return (EMAP_EFAULT);

        EMAP_PDB_INITIALIZED(pdb);
        EMAP_PDB_LOADED(pdb);

        /*
         * prepare pointer arrays
         */
#ifndef NDEBUG
        fprintf(stderr, "DEBUG: Allocating %zu MiB of memory for the `psort' array\n",
                (sizeof(emap_point_t *) * pdb->count * pdb->arity)/(1024*1024));
#endif
        pdb->psort = alloc_array(emap_point_t *, pdb->count * pdb->arity);

        if (pdb->psort == NULL) {
#ifndef NDEBUG
                fprintf(stderr, "DEBUG: Allocation failed.\n");
#endif
                return (EMAP_ENOMEM);
        }

        /*
         * sort by the dependent variable first, since it would invalidate the
         * pointer array otherwise.
         */
#ifndef NDEBUG
        fprintf(stderr, "DEBUG: Sorting `point' array %p, %zu items, item size is %zu (arity is %u)\n",
                pdb->point, pdb->count, EMAP_POINT_SIZE(pdb->arity), pdb->arity);
#endif
        qsort(&pdb->point[0], pdb->count, EMAP_POINT_SIZE(pdb->arity),
              (int(*)(const void *, const void *))_pointcmp_stage1);

        /*
         * prepare psort arrays for sorting
         */
        for (i = 0; i < pdb->count; ++i) /* don't parallelize */
                pdb->psort[i] = emap_pointp(pdb, i);

        for (i = 1; i < pdb->arity; ++i)
                memcpy(pdb->psort + (i * pdb->count), pdb->psort, sizeof(emap_point_t *) * pdb->count);

        assert(sizeof(uint8_t) <= sizeof(void *));

        /*
         * sort the psort arrays
         */
#pragma omp parallel for
        for (i8 = 0; i8 < pdb->arity; ++i8) {
#ifndef NDEBUG
                fprintf(stderr, "DEBUG: Sorting `psort' array %p (x = %u)\n", pdb->psort + (i8 * pdb->count), i8);
#endif
                qsort_r(pdb->psort + (i8 * pdb->count), pdb->count, sizeof(emap_point_t *),
                        (int(*)(const void *, const void *, void *))_pointcmp_stage2, (void *)(i8));
        }

        pdb->flags |= EMAP_PDBF_SORT;

#ifndef NDEBUG
        {
                emap_float prev;
                /*
                 * Check that everything is sorted properly
                 */
                prev = emap_pointp(pdb, 0)->y;

                for (i = 1; i < pdb->count; ++i) {
                        assert(prev <= emap_pointp(pdb, i)->y);
                        prev = emap_pointp(pdb, i)->y;
                }

                for (x = 0; x < pdb->arity; ++x) {
                        prev = emap_psortp(pdb, x, 0)->x[x];

                        for (i = 1; i < pdb->count; ++i) {
                                assert(prev <= emap_psortp(pdb, x, i)->x[x]);
                                prev = emap_psortp(pdb, x, i)->x[x];
                        }
                }
        }
#endif
        /*
         * Walk thru all psort arrays and get the key segments (index of each variable)
         */
        pdb->keymax = alloc_array(uint32_t, pdb->arity);

#pragma omp parallel for private(i, k) /* `x' is private by default */
        for (x = 0; x < pdb->arity; ++x) {
                emap_psortp(pdb, x, 0)->ptkey[x] = k = 0;

                for(i = 1; i < pdb->count; ++i) {
                        if (emap_psortp(pdb, x, i - 1)->x[x] < emap_psortp(pdb, x, i)->x[x])
                                ++k;
                        else if (emap_psortp(pdb, x, i - 1)->x[x] > emap_psortp(pdb, x, i)->x[x])
                                ++k;

                        emap_psortp(pdb, x, i)->ptkey[x] = k;
                }

                pdb->keymax[x] = k;
#ifndef NDEBUG
                fprintf(stderr, "DEBUG: highest index for psort(%p, %d) is %d\n", pdb, x, k);
#endif
        }

        if (pdb->pindex != NULL)
                rbt_u32mem_free(pdb->pindex);

        pdb->pindex = rbt_u32mem_new();

        for (i = 0; i < pdb->count; ++i) {
                if (rbt_u32mem_add(pdb->pindex,
                                   emap_pointp(pdb, i)->ptkey, pdb->arity, emap_pointp(pdb, i)) != 0)
                {
                        int a;

                        fprintf(stderr, "DEBUG: rbt_u32mem_add(%p, %p, %u, %p) failed: i=%u\n",
                                pdb->pindex, emap_pointp(pdb, i)->ptkey, pdb->arity, emap_pointp(pdb, i), i);
                        fprintf(stderr, "ERROR: there seems to be more than one occurence of the independent variable vector ");

                        fprintf(stderr, "(");
                        for (a = 0; a < pdb->arity - 1; ++a)
                                fprintf(stderr, "%"EMAP_FLTFMT" ", emap_pointp(pdb, i)->x[a]);
                        fprintf(stderr, "%"EMAP_FLTFMT")\n", emap_pointp(pdb, i)->x[a]);

                        return (EMAP_EUNKNOWN);
                }
        }

        return (EMAP_SUCCESS);
}

void emap_pointdb_apply(emap_pointdb_t *pdb, void (*fn)(emap_pointdb_t *, emap_point_t *))
{
        register size_t i;

        EMAP_PDB_INITIALIZED(pdb);
        EMAP_PDB_LOADED(pdb);
        EMAP_PDB_SORTED(pdb);

#pragma omp parallel for
        for (i = 0; i < pdb->count; ++i)
                fn(pdb, emap_pointp(pdb, i));

        return;
 }

void emap_pointdb_apply_r(emap_pointdb_t *pdb, void (*fn)(emap_pointdb_t *, emap_point_t *, void *), void *fnarg)
{
        register size_t i;

        EMAP_PDB_INITIALIZED(pdb);
        EMAP_PDB_LOADED(pdb);
        EMAP_PDB_SORTED(pdb);

#pragma omp parallel for
        for (i = 0; i < pdb->count; ++i)
                fn(pdb, emap_pointp(pdb, i), fnarg);

        return;
}

static uint64_t three_to_n(int n)
{
        uint64_t r = 1;

        while (n > 0) {
                r *= 3;
                --n;
        }

        return (r);
}

/**
 * Apply a predicate function `fn' on the neighbors of `p'.
 * @return true if the predicate holds for all neighbors
 */
uint32_t emap_pointnb_applyp(emap_pointdb_t *pdb, emap_point_t *p, uint32_t (*fn)(emap_pointdb_t *, emap_point_t *, emap_point_t *))
{
        /*
         * Number of neighbors should be 3^n - 1
         */
        register uint64_t  n, n_max;
        register uint32_t  res = UINT32_MAX; /* set all bits to 1 */
        register uint32_t *d_key, *n_key;
        register int i;
        emap_point_t *n_point = NULL;
        int32_t d[3] = { 0, -1, 1 };

#ifndef NDEBUG
        EMAP_PDB_INITIALIZED(pdb);
        EMAP_PDB_LOADED(pdb);
        EMAP_PDB_SORTED(pdb);
#endif
        d_key = alloc_array(uint32_t, pdb->arity);
        n_key = alloc_array(uint32_t, pdb->arity);
        n_max = three_to_n(pdb->arity);

        /* Initialize neighbor key delta */
        for (i = 1; i < pdb->arity; ++i)
                d_key[i] = 0;

        /* Skip (0, 0, ..., 0) */
        d_key[0] = 1;
        n        = 1;
        goto __first;

        for (;;) {
                while (d_key[0] < 3) {
                __first:
                        /*
                         * p_key + d_key = n_key
                         */
                        for (i = 0; i < pdb->arity; ++i)
                                n_key[i] = p->ptkey[i] + d[d_key[i]];

                        /*
                         * look for the point with key `n_key' in the point index
                         */
                        if (__predict(rbt_u32mem_get(pdb->pindex, n_key, pdb->arity, (void **)&n_point) != 0, 0)) {
#if 0
                                int a;
                                fprintf(stderr, "ERROR: can't find neighboring point: p=(");

                                /* p */
                                for (a = 0; a < pdb->arity - 1; ++a)
                                        fprintf(stderr, "%u ", p->ptkey[a]);
                                fprintf(stderr, "%u), d=(", p->ptkey[a]);

                                /* n */
                                for (a = 0; a < pdb->arity - 1; ++a)
                                        fprintf(stderr, "%d ", d[d_key[a]]);
                                fprintf(stderr, "%d)\n", d[d_key[a]]);
#endif
                                res = 0;
                                goto __ret;
                        }

#if 0
                        if (p == n_point) {
                                int a;
                                fprintf(stderr, "DEBUG: p == n_point, d=(");
                                /* n */
                                for (a = 0; a < pdb->arity - 1; ++a)
                                        fprintf(stderr, "%u ", d[d_key[a]]);
                                fprintf(stderr, "%u)\n", d[d_key[a]]);
                        }

                        assert(p != n_point);
#endif
                        res = res & fn(pdb, p, n_point);

                        if (res == 0 || ++n == n_max) {
//#ifndef NDEBUG
//                                fprintf(stderr, "DEBUG: res=%s, n = %"PRIu64"\n", res ? "true" : "false", n);
//#endif
                                goto __ret;
                        }

                        d_key[0] += 1;
                }

                i = 0;

                do {
                        d_key[i] = 0;
                        assert(i < pdb->arity);

                        ++i;
                        d_key[i] += 1;
#if 0
                        int a;
                        fprintf(stderr, "DEBUG: d=(");
                        /* n */
                        for (a = 0; a < pdb->arity - 1; ++a)
                                fprintf(stderr, "%d ", d[d_key[a]]);
                        fprintf(stderr, "%d)\n", d[d_key[a]]);
#endif
                } while(d_key[i] == 3);
        }

__ret:
        free(n_key);
        free(d_key);

        return (res);
}

int emap_point_keydistance(const emap_point_t *a, const emap_point_t *b, emap_pointdb_t *pdb)
{
        register int i, r = 0;

        for (i = 0; i < pdb->arity; ++i)
                r += abs(a->ptkey[i] - b->ptkey[i]);

        return r;
}

bool emap_point_keyneighbor(const emap_point_t *a, const emap_point_t *b, emap_pointdb_t *pdb)
{
        register int  i;
        register bool r = true;

        for (i = 0; i < pdb->arity; ++i)
                r &= abs(a->ptkey[i] - b->ptkey[i]) <= 1 ? true : false;

        return r;
}
