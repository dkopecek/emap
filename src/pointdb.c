#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <errno.h>
#include "emap.h"
#include "pointdb.h"

emap_pointdb_t *emap_pointdb_init(emap_pointdb_t *pdb)
{
        if (pdb == NULL)
                return (NULL);

        pdb->flags = EMAP_PDBF_INIT;
        pdb->arity = 0;
        pdb->count = 0;
        pdb->point = NULL;
        pdb->psort = NULL;

        return (pdb);
}

void emap_pointdb_free(emap_pointdb_t *pdb)
{
        if (pdb->point != NULL)
                free(pdb->point);

        if (pdb->flags & EMAP_PDBF_FREE)
                free(pdb);
}

int emap_pointdb_load(emap_pointdb_t *pdb, const char *path, uint32_t y_n, const char *c_chars)
{
        fpos_t fl_pos;
        FILE  *fp;
        char   line_buffer[EMAP_LINEBUFFER_SIZE], *row, **tok, *toksave;
        register int i, ts, lines, gi, xi, pi;

#if defined(EMAP_POINT_DOUBLE)
        double res;
#else
        float  res;
#endif

        if (pdb == NULL || path == NULL)
                return (EMAP_EFAULT);

        EMAP_PDB_INITIALIZED(pdb);

        fp = fopen(path, "r");

        if (fp == NULL)
                return (EMAP_ESYSTEM);

        /*
         * Read in the first data row and determine parameters for subsequent reads
         */
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
        tok = malloc(sizeof(char *) * ts);
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

        pdb->arity = ts - 1;
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

        if (y_n == 0)
                y_n = ts - 1;
        else
                --y_n;

        pi = 0;

        do {
                /*
                 * convert
                 */
                for (gi = 0, xi = 0; xi < pdb->arity; ++xi) {
                        toksave = NULL;
                        errno   = 0;

#if defined(EMAP_POINT_DOUBLE)
# define STR2FLT(n, p) strtod((n), (p))
#else
# define STR2FLT(n, p) strtof((n), (p))
#endif
                        /* convert string token to float/double */
                        if (gi == y_n)
                                pdb->point[pi].y     = res = STR2FLT(tok[gi], &toksave);
                        else
                                pdb->point[pi].x[xi] = res = STR2FLT(tok[gi], &toksave);

                        /* error check */
                        if ((res == 0 && toksave == tok[gi]) || errno == ERANGE) {
#ifndef NDEBUG
                                fprintf(stderr, "DEBUG: str2flt failed: xi=%u, gi=%u, tok=%s\n", xi, gi, tok[gi]);
#endif
                                /* TODO: free resources */
                                return (EMAP_EINVAL);
                        }
                }

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
        emap_float d = a->y - b->y;

        if (emap_float_abs(d) <= EMAP_FLTCMP_DELTA)
                return 0;
        if (d > 0)
                return 1;
        else
                return -1;
}

static int _pointcmp_stage2_x = -1;

static int _pointcmp_stage2(const emap_point_t **a, const emap_point_t **b)
{
        emap_float d = (*a)->x[_pointcmp_stage2_x] - (*b)->x[_pointcmp_stage2_x];

        assert(_pointcmp_stage2_x != -1);

        if (emap_float_abs(d) <= EMAP_FLTCMP_DELTA)
                return 0;
        if (d > 0)
                return 1;
        else
                return -1;
}

int emap_pointdb_sort(emap_pointdb_t *pdb)
{
        register int i;

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
        pdb->psort = malloc((sizeof(emap_point_t *) * pdb->count) * pdb->arity);

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
        qsort(pdb->point, pdb->count, EMAP_POINT_SIZE(pdb->arity),
              (int(*)(const void *, const void *))_pointcmp_stage1);

        /*
         * prepare psort arrays for sorting
         */
        for (i = 0; i < pdb->count; ++i)
                pdb->psort[i] = pdb->point + i;

        for (i = 1; i < pdb->arity; ++i)
                memcpy(pdb->psort + (i * pdb->count), pdb->psort, sizeof(emap_point_t *) * pdb->count);

        /*
         * sort the psort arrays
         */
        for (i = 0; i < pdb->arity; ++i) {
                _pointcmp_stage2_x = i;
                qsort(pdb->psort + (i * pdb->count), pdb->count, sizeof(emap_point_t *),
                      (int(*)(const void *, const void *))_pointcmp_stage2);
        }

        _pointcmp_stage2_x = -1;
        pdb->flags |= EMAP_PDBF_SORT;

        return (EMAP_SUCCESS);
}
