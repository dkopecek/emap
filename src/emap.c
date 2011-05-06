#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libgen.h>
#include "emap.h"
#include "pointdb.h"

#define EMAP_SHORT_OPTIONS "c:e:o:hv"

static bool local_minimum_p(emap_pointdb_t *pdb, emap_point_t *p, emap_point_t *n);
static void collect_minima(emap_pointdb_t *pdb, emap_point_t *p);

static void fprintu(FILE *stream, char *arg0)
{
        fprintf(stream, "\n");
        fprintf(stream, " Usage: %s [OPTIONS] <input>\n", basename(arg0));
        fprintf(stream, "\n");
        fprintf(stream, "\t-o <file>      Output file or `-' for stdout (default is stdout).\n\n");
        fprintf(stream, "\t-e <n>         Column number representing the dependent variable:\n");
        fprintf(stream, "\t                  0 ... last column (default)\n");
        fprintf(stream, "\t                  n ... n-th column (max. %u)\n\n", (1<<16) - 1);
        fprintf(stream, "\t-c <str>       Comment characters (lines starting wich such\n");
        fprintf(stream, "\t               character will be skipped; default is \"%s\").\n\n", EMAP_COMMENT_CHARS);
        fprintf(stream, "\t-v             Be verbose.\n");
        fprintf(stream, "\t-h             This help.\n");
        fprintf(stream, "\n");
}

int opt_verbose = 0;

int main(int argc, char *argv[])
{
        char *arg0, *path_in, *path_out;
        int   opt, ret;
        int   opt_yn = 0;
        char *opt_comment = strdup(EMAP_COMMENT_CHARS);

        emap_pointdb_t pdb;

        arg0     = argv[0];
        path_in  = NULL;
        path_out = NULL;

        while ((opt = getopt(argc, argv, EMAP_SHORT_OPTIONS)) != -1) {
                switch(opt) {
                case 'o':
                        path_out = strdup(optarg);
                        break;
                case 'e':
                        opt_yn = atoi(optarg);

                        if (opt_yn < 0 || opt_yn >= 1<<16) {
                                fprintf(stderr,
                                        "Invalid value for option `-e': only values from the interval <%u,%u> are allowed\n", 0, (1<<16) - 1);
                                fprintu(stderr, arg0);
                                return (EXIT_FAILURE);
                        }

                        break;
                case 'c':
                        opt_comment = strdup(optarg);
                        break;
                case 'h':
                        fprintu(stdout, arg0);
                        return (EXIT_SUCCESS);
                case 'v':
                        opt_verbose = 1;
                        break;
                default:
                        fprintf(stderr, "Unknown option `%c'\n", optopt);
                        fprintu(stderr, arg0);
                        return (EXIT_FAILURE);
                }
        }

        argc -= optind;
        argv += optind;

        if (argc != 1) {
                fprintf(stderr, "Expecting exactly one input file.\n");
                fprintu(stderr, arg0);
                return (EXIT_FAILURE);
        }

        path_in = strdup(argv[0]);

        emap_pointdb_init(&pdb);
        ret = emap_pointdb_load(&pdb, path_in, opt_yn, opt_comment);

        if (ret != EMAP_SUCCESS) {
                fprintf(stderr, "Unable to load data points from \"%s\": ret=%d\n", path_in, ret);
                return (EXIT_FAILURE);
        }

        ret = emap_pointdb_index(&pdb);

        if (ret != EMAP_SUCCESS) {
                fprintf(stderr, "Unable to index the data points from \"%s\": ret=%d\n", path_in, ret);
                return (EXIT_FAILURE);
        }

        emap_pointdb_apply(&pdb, collect_minima);
        emap_pointdb_free(&pdb);

        if (path_in != NULL)
                free(path_in);
        if (path_out != NULL)
                free(path_out);

        free(opt_comment);

        return (EXIT_SUCCESS);
}

static bool local_minimum_p(emap_pointdb_t *pdb, emap_point_t *p, emap_point_t *n)
{
        return (p->y <= n->y ? true : false);
}

static void collect_minima(emap_pointdb_t *pdb, emap_point_t *p)
{
        bool local_minimum;

        local_minimum = emap_pointnb_applyp(pdb, p, local_minimum_p);
#ifndef NDEBUG
        if (local_minimum) {
                int n;
                fprintf(stderr, "DEBUG: found local minimum y=%f at (", p->y);
                for (n = 0; n < pdb->arity; ++n)
                        fprintf(stderr, "%f ", p->x[n]);
                fprintf(stderr, ")\n");
        }
#endif
}
