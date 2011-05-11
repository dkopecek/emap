#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libgen.h>
#include "emap.h"
#include "helpers.h"
#include "pointdb.h"
#include "POI.h"
#include "PES.h"
#include "DG.h"

#define EMAP_SHORT_OPTIONS "c:y:o:E:hv"

static void fprintu(FILE *stream, char *arg0)
{
        fprintf(stream, "\n");
        fprintf(stream, " Usage: %s [OPTIONS] <input>\n", basename(arg0));
        fprintf(stream, "\n");
        fprintf(stream, "\t-o <file>      Output file or `-' for stdout (default is stdout).\n\n");
        fprintf(stream, "\t-y <n>         Column number representing the dependent variable:\n");
        fprintf(stream, "\t                  0 ... last column (default)\n");
        fprintf(stream, "\t                  n ... n-th column (max. %u)\n\n", (1<<16) - 1);
        fprintf(stream, "\t-E <n>[,...]   Exclude n-th column. Multiple columns may be specified.\n\n");
        fprintf(stream, "\t-c <str>       Comment characters (lines starting wich such characters\n");
        fprintf(stream, "\t               will be skipped; default is \"%s\").\n\n", EMAP_COMMENT_CHARS);
        fprintf(stream, "\t-v             Be verbose.\n");
        fprintf(stream, "\t-h             This help.\n");
        fprintf(stream, "\n");
}

int opt_verbose = 0;

int main(int argc, char *argv[])
{
        char *arg0, *path_in, *path_out, *tok;
        int   opt, ret;
        int   opt_yn = 0;
        char *opt_comment = strdup(EMAP_COMMENT_CHARS);

        emap_pointdb_t pdb;
        POIdb_t POIdb;
        DG_t *dg;
        emap_surface_t *es;

        uint32_t *skip_x_n = NULL, x_n;
        size_t    skip_n   = 0;

        arg0     = argv[0];
        path_in  = NULL;
        path_out = NULL;

        while ((opt = getopt(argc, argv, EMAP_SHORT_OPTIONS)) != -1) {
                switch(opt) {
                case 'o':
                        path_out = strdup(optarg);
                        break;
                case 'y':
                        opt_yn = atoi(optarg);

                        if (opt_yn < 0 || opt_yn >= 1<<16) {
                                fprintf(stderr,
                                        "Invalid value for option `-y': only values from the interval <%u,%u> are allowed\n", 0, (1<<16) - 1);
                                fprintu(stderr, arg0);
                                return (EXIT_FAILURE);
                        }

                        break;
                case 'E':
                        if (skip_x_n != NULL)
                                free(skip_x_n);

                        skip_n   = 0;
                        skip_x_n = NULL;

                        tok = strtok(optarg, ",");

                        do {
                                x_n = atoi(tok);

                                if (x_n <= 0) {
                                        fprintf(stderr,
                                                "Invalid value for option `-E': only integers greater than zero are allowed\n");
                                        fprintu(stderr, arg0);

                                        if (skip_x_n != NULL)
                                                free(skip_x_n);

                                        return (EXIT_FAILURE);
                                }

                                skip_x_n = realloc(skip_x_n, sizeof(uint32_t) * ++skip_n);
                                skip_x_n[skip_n - 1] = x_n - 1;
                        } while((tok = strtok(NULL, ",")) != NULL);

                        break;
                case 'c':
                        opt_comment = strdup(optarg);
                        break;
                case 'h':
                        fprintu(stdout, arg0);
                        return (EXIT_SUCCESS);
                case 'v':
                        opt_verbose = 1;
                        setbuf(stdout, NULL);
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

        if (opt_verbose) {
#ifdef _OPENMP
                printf("[i] Compiled with OpenMP support\n");
#endif
                printf("[i] Loading data points from \"%s\"... ", path_in);
        }

        ret = emap_pointdb_load(&pdb, path_in, opt_yn, skip_x_n, skip_n, opt_comment);

        if (ret != EMAP_SUCCESS) {
                if (opt_verbose)
                        fprintf(stdout, "FAILED\n");
                fprintf(stderr, "Unable to load data points from \"%s\": ret=%d\n", path_in, ret);
                return (EXIT_FAILURE);
        }

        if (opt_verbose) {
                printf("OK\n");
                printf("[i] Loaded %zu data points, # of independent variables is %u\n", pdb.count, pdb.arity);
                printf("[i] Generating index... ");
        }

        ret = emap_pointdb_index(&pdb);

        if (ret != EMAP_SUCCESS) {
                if (opt_verbose)
                        fprintf(stdout, "FAIL\n");
                fprintf(stderr, "Unable to index the data points from \"%s\": ret=%d\n", path_in, ret);
                return (EXIT_FAILURE);
        }

        if (opt_verbose) {
                int a;
                printf("OK\n");
                printf("[i] Key component boundaries:\n");
                for (a = 0; a < pdb.arity; ++a)
                        printf("    x[%u] ... <0,%u>\n", a, pdb.keymax[a]);
        }

#ifdef _OPENMP
        if (pthread_mutex_init(&POIdb.mutex, NULL) != 0) {
                fprintf(stderr, "ERROR: pthread_mutex_init failed\n");
                return (EXIT_FAILURE);
        }
#endif
        POIdb.cntinc = sysconf(_SC_PAGESIZE) / sizeof(emap_point_t *);

        POIdb.palloc = POIdb.cntinc;
        POIdb.pcount = 0;
        POIdb.points = alloc_array(emap_point_t *, POIdb.palloc);

        POIdb.minima = NULL;
        POIdb.mcount = 0;

        POIdb.transp = NULL;
        POIdb.tcount = 0;

        emap_pointdb_apply_r(&pdb, collect_POI, &POIdb);

        if (opt_verbose) {
                printf("[i] Found %zu POIs\n"
                       "       local minima: %zu\n"
                       "  transition points: %zu\n"
                       "            overlap: %s\n",
                       POIdb.pcount, POIdb.mcount, POIdb.tcount,
                       POIdb.tcount + POIdb.mcount > POIdb.pcount ? "yes" : "no");
        }

        es = POI_postprocess(&pdb, &POIdb);

        if (es == NULL) {
                fprintf(stderr, "ERROR: post-processing of POIs failed\n");
                return (EXIT_FAILURE);
        }

        dg = DG_create(es);

        if (dg == NULL) {
                fprintf(stderr, "ERROR: unable to create a disconnectivity graph from the POI set\n");
                return (EXIT_FAILURE);
        }

        if (DG_write(dg, path_out)) {
                fprintf(stderr, "ERROR: unable to write the disconnectivity graph to \"%s\"\n", path_out);
                return (EXIT_FAILURE);
        }

        /*
         * Cleanup & exit
         */
        emap_pointdb_free(&pdb);

#ifdef _OPENMP
        if (pthread_mutex_destroy(&POIdb.mutex) != 0) {
                fprintf(stderr, "ERROR: pthread_mutex_destroy failed\n");
                return (EXIT_FAILURE);
        }
#endif
        if (path_in != NULL)
                free(path_in);
        if (path_out != NULL)
                free(path_out);

        free(opt_comment);

        return (EXIT_SUCCESS);
}
