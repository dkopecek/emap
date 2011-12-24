#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <inttypes.h>
#include <libgen.h>
#include "emap.h"
#include "helpers.h"
#include "pointdb.h"
#include "POI.h"
#include "PES.h"
#include "DG.h"

#define EMAP_SHORT_OPTIONS "c:y:s:o:O:E:t:T:m:nhv"

static void fprintu(FILE *stream, char *arg0)
{
        fprintf(stream, "\n");
        fprintf(stream, " Usage: %s [OPTIONS] <input>\n", basename(arg0));
        fprintf(stream, "\n");
        fprintf(stream, "\t-o <file>      Output file (dot/graphviz) or `-' for stdout (default is stdout).\n");
	fprintf(stream, "\t-O <file>      Output file (emap) or `-' for stdout (default is stdout).\n\n");
        fprintf(stream, "\t-y <n>         Column number representing the dependent variable:\n");
        fprintf(stream, "\t                  0 ... last column (default)\n");
        fprintf(stream, "\t                  n ... n-th column (max. %u)\n\n", (1<<16) - 1);
	fprintf(stream, "\t-T <str>       Transform the dependent variables using a function:\n");
	fprintf(stream, "\t                  log   ... logarithm to the base e\n");
	fprintf(stream, "\t                  log2  ... logarithm to the base 2\n");
	fprintf(stream, "\t                  log10 ... logarithm to the base 10\n\n");
        fprintf(stream, "\t-s <n>         Height of the energy band used in building the DG.\n");
        fprintf(stream, "\t               The value should be a floating point number larger\n");
        fprintf(stream, "\t               than 0.\n\n");
        fprintf(stream, "\t-E <n>[,...]   Exclude n-th column. Multiple columns may be specified.\n\n");
        fprintf(stream, "\t-c <str>       Comment characters (lines starting wich such characters\n");
        fprintf(stream, "\t               will be skipped; default is \"%s\").\n\n", EMAP_COMMENT_CHARS);
        fprintf(stream, "\t-m <file>      Dump found local minima into `file'.\n");
        fprintf(stream, "\t-t <file>      Dump transition points into `file'.\n");
        fprintf(stream, "\t-n             Skip superbasin analysis and, therefore, building of the DG.\n");
        fprintf(stream, "\t-v             Be verbose.\n");
        fprintf(stream, "\t-h             This help.\n");
        fprintf(stream, "\n");
}

int opt_verbose = 0;

#include <time.h>
#define pI(fmt, ...) printf("[i: %.2fs] "fmt, ((double)clock()/(double)(CLOCKS_PER_SEC)), ##__VA_ARGS__)

int main(int argc, char *argv[])
{
        char *arg0, *path_in, *path_out, *tok, *path_out_emap;
        int   opt, ret;
        int   opt_yn = 0;
        char *opt_comment = strdup(EMAP_COMMENT_CHARS);

        emap_pointdb_t pdb;
        POIdb_t POIdb;
        DG_t *dg;
        emap_surface_t *es;

        uint32_t *skip_x_n = NULL;
	int32_t x_n;
        size_t    skip_n   = 0;

        emap_float dE = 0.0;
        bool skip_DG = false;
        char *dump_mfile = NULL, *dump_tfile = NULL;
	enum emap_transfn opt_ytransform = EMAP_TRANSFORM_NONE;
	char *opt_ytransform_str = "none";

        arg0     = argv[0];
        path_in  = NULL;
        path_out = path_out_emap = NULL;

        while ((opt = getopt(argc, argv, EMAP_SHORT_OPTIONS)) != -1) {
                switch(opt) {
                case 'o':
                        path_out = strdup(optarg);
                        break;
		case 'O':
			path_out_emap = strdup(optarg);
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
                case 's':
                        dE = emap_strtoflt(optarg, NULL);

                        if (dE <= 0) {
                                fprintf(stderr, "Invalid value for option `-s': only a value larger than 0 is allowed\n");
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
                case 'n':
                        skip_DG = true;
                        break;
                case 'm':
                        dump_mfile = optarg;
                        break;
                case 't':
                        dump_tfile = optarg;
                        break;
                case 'T':
			if (strcmp(optarg, "log") == 0)
				opt_ytransform = EMAP_TRANSFORM_LOG;
			else if (strcmp(optarg, "log2") == 0)
				opt_ytransform = EMAP_TRANSFORM_LOG2;
			else if (strcmp(optarg, "log10") == 0)
				opt_ytransform = EMAP_TRANSFORM_LOG10;
			else {
				fprintf(stderr, "Unsupported transformation: %s\n", optarg);
				fprintu(stderr, arg0);
				return (EXIT_FAILURE);
			}

			opt_ytransform_str = strdup(optarg);
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

        if (dE <= 0.0 && !skip_DG) {
                fprintf(stderr, "Required option `-s' not specified!\n");
                fprintu(stderr, arg0);
                return (EXIT_FAILURE);
        }

        path_in = strdup(argv[0]);

        emap_pointdb_init(&pdb);

        if (opt_verbose) {
#ifdef _OPENMP
                pI("Compiled with OpenMP support\n");
#endif
                pI("Loading data points from \"%s\"...\n", path_in);

		if (opt_ytransform != EMAP_TRANSFORM_NONE) {
			pI("Note: Using %s transformation for the dependent variable\n", opt_ytransform_str);
			free(opt_ytransform_str);
		}
        }

	pdb.y_trans = opt_ytransform;
        ret = emap_pointdb_load(&pdb, path_in, opt_yn, skip_x_n, skip_n, opt_comment);

        if (ret != EMAP_SUCCESS) {
                if (opt_verbose)
                        fprintf(stdout, "FAILED\n");
                fprintf(stderr, "Unable to load data points from \"%s\": ret=%d\n", path_in, ret);
                return (EXIT_FAILURE);
        }

        if (opt_verbose) {
                pI("OK\n");
                pI("Loaded %zu data points, # of independent variables is %u\n", pdb.count, pdb.arity);
                pI("Generating index...\n");
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
                pI("OK\n");
                pI("Key component boundaries:\n");
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

        POIdb.transp  = NULL;
        POIdb.tcount  = 0;

        emap_pointdb_apply_r(&pdb, collect_POI, &POIdb);

        if (opt_verbose) {
                pI("Found %zu POIs\n"
                   "       local minima candidates: %zu\n"
                   "  transition points candidates: %zu\n"
                   "                       overlap: %s\n",
		   POIdb.pcount, POIdb.mcount, POIdb.tcount,
		   POIdb.tcount + POIdb.mcount > POIdb.pcount ? "yes" : "no");
                pI("Preparing an abstract surface for superbasin analysis...\n");
        }

        es = POI_postprocess(&pdb, &POIdb, opt_verbose);

        if (es == NULL) {
                if (opt_verbose)
                        printf("FAILED\n");
                fprintf(stderr, "ERROR: post-processing of POIs failed\n");
                return (EXIT_FAILURE);
        }

        if (opt_verbose) {
                pI("OK\n");
                pI("The abstract surface has a total of %zu points\n"
		   "       local minima: %zu\n"
		   "  transition points: %zu\n",
		   es->mcount + es->tcount, es->mcount, es->tcount);
        }

        if (dump_tfile != NULL) {
                FILE *fp;
                size_t i, a;

                fp = fopen(dump_tfile, "w");

                /* write: x0 x1 ... xn y */
                for (i = 0; i < es->tcount; ++i) {
                        fprintf(fp, "%"PRIu32" ", es->tpoint[i]->cmaximum->line);

                        for (a = 0; a < pdb.arity; ++a)
                                fprintf(fp, "%"EMAP_FLTFMT" ", es->tpoint[i]->cmaximum->x[a]);
                        fprintf(fp, "%"EMAP_FLTFMT"\n", es->tpoint[i]->cmaximum->y);
                }

                fclose(fp);
        }

        if (dump_mfile != NULL) {
                FILE *fp;
                size_t i, a;

                fp = fopen(dump_mfile, "w");

                for (i = 0; i < es->mcount; ++i) {
                        fprintf(fp, "%"PRIu32" ", es->mpoint[i]->cminimum->line);

                        for (a = 0; a < pdb.arity; ++a)
                                fprintf(fp, "%"EMAP_FLTFMT" ", es->mpoint[i]->cminimum->x[a]);
                        fprintf(fp, "%"EMAP_FLTFMT"\n", es->mpoint[i]->cminimum->y);
                }

                fclose(fp);
        }

        if (!skip_DG) {
                if (opt_verbose)
                        pI("Building the DG with dE set to %"EMAP_FLTFMT"... \n", dE);

                dg = DG_create(&pdb, es, dE, opt_verbose ? true : false);

                if (dg == NULL) {
                        fprintf(stderr, "ERROR: unable to create a disconnectivity graph from the abstract surface\n");
                        return (EXIT_FAILURE);
                }

                if (DG_write(dg, path_out)) {
                        fprintf(stderr, "ERROR: unable to write the disconnectivity graph to \"%s\"\n", path_out);
                        return (EXIT_FAILURE);
                }

                if (DG_write_emap(dg, path_out_emap)) {
                        fprintf(stderr, "ERROR: unable to write the disconnectivity graph to \"%s\"\n", path_out_emap);
                        return (EXIT_FAILURE);
                }
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
        if (path_out_emap != NULL)
                free(path_out_emap);

        free(opt_comment);

	if (opt_verbose)
		pI("\nDone.\n");

        return (EXIT_SUCCESS);
}
