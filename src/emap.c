#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <libgen.h>

#define EMAP_SHORT_OPTIONS "o:e:hv"

static void fprintu(FILE *stream, char *arg0)
{
        fprintf(stream, "\n");
        fprintf(stream, " Usage: %s [OPTIONS] <input>\n", basename(arg0));
        fprintf(stream, "\n");
        fprintf(stream, "\t-o <file>      Output file or `-' for stdout (default is stdout).\n");
        fprintf(stream, "\t-e <n>         Column number representing the dependent variable:\n");
        fprintf(stream, "\t                  0 ... last column (default)\n");
        fprintf(stream, "\t                  n ... n-th column\n");
        fprintf(stream, "\t-v             Be verbose.\n");
        fprintf(stream, "\t-h             This help.\n");
        fprintf(stream, "\n");
}

int main(int argc, char *argv[])
{
        char *arg0, *path_in, *path_out;
        int   opt;

        arg0     = argv[0];
        path_in  = NULL;
        path_out = NULL;

        while ((opt = getopt(argc, argv, EMAP_SHORT_OPTIONS)) != -1) {
                switch(opt) {
                case 'o':
                        break;
                case 'h':
                        fprintu(stdout, arg0);
                        return (EXIT_SUCCESS);
                case 'v':
                        break;
                default:
                        fprintf(stderr, "Unknown option `%c'\n", optopt);
                        fprintu(stderr, arg0);
                        return (EXIT_FAILURE);
                }
        }

        argc -= optind;
        argv -= optind;

        if (argc != 1) {
                fprintf(stderr, "Expecting exactly one input file.\n");
                fprintu(stderr, arg0);
                return (EXIT_FAILURE);
        }

        path_in = argv[0];

        return (EXIT_SUCCESS);
}
