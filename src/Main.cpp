#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>



static void usage(FILE *fp);
//const char* createCommandLine(convert_args_t *args, const char *optionName); 
int main_m3vcfconcat(int argc, char *argv[]);
int main_m3vcfcompress(int argc, char *argv[]);
int main_m3vcfconvert(int argc, char *argv[]);


typedef struct
{
    int (*func)(int, char*[]);
    const char *alias, *help;
}
cmd_t;

static cmd_t cmds[] =
{
    { .func  = NULL,
      .alias = "M3VCF manipulation",
      .help  = NULL
    },

    { .func  = main_m3vcfcompress, //main_vcfannotate,
      .alias = "compress",
      .help  = "compress VCF file to M3VCF",
    },
    { .func  = main_m3vcfconcat, //main_vcfconcat,
      .alias = "concat",
      .help  = "concatenate M3VCF files from the same set of samples"
    },
    { .func  = main_m3vcfconvert, //main_vcfconvert,
      .alias = "convert",
      .help  = "convert M3VCF files to different formats and back"
    },
    { .func  = NULL,
      .alias = NULL,
      .help  = NULL
    }

};

const char *m3vcftools_version(void)
{
    return VERSION;
}

int main(int argc, char **argv)
{
    if (argc < 2) { usage(stderr); return 1; }
    if (strcmp(argv[1], "version") == 0 || strcmp(argv[1], "--version") == 0 || strcmp(argv[1], "-v") == 0) {
        printf("m3vcftools %s\n", m3vcftools_version());
        printf("This is free software: you are free to change and redistribute it.\nThere is NO WARRANTY, to the extent permitted by law.\n");
        return 0;
    }
    else if (strcmp(argv[1], "--version-only") == 0) {
        printf("%s\n", m3vcftools_version());
        return 0;
    }
    else if (strcmp(argv[1], "help") == 0 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) {
        if (argc == 2) { usage(stdout); return 0; }
        // Otherwise change "m3vcftools help COMMAND [...]" to "m3vcftools COMMAND";
        // main_xyz() functions by convention display the subcommand's usage
        // when invoked without any arguments.
        argv++;
        argc = 2;
    }
    int i = 0;
    while (cmds[i].alias)
    {
        if (cmds[i].func && strcmp(argv[1],cmds[i].alias)==0)
        {
            return cmds[i].func(argc-1,argv+1);
        }
        i++;
    }
    fprintf(stderr, "[ERROR:] unrecognized command '%s'\n", argv[1]);
    return 1;
}





static void usage(FILE *fp)
{
    fprintf(fp, "\n");
    fprintf(fp, " -------------------------------------------------------------------------------- \n");
	fprintf(fp, "                  m3vcftools - A Tool for Manipulating M3VCF Files\n");
	fprintf(fp, " --------------------------------------------------------------------------------\n");
    fprintf(fp, "\n (c) 2017 - Sayantan Das, Goncalo Abecasis \n\n");


    fprintf(fp,  " Version: %s\n", m3vcftools_version());
    fprintf(fp,  " Built  : %s\n", DATE);
    fprintf(fp,  " User   : %s\n", USER);
    fprintf(fp,  " URL    : http://genome.sph.umich.edu/wiki/m3vcftools\n");


    fprintf(fp,"\n Usage  : m3vcftools [--version|--version-only] [--help] <command> <argument>\n");
    fprintf(fp, "\n");
    fprintf(fp, " Commands:\n");

    int i = 0;
    const char *sep = NULL;
    while (cmds[i].alias)
    {
        if ( !cmds[i].func ) sep = cmds[i].alias;
        if ( sep )
        {
            fprintf(fp, "\n -- %s\n", sep);
            sep = NULL;
        }
        if ( cmds[i].func && cmds[i].help[0]!='-' ) fprintf(fp, "    %-12s %s\n", cmds[i].alias, cmds[i].help);
        i++;
    }
    fprintf(fp,"\n");
    fprintf(fp,"\n");
}

