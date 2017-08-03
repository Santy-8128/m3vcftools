

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <iostream>
#include "m3vcfHeader.h"
#include "m3vcfBlock.h"
#include "m3vcfFileWriter.h"
#define VCF 2
#define ZIPVCF 1
#define M3VCF 4
#define ZIPM3VCF 3



using namespace std;


struct args_t
{
//    bcf_srs_t *files;
//    htsFile *out_fh;
    int output_type, n_threads, record_cmd_line;
    m3vcfHeader out_hdr;
    m3vcfFileWriter <m3vcfHeader> outFile;


//    int *seen_seq;
//
//    // phasing
//    int *start_pos, start_tid, ifname;
//    int *swap_phase, nswap, *nmatch, *nmism;
//    bcf1_t **buf;
//    int nbuf, mbuf, prev_chr, min_PQ, prev_pos_check;
//    int32_t *GTa, *GTb, mGTa, mGTb, *phase_qual, *phase_set;
//
    vector<int> start_pos;
    char **argv, **fnames;
    const char *output_fname, *file_list; //, , , , *remove_dups, *regions_list;
    int phased_concat, nfnames, argc; //, , allow_overlaps, , regions_is_file;
//    int compact_PS, phase_set_changed, naive_concat;
};



void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(-1);
}

static void init_data(args_t *args)
{
    if ( args->phased_concat )
    {
        args->start_pos.resize(args->nfnames,0);
    }

    // Gather and Merge Header Information
    for (int i=0; i<args->nfnames; i++)
    {
        IFILE m3vcfxStream = ifopen(args->fnames[i], "r"); if ( !m3vcfxStream ) error("Failed to open: %s\n", args->fnames[i]);
        m3vcfHeader myHeader;
        myHeader.read(m3vcfxStream);


        if (!args->out_hdr.checkMergeHeader(myHeader) )
            error("Different sample information in %s.\n", args->fnames[i]);
        else
            args->out_hdr.mergeHeader(myHeader);

        if ( args->phased_concat )
        {
             m3vcfBlock CurrentBlock;
             CurrentBlock.read(m3vcfxStream,myHeader,true);
             args->start_pos[i] = CurrentBlock.getStartBasePosition();
        }
        ifclose(m3vcfxStream);
    }


    if (args->record_cmd_line) args->out_hdr.appendMetaLine("m3vcftools_concat");
    args->outFile.open(args->output_fname,args->out_hdr);

    args->outFile.writeHeader(args->out_hdr);

    IFILE firstFile = NULL, secondFile = NULL;
    m3vcfBlock firstFileBlock, secondFileBlock;
    m3vcfRecord firstFileRecord, secondFileRecord;
    m3vcfHeader myHeader;

    int i=0;
    firstFile = ifopen(args->fnames[i], "r"); if ( !firstFile ) error("Failed to open: %s\n", args->fnames[i]);
    myHeader.read(firstFile);

    while(firstFileBlock.read(firstFile,args->out_hdr) && firstFileBlock.getEndBasePosition()<=args->start_pos[1])
    {
        args->outFile.writeBlock(firstFileBlock);
        while(!firstFileBlock.isBlockFinished())
        {
//            abort();
            firstFileRecord.read(firstFile,firstFileBlock);
            args->outFile.writeRecord(firstFileRecord);
        }
    }




//    CurrentBlock.write()
//
//
//       m3vcfBlock CurrentBlock;
//             CurrentBlock.read(m3vcfxStream,myHeader,true);
//       m3vcfBlock CurrentBlock;
//             CurrentBlock.read(m3vcfxStream,myHeader,true);


//    IFILE m3vcfxStream = ifopen(args->fnames[0], "r");
//
//    m3vcfHeader myHeader;
//    m3vcfBlock CurrentBlock;
//
//    myHeader.read(m3vcfxStream);
//    while(CurrentBlock.read(m3vcfxStream,myHeader))
//    {
//
//         cout<<" My = "<<CurrentBlock.getChromStr()<<" == THIS" <<endl;
//        cout<<" My = "<<CurrentBlock.getStartBasePosition()<<" == THIS" <<endl;
//        cout<<" My = "<<CurrentBlock.getEndBasePosition()<<" == THIS" <<endl;
//        cout<<" My = "<<CurrentBlock.getNumMarkers()<<" == THIS" <<endl;
//        cout<<" My = "<<CurrentBlock.getNumUniqueReps()<<" == THIS" <<endl<<endl;
//
////abort();
//    }
//



}


static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, " About:   Concatenate or combine M3VCF files. All source files must have the same sample\n");
    fprintf(stderr, "          columns appearing in the same order. The program can be used, for example, to\n");
    fprintf(stderr, "          ligate haplotypes \n");
    fprintf(stderr, " Usage:   bcftools concat [options] <A.vcf.gz> [<B.vcf.gz> [...]]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, " Options:\n");
    fprintf(stderr, "   -f, --file-list <file>         Read the list of files from a file.\n");
    fprintf(stderr, "   -l, --ligate                   Ligate phased VCFs by matching phase at overlapping haplotypes\n");
    fprintf(stderr, "       --no-version               do not append version and command line to the header\n");
    fprintf(stderr, "   -o, --output <file>            Write output to a file [standard output]\n");
    fprintf(stderr, "   -O, --output-type <z|v|m|M>    z: compressed VCF,   v: uncompressed VCF, \n");
    fprintf(stderr, "                                  m: compressed M3VCF, M: uncompressed M3VCF [M]\n");
    fprintf(stderr, "       --threads <int>            Number of extra output compression threads [0]\n");
    fprintf(stderr, "\n");
    exit(1);
}




int main_m3vcfconcat(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->output_fname = "-";
    args->output_type = M3VCF;
    args->n_threads = 0;
    args->record_cmd_line = 1;

    static struct option loptions[] =
    {

        {"ligate",no_argument,NULL,'l'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"file-list",required_argument,NULL,'f'},
        {"no-version",no_argument,NULL,8},
        {NULL,0,NULL,0}
    };

    while ((c = getopt_long(argc, argv, "lo:O:f:",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'l': args->phased_concat = 1; break;
            case 'f': args->file_list = optarg; break;
            case 'o': args->output_fname = optarg; break;
            case 'O':
                switch (optarg[0]) {
                    case 'v': args->output_type = VCF; break;
                    case 'z': args->output_type = ZIPVCF; break;
                    case 'm': args->output_type = ZIPM3VCF; break;
                    case 'M': args->output_type = M3VCF; break;
                    default: error("The output type \"%s\" not recognised\n", optarg);
                };
                break;
            case  9 : args->n_threads = strtol(optarg, 0, 0); break;
            case  8 : args->record_cmd_line = 0; break;
            case '?': usage(args); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    while ( optind<argc )
    {
        args->nfnames++;
        args->fnames = (char **)realloc(args->fnames,sizeof(char*)*args->nfnames);
        args->fnames[args->nfnames-1] = strdup(argv[optind]);
        optind++;
    }
    if ( !args->nfnames ) usage(args);

    init_data(args);
    cout<<" WOOOWWWWWWWWWWWW \n";
    return 0;
}

