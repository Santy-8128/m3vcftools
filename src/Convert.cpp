

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <iostream>
#include "VcfFileReader.h"
#include "m3vcfHeader.h"
#include "m3vcfBlock.h"
#include "m3vcfFileWriter.h"
#include "Unique.h"
#define VCF 2
#define ZIPVCF 1
#define M3VCF 4
#define ZIPM3VCF 3
using namespace std;


typedef struct _args_t
{

    //VCF File Variables
    VcfFileReader vcfFile;
    VcfHeader myVcfHeader;
    VcfRecord myVcfRecord;
    VcfSubsetSamples MySubsetSamples;

    // M3VCF File Variables

    m3vcfHeader out_hdr;
    m3vcfFileWriter <m3vcfHeader> outFile;
    vector<m3vcfRecord> myM3vcfRecordList;
    m3vcfBlock myM3vcfBlock;
    bool contiguous, keepInfo;
    int numHaplotypes;
    
    // Haplotype Data
    vector<string> Haplotypes;int bufferSize;

    //Common File and Argument Variables
    char **argv;
    const char *output_fname, *fname;
    int argc, output_type, record_cmd_line;
    char *sample_include_list, *sample_exlcude_list;

}
args_t;

static void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(-1);
}

         


static void compress_Data(args_t *args)
{
    
    IFILE firstFile = NULL, secondFile = NULL;
    m3vcfBlock firstFileBlock, secondFileBlock;
    m3vcfRecord firstFileRecord, secondFileRecord;
    m3vcfHeader myHeader;

    int i=0;
    firstFile = ifopen(args->fname, "r"); if ( !firstFile ) error("Failed to open: %s\n", args->fname);
    myHeader.read(firstFile);
////    args->out_hdr.mergeHeader(myHeader);
//      
    args->outFile.open(args->output_fname,myHeader);

//    while(firstFileBlock.read(firstFile,myHeader))
//    {
//        args->outFile.writeBlock(firstFileBlock);
//        while(!firstFileBlock.isBlockFinished())
//        {
////            abort();
//            firstFileRecord.read(firstFile,firstFileBlock);
//            args->outFile.writeRecord(firstFileRecord);
//        }
//    }
    
    

}


static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, " About:   Convert M3VCF file to VCF file  \n");
    fprintf(stderr, " Usage:   m3vcftools compress [options] <A.m3vcf.gz> \n");
    fprintf(stderr, "\n");
    fprintf(stderr, " Options:\n");
    fprintf(stderr, "       --no-version               do not append version and command line to the header\n");
//    fprintf(stderr, "   -S, --samples-file [^]<file>   file of samples to include (or exclude with \"^\" prefix)\n");
//    fprintf(stderr, "   -k, --keep-info                Copy INFO column to M3VCF file \n");
//    fprintf(stderr, "   -b, --buffer <int>             Number of variants to compress at a time [1000]\n");
    fprintf(stderr, "   -o, --output <file>            Write output to a file [standard output]\n");
    fprintf(stderr, "   -O, --output-type <z|v>        z: compressed VCF, v: uncompressed VCF [v] \n");
 // fprintf(stderr, "   -n, --non-contiguous           Consecutive blocks do NOT overlap at the boundary \n");
   fprintf(stderr, "\n");
    exit(1);
}



int main_m3vcfconvert(int argc, char *argv[])
{
    int c;
    args_t* args = new args_t();
    args->argc    = argc; args->argv = argv;
    args->output_fname = "-";
    args->output_type = M3VCF;
    args->bufferSize = 1000;
    args->record_cmd_line = 1;
    args->keepInfo = false;
    args->contiguous = true;
    args->sample_include_list = NULL;
    args->sample_exlcude_list = NULL;
    
    static struct option loptions[] =
    {
//        {"samples-file",required_argument,NULL,'S'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"no-version",no_argument,NULL,8},
        {NULL,0,NULL,0}
    };

    while ((c = getopt_long(argc, argv, "b:o:O:kS:",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'o': args->output_fname = optarg; break;
            case 'S': args->sample_include_list = optarg; break;
            case 'O':
                switch (optarg[0]) {
                    case 'm': args->output_type = ZIPM3VCF; break;
                    case 'M': args->output_type = M3VCF; break;
                    default: error("The output type \"%s\" not recognised\n", optarg);
                };
                break;
            case  8 : args->record_cmd_line = 0; break;
            case '?': usage(args); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }


    args->fname = NULL;
    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->fname = "-";  // reading from stdin
        else usage(args);
    }
    else args->fname = argv[optind];
 
    if(args->sample_include_list!=NULL)
    {
        if(args->sample_include_list[0]=='^')
        {
            args->sample_exlcude_list=&args->sample_include_list[1];
            args->sample_include_list=NULL;
        }
        
    }
    compress_Data(args);
    return 0;
}

