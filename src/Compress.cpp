

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
    int ThisIsTheFirstMarker;
    
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

static void CompressAndFlushChunk(args_t *args)
{
    Compressor <vector<string> > ThisChunkCompressor;
    ThisChunkCompressor.CompressChunk(args->Haplotypes);

    int markerIndex=0;
    while(markerIndex<args->Haplotypes[0].length())
    {
        if(ThisChunkCompressor.NewBlockReady())
        {
            if(markerIndex>0)
                markerIndex--;
            ThisChunkCompressor.GetBlockHeader(args->Haplotypes,args->myM3vcfBlock);
            // Copy ploidy information from last read VcfRecord
            
//            abort();
            args->myM3vcfRecordList[markerIndex].copyStartInfotoBlock(args->myM3vcfBlock);
            args->myM3vcfRecordList[ThisChunkCompressor.getBlockHeaderEndPosition()].copyEndInfotoBlock(args->myM3vcfBlock);
            args->outFile.writeBlock(args->myM3vcfBlock);
        }


        ThisChunkCompressor.GetM3vcfRecord(args->myM3vcfRecordList[markerIndex],args->Haplotypes);
        args->outFile.writeRecord(args->myM3vcfRecordList[markerIndex]);

        markerIndex++;
    }
}


static int getNumHaplotypes(args_t *args)
{
    int Count = 0;
    for (int i = 0; i<(args->myVcfRecord.getNumSamples()); i++) Count+=args->myVcfRecord.getNumGTs(i);
    return Count;
}


static void AnalyseHeader(args_t *args)
{
    
    if (!args->vcfFile.open(args->fname, args->myVcfHeader))error("Failed to open: %s\n", args->fname);
    if(!args->MySubsetSamples.init(args->myVcfHeader, args->sample_include_list, 
                                   NULL, args->sample_exlcude_list)) error("Failed to read subset samples file ...\n");   
     
    args->out_hdr.copyHeader(args->myVcfHeader);
    args->outFile.open(args->output_fname, args->out_hdr, args->output_type==4? InputFile::UNCOMPRESSED : InputFile::GZIP);   
}


static void readThisVcfRecord(args_t *args)
{
    int haplotype_index=0;
    for (int i = 0; i<(args->myVcfHeader.getNumSamples()); i++)
    {
        for (int j = 0; j<args->myVcfRecord.getNumGTs(i); j++)
        {
            if(args->myVcfRecord.getGT(i, j)>=0)
                args->Haplotypes[haplotype_index]+= (char)('0'+args->myVcfRecord.getGT(i, j));
            else
            {
                error("Missing value not supported in input VCF file \n");
            }
            haplotype_index++;
        }
    }   
}     
        
        
static void GetNextChunkRead(args_t *args)
{
    
    // If contiguous boundaries are wanted then put the last variant information back
    // before starting to read the next set of variants.
    if(args->contiguous)
    {
        int TotLength = args->Haplotypes[0].length();
        for (int i = 0; i < args->numHaplotypes; i++)
        {
            char tempVal = args->Haplotypes[i][TotLength-1];
            args->Haplotypes[i].clear();
            args->Haplotypes[i]+=tempVal;
        }
        args->myM3vcfRecordList[0] =  args->myM3vcfRecordList[TotLength-1];
    }
    // If not then just clear the haplotype data and start reading next chunk
    else
    {
        for (int i = 0; i < args->numHaplotypes; i++)
        {
            args->Haplotypes[i].clear();
        }
    }

    
}


static void copyVcfRecordToM3vcfRecordList(args_t *args)
{
    args->myM3vcfRecordList[args->Haplotypes[0].length()-1].copyFromVcfRecord(args->myVcfRecord);
    if(!args->keepInfo)
        args->myM3vcfRecordList[args->Haplotypes[0].length()-1].clearInfo(); 
}


static void InitializeHaplotypeData(args_t *args)
{
    args->numHaplotypes=getNumHaplotypes(args);
    args->Haplotypes.resize(args->numHaplotypes);
}
           
            
static void compress_Data(args_t *args)
{

    // Gather Header Information
    AnalyseHeader(args);
 
    args->myM3vcfRecordList.resize(args->bufferSize);
    
    args->ThisIsTheFirstMarker = 0;
    while (args->vcfFile.readRecord(args->myVcfRecord, &args->MySubsetSamples))
    {
        // If this is the first marker being read, initialize haplotype
        // data based on ploidy information at first marker
        if(args->ThisIsTheFirstMarker++==0) {args->myM3vcfBlock.updatePloidy(args->myVcfRecord); InitializeHaplotypeData(args);}
      
        // After that Read Haplotype Data till bufferSize markers
              
        // Read the Record from VCF File and append into Haplotypes
        readThisVcfRecord(args);
            
	// Copy variant info from this VcfRecord and store them in the 
        // vector myM3vcfRecordList
        copyVcfRecordToM3vcfRecordList(args);
        
        // If bufferSize variants have been read, start compressing them
        if (args->Haplotypes[0].length() >= args->bufferSize)
        {
            // Call function that does the entire compression and streams the M3VCF format out
            CompressAndFlushChunk(args);
            GetNextChunkRead(args);
        }
    }

    // If more than one variant is left compress whatever is left behind
    if(args->Haplotypes[0].length()>1)
        CompressAndFlushChunk(args);
    
}


static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, " About:   Compress VCF files to M3VCF files.  \n");
    fprintf(stderr, " Usage:   m3vcftools compress [options] <A.vcf.gz> \n");
    fprintf(stderr, "\n");
    fprintf(stderr, " Options:\n");
    fprintf(stderr, "       --no-version               do not append version and command line to the header\n");
    fprintf(stderr, "   -S, --samples-file [^]<file>   file of samples to include (or exclude with \"^\" prefix)\n");
    fprintf(stderr, "   -k, --keep-info                Copy INFO column to M3VCF file \n");
    fprintf(stderr, "   -b, --buffer <int>             Number of variants to compress at a time [1000]\n");
    fprintf(stderr, "   -o, --output <file>            Write output to a file [standard output]\n");
    fprintf(stderr, "   -O, --output-type <m|M>        m: compressed M3VCF, M: uncompressed M3VCF [M] \n");
 // fprintf(stderr, "   -n, --non-contiguous           Consecutive blocks do NOT overlap at the boundary \n");
   fprintf(stderr, "\n");
    exit(1);
}



int main_m3vcfcompress(int argc, char *argv[])
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
        {"samples-file",required_argument,NULL,'S'},
        {"non-contiguous",no_argument,NULL,'n'},
        {"keep-info",no_argument,NULL,'k'},
        {"buffer",required_argument,NULL,'b'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"no-version",no_argument,NULL,8},
        {NULL,0,NULL,0}
    };

    while ((c = getopt_long(argc, argv, "b:o:O:kS:",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'b': args->bufferSize = strtol(optarg, 0, 0); break;
            case 'o': args->output_fname = optarg; break;
            case 'n': args->contiguous = false; break;
            case 'k': args->keepInfo = true; break;
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
    if ( args->bufferSize<10 )
    {
        error("Invalid buffer (must be greater than 10): %d\n", args->bufferSize);
    }

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

