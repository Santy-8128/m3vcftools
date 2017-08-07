

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <iostream>
#include "m3vcfHeader.h"
#include "m3vcfBlockHeader.h"
#include "m3vcfFileWriter.h"
#include "m3vcfBlock.h"
#define MAX_BASE_POSITION 999999999
#define VCF 2
#define ZIPVCF 1
#define M3VCF 4
#define ZIPM3VCF 3



using namespace std;


struct concat_args_t
{
//    bcf_srs_t *files;
//    htsFile *out_fh;
    int output_type, n_threads, record_cmd_line;
    m3vcfHeader out_hdr;
    m3vcfFileWriter <m3vcfHeader> outFile;


    m3vcfBlockHeader firstFileBlockHeader, secondFileBlockHeader;
    m3vcfRecord firstFileRecord, secondFileRecord;
    m3vcfHeader myCommonHeader;
    vector<m3vcfBlock*> FirstFileBlocks, SecondFileBlocks;


    // File Handling Variables
    IFILE firstFile, secondFile;



    // Variables for frameworking
    vector<int> start_pos;


    //Common File and Argument Variables
    char **argv, **fnames;
    const char *output_fname, *file_list;
    int phased_concat, nfnames, argc;

};



void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(-1);
}

static void Initialize(concat_args_t *args)
{
    if ( args->phased_concat )
    {
        args->start_pos.resize(args->nfnames+1,0);
    }
    args->firstFile = NULL, args->secondFile = NULL;

}


static void AnalyseHeader(concat_args_t *args)
{
    // Gather and Merge Header Information
    int i;
    for (i=0; i<args->nfnames; i++)
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
            m3vcfBlockHeader CurrentBlock;
            CurrentBlock.read(m3vcfxStream,myHeader,true);
            args->start_pos[i] = CurrentBlock.getStartBasePosition();
        }
        ifclose(m3vcfxStream);
    }
    args->start_pos[i]=MAX_BASE_POSITION;

    if (args->record_cmd_line) args->out_hdr.appendMetaLine("m3vcftools_concat");
    args->outFile.open(args->output_fname,args->out_hdr);
}


static void FlushFirstFileNonOverlappingPart(concat_args_t *args)
{
    args->firstFile = ifopen(args->fnames[0], "r"); if ( !args->firstFile ) error("Failed to open: %s\n", args->fnames[0]);
    args->myCommonHeader.read(args->firstFile);

    while(args->firstFileBlockHeader.read(args->firstFile,args->out_hdr) && args->firstFileBlockHeader.getEndBasePosition()<=args->start_pos[1])
    {
        args->outFile.writeBlock(args->firstFileBlockHeader);
        while(!args->firstFileBlockHeader.isBlockFinished())
        {
            args->firstFileRecord.read(args->firstFile,args->firstFileBlockHeader);
            args->outFile.writeRecord(args->firstFileRecord);
        }
    }

}


static void FlushFirstFileBlockNonOverlappingPart(int chunkNo, concat_args_t *args)
{
    args->outFile.writeBlock(args->firstFileBlockHeader);
    args->firstFileRecord.read(args->firstFile,args->firstFileBlockHeader);
    args->outFile.writeRecord(args->firstFileRecord);

    while(!args->firstFileBlockHeader.isBlockFinished() && args->firstFileRecord.getBasePosition() < args->secondFileRecord.getBasePosition())
    {
        args->firstFileRecord.read(args->firstFile,args->firstFileBlockHeader);
        args->outFile.writeRecord(args->firstFileRecord);

    }

    if(args->firstFileRecord.IsMatching(args->secondFileRecord) == 0)
        error("[ERROR:] Variant [%s] present only in one chunk [%s]\n", args->firstFileRecord.PrintVariant().c_str(),args->fnames[chunkNo]);

}



static void InitSecondFile(int chunkNo, concat_args_t *args)
{
    args->secondFile = ifopen(args->fnames[chunkNo], "r"); if ( !args->secondFile ) error("Failed to open: %s\n", args->fnames[chunkNo]);
    args->myCommonHeader.read(args->secondFile);
}



static void CheckIfOverlaps(int chunkNo, concat_args_t *args)
{
    if(args->firstFileBlockHeader.isBlockFinished())
    {
        if(!args->firstFileBlockHeader.read(args->firstFile,args->out_hdr))
            error("[ERROR:] Only one variant overlaps between chunk [%s] and [%s]\n", args->fnames[chunkNo-1], args->fnames[chunkNo] );
        else
            args->outFile.writeBlock(args->firstFileBlockHeader);
    }

    if(args->secondFileBlockHeader.isBlockFinished())
        if(!args->secondFileBlockHeader.read(args->secondFile,args->out_hdr))
            error("[ERROR:] No overlap between chunk [%s] and [%s]\n", args->fnames[chunkNo-1], args->fnames[chunkNo] );

}

static void ReadOverlappingRegion(int chunkNo, concat_args_t *args)
{
    while(!args->firstFileBlockHeader.isBlockFinished() && !args->secondFileBlockHeader.isBlockFinished() )
    {
        args->firstFileRecord.read(args->firstFile,args->firstFileBlockHeader);
        args->outFile.writeRecord(args->firstFileRecord);
        args->secondFileRecord.read(args->secondFile,args->secondFileBlockHeader);

        if(args->firstFileBlockHeader.isBlockFinished())
        {
            if(!args->firstFileBlockHeader.read(args->firstFile,args->out_hdr))
                break;
            else
                args->outFile.writeBlock(args->firstFileBlockHeader);
        }

        if(args->secondFileBlockHeader.isBlockFinished())
            if(!args->secondFileBlockHeader.read(args->secondFile,args->out_hdr))
                error("[ERROR:] Chunk [%s] ends before the chunk preceding it [%s] ", args->fnames[chunkNo-1], args->fnames[chunkNo] );
    }
}


static void SaveFirstFileOverlappingRegion(int chunkNo, concat_args_t *args)
{
    args->FirstFileBlocks.clear();
    args->FirstFileBlocks.resize(1);
    args->FirstFileBlocks[0]=new m3vcfBlock;
    args->FirstFileBlocks[0]->CopyBlockHeader(args->firstFileBlockHeader);

    do
    {
        while(!args->FirstFileBlocks.back()->isBlockFinished())
        {
            (args->FirstFileBlocks.back())->readRecord(args->firstFile);
        }

        args->FirstFileBlocks.resize(args->FirstFileBlocks.size()+1);
        args->FirstFileBlocks.back()=new m3vcfBlock;
    }while((args->FirstFileBlocks.back()->readHeader(args->firstFile,args->out_hdr)));

    args->FirstFileBlocks.pop_back();
    ifclose(args->firstFile);
}



static void SaveSecondFileOverlappingRegion(int chunkNo, concat_args_t *args)
{

    int PrevLastPosition = args->FirstFileBlocks.back()->getEndBasePosition();

    InitSecondFile(chunkNo, args);

    args->SecondFileBlocks.clear();
    args->SecondFileBlocks.resize(1);
    args->SecondFileBlocks[0]=new m3vcfBlock;


    while((args->SecondFileBlocks.back()->readHeader(args->secondFile,args->out_hdr))
          && args->SecondFileBlocks.back()->getStartBasePosition()<PrevLastPosition)
    {
        while(!args->SecondFileBlocks.back()->isBlockFinished())
        {
            (args->SecondFileBlocks.back())->readRecord(args->secondFile);
        }

        args->SecondFileBlocks.resize(args->SecondFileBlocks.size()+1);
        args->SecondFileBlocks.back()=new m3vcfBlock;
    }


}

static void FlushSecondFileBlockNonOverlappingPart(int chunkNo, concat_args_t *args)
{
    args->SecondFileBlocks.back()->CopyToBlockHeader(args->secondFileBlockHeader);

    do
    {
        args->outFile.writeBlock(args->secondFileBlockHeader);
        while(!args->secondFileBlockHeader.isBlockFinished())
        {
            args->secondFileRecord.read(args->secondFile,args->secondFileBlockHeader);
            args->outFile.writeRecord(args->secondFileRecord);
        }
    }
    while((args->secondFileBlockHeader.read(args->secondFile,args->out_hdr))
          && args->secondFileBlockHeader.getEndBasePosition()<=args->start_pos[chunkNo+1]);


}


static void FlushOverlappingPart(int chunkNo, concat_args_t *args)
{
    int FirstBlockHeaderIndex=0, FirstRecordIndex=0;

    while(args->FirstFileBlocks[0]->getM3vcfRecord(FirstRecordIndex)->getBasePosition() < args->start_pos[chunkNo]) {FirstRecordIndex++;}


    for(int i=0; i<(int)args->SecondFileBlocks.size()-1; i++)
    {

        for(int j=0; j< args->SecondFileBlocks[i]->getNumMarkers(); j++)
        {

            if(args->SecondFileBlocks[i]->getM3vcfRecord(j)->IsMatching(*args->FirstFileBlocks[FirstBlockHeaderIndex]->getM3vcfRecord(FirstRecordIndex))!=1)
                error("[ERROR:] Variant [%s] present only in one chunk [%s]\n",
                      args->SecondFileBlocks[i]->getM3vcfRecord(j)->getVariantID().c_str(),
                      args->fnames[chunkNo-1]);

            FirstRecordIndex++;
            if(FirstRecordIndex>=args->FirstFileBlocks[FirstBlockHeaderIndex]->getNumMarkers())
            {
                if(FirstBlockHeaderIndex==args->FirstFileBlocks.size()-1) break;
                FirstRecordIndex=0;
                FirstBlockHeaderIndex++;
                j--;
            }

            int pp=1;
        }
        FirstRecordIndex--;

    }



    int h=1;
}

static void ConcatenateThisChunkToPrevious(int chunkNo, concat_args_t *args)
{
    // Read the overlapping part from first file
    SaveFirstFileOverlappingRegion(chunkNo, args);

    // Read the overlapping part from second file
    SaveSecondFileOverlappingRegion(chunkNo, args);


    //Flush Overlapping Part
    FlushOverlappingPart(chunkNo, args);


    //Flush Remaining non-overlapping part of second file
    FlushSecondFileBlockNonOverlappingPart(chunkNo, args);

    args->firstFileBlockHeader=args->secondFileBlockHeader;
    args->firstFile=args->secondFile;

    // Open and Initialize the Second File and read Second Block and Record


//    InitSecondFile(chunkNo, args);
//
//    // First Flush the First File Block until common marker
//    FlushFirstFileBlockNonOverlappingPart(chunkNo, args);
//
//    // Check Overlap Eligibility
//    CheckIfOverlaps(chunkNo, args);
//
//    // Once Overlaps are verified, read both blocks and records together till first file ends
//    // Throw error if first file ends before second file
//    ReadOverlappingRegion(chunkNo, args);
//

    // Next, read the second file, till the start of the next file
//    do
//    {
//        args->outFile.writeBlock(args->secondFileBlockHeader);
//        while(!args->secondFileBlockHeader.isBlockFinished())
//        {
//            args->secondFileRecord.read(args->secondFile,args->secondFileBlockHeader);
//            args->outFile.writeRecord(args->secondFileRecord);
//        }
//    }while(args->secondFileBlockHeader.read(args->secondFile,args->out_hdr)
//           && args->secondFileBlockHeader.getEndBasePosition()<=args->start_pos[chunkNo+1]);
//






//    while(args->firstFileBlockHeader.read(args->firstFile,args->out_hdr) && args->firstFileBlockHeader.getEndBasePosition()<=args->start_pos[1])
//    {
//        args->outFile.writeBlock(args->firstFileBlockHeader);
//        while(!args->firstFileBlockHeader.isBlockFinished())
//        {
//            args->firstFileRecord.read(args->firstFile,args->firstFileBlockHeader);
//            args->outFile.writeRecord(args->firstFileRecord);
//        }
//    }
//

}

static void PerformPhasedConcat(concat_args_t *args)
{

    FlushFirstFileNonOverlappingPart(args);

    for (int i=1; i<args->nfnames; i++)
    {
        ConcatenateThisChunkToPrevious(i,args);
    }
}


static void init_data(concat_args_t *args)
{
    Initialize(args);
    AnalyseHeader(args);

    if(args->phased_concat)
    {
        PerformPhasedConcat(args);
    }




//    CurrentBlock.write()
//
//
//       m3vcfBlockHeader CurrentBlock;
//             CurrentBlock.read(m3vcfxStream,myHeader,true);
//       m3vcfBlockHeader CurrentBlock;
//             CurrentBlock.read(m3vcfxStream,myHeader,true);


//    IFILE m3vcfxStream = ifopen(args->fnames[0], "r");
//
//    m3vcfHeader myHeader;
//    m3vcfBlockHeader CurrentBlock;
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


static void usage(concat_args_t *args)
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
    concat_args_t* args = new concat_args_t();

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

