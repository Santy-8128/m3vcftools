

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
    m3vcfHeader out_hdr;
    m3vcfFileWriter <m3vcfHeader> outFile;
    m3vcfBlockHeader firstFileBlockHeader, secondFileBlockHeader;
    m3vcfRecord firstFileRecord, secondFileRecord;
    m3vcfHeader myCommonHeader;
    vector<m3vcfBlock*> FirstFileBlocks, SecondFileBlocks;
    m3vcfBlock SecondFileBlock;

    // File Handling Variables
    IFILE firstFile, secondFile;

    // Hammond Distance 
    vector<int> HammondDist;
    vector<int> SamplesToBeSwapped;
    int NoOverlappingMarkers;
    int FirBPos, FirRPos, SecBPos, SecRPos;
    
    // Variables for frameworking
    vector<int> start_pos;


    //Common File and Argument Variables
    char **argv;
    vector<char*> fnames;
    const char *output_fname, *file_list;
    int phased_concat, nfnames, argc;
    int output_type, record_cmd_line;
 
};



void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(-1);
}


static const char* createCommandLine(concat_args_t &args, const char *optionName)
{
    string str = "##m3vcftools_Command=" + (string)optionName;
    for (int i=1; i<args.argc; i++){str+=" "; str += +args.argv[i];}
    time_t tt;
    str+="; Date=";
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    str+= asctime (timeinfo)  ;
    return  str.substr(0,str.size()-1).c_str();
}



static void Initialize(concat_args_t &args)
{
    if ( args.phased_concat == 1 )
    {
        args.start_pos.clear();
        args.start_pos.resize(args.nfnames+1,0);
    }
    args.firstFile = NULL, args.secondFile = NULL;

}


static void AnalyseHeader(concat_args_t &args)
{
    // Gather and Merge Header Information
    int i;
    for (i=0; i<args.nfnames; i++)
    {
        IFILE m3vcfxStream = ifopen(args.fnames[i], "r"); if ( !m3vcfxStream ) error("[ERROR:] Failed to open: %s\n", args.fnames[i]);
        m3vcfHeader myHeader;
        myHeader.read(m3vcfxStream);

        // Check if samples are same and in same order in different files
        if (!args.out_hdr.checkMergeHeader(myHeader) )
            error("[ERROR:] Different sample information in %s.\n", args.fnames[i]);
        else
            args.out_hdr.mergeHeader(myHeader);

        // If phased concat then store first position to get framework
        if ( args.phased_concat == 1 )
        {
            m3vcfBlockHeader CurrentBlock;
            CurrentBlock.read(m3vcfxStream,myHeader,true);
            args.start_pos[i] = CurrentBlock.getStartBasePosition();
        }
        ifclose(m3vcfxStream);
    }
    
    if ( args.phased_concat == 1 )
    {
        args.start_pos[i]=MAX_BASE_POSITION;
        for (i=1; i<args.nfnames; i++)
        {
            if(args.start_pos[i]<=args.start_pos[i-1]) 
            error("[ERROR:] Chunk [%s] and [%s] are in the wrong order\n", args.fnames[i-1], args.fnames[i]);
        }
    
        
    }
    if (args.record_cmd_line==1) args.out_hdr.appendMetaLine(createCommandLine(args,"concat"));
    args.outFile.open(args.output_fname, args.out_hdr, args.output_type==4? InputFile::UNCOMPRESSED : InputFile::GZIP);   
}


static void FlushFirstFileNonOverlappingPart(concat_args_t &args)
{
    args.firstFile = ifopen(args.fnames[0], "r"); if ( !args.firstFile ) error("[ERROR:] Failed to open: %s\n", args.fnames[0]);
    args.myCommonHeader.read(args.firstFile);
    bool ReadStatus = args.firstFileBlockHeader.read(args.firstFile,args.out_hdr);
    
    // Use ReadStatus to see if the end of firstFile is reached before the start
    // position of the second file. This is ojectively the first file, so does
    // NOT need any phase swapping
    while(ReadStatus && args.firstFileBlockHeader.getEndBasePosition()<args.start_pos[1])
    {
        args.outFile.writeBlock(args.firstFileBlockHeader);
        while(!args.firstFileBlockHeader.isBlockFinished())
        {
            args.firstFileRecord.read(args.firstFile,args.firstFileBlockHeader);
            args.outFile.writeRecord(args.firstFileRecord);
        }
        ReadStatus = args.firstFileBlockHeader.read(args.firstFile,args.out_hdr);
    }
    
    if(ReadStatus==false)
        error("[ERROR:] Consecutive chunks [%s] and [%s] do NOT overlap\n", args.fnames[0], args.fnames[1]);

}

static void InitSecondFile(int chunkNo, concat_args_t &args)
{
    args.secondFile = ifopen(args.fnames[chunkNo], "r"); if ( !args.secondFile ) error("[ERROR:] Failed to open: %s\n", args.fnames[chunkNo]);
    args.myCommonHeader.read(args.secondFile);
}


static void SaveFirstFileOverlappingRegion(int chunkNo, concat_args_t &args)
{
    args.FirstFileBlocks.clear();
    args.FirstFileBlocks.resize(1);
    args.FirstFileBlocks[0]=new m3vcfBlock;
    args.FirstFileBlocks[0]->CopyBlockHeader(args.firstFileBlockHeader, args.myCommonHeader);

    do
    {
        // We should swap the phase as we are reading the data.
        // This way we ensure that the new read data of the overlapping part
        // is already phased-aligned correctly. This way we check for phase-
        // alignment with the next file in the proper way
        args.FirstFileBlocks.back()->swapPhase(args.SamplesToBeSwapped);
        while(!args.FirstFileBlocks.back()->isBlockFinished())
        {
            (args.FirstFileBlocks.back())->readRecord(args.firstFile);
        }

        args.FirstFileBlocks.resize(args.FirstFileBlocks.size()+1);
        args.FirstFileBlocks.back()=new m3vcfBlock;
    }while((args.FirstFileBlocks.back()->readHeader(args.firstFile,args.out_hdr)));
    
    // Remove the last information from FirstFileBlocks since it is empty
    args.FirstFileBlocks.pop_back();
    ifclose(args.firstFile);
}

static void SaveSecondFileOverlappingRegion(int chunkNo, concat_args_t &args)
{
    // Get the Last position of the first file chunk (that would determine 
    // the overlapping region)
    int PrevLastPosition = args.FirstFileBlocks.back()->getEndBasePosition();

    InitSecondFile(chunkNo, args);

    args.SecondFileBlocks.clear();
    args.SecondFileBlocks.resize(1);
    args.SecondFileBlocks[0]=new m3vcfBlock;
    args.NoOverlappingMarkers = 0;
    bool ReadStatus=args.SecondFileBlocks.back()->readHeader(args.secondFile,args.out_hdr);
    
    while(ReadStatus && args.SecondFileBlocks.back()->getStartBasePosition()<PrevLastPosition)
    {
        // Here we do NOT do any phase alignment since this is the new file
        // overlapping part that is being read.
        while(!args.SecondFileBlocks.back()->isBlockFinished())
        {
            (args.SecondFileBlocks.back())->readRecord(args.secondFile);
            args.NoOverlappingMarkers++;
        }

        args.NoOverlappingMarkers--;
        args.SecondFileBlocks.resize(args.SecondFileBlocks.size()+1);
        args.SecondFileBlocks.back()=new m3vcfBlock;
        ReadStatus=args.SecondFileBlocks.back()->readHeader(args.secondFile,args.out_hdr);
    
    }
    if(ReadStatus==false)
        error("[ERROR:] Chunk [%s] is contained inside chunk [%s]\n", args.fnames[chunkNo-1], args.fnames[chunkNo]);

    if(args.SecondFileBlocks.size()==1 && args.SecondFileBlocks.back()->getStartBasePosition()>PrevLastPosition)
    error("[ERROR:] Consecutive chunks [%s] and [%s] do NOT overlap\n", args.fnames[chunkNo-1], args.fnames[chunkNo]);
     


}

static void FlushSecondFileBlockNonOverlappingPart(int chunkNo, concat_args_t &args)
{
    
    args.SecondFileBlocks.back()->CopyToBlockHeader(args.SecondFileBlock);
    bool ReadStatus=true;
      
    do
    {
        int MarkerCounter = 0;
        // Swap the Phase as you print out the blocks for the non-overlapping part
        // of the second file
        args.SecondFileBlock.swapPhase(args.SamplesToBeSwapped);
        args.SecondFileBlock.writeHeader(args.outFile);
        while(!args.SecondFileBlock.isBlockFinished())
        {
            args.SecondFileBlock.readRecord(args.secondFile);
            args.SecondFileBlock.writeRecord(args.outFile, MarkerCounter++);
        }
        ReadStatus=args.SecondFileBlock.readHeader(args.secondFile,args.out_hdr);
    }
    while(ReadStatus && args.SecondFileBlock.getEndBasePosition()<args.start_pos[chunkNo+1]);

    // If end of second file is reached before the beginning of the next and if this is not
    // the last pair of chunks, then throw and error for non-overlapping
    if(ReadStatus==false && chunkNo<args.nfnames-1)
        error("[ERROR:] Consecutive chunks [%s] and [%s] do NOT overlap\n", args.fnames[chunkNo], args.fnames[chunkNo+1]);

    // Swap First File and Second File and their Block Header 
    args.SecondFileBlock.CopyToBlockHeader(args.firstFileBlockHeader);
    args.firstFile=args.secondFile;

}

   
   
static void UpdateHammondDistance(concat_args_t &args, 
        int FirstHeaderIndex, int FirstRecordIndex,
        int SecondHeaderIndex, int SecondRecordIndex)
{
    
    args.FirstFileBlocks[FirstHeaderIndex]->getM3vcfRecord(FirstRecordIndex)->Deserialize();
    args.SecondFileBlocks[SecondHeaderIndex]->getM3vcfRecord(SecondRecordIndex)->Deserialize();
    
    int HaploIndex = 0;
    for(int i=0; i<args.myCommonHeader.getNumSamples(); i++)
    {
        
        if(args.FirstFileBlocks[FirstHeaderIndex]->getSamplePloidy(i)==2)
        {
            AlleleType L1 = args.FirstFileBlocks[FirstHeaderIndex]->getAllele(FirstRecordIndex, HaploIndex);
            AlleleType R1 = args.SecondFileBlocks[SecondHeaderIndex]->getAllele(SecondRecordIndex, HaploIndex);
            AlleleType R2 = args.SecondFileBlocks[SecondHeaderIndex]->getAllele(SecondRecordIndex, HaploIndex+1);
            
            if(L1!=R1)
                args.HammondDist[2*i]++;
            if(L1!=R2)
                args.HammondDist[2*i+1]++;
        }
        HaploIndex+=args.FirstFileBlocks[FirstHeaderIndex]->getSamplePloidy(i);
    }
   
}

static void FlushOverlappingPart(int chunkNo, concat_args_t &args)
{
    
    args.FirstFileBlocks[args.FirBPos]->DeleteVariantsFrom(args.FirRPos+1);
    for(int i=0; i<=args.FirBPos; i++)
    {
        args.FirstFileBlocks[i]->write(args.outFile);
    }
    args.SecondFileBlocks[args.SecBPos]->DeleteVariantsTill(args.SecRPos-1);
    for(int i=args.SecBPos; i<(int)args.SecondFileBlocks.size()-1; i++)
    {
        args.SecondFileBlocks[i]->swapPhase(args.SamplesToBeSwapped);
        args.SecondFileBlocks[i]->write(args.outFile);
    }
    
}


static void CreateSamplesToBeSwappedList(int chunkNo, concat_args_t &args)
{
    args.SamplesToBeSwapped.clear();
    for(int i=0; i<args.myCommonHeader.getNumSamples(); i++)
    {
        if(args.SecondFileBlocks[0]->getSamplePloidy(i)==2)
        {
            if(args.HammondDist[2*i+1] < args.HammondDist[2*i])
                args.SamplesToBeSwapped.push_back(i);
        }
    }
    
}


static void AnalyseAndFlushOverlappingPart(int chunkNo, concat_args_t &args)
{
    int FirstBlockHeaderIndex=0, FirstRecordIndex=0;
    while(args.FirstFileBlocks[0]->getM3vcfRecord(FirstRecordIndex)->getBasePosition() < args.start_pos[chunkNo]) 
    {
        FirstRecordIndex++;
    }
    args.HammondDist.clear(); args.HammondDist.resize(2*args.myCommonHeader.getNumSamples(),0);
    
    int NoMarkersRead = 0;
    for(int i=0; i<(int)args.SecondFileBlocks.size()-1; i++)
    {

        for(int j=0; j< args.SecondFileBlocks[i]->getNumMarkers(); j++)
        {

            if(args.SecondFileBlocks[i]->getM3vcfRecord(j)->IsMatching( *(args.FirstFileBlocks[FirstBlockHeaderIndex]->getM3vcfRecord(FirstRecordIndex) )  ) !=1)
                error("[ERROR:] Variant [%s] present only in one chunk [%s]\n",
                      args.SecondFileBlocks[i]->getM3vcfRecord(j)->getVariantID().c_str(),
                      args.fnames[chunkNo]);

            UpdateHammondDistance(  args, FirstBlockHeaderIndex, FirstRecordIndex, i, j);
            NoMarkersRead++;
            if(NoMarkersRead==args.NoOverlappingMarkers/2)
            {
                args.FirBPos=FirstBlockHeaderIndex;
                args.FirRPos=FirstRecordIndex;
                args.SecBPos=i;
                args.SecRPos=j;
            }
                
            
            FirstRecordIndex++;
            if(FirstRecordIndex>=args.FirstFileBlocks[FirstBlockHeaderIndex]->getNumMarkers())
            {
                if(FirstBlockHeaderIndex==args.FirstFileBlocks.size()-1) break;
                FirstRecordIndex=0;
                FirstBlockHeaderIndex++;
                j--;
            }
        }
        FirstRecordIndex--;

    }

    CreateSamplesToBeSwappedList(chunkNo, args);
    FlushOverlappingPart(chunkNo, args);
  
   
}

static void ConcatenateThisChunkToPrevious(int chunkNo, concat_args_t &args)
{
    // Read the overlapping part from first file
    SaveFirstFileOverlappingRegion(chunkNo, args);

    // Read the overlapping part from second file
    SaveSecondFileOverlappingRegion(chunkNo, args);

    //Flush Overlapping Part
    AnalyseAndFlushOverlappingPart(chunkNo, args);

    //Flush Remaining non-overlapping part of second file
    FlushSecondFileBlockNonOverlappingPart(chunkNo, args);
}

static void PerformPhasedConcat(concat_args_t &args)
{
    FlushFirstFileNonOverlappingPart(args);
    for (int i=1; i<args.nfnames; i++)
    {
        ConcatenateThisChunkToPrevious(i,args);
    }
}


static void concat_Data(concat_args_t &args)
{
    Initialize(args);
    AnalyseHeader(args);

    if(args.phased_concat == 1)
    {
        PerformPhasedConcat(args);
    }
    else
    {
        error("[ERROR:] Only ligation is supported currently \n");
    }

}


static void usage(concat_args_t &args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, " About:   Concatenate or combine M3VCF files. All source files must have the same sample\n");
    fprintf(stderr, "          columns appearing in the same order. The program can be used, for example, to\n");
    fprintf(stderr, "          ligate haplotypes \n");
    fprintf(stderr, " Note :   STDIN(-) cannot be used as in input file for concatenation\n");
    fprintf(stderr, " Usage:   m3vcftools concat [options] <A.m3vcf.gz> [<B.m3vcf.gz> [...]]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, " Options:\n");
    fprintf(stderr, "   -f, --file-list <file>         Read the list of files from a file.\n");
    fprintf(stderr, "   -l, --ligate                   Ligate phased VCFs by matching phase at overlapping haplotypes\n");
    fprintf(stderr, "       --no-version               do not append version and command line to the header\n");
    fprintf(stderr, "   -o, --output <file>            Write output to a file [standard output]\n");
    fprintf(stderr, "   -O, --output-type <m|M>        m: compressed M3VCF, M: uncompressed M3VCF [m]\n");
    fprintf(stderr, "\n");
    exit(1);
}




int main_m3vcfconcat(int argc, char *argv[])
{
    int c;
    concat_args_t args;

    args.argc    = argc; args.argv = argv;
    args.output_fname = "-";
    args.output_type = ZIPM3VCF;
    args.record_cmd_line = 1;
    args.nfnames=0;
    args.fnames.clear(); 
    args.file_list=NULL;
    args.phased_concat = 0;
    
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
            case 'l': args.phased_concat = 1; break;
            case 'f': args.file_list = optarg; break;
            case 'o': args.output_fname = optarg; break;
            case 'O':
                switch (optarg[0]) {
                    case 'm': args.output_type = ZIPM3VCF; break;
                    case 'M': args.output_type = M3VCF; break;
                    default: error("[ERROR:] The output type \"%s\" not recognized\n", optarg);
                };
                break;
            case  8 : args.record_cmd_line = 0; break;
            case '?': usage(args); break;
            default: error("[ERROR:] Unknown argument: %s\n", optarg);
        }
    }
    
    while ( optind<argc )
    {
        args.nfnames++;
        args.fnames.resize(args.nfnames);
        args.fnames[args.nfnames-1] = strdup(argv[optind]);
        optind++;
    }
    if ( args.file_list )
    {
        if ( args.nfnames>0 ) error("[ERROR:] Cannot combine -f with file names on command line.\n");
        
        IFILE myfile = ifopen(args.file_list, "r");
        if ( !myfile ) error("[ERROR:] Could not read the file: %s\n", args.file_list);
        string line;
        while (myfile->readLine(line)!=-1)
        {
            args.nfnames++;
            args.fnames.resize(args.nfnames);
            args.fnames[args.nfnames-1] = strdup(line.c_str());
            line.clear();  
        }
      }
    if ( args.nfnames == 0 ) usage(args);

    concat_Data(args);
    return 0;
}

