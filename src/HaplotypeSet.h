#ifndef HAPLOTYPESET_H_INCLUDED
#define HAPLOTYPESET_H_INCLUDED

#include "StringBasics.h"
#include "VcfFileReader.h"
#include "Unique.h"
//#include "pFile.h"
#include <fstream>
#include <functional>
#include <sstream>
#include <algorithm>


using namespace std;



class HaplotypeSet
{

	public:

		// Basic Properties

		int         numHaplotypes,numSamples,OrignumHaplotypes,OrigNumSamples, NoBlocks;
		int         numMarkers;
		bool        unphasedOutput,vcfType,m3vcfxType;
		vector<string> HeaderInfo;
        float m3vcfVersion;


		// Variant Properties

        vector<string> markerName;
		vector<variant> VariantList;
        string finChromosome;
        bool PseudoAutosomal;



		// Sample Properties
		vector<string> individualName,maleindividualName,femaleindividualName,HaploName;
		vector<int> SampleNoHaplotypes;
        vector<bool> SampleInputIndex,HapInputIndex;


		vector<int>        optEndPoints;
		vector<int>        ScaffoldIndex;
		vector<int>        UnScaffoldIndex;
		vector<ReducedHaplotypeInfo> ReducedStructureInfo;
		vector<vector<ReducedHaplotypeInfo> > ReducedStructureInfoBuffer;



		vector<vector<bool> >     haplotypesUnscaffolded;
		vector<vector<bool> >     MissingSampleUnscaffolded;
        String FileName;


		vector<vector<float> > Dosage;
		vector<vector<bool> > ImputedAlleles;
		int PrintStartIndex,PrintEndIndex;
        bool GT,DS,GP;
        bool onlyCompress,filter;
        bool EstimateNcompress;
        bool AllMaleTarget;
        bool Duplicates;
        int transFactor, cisFactor;
        string MyChromosome;

        vector<double>      Recom,Error;
		vector<char> refAlleleList;
		vector<bool> major, minor;



		vector<bool> missing, MarkerIndices;
		bool allowMissing, machType;
        bool Check;



        HaplotypeSet()
		{
			numHaplotypes = 0;
			numMarkers = 0;
            vcfType = false;




			ReducedStructureInfo.clear();
			individualName.clear();
			SampleNoHaplotypes.clear();


			markerName.clear();
			refAlleleList.clear();
			major.clear();
			minor.clear();
			missing.clear();
			MarkerIndices.clear();
			allowMissing = true;
            m3vcfxType=false;
			machType=false;
			PrintStartIndex=0;
			PrintEndIndex=0;
            ReadDUPLICATES=false;

		}






                bool            readm3vcfFile                               (String filename);
                bool            LoadInputHaplotypes                         (String filename);
                string          DetectInputFileType                         (String filename);
                bool            readvcfFile                                 (String filename);
                void            printm3VcfErr                               (String filename);



                void            GetInfo                             (variant &tempVariant,char *pch3);
                bool            FilterIndicator                       (VcfRecord &ThisRecord);
                bool            FilterIndicator                       (ReducedHaplotypeInfo &ThisRecord);
                bool            FilterIndicator                       (variant &ThisRecord);

                void            ReadThisBlock                               (int blockIndex,
                                                                            int tempRepCount,
                                                                            ReducedHaplotypeInfo &tempBlock,
                                                                             int tempVarCount,IFILE m3vcfxStream,
                                                                              int &NoMarkersImported, String filename);

                bool            ReadBlockHeader                             (ReducedHaplotypeInfo &tempBlocktoCheck,
                                                                              string &line, int blockIndex,
                                                                              int &tempVarCount, int &tempRepCount,
                                                                              String filename);



                bool            UpdateInputSampleList                       (String filename);
                bool            CreateSampleNoHaplotypes                (int &OrignumHaplotypes,String filename);

                void            GetHeader(VcfHeader& Header);


                string readFirstFileAfterHeader(IFILE m3vcfxStream);
                bool ReportSampleCount(string &line);
                void            readFile                                    (vector<string> &Vector, String FileName);












        bool    getScaffoldedHaplotype                      (int sample,int marker);
        bool    getMissingScaffoldedHaplotype                      (int sample,int marker);


        void    Create                                      (vector<bool> &tempHaplotype);
        void    calculateFreq                               ();
		void    CalculateFreq                               ();
		bool    LoadSnpList                                 (String filename);
		bool    LoadMachHaplotypes                          (String filename, String targetSnpfile, vector<string> &refSnpList);
		bool    LoadMachHaplotypes                          (String filename, String targetSnpfile);
        char    convertAlleles                              (string markerId, string indivId, const char *alleles, string refAlleleString, string altAlleleString);
		bool    LoadHaplotypes                              (String filename, String snpNames, int maxIndiv, int maxMarker,bool optmx,String CNO,int START,int END);
//		bool    readvcfFile                          (String filename, int maxIndiv, int maxMarker,String CNO,int START,int END,int WINDOW,bool rsid,bool compressOnly,bool filter);

        bool    LoadTargetHaplotypes                        (String filename, String targetSnpfile, vector<string> &refSnpList, HaplotypeSet &rHap);
		bool    LoadVcfTargetHaplotypes                     (String filename, String snpNames, vector<string> &refSnpList,HaplotypeSet &rHap);
        void    convertToReducedStructure                   (vector<int> &optStructure);
        void    writem3vcfFile                              (String &filename,bool &gzip);

        void    reconstructHaplotype                        (vector<bool> &reHaplotypes,int &index);
        void    SaveDosageForVcfOutput                      (int hapID,vector<float> dose,vector<bool> impAlleles);
        void    SaveDosageForVcfOutputSampleWise            (int SamID,string &SampleName, vector<float> &dose1,vector<float> &dose2,vector<bool> &impAlleles1,vector<bool> &impAlleles2);
        void    InitializeDosageForVcfOutput                (int NHaps,int NMarkers);
        void    InitializePartialDosageForVcfOutput         (int NHaps,int NMarkers, vector<bool> &Format);
        void    InitializePartialDosageForVcfOutputMaleSamples    (int NHaps,int NMarkers, vector<bool> &Format);
        void    PrintDosageForVcfOutputForID                (IFILE vcfdose, int MarkerIndex,bool majorIsReference,char refAllele);
        void    PrintPartialDosageForVcfOutputForID         (IFILE vcfdose, int MarkerIndex,bool majorIsReference,char refAllele);
        bool    BasicCheckForTargetHaplotypes               (String filename);
        string  DetectTargetFileType                        (String filename);
        void    SaveDosageForVcfOutputSampleWiseChrX        (int SamID,string &SampleName, vector<float> &dose1,
                                                            vector<bool> &impAlleles1);
        bool    CheckValidChrom                             (string chr);

        void PrintDosageForVcfOutputForIDMaleSamples(IFILE vcfdose, int MarkerIndex,bool majorIsReference,char refAllele);

        void updateCoeffs(int trans,int cis)
        {
            transFactor = trans;
            cisFactor = cis;

        }

        bool FilterPASS,FilterMaxIndiv,ReadDUPLICATES;
        string KeepFiltered,FilterCHR,FilterNOTCHR,RemoveFiltered,FilterPOSITIONS,FilterEXCLUDEPOSITIONS;
        int NoDUPLICATES,NoMULTIALLELIC,NoFAIL,CPU,FilterSTART,FilterEND;


        String KeepSample,RemoveSample;
        int MaxIndiv;

        void    InitializeSiteFilterParameters(String CNO,String NOTCNO,int START,int END,
                                           String &positions, String &excludepositions,

                                           bool filter
                                                  ,String &keepFiltered,
                                                  String &removeFiltered,
                                                  bool readDuplicates,int cpus, bool check)

        {
            FilterPASS=filter;

            FilterCHR=CNO.c_str();
            FilterNOTCHR=NOTCNO.c_str();
            FilterSTART=START;
            FilterEND=END;
            FilterPOSITIONS=positions.c_str();
            FilterEXCLUDEPOSITIONS=excludepositions.c_str();


            FilterMaxIndiv=0;
            CPU=cpus;
            NoDUPLICATES=0;
            NoMULTIALLELIC=0;
            NoFAIL=0;
            KeepFiltered=keepFiltered.c_str();
            RemoveFiltered=removeFiltered.c_str();

            ReadDUPLICATES=readDuplicates;
            CPU=cpus;
            Check=check;


        }



        void    InitializeSampleFilterParameters(String keepSample,String removeSample,int max_Indiv)
        {
            KeepSample=keepSample.c_str();
            RemoveSample=removeSample.c_str();
            MaxIndiv=max_Indiv;

        }










        bool Recode;

        void    InitializeOutputParameters(bool recode)
        {
            Recode=recode;
        }


};


#endif // HAPLOTYPESET_H_INCLUDED
