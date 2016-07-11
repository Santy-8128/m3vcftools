

#include <stdio.h>
#include <iostream>
#include <ctime>
#include <unistd.h>
#include <fstream>
#include <ostream>
#include <omp.h>


#include "Parameters.h"
#include "StringBasics.h"
#include "HaplotypeSet.h"
#include "Unique.h"
#include "Analyse.h"
int transFactor = 3;
int cisFactor = 2;


using namespace std;
void m3vcftoolsVersion();
void helpFile();

int main(int argc, char ** argv)
{
	// Parameter Options

    String infile = "";
	String haps = "", snps = "",removeSam="";

	// FILTER FLAG FILTERING
	String keepFiltered="",removeFiltered="";
	bool passOnly = false;

	// POSITION FILTERING
	String positions="",excludePositions="",chr="",notChr="", hapr2PosFile="";
	int start=0, end=0;

    // OUTPUT OPTIONS
    bool recode=false,uncode=false,freq=false,counts=false,hapr2=false;


	String outfile = "m3vcftools.Output";
	String format = "GT,DS";
	String recFile = "", errFile = "",golden="";
	int cpus=1,  max_Indiv = 0, max_marker = 0;
    vector<bool> formatVector(3,false);
    #ifdef _OPENMP
    cpus=5;
    #endif


    int ldWindow=-999999999, ldWindowMin=-999999999, ldWindowBp=-999999999, ldWindowBpMin=-999999999;
	double minHapR2=0.0;
	bool check = false, log = false, duplicates=false, unphasedOutput=false, phased = false, doseOutput = false, vcfOutput = true, gzip = true, nobgzip = false;
	bool processReference=false,updateModel=false, onlyRefMarkers=false, help = false, params = false;
    String MyChromosome="";

    String keepSample="",removeSample="";


	ParameterList inputParameters;
	PhoneHome::allThinning = 50;

	BEGIN_LONG_PARAMETERS(longParameterList)
		LONG_PARAMETER_GROUP("MANDATORY OPTIONS")
		LONG_STRINGPARAMETER("in", &infile)

    BEGIN_LEGACY_PARAMETERS()

    	LONG_PARAMETER_GROUP("INPUT FILE OPTIONS")
		LONG_PARAMETER("check", &check)
		LONG_STRINGPARAMETER("vcf", &infile)

		LONG_PARAMETER_GROUP("OUTPUT FILE OPTIONS")
		LONG_STRINGPARAMETER("out", &outfile)
        LONG_PARAMETER("recode", &recode)
        LONG_PARAMETER("uncode", &uncode)




		LONG_PARAMETER_GROUP("FILTER FLAG FILTERING")
		LONG_PARAMETER("remove-filtered-all", &passOnly)
        LONG_STRINGPARAMETER("keep-filtered", &keepFiltered)
        LONG_STRINGPARAMETER("remove-filtered", &removeFiltered)


        LONG_PARAMETER_GROUP("SAMPLE FILTERING")
        LONG_STRINGPARAMETER("keep", &keepSample)
        LONG_STRINGPARAMETER("remove", &removeSample)
        LONG_INTPARAMETER("max-indv", &max_Indiv)


        LONG_PARAMETER_GROUP("POSITION FILTERING")
		LONG_STRINGPARAMETER("chr", &chr)
        LONG_STRINGPARAMETER("not-chr", &notChr)
        LONG_INTPARAMETER("from-bp", &start)
        LONG_INTPARAMETER("to-bp", &end)
        LONG_STRINGPARAMETER("positions ", &positions)
        LONG_STRINGPARAMETER("exclude-positions", &excludePositions)



		LONG_PARAMETER_GROUP("OUTPUT ALLELE STATISTICS")
		LONG_PARAMETER("freq", &freq)
        LONG_PARAMETER("counts", &counts)


		LONG_PARAMETER_GROUP("OUTPUT LD STATISTICS")
        LONG_PARAMETER("hap-r2", &hapr2)
        LONG_STRINGPARAMETER("hap-r2-positions ", &hapr2PosFile)
        LONG_INTPARAMETER("ld-window", &ldWindow)
        LONG_INTPARAMETER("ld-window-min", &ldWindowMin)
        LONG_INTPARAMETER("ld-window-bp", &ldWindowBp)
        LONG_INTPARAMETER("ld-window-bp-min", &ldWindowBpMin)
        LONG_DOUBLEPARAMETER("min-r2", &minHapR2)







//		LONG_STRINGPARAMETER("snps", &snps)
		LONG_PARAMETER_GROUP("Output Parameters")
		LONG_PARAMETER("processReference", &processReference)
		LONG_PARAMETER("updateModel", &updateModel)
		LONG_PARAMETER("nobgzip", &nobgzip)
		LONG_PARAMETER("vcfOutput", &vcfOutput)
		LONG_PARAMETER("doseOutput", &doseOutput)
		LONG_PARAMETER("hapOutput", &phased)
		LONG_STRINGPARAMETER("format", &format)
		LONG_PARAMETER_GROUP("Other Parameters")
		LONG_PARAMETER("log", &log)
		LONG_PARAMETER("help", &help)
		LONG_INTPARAMETER("cpus", &cpus)
		LONG_PARAMETER("params", &params)
		LONG_PHONEHOME(VERSION)
		LONG_STRINGPARAMETER("MyChromosome", &MyChromosome)
		LONG_PARAMETER("onlyRefMarkers", &onlyRefMarkers)
		LONG_INTPARAMETER("transFactor", &transFactor)
		LONG_INTPARAMETER("cisFactor", &cisFactor)
//		LONG_INTPARAMETER("sample", &max_indiv)
		LONG_INTPARAMETER("marker", &max_marker)
		LONG_PARAMETER("duplicates", &duplicates)
		LONG_PARAMETER("unphasedOutput", &unphasedOutput)
		END_LONG_PARAMETERS();


	inputParameters.Add(new LongParameters(" Command Line Options: ",longParameterList));

    String compStatus;
	inputParameters.Read(argc, &(argv[0]));

    FILE *LogFile=NULL;
    if(log)
        LogFile=freopen(outfile+".logfile","w",stdout);
    dup2(fileno(stdout), fileno(stderr));


    m3vcftoolsVersion();
	if (help)
	{
		helpFile();
		return(-1);
	}

	inputParameters.Status();

    #ifdef _OPENMP
        omp_set_num_threads(cpus);
    #else
        cpus=1;
    #endif


    if(nobgzip)
        gzip=false;

    cout<<endl<<endl;
	if (infile == "" )
	{
		cout<< " Missing \"--in\"(or \"--vcf\"), a required parameter.\n";
		cout<< " Type \"--help\" for more help.\n\n";
		compStatus="Command.Line.Error";
		PhoneHome::completionStatus(compStatus.c_str());
		return(-1);
	}





    HaplotypeSet InputFile;

    InputFile.InitializeSiteFilterParameters(chr,notChr,start,end,positions,excludePositions,

                                         passOnly,keepFiltered,removeFiltered,

                                         duplicates,cpus,check);

    InputFile.InitializeSampleFilterParameters(keepSample, removeSample, max_Indiv);

    InputFile.InitializeOutputParameters(recode);



	int start_time = time(0);
	int time_prev = start_time;

    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                              INPUT HAPLOTYPE FILE                           "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;



    InputFile.updateCoeffs(transFactor,cisFactor);


	if (!InputFile.LoadInputHaplotypes(infile))
	{
		cout << "\n Program Exiting ... \n\n";
		compStatus="Input.File.Load.Error";
        PhoneHome::completionStatus(compStatus.c_str());
        return(-1);
	}


    cout<<endl;


	int time_load = time(0) - time_prev;
	time_prev = time(0);



    cout << "\n Time taken to load reference haplotype set = " << time_load << " seconds."<<endl<<endl;



    Analyse ThisData(gzip,outfile);

    if(recode)
        ThisData.writem3vcfFile(InputFile);

//    if(uncode)
//        ThisData.writeVcfFile(InputFile);







    if(freq)
        {
            ThisData.CalculateFreq(InputFile);
            ThisData.PrintFreqOutput(InputFile);
        }

    if(counts)
        {
            ThisData.CalculatCount(InputFile);
            ThisData.PrintCountOutput(InputFile);
        }

    if(hapr2)
        {
            ThisData.InitializeLDParameters(hapr2PosFile, ldWindow, ldWindowMin, ldWindowBp, ldWindowBpMin,minHapR2);


//cout<<" UMMM "<<ThisData.LdWindow<<endl;

            ThisData.CalculatHapR2(InputFile);
        }






    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                                END OF PROGRAM                                 "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;

    time_load = time(0) - time_prev;
	int time_tot = time(0) - start_time;

    cout << "\n Program Successfully Implemented... \n ";


	printf("\n Total Run completed in %d hours, %d mins, %d seconds.\n",
		time_tot / 3600, (time_tot % 3600) / 60, time_tot % 60);

    cout<<"\n Thank You for using m3vcftools !!! "<<endl<<endl;

    if(log)
        fclose (LogFile);


    compStatus="Success";
    PhoneHome::completionStatus(compStatus.c_str());

	return 0;

}




void m3vcftoolsVersion()
{
	printf("\n\n -------------------------------------------------------------------------------- \n");
	printf("            m3vcftools - A Tool for Fast Querying on M3VCF Files\n");
	printf(" --------------------------------------------------------------------------------\n");
    printf(" (c) 2015 - Sayantan Das, Goncalo Abecasis \n");
//	printf(" Version	: Undocumented Release\n");
//	printf(" Built		: sayantan\n\n");
	cout<<"\n Version: " << VERSION<< ";\n Built: " << DATE << " by " << USER << std::endl;
    printf("\n URL = http://genome.sph.umich.edu/wiki/m3vcftools\n");



}

void helpFile()
{




    printf("\n\n\t  m3vcftools is a lower memory and more computationally efficient implementation of \"minimac2\".\n");


    printf("\t It is an algorithm for genotypic imputation that works on phased genotypes (say from MaCH).\n");
    printf("\t m3vcftools is designed to handle very large reference panels in a more computationally efficient \n");
    printf("\t way with no loss of accuracy. This algorithm analyzes only the unique sets of haplotypes in \n");
    printf("\t small genomic segments, thereby saving on time-complexity, computational memory but no loss\n");
    printf("\t in degree of accuracy.\n");

printf("\n\n ----------------------------------------------------------------------------------------- \n");
	printf("                            m3vcftools - List of Usage Options \n");
	printf(" -----------------------------------------------------------------------------------------\n\n");

    printf(" --------- Reference Haplotypes --------- \n");
    printf("\n     --infile filename   : VCF file or M3VCF file containing haplotype data for reference panel.\n");
    printf("             --passOnly   : This option only imports variants with FILTER = PASS.\n");



  printf("\n --------- GWAS Haplotypes --------- \n");
  printf("\n        --haps filename   : File containing haplotype data for target (gwas) samples. Must be VCF \n");
    printf("                            file. Zipped versions allowed.\n");

  printf("\n --------- Output Parameters --------- \n");
  printf("\n        --prefix output   : Prefix for all output files generated. By default: [m3vcftools.Output]\n");
    printf("     --processReference   : This option will only convert an input VCF file to M3VCF format\n");
    printf("                            (maybe for a later run of imputation). If this option is ON, \n");
    printf("                            no imputation would be performed.\n");
    printf("              --nobgzip   : If ON, output files will NOT be gzipped.\n");
    printf("          --updateModel   : If ON, the GWAS panel will also be used to update the parameter \n");
    printf("                            estimates (if and when estimates are found in M3VCF files)\n");
    printf("           --doseOutput   : If ON, imputed data will be output as MaCH dosage file [Default: OFF].\n");
    printf("            --hapOutput   : If ON, phased imputed data will be output as well [Default: OFF]. \n");
    printf("               --format   : Specifies which fields to output for the FORMAT field in output \n");
    printf("                            VCF file. Available handles: GT,DS,GP [Default: GT,DS].\n");



  printf("\n --------- Subset Parameters --------- \n");
  printf("\n               --chr 22   : Chromosome number for which we will carry out imputation.\n");
    printf("         --start 100000   : Start position for imputation by chunking.\n");
    printf("           --end 200000   : End position for imputation by chunking. \n");
    printf("         --window 20000   : Length of buffer region on either side of --start and --end.\n");


  printf("\n --------- Estimation Parameters --------- \n");
  printf("\n             --rounds 5   : Rounds of optimization for model parameters, which describe population \n");
    printf("                            recombination rates and per SNP error rates. By default = 5.\n");
    printf("           --states 200   : Maximum number of reference (or target) haplotypes to be examined  \n");
    printf("                            during parameter optimization. By default = 200.\n");
    printf("               --cpus 5   : Number of cpus for parallel computing. Works only with m3vcftools-omp.\n\n");

  printf("\n Please visit <http://genome.sph.umich.edu/wiki/m3vcftools> for detailed documentation ...\n\n");
    cout<<endl;



	return;
}


