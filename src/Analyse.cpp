#include "Analyse.h"
#include <omp.h>

void Analyse::writem3vcfFile(HaplotypeSet &ThisData)
{


    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                              OUTPUT PROCESSING                            "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;



    std::cout << "\n Writing compresses haplotype information to [.m3vcf] file  : " <<OutputFileName+".m3vcf" + (BgZip ? ".gz" : "")<< endl<<endl;

    m3vcfVersion=ThisData.m3vcfVersion;

    IFILE m3vcffile = ifopen(OutputFileName + ".m3vcf" + (BgZip ? ".gz" : ""), "wb",(BgZip ? InputFile::BGZF : InputFile::UNCOMPRESSED));
    ifprintf(m3vcffile, "##fileformat=M3VCF\n");
    ifprintf(m3vcffile, "##version=%2.1f\n",m3vcfVersion);
    ifprintf(m3vcffile, "##source=m3vcftools\n");
    ifprintf(m3vcffile, "##compression=block\n");
    ifprintf(m3vcffile, "##n_blocks=%d\n",ThisData.ReducedStructureInfo.size());
    ifprintf(m3vcffile, "##n_haps=%d\n",ThisData.numHaplotypes);
    ifprintf(m3vcffile, "##n_markers=%d\n",ThisData.numMarkers);
    if(ThisData.finChromosome=="X")
        ifprintf(m3vcffile, "##chrxRegion=%s\n",ThisData.PseudoAutosomal?"PseudoAutosomal":"NonPseudoAutosomal");
    ifprintf(m3vcffile, "##<Note=This is NOT a VCF File and cannot be read by vcftools>\n");
    ifprintf(m3vcffile, "#CHROM\t");
    if(m3vcfVersion<2)
    {
       ifprintf(m3vcffile, "POS\t");
    }
    else
    {
        ifprintf(m3vcffile, "START\t");
        ifprintf(m3vcffile, "END\t");
    }
    ifprintf(m3vcffile, "ID\t");
    ifprintf(m3vcffile, "REF\t");
    ifprintf(m3vcffile, "ALT\t");
    ifprintf(m3vcffile, "QUAL\t");
    ifprintf(m3vcffile, "FILTER\t");
    ifprintf(m3vcffile, "INFO\t");
    ifprintf(m3vcffile, "FORMAT");
    int i,j,k;



    for(i=0;i<(int)ThisData.individualName.size();i++)
    {
        ifprintf(m3vcffile, "\t%s_HAP_1",ThisData.individualName[i].c_str());

        if(ThisData.finChromosome!="X")
            ifprintf(m3vcffile, "\t%s_HAP_2",ThisData.individualName[i].c_str());
        else if(ThisData.SampleNoHaplotypes[i]==2)
            ifprintf(m3vcffile, "\t%s_HAP_2",ThisData.individualName[i].c_str());
    }
    ifprintf(m3vcffile, "\n");

    int length=(int)ThisData.ReducedStructureInfo.size();
    string cno;

    for(i=0;i<length;i++)
    {

        cno=ThisData.VariantList[ThisData.ReducedStructureInfo[i].startIndex].chr;
//        int startBp=ThisData.VariantList[ThisData.ReducedStructureInfo[i].startIndex].bp;
//        int endBp=ThisData.VariantList[ThisData.ReducedStructureInfo[i].endIndex].bp;
        int nvariants=ThisData.ReducedStructureInfo[i].endIndex-ThisData.ReducedStructureInfo[i].startIndex+1;
        int reps=ThisData.ReducedStructureInfo[i].uniqueCardinality.size();


        ifprintf(m3vcffile, "%s\t",cno.c_str());

        if(m3vcfVersion<2)
        {
            ifprintf(m3vcffile, "%d-%d\t",ThisData.VariantList[ThisData.ReducedStructureInfo[i].startIndex].bp
                     ,ThisData.VariantList [ThisData.ReducedStructureInfo[i].endIndex].bp);
        }
        else
        {
            ifprintf(m3vcffile, "%d\t",ThisData.VariantList[ThisData.ReducedStructureInfo[i].startIndex].bp);
            ifprintf(m3vcffile, "%d\t",ThisData.VariantList [ThisData.ReducedStructureInfo[i].endIndex].bp);
        }

        ifprintf(m3vcffile, "<BLOCK:%d-%d>\t.\t.\t.\t.\t",ThisData.ReducedStructureInfo[i].startIndex,ThisData.ReducedStructureInfo[i].endIndex);

        ifprintf(m3vcffile, "B%d;VARIANTS=%d;REPS=%d\t.",i+1,nvariants,reps);


        for(j=0;j<ThisData.numHaplotypes;j++)
            ifprintf(m3vcffile, "\t%d",ThisData.ReducedStructureInfo[i].uniqueIndexMap[j]);

        ifprintf(m3vcffile, "\n");

        for(j=0;j<nvariants;j++)
        {
            ifprintf(m3vcffile, "%s\t",cno.c_str());
            int VarPosition=j+ThisData.ReducedStructureInfo[i].startIndex;
            variant *ThisVariant=&ThisData.VariantList[VarPosition];

            ifprintf(m3vcffile, "%d\t",(*ThisVariant).bp);
            if(m3vcfVersion>=2.0)
                ifprintf(m3vcffile, "%d\t",(*ThisVariant).bp);
            ifprintf(m3vcffile, "%s\t",(*ThisVariant).name.c_str());
            ifprintf(m3vcffile, "%s\t%s\t%.0f\t%s\t",(*ThisVariant).refAlleleString.c_str(),
                     (*ThisVariant).altAlleleString.c_str(),
                     (*ThisVariant).Qual,
                     (*ThisVariant).Filter.c_str());


            ifprintf(m3vcffile, "B%d.M%d",i+1,j+1);

            if((*ThisVariant).FullInfo!="")
            {
                ifprintf(m3vcffile, ";%s",(*ThisVariant).FullInfo.c_str());
            }

            ifprintf(m3vcffile, "\t");

            for(k=0;k<reps;k++)
            {
                ifprintf(m3vcffile,"%d",!ThisData.ReducedStructureInfo[i].uniqueHaps[k][j]? 0:1);
//                int tempVal=ThisData.ReducedStructureInfo[i].uniqueIntHaps[k][j];
//
//                if(tempVal>=0)
//                    ifprintf(m3vcffile,"%d",ThisData.ReducedStructureInfo[i].uniqueIntHaps[k][j]);
//                else
//                    ifprintf(m3vcffile,"-");





//                uniqueHaps
            }
            ifprintf(m3vcffile, "\n");
        }
    }

    std::cout << " Successfully written file ... "<<endl;
    ifclose(m3vcffile);

}





void Analyse::CalculateFreq(HaplotypeSet &ThisData)
{
    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                              FREQUENCY PROCESSING                            "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;

    std::cout << "\n Calculating Allele Frequency ... " <<endl<<endl;


    AlleleFreq.clear();
    AlleleFreq.resize(ThisData.numMarkers, 0.0);

//    int i,j,k;
    #ifdef _OPENMP
        omp_set_num_threads(10);
    #endif

    #pragma omp parallel for
    for(int k=0;k<(int)ThisData.ReducedStructureInfo.size();k++)
    {
//        cout<<k<<endl;
        ReducedHaplotypeInfo *ThisInfo=&ThisData.ReducedStructureInfo[k];

        for (int i = 0; i<(int)ThisInfo->uniqueCardinality.size(); i++)
        {
            int j;
            for(j=ThisInfo->startIndex;j<ThisInfo->endIndex;j++)
            {
                if(ThisInfo->uniqueHaps[i][j-ThisInfo->startIndex])
                    {
                        AlleleFreq[j]+=ThisInfo->uniqueCardinality[i];
                    }
            }

            if(k==(int)ThisData.ReducedStructureInfo.size()-1)
                if(ThisInfo->uniqueHaps[i][j-ThisInfo->startIndex])
                    AlleleFreq[j]+=ThisInfo->uniqueCardinality[i];
        }
    }

//	major.resize(ThisData.numMarkers, false);
//	minor.resize(ThisData.numMarkers, 0);

	for (int i = 0; i<ThisData.numMarkers; i++)
	{

		AlleleFreq[i] /= (double)ThisData.numHaplotypes;

//		if (AlleleFreq[i]>0.5)
//        {
//            major[i] = true;
//            ThisData.VariantList[i].MinAlleleString=ThisData.VariantList[i].refAlleleString;
//            ThisData.VariantList[i].MajAlleleString=ThisData.VariantList[i].altAlleleString;
//        }
//        else
//        {
//            ThisData.VariantList[i].MinAlleleString=ThisData.VariantList[i].altAlleleString;
//            ThisData.VariantList[i].MajAlleleString=ThisData.VariantList[i].refAlleleString;
//        }
	}

}
void Analyse::PrintFreqOutput(HaplotypeSet &ThisData)
{

    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                              OUTPUT PROCESSING                            "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;


    std::cout << "\n Writing Frequency to [.frq] file  : " <<OutputFileName+".frq"<< endl<<endl;


    IFILE freqFile = ifopen(OutputFileName + ".frq" , "wb");


    ifprintf(freqFile, "#CHROM\t");
    ifprintf(freqFile, "POS\t");
    ifprintf(freqFile, "N_ALLELES\t");
    ifprintf(freqFile, "N_CHR\t");
    ifprintf(freqFile, "{ALLELE:FREQ}\n");

    for(int i=0;i<ThisData.numMarkers;i++)
    {

        ifprintf(freqFile, "%s\t",ThisData.VariantList[i].chr.c_str());
        ifprintf(freqFile, "%d\t",ThisData.VariantList[i].bp);
        ifprintf(freqFile, "2\t");
        ifprintf(freqFile, "2\t",ThisData.numHaplotypes);
        ifprintf(freqFile, "%s:%f\t",ThisData.VariantList[i].refAlleleString.c_str(),1-AlleleFreq[i]);
        ifprintf(freqFile, "%s:%f",ThisData.VariantList[i].altAlleleString.c_str(),AlleleFreq[i]);
        ifprintf(freqFile, "\n");

    }
    ifclose(freqFile);

}




void Analyse::CalculatCount(HaplotypeSet &ThisData)
{
    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                           ALLELE COUNT PROCESSING                            "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;

    std::cout << "\n Calculating Allele Count ... " <<endl<<endl;


    AlleleCount.clear();
    AlleleCount.resize(ThisData.numMarkers, 0);

    int i,j,k;
    for(k=0;k<(int)ThisData.ReducedStructureInfo.size();k++)
    {
        ReducedHaplotypeInfo *ThisInfo=&ThisData.ReducedStructureInfo[k];

        for (i = 0; i<(int)ThisInfo->uniqueCardinality.size(); i++)
        {
            for(j=ThisInfo->startIndex;j<ThisInfo->endIndex;j++)
            {
                if(ThisInfo->uniqueHaps[i][j-ThisInfo->startIndex])
                    {
                        AlleleCount[j]+=ThisInfo->uniqueCardinality[i];
                    }
            }

            if(k==(int)ThisData.ReducedStructureInfo.size()-1)
                if(ThisInfo->uniqueHaps[i][j-ThisInfo->startIndex])
                    AlleleCount[j]+=ThisInfo->uniqueCardinality[i];
        }
    }

}
void Analyse::PrintCountOutput(HaplotypeSet &ThisData)
{

    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                              OUTPUT PROCESSING                            "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;


    std::cout << "\n Writing Allele Count to [.frq.count] file  : " <<OutputFileName+".frq.count"<< endl<<endl;


    IFILE freqFile = ifopen(OutputFileName + ".frq.count" , "wb");


    ifprintf(freqFile, "#CHROM\t");
    ifprintf(freqFile, "POS\t");
    ifprintf(freqFile, "N_ALLELES\t");
    ifprintf(freqFile, "N_CHR\t");
    ifprintf(freqFile, "{ALLELE:COUNT}\n");

    for(int i=0;i<ThisData.numMarkers;i++)
    {

        ifprintf(freqFile, "%s\t",ThisData.VariantList[i].chr.c_str());
        ifprintf(freqFile, "%d\t",ThisData.VariantList[i].bp);
        ifprintf(freqFile, "2\t");
        ifprintf(freqFile, "2\t",ThisData.numHaplotypes);
        ifprintf(freqFile, "%s:%d\t",ThisData.VariantList[i].refAlleleString.c_str(),ThisData.numHaplotypes-AlleleCount[i]);
        ifprintf(freqFile, "%s:%d",ThisData.VariantList[i].altAlleleString.c_str(),AlleleCount[i]);
        ifprintf(freqFile, "\n");

    }
    ifclose(freqFile);

}




double Analyse::CalculateDPrime(double &ThisD, double &iFreq, double &jFreq)
{

    double Min,Y;

    if(ThisD<0)
    {
        Min=( iFreq*jFreq );
        Y=( (1-iFreq) * (1-jFreq) );

        if(Y<Min)
            Min=Y;
    }
    else
    {
        Min=( (iFreq) * (1-jFreq) );
        Y=( (1-iFreq) * (jFreq) );

        if(Y<Min)
            Min=Y;
    }

    return(ThisD/Min);

}
void Analyse::CalculateHapCrossBlock(vector<vector<double> > &HapR2, HaplotypeSet &ThisData, int Block1, int Block2)
{

//    HapR2.clear();
    ReducedHaplotypeInfo *ThisInfo1=&ThisData.ReducedStructureInfo[Block1];
    ReducedHaplotypeInfo *ThisInfo2=&ThisData.ReducedStructureInfo[Block2];
    int LastElement= (Block2==(int)ThisData.ReducedStructureInfo.size()-1?1:0);



    int ThisSize1=ThisInfo1->endIndex-ThisInfo1->startIndex;
    int ThisSize2=ThisInfo2->endIndex-ThisInfo2->startIndex+LastElement;

    int minMarkerinBetween = ThisInfo2->startIndex - ThisInfo1->endIndex+1;
    int minBPinBetween = ThisData.VariantList[ThisInfo2->startIndex].bp - ThisData.VariantList[ThisInfo1->endIndex-1].bp;

    int maxMarkerinBetween = ThisInfo2->endIndex - 1 + LastElement - ThisInfo1->startIndex;
    int maxBPinBetween = ThisData.VariantList[ThisInfo2->endIndex-1+LastElement].bp - ThisData.VariantList[ThisInfo1->startIndex].bp;

    if(LdWindow!=-999999999)
    {
        if(minMarkerinBetween>LdWindow)
            return;
    }

    if(LdWindowBp!=-999999999)
    {
        if(minBPinBetween>LdWindowBp)
            return;
    }

    if(LdWindowMin!=-999999999)
    {
        if(maxMarkerinBetween<LdWindowMin)
            return;
    }
    if(LdWindowBpMin!=-999999999)
    {
        if(maxBPinBetween<LdWindowBpMin)
            return;
    }




    std::map<string,int> HashCrossMap;
    HapR2.resize(ThisSize1);
    for(int i=0;i<ThisSize1;i++)
    {
        HapR2[i].resize(ThisSize2,0.0);
    }

     for(int i=0; i<ThisData.numHaplotypes; i++)
     {
        stringstream strs1,strs2;
        string tempString3;

        strs1<<ThisInfo1->uniqueIndexMap[i];
        strs2<<ThisInfo2->uniqueIndexMap[i];

        tempString3=(string)(strs1.str())+"_"+(string)(strs2.str());
        HashCrossMap[tempString3]++;

     }


    for (std::map<string,int>::iterator it=HashCrossMap.begin(); it!=HashCrossMap.end(); ++it)
    {
        string tempString= it->first;
        string delimiter = "_";

        int pos=tempString.find(delimiter);

        int i=atoi(tempString.substr(0,pos).c_str());
        tempString.erase(0, pos + delimiter.length());

        pos=tempString.find(delimiter);
        int j=atoi(tempString.substr(0,pos).c_str());

//            std::cout << i<<" _ "<<j << " => " << it->second << '\n';

        vector<bool> *UniqueHap1=&ThisInfo1->uniqueHaps[i];
        vector<bool> *UniqueHap2=&ThisInfo2->uniqueHaps[j];

        for(int marker1=ThisInfo1->startIndex;marker1<ThisInfo1->endIndex;marker1++)
        {
            for(int marker2=ThisInfo2->startIndex;marker2<ThisInfo2->endIndex+LastElement;marker2++)
            {

                if((*UniqueHap1)[marker1-ThisInfo1->startIndex] && (*UniqueHap2)[marker2-ThisInfo2->startIndex] )
                {
                    HapR2[marker1-ThisInfo1->startIndex][marker2-ThisInfo2->startIndex]+=it->second ;
                }
            }
        }
    }

//    PrintHapR2Output(ThisData,Block1,Block2);

}
void Analyse::PrintHapR2Output(vector<vector<double> > &HapR2, HaplotypeSet &ThisData, int Block1, int Block2)
{

    IFILE LDFile = ifopen(OutputFileName + ".hap.ld" , "a");
    ReducedHaplotypeInfo *ThisInfo1=&ThisData.ReducedStructureInfo[Block1];
    ReducedHaplotypeInfo *ThisInfo2=&ThisData.ReducedStructureInfo[Block2];

    int LastElement1= (Block1==(int)ThisData.ReducedStructureInfo.size()-1?1:0);
    int LastElement2= (Block2==(int)ThisData.ReducedStructureInfo.size()-1?1:0);


    int ThisSize1=ThisInfo1->endIndex-ThisInfo1->startIndex+LastElement1;;
    int ThisSize2=ThisInfo2->endIndex-ThisInfo2->startIndex+LastElement2;


    for(int i=0;i<ThisSize1;i++)
    {
        int StartValue=(Block1==Block2? i+1: 0 );

        double iFreq=AlleleFreq[i+ThisInfo1->startIndex];

        for(int j=StartValue;j<ThisSize2;j++)
        {


            int NoMarkerinBetween = j+ThisInfo2->startIndex - i-ThisInfo1->startIndex ;
            int NoBPinBetween = ThisData.VariantList[j+ThisInfo2->startIndex].bp - ThisData.VariantList[i+ThisInfo1->startIndex].bp ;


//cout<<NoMarkerinBetween<<" " <<LdWindow<<endl;
            if(LdWindow!=-999999999)
            {
                if(NoMarkerinBetween>LdWindow)
                    continue;
            }

//cout<<ThisSize1<<" "<<ThisSize2<<endl;


            if(LdWindowBp!=-999999999)
            {
                if(NoBPinBetween>LdWindowBp)
                    continue;
            }


            if(LdWindowMin!=-999999999)
            {
                if(NoMarkerinBetween<LdWindowMin)
                    continue;
            }

            if(LdWindowBpMin!=-999999999)
            {
                if(NoBPinBetween<LdWindowBpMin)
                    continue;
            }


            double jFreq=AlleleFreq[j+ThisInfo2->startIndex];
            double ThisD=HapR2[i][j];
            ThisD /= (double)ThisData.numHaplotypes;
            ThisD -= (iFreq*jFreq);
            double ThisR2=(ThisD*ThisD)/( (iFreq)*(1-iFreq)*(jFreq)*(1-jFreq));

            if(MinHapR2<=0.0 || ThisR2>=MinHapR2)
            {

            ifprintf(LDFile, "%s\t",ThisData.VariantList[i+ThisInfo1->startIndex].chr.c_str());
                ifprintf(LDFile, "%d\t",ThisData.VariantList[i+ThisInfo1->startIndex].bp);
                ifprintf(LDFile, "%d\t",ThisData.VariantList[j+ThisInfo2->startIndex].bp);
                ifprintf(LDFile, "%d\t",ThisData.numHaplotypes);

                ifprintf(LDFile, "%f\t",ThisR2);
                ifprintf(LDFile, "%f\t",ThisD);
                ifprintf(LDFile, "%f\t",CalculateDPrime(ThisD,iFreq,jFreq));
                ifprintf(LDFile, "\n");
            }
        }
    }
    ifclose(LDFile);

}
void Analyse::CalculateHapSelfBlock(vector<vector<double> > &HapR2, HaplotypeSet &ThisData, int Block1)
{

//    HapR2.clear();
    ReducedHaplotypeInfo *ThisInfo=&ThisData.ReducedStructureInfo[Block1];


    int LastElement= (Block1==(int)ThisData.ReducedStructureInfo.size()-1?1:0);


    int ThisSize=ThisInfo->endIndex-ThisInfo->startIndex+LastElement;

    HapR2.resize(ThisSize);
    for(int i=0;i<ThisSize;i++)
    {
        HapR2[i].resize(ThisSize,0.0);
    }


    for (int i = 0; i<(int)ThisInfo->uniqueCardinality.size(); i++)
    {
        for(int marker1=ThisInfo->startIndex;marker1<ThisInfo->endIndex+LastElement;marker1++)
        {
            for(int marker2=marker1+1;marker2<ThisInfo->endIndex+LastElement;marker2++)
            {
                if(ThisInfo->uniqueHaps[i][marker1-ThisInfo->startIndex] && ThisInfo->uniqueHaps[i][marker2-ThisInfo->startIndex] )
                {
                    HapR2[marker1-ThisInfo->startIndex][marker2-ThisInfo->startIndex]+=ThisInfo->uniqueCardinality[i];
                }
            }

        }

    }

//    PrintHapR2Output(ThisData,Block1,Block1);

}
void Analyse::CalculatHapR2(HaplotypeSet &ThisData)
{
    CalculateFreq(ThisData);

    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                            HAPLOTYPE R2 PROCESSING                            "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;



    std::cout << "\n Calculating Haplotype R2 ... " <<endl;
    std::cout << "\n Writing Haplotype R2 to [.hap.ld] file  : " <<OutputFileName+".hap.ld"<< endl<<endl;


    IFILE freqFile = ifopen(OutputFileName + ".hap.ld" , "wb");
    ifprintf(freqFile, "#CHROM\t");
    ifprintf(freqFile, "POS1\t");
    ifprintf(freqFile, "POS2\t");
    ifprintf(freqFile, "N_CHR\t");

    ifprintf(freqFile, "R^2\t");
    ifprintf(freqFile, "D\t");
    ifprintf(freqFile, "Dprime\n");
    ifclose(freqFile);


    int i,k;
//    for(k=0;k<(int)ThisData.ReducedStructureInfo.size();k++)
//    {
//        CalculateHapSelfBlock(ThisData,k);
//    }

    #ifdef _OPENMP
        omp_set_num_threads(10);
    #endif

    #pragma omp parallel for
    for(i=0;i<(int)ThisData.ReducedStructureInfo.size();i++)
        {

         #pragma omp parallel for

         for(k=i;k<(int)ThisData.ReducedStructureInfo.size();k++)
        {
//            cout<<i<<" "<< k<<endl;


            vector<vector<double> > HapR2;
            HapR2.clear();


            if(k==i)
                CalculateHapSelfBlock(HapR2,ThisData,k);
            else
                CalculateHapCrossBlock(HapR2,ThisData,i,k);


            #pragma omp critical
            PrintHapR2Output(HapR2,ThisData,i,k);
        }
        }

}




