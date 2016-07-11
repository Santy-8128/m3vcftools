#include "HaplotypeSet.h"



void HaplotypeSet::GetInfo(variant &tempVariant,char *pch3)
{

    tempVariant.FullInfo="";

    char *end_str1;
    char  *pch  = strtok_r (pch3,";", &end_str1);
    pch  = strtok_r (NULL,";", &end_str1);
    if(pch!=NULL)
    {

        tempVariant.FullInfo+=(string)pch;
        pch  = strtok_r (NULL,";", &end_str1);
        while(pch!=NULL)
        {

            tempVariant.FullInfo+=(";"+(string)pch);
//            cout<<tempVariant.FullInfo<<endl;
            pch= strtok_r (NULL, ";", &end_str1);
        }


    }

}
bool HaplotypeSet::ReadBlockHeader(ReducedHaplotypeInfo &tempBlocktoCheck, string &line, int blockIndex,
                                    int &tempVarCount, int &tempRepCount, String filename)
{
    string tempString2,tempString3;

    char *end_str_new;
    char * pch;


    pch = strtok_r ((char*)line.c_str(),"\t", &end_str_new);
    tempBlocktoCheck.chr=string(pch);

    pch = strtok_r (NULL,"\t", &end_str_new);
    string tempBlockPos=string(pch);

    char *end_str_new1;
    char *pch1;

    if( m3vcfVersion>=2.0)
    {
        tempBlocktoCheck.startIndex=atoi(pch);
        pch = strtok_r (NULL,"\t", &end_str_new);
        tempBlocktoCheck.endIndex=atoi(pch);
    }
    else
    {

        pch1 = strtok_r ((char*)tempBlockPos.c_str(),"-", &end_str_new1);
        tempBlocktoCheck.startIndex=atoi(pch1);


        pch1 = strtok_r (NULL,"-", &end_str_new1);
        tempBlocktoCheck.endIndex=atoi(pch1);
    }


    bool flag=1-FilterIndicator(tempBlocktoCheck);
//abort();


    pch = strtok_r (NULL,"\t", &end_str_new);
    string blockName=string(pch);

    pch = strtok_r (NULL,"\t", &end_str_new);
    pch = strtok_r (NULL,"\t", &end_str_new);
    pch = strtok_r (NULL,"\t", &end_str_new);
    pch = strtok_r (NULL,"\t", &end_str_new);

    pch = strtok_r (NULL,"\t", &end_str_new);
    tempString2=string(pch);

    char *pch2,*end_str_new2;
    pch2 = strtok_r ((char*)tempString2.c_str(),";", &end_str_new2);

    if(Check)
    {
        stringstream strs;
        strs<<(blockIndex+1);

        tempString3="B"+ (string)(strs.str());

        if(tempString3.compare((string)pch2)!=0)
        {
            cout<<endl<<" Error in INFO column (Block Identifier) for block : "<<blockName <<endl;
            cout<<" Block Identifier should be : "<< tempString3<<" but is : "<<pch2<<endl<<endl;
            printm3VcfErr(filename);
        }

    }



    tempString3 = (string)strtok_r (NULL,";", &end_str_new2);
    char *pch3,*end_str_new3;
    pch3 = strtok_r ((char*)tempString3.c_str(),"=", &end_str_new3);
    pch3 = strtok_r (NULL,"=", &end_str_new3);
    tempVarCount=atoi(pch3);

    if(flag)
        return true;



    tempString3 = (string)strtok_r (NULL,";", &end_str_new2);
    end_str_new3=NULL;
    pch3 = strtok_r ((char*)tempString3.c_str(),"=", &end_str_new3);
    pch3 = strtok_r (NULL,"=", &end_str_new3);
    tempRepCount=atoi(pch3);


    tempBlocktoCheck.uniqueCardinality.resize(tempRepCount,0.0);
//    tempBlocktoCheck.InvuniqueCardinality.resize(tempRepCount,0.0);

    tempBlocktoCheck.uniqueIndexMap.resize(numHaplotypes);

    pch = strtok_r (NULL,"\t", &end_str_new);
    pch = strtok_r (NULL,"\t", &end_str_new);



    int check=0,inputcheck=0;
    while(pch!=NULL)
    {
        if(HapInputIndex[check++])

        {
            int tempval=atoi(pch);
            tempBlocktoCheck.uniqueIndexMap[inputcheck++]=tempval;
            tempBlocktoCheck.uniqueCardinality[tempval]++;
        }
          pch = strtok_r (NULL,"\t", &end_str_new);

    }

    if(check!=OrignumHaplotypes)
    {
        cout<<endl<<" Error in Data consistency !!! "<<endl;
        cout<<" Number of Haplotypes should be : "<< OrignumHaplotypes<<
        ", but "<<check<<" indices found in block"<<blockName << endl<<endl;
        printm3VcfErr(filename);
    }

//    for (int i = 0; i < tempRepCount; i++)
//    {
//        tempBlocktoCheck.InvuniqueCardinality[i]=1.0/(float)tempBlocktoCheck.uniqueCardinality[i];
//    }


    return false;

}
void HaplotypeSet::ReadThisBlock(int blockIndex, int tempRepCount, ReducedHaplotypeInfo &tempBlock,
                                 int tempVarCount,IFILE m3vcfxStream, int &NoMarkersImported,
                                 String filename)
{

    string line, tempString;
    int blockEnterFlag=0;

//abort();

    vector<bool> UniqueHapInputIndex(tempRepCount,false);

    int newRepCount=0;
    vector<int> newCardinality(newRepCount);
    vector<int> newIndexMap(numHaplotypes);

    vector<int> mapper(tempRepCount,-9);

    for(int index=0;index<tempRepCount;index++)
    {
        if(tempBlock.uniqueCardinality[index]>0)
        {

            mapper[index]=newRepCount;
            UniqueHapInputIndex[index]=true;
            newRepCount++;
            newCardinality.push_back(tempBlock.uniqueCardinality[index]);
        }
    }
    for(int index=0;index<numHaplotypes;index++)
    {

        int temp=mapper[tempBlock.uniqueIndexMap[index]];

        assert(temp!=-9);

        tempBlock.uniqueIndexMap[index]=temp;

    }



    tempBlock.uniqueHaps.resize(newRepCount);
    tempBlock.uniqueCardinality=newCardinality;


    for(int tempIndex=0;tempIndex<tempVarCount;tempIndex++)
    {


        int flag=0;
        line.clear();
        m3vcfxStream->readLine(line);



        string tempString2,tempString3;
        char  *end_str_new3;
        char * pch3;

        pch3 = strtok_r ((char*)line.c_str(),"\t", &end_str_new3);
        string tempChr=(string)pch3;


        pch3 = strtok_r (NULL,"\t", &end_str_new3);
        int tempPos=atoi(pch3);

        if( m3vcfVersion>=2.0)
        {
            pch3 = strtok_r (NULL,"\t", &end_str_new3);
        }

        pch3 = strtok_r (NULL,"\t", &end_str_new3);
        string tempName=(string)(pch3);


        variant tempVariant;

        tempVariant.assignValues(tempName,tempChr,tempPos);

        pch3 = strtok_r (NULL,"\t", &end_str_new3);
        string tempString98=(string)(pch3);

        pch3 = strtok_r (NULL,"\t", &end_str_new3);
        string tempString99=(string)(pch3);

        tempVariant.assignRefAlt(tempString98,tempString99);


        pch3 = strtok_r (NULL,"\t", &end_str_new3);
        tempVariant.Qual=atof(pch3);

        pch3 = strtok_r (NULL,"\t", &end_str_new3);
        tempVariant.Filter=(string)(pch3);

        pch3 = strtok_r (NULL,"\t", &end_str_new3);

        GetInfo(tempVariant,pch3);

//        tempVariant.FullInfo=(string)(pch3);


        if(Check)
        {

            tempString=(string)(pch3);
            char *pch4, *end_str_new4;
            pch4=strtok_r ((char*)tempString.c_str(),";", &end_str_new4);
            stringstream strs1,strs2;
            strs1<<(blockIndex+1);
            strs2<<(tempIndex+1);
            tempString3="B"+ (string)(strs1.str())+ ".M"+ (string)(strs2.str());
            if(tempString3.compare((string)pch4)!=0)
            {
                cout<<endl<<" Error in INFO column (Block Identifier) for variant : "<<tempName <<endl;
                cout<<" Block Identifier should be : "<< tempString3<<" but is : "<<pch4<<endl<<endl;
                printm3VcfErr(filename);
            }

        }




        flag=1-FilterIndicator(tempVariant);

        if(tempIndex==(tempVarCount-1) && NoMarkersImported==0)
            flag=1;


        if(flag==1)
            continue;



        if(tempIndex<(tempVarCount-1) || blockIndex==(NoBlocks-1))
        {

            stringstream strs1;
            strs1<<(tempVariant.bp);

            VariantList.push_back(tempVariant);


			markerName.push_back( tempVariant.chr+":"+(string)strs1.str());

            NoMarkersImported++;
        }


        if(blockEnterFlag==0)
        {
            tempBlock.startIndex=NoMarkersImported-1;
            blockEnterFlag=1;
        }

        tempBlock.endIndex=NoMarkersImported;
        if(blockIndex==(NoBlocks-1))
            tempBlock.endIndex--;



//if(NoMarkersImported==5000)
//    abort();


        pch3 = strtok_r (NULL,"\t", &end_str_new3);
        tempString=(string)(pch3);




        int inputIndex=0;
        for(int index=0;index<tempRepCount;index++)
        {

            if(UniqueHapInputIndex[index])
            {
                char t=tempString[index];
                tempBlock.uniqueHaps[inputIndex++].push_back(t=='0'? false:true);
            }

        }




    }

//abort();


}
string HaplotypeSet::readFirstFileAfterHeader(IFILE m3vcfxStream)
{
    string line,tempString;
    m3vcfxStream->readLine(line);

    if(line.compare("##fileformat=M3VCF")!=0 && line.compare("##fileformat=OPTM")!=0)
    {
        cout<<" Incorrect Header Information : "<<line<<endl;
        cout<<" Header line should be : ##fileformat=M3VCF "<<endl<<endl;
        printm3VcfErr(FileName);
    }

    bool Header=true;
    char * pch;
    char *end_str1;

    while(Header)
    {


        line.clear();
        m3vcfxStream->readLine(line);
        if(line.substr(0,2).compare("##")==0)
            Header=true;
        else
            break;
        tempString = (string) strtok_r ((char*)line.c_str(),"=", &end_str1);
        pch = strtok_r (NULL, "=", &end_str1);

        if(tempString.compare("##n_blocks")==0)
        {
            NoBlocks=atoi(pch);
            continue;
        }
        else if(tempString.compare("##version")==0)
        {
            m3vcfVersion=atof(pch);
            continue;
        }
        else if(tempString.compare("##n_haps")==0)
        {
            numHaplotypes=atoi(pch);
            continue;
        }
        else if(tempString.compare("##n_markers")==0)
        {
            numMarkers=atoi(pch);
            continue;

        }
        else if(tempString.compare("##chrxRegion")==0)
        {
            finChromosome="X";
            tempString = (string) pch;
            if(tempString.compare("NonPseudoAutosomal")==0)
                PseudoAutosomal=false;
            else if(tempString.compare("PseudoAutosomal")==0)
                PseudoAutosomal=true;
            else
            {
                cout << "\n Inconsistent Tag for Chr X. "<<endl;
                cout << " Please check the file properly..\n";
                cout << " Program Aborting ... "<<endl;
                abort();

            }
            continue;
        }
    }


    std::cout << " Number of Haplotypes in M3VCF File                  : " << numHaplotypes << endl;
    std::cout << " Number of Markers in M3VCF File                     : " << numMarkers << endl;


    return line;
}
bool HaplotypeSet::ReportSampleCount(string &line)
{

    char * pch;
    char *end_str2;

    int colCount=0;

    pch = strtok_r ((char*)line.c_str(),"\t", &end_str2);

    int lastIndex=-1;
    int MaleCount=0,FemaleCount=0;

    string tempString2;

    int SkipCol=m3vcfVersion>=2.0?10:9;

    while(pch!=NULL)
    {
        colCount++;

        if(colCount>SkipCol)
        {

            tempString2=string(pch);
            HaploName.push_back(tempString2.substr(0,tempString2.size()-6));

            int HapNo=atoi(tempString2.substr(tempString2.size()-1,tempString2.size()).c_str());

            if(HapNo==1)
            {
                individualName.push_back(tempString2.substr(0,tempString2.size()-6));
                SampleNoHaplotypes.push_back(1);
                lastIndex++;
                MaleCount++;
            }
            else if (HapNo==2)
            {
                   SampleNoHaplotypes[lastIndex]=2;
                   FemaleCount++;

            }
            else
            {
                cout << "\n Inconsistent HAPLO NO for Chr X = "<<HapNo<<endl;
                cout << " Please check the file properly..\n";
                cout << " Program Aborting ... "<<endl;
                abort();
            }

        }

        pch = strtok_r (NULL,"\t", &end_str2);
    }

    MaleCount=MaleCount-FemaleCount;
    numSamples=(int) individualName.size();
    OrigNumSamples=numSamples;
    OrignumHaplotypes=numHaplotypes;

    if((int)HaploName.size()!=numHaplotypes)
    {
        cout<<endl<<" Error in Data consistency !!! "<<endl<<endl;
        cout<<" Number of Haplotypes should be : "<< numHaplotypes<<", but "
        <<HaploName.size()<<" Haplotypes found in header row !!!"<< endl<<endl;
        return false;
    }


    std::cout << " Number of Samples in M3VCF File                     : " << numSamples << endl;



    if(finChromosome=="X")
    {


        cout<<"\n Chromosome X Detected !!!";
        cout<<endl;

        if(!PseudoAutosomal)
        {
            std::cout << " Number of MALE Samples (Haplotypes)                 : " << MaleCount << " ("<<MaleCount <<") "<<endl;
            std::cout << " Number of FEMALE Samples (Haplotypes)               : " << FemaleCount<< " ("<<FemaleCount*2 <<") "<<endl;
        }
        else
        {
            std::cout << "\n All " << FemaleCount<<" samples have two alleles on Chromosome X !!! "<< endl;
            cout<<" Cannot determine number of male/female samples from Pseudo-Autosomal Region ..."<<endl;
        }
    }
    return true;

}
bool HaplotypeSet::readm3vcfFile(String filename)
{
    string line,tempString,tempString2,tempString3,tempName;
    variant tempVariant;
    variant tempVariant2;
    int NoMarkersImported=0;
    int InitialNMarkers=0,blockIndex,tempVarCount,tempRepCount;
    FileName=filename;

    cout<<"\n Loading Input Haplotype Set from M3VCF files : "<<filename<<endl<<endl;
    optEndPoints.clear();
    individualName.clear();
    SampleNoHaplotypes.clear();
    HaploName.clear();
    finChromosome=".";


    IFILE m3vcfxStream = ifopen(filename, "r");

    if(m3vcfxStream)
    {

        line=readFirstFileAfterHeader(m3vcfxStream);

        InitialNMarkers=numMarkers;

        if(!ReportSampleCount(line))
            printm3VcfErr(filename);

        cout<<"\n After sample filtering ..."<<endl<<endl;
        std::cout << " Number of Markers to be read                        : " << InitialNMarkers<< endl;
        UpdateInputSampleList(filename);


        cout<<"\n Reading  "<<numHaplotypes<< " haplotypes from data ..."<<endl<<endl;


        for(blockIndex=0;blockIndex<NoBlocks;blockIndex++)
        {


            if (blockIndex % 1000 == 0)
			{
			    printf("  Loading Block %d out of %d blocks to be loaded... [%.1f%%] "
                , blockIndex + 1, NoBlocks, 100*(double)(blockIndex + 1)/(double)NoBlocks);
                cout<<endl;
            }

            line.clear();
            m3vcfxStream->readLine(line);
            ReducedHaplotypeInfo tempBlock;
            tempBlock.endIndex=0;
            tempBlock.startIndex=0;

            if(ReadBlockHeader(tempBlock, line, blockIndex, tempVarCount, tempRepCount, filename))
            {
                for(int tempIndex=0;tempIndex<tempVarCount;tempIndex++)
                    m3vcfxStream->discardLine();

                continue;
            }

            tempBlock.endIndex=0;
            tempBlock.startIndex=0;


            ReadThisBlock(blockIndex, tempRepCount,
                          tempBlock, tempVarCount,
                           m3vcfxStream, NoMarkersImported, filename);


            if(tempBlock.endIndex>tempBlock.startIndex)
            {
                optEndPoints.push_back(tempBlock.startIndex);
                finChromosome=tempBlock.chr;

                ReducedStructureInfo.push_back(tempBlock);
            }

        }

        if(ReducedStructureInfo.size()>0)
            optEndPoints.push_back(ReducedStructureInfo[ReducedStructureInfo.size()-1].endIndex);

    }
    else
    {
        cout<<" Following M3VCF File Not Available : "<<filename<<endl;
        cout<<" Program Exiting ... "<<endl<<endl;
        return false;
    }


    numMarkers=VariantList.size();

    cout<<endl<<" Reference Haplotype information succesfully recorded. "<<endl;

	std::cout << "\n Number of Markers Recorded                          : " << numMarkers << endl;
    std::cout << " Number of Haplotypes Recorded                       : " << numHaplotypes << endl;

    ifclose(m3vcfxStream);

    if(numMarkers<2)
    {
        cout << "\n None/Single marker left after filtering from Input File : "<<filename<<endl;
		cout << " Please check the file or the filtering options properly ...\n";
		cout << " Program Aborting ... "<<endl;
		return false;
    }



    std::cout << "\n Haplotype Set successfully loaded from M3VCF File     : " << filename << endl;

    return true;




}
void HaplotypeSet::printm3VcfErr(String filename)
{
    cout<<" Error in M3VCF File !!! "<<endl;
    cout<<" Please re-construct the following [.m3vcf] file using m3vcftools and try again ..."<<endl;
    cout<<" [ "<< filename<<" ] "<<endl;
    cout<<" Contact author if problem still persists : sayantan@umich.edu "<<endl;
    cout<<" Program Exiting ..."<<endl<<endl;
    abort();
}






string HaplotypeSet::DetectInputFileType(String filename)
{
    IFILE fileStream = ifopen(filename, "r");
    string line;
    if(fileStream)
    {
        fileStream->readLine(line);
        if(line.length()<1)
            return "Invalid";
        string tempString;
        tempString=(line.substr(0,17));

        char temp[tempString.length() + 1];
        std::strcpy(temp,tempString.c_str());
        for (char *iter = temp; *iter != '\0'; ++iter)
        {
           *iter = std::tolower(*iter);
        }
        if(((string)temp).compare("##fileformat=m3vc")==0)
        {
            return "m3vcf";
        }
        else if(((string)temp).compare("##fileformat=vcfv")==0)
        {
            return "vcf";
        }
        else
            return "Invalid";

    }
    else
    {
        return "NA";
    }

    ifclose(fileStream);

    return "NA";
}
bool HaplotypeSet::LoadInputHaplotypes(String filename)
{


    string FileType=DetectInputFileType(filename);

    if(FileType.compare("m3vcf")==0)
    {
        cout<<"\n Input Format = M3VCF (Minimac3 VCF File) "<<endl;
        vcfType=false;
        m3vcfxType=true;
        return readm3vcfFile(filename);
    }
    else if(FileType.compare("NA")==0)
    {
        cout<<"\n Following File File Not Available : "<<filename<<endl;
        return false;
    }
    else if(FileType.compare("Invalid")==0)
    {

        cout << "\n Input File provided by \"--in\" must be a VCF or M3VCF file !!! \n";
        cout << " Please check the following file : "<<filename<<endl;
        return false;
    }
    cout<<"\n Input Format = VCF (Variant Call Format) "<<endl;
    vcfType=true;
    m3vcfxType=false;
    return readvcfFile(filename);

}



bool HaplotypeSet::FilterIndicator(variant &ThisRecord)
{

//    string filterValue=ThisRecord.Filter;
//    if (filterValue.compare("PASS") != 0)
//    {
//        NoFAIL++;
//        if(FilterPASS)
//            return false;
//    }
//    if (KeepFiltered!="" && filterValue.compare(KeepFiltered) != 0)
//    {
//        return false;
//    }
//    if (RemoveFiltered!="" && filterValue.compare(RemoveFiltered) == 0)
//    {
//        return false;
//    }


    if(FilterCHR!="")
    {
        if(ThisRecord.chr.compare(FilterCHR.c_str())!=0)
            return false;
        else
        {
            if(FilterEND>0)
            {
                if(ThisRecord.bp>FilterEND || ThisRecord.bp<FilterSTART)
                    return false;
            }
            else
                if(ThisRecord.bp<FilterSTART)
                    return false;
        }

    }

    return true;
}
bool HaplotypeSet::FilterIndicator(VcfRecord &ThisRecord)
{

    if (ThisRecord.getNumAlts()>1)
    {
        NoMULTIALLELIC++;
//        return false;
    }


    // CHECK FOR FILTER VALUES
    {
    string filterValue=ThisRecord.getFilter().getString(0);
    if (filterValue.compare("PASS") != 0)
    {
        NoFAIL++;
        if(FilterPASS)
            return false;
    }
    if (KeepFiltered!="" && filterValue.compare(KeepFiltered) != 0)
    {
        return false;
    }
    if (RemoveFiltered!="" && filterValue.compare(RemoveFiltered) == 0)
    {
        return false;
    }
    }


    string refAllele = ThisRecord.getRefStr();
    string altAllele = ThisRecord.getAltStr();
    string cno=ThisRecord.getChromStr();
    int bp=ThisRecord.get1BasedPosition();

    if(FilterCHR!="")
    {
        if(cno.compare(FilterCHR.c_str())!=0)
            return false;
        else
        {
            if(FilterEND>0)
            {
                if(bp>FilterEND || bp<FilterSTART)
                    return false;
            }
            else
                if(bp<FilterSTART)
                    return false;
        }

    }

    return true;
}
bool HaplotypeSet::FilterIndicator(ReducedHaplotypeInfo &ThisRecord)
{
    if(FilterCHR!="")
    {
        if(ThisRecord.chr.compare(FilterCHR.c_str())!=0)
            return false;
        else
        {
            if(FilterEND>0)
            {
                if(ThisRecord.startIndex>FilterEND)
                    return false;
            }
            else
                if(ThisRecord.endIndex<FilterSTART)
                    return false;
        }

    }

    return true;

}




bool HaplotypeSet::CreateSampleNoHaplotypes(int &OrignumHaplotypes, String filename)
{

    VcfFileReader inFile;
	VcfHeader header;
	VcfRecord record;

    SampleNoHaplotypes.clear();
    int maleCount=0,femaleCount=0;
    numHaplotypes=0;


    if(finChromosome!="X")
    {
        SampleNoHaplotypes.resize(numSamples,2);
    }
    else
    {
        cout<<"\n Chromosome X Detected !!! \n";

        PseudoAutosomal=false;
        inFile.open(filename, header);
        inFile.setSiteOnly(false);
        inFile.readRecord(record);

        for (int i = 0; i<(numSamples); i++)
        {
            if(record.getNumGTs(i)==0)
            {
                std::cout << "\n Empty Value for Individual : " << individualName[i] << " at Marker : " << VariantList[0].name << endl;
                std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
                return false;
            }
            else
            {
                if(record.getNumGTs(i)==1)
                    maleCount++;
                else if(record.getNumGTs(i)==2)
                    femaleCount++;
                else
                {
                    std::cout << "\n Abnormal Ploidy for for Individual : " << individualName[i] << " at Marker : " << VariantList[0].name << endl;
                    std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
                    return false;
                }

                SampleNoHaplotypes.push_back(record.getNumGTs(i));
                OrignumHaplotypes+=record.getNumGTs(i);
            }

        }
        inFile.close();

        if(maleCount>0)
        {
            std::cout << " Number of MALE Samples (Haplotypes)                 : " << maleCount << " ("<<maleCount <<") "<<endl;
            std::cout << " Number of FEMALE Samples (Haplotypes)               : " << femaleCount<< " ("<<femaleCount*2 <<") "<<endl;
        }
        else
        {
            std::cout << "\n All " << numSamples<<" samples have two alleles on Chromosome X !!! "<< endl;
            cout<<" Cannot determine number of male/female samples from Pseudo-Autosomal Region ..."<<endl;
            PseudoAutosomal=true;
        }

    }
    return true;
}



void HaplotypeSet::GetHeader(VcfHeader& Header)
{

    HeaderInfo.clear();
    HeaderInfo.resize(Header.getNumMetaLines());

    for(int i=0;i<(int)HeaderInfo.size();i++)
    {
        HeaderInfo[i]=Header.getMetaLine(i);
    }
}


bool HaplotypeSet::readvcfFile(String filename)
{

    VcfFileReader inFile;
	VcfHeader header;
	VcfRecord record;
    m3vcfVersion=1.0;


	if (!inFile.open(filename, header))
	{
		cout << "\n Program could NOT open file : " << filename << endl;
		return false;
	}

    std::cout << "\n Loading Input Haplotype Set from VCF File       : " << filename << endl;

    GetHeader(header);

    PrintStartIndex=0;
    PrintEndIndex=0;
    int     numReadRecords = 0;
	int numtoBeWrittenRecords = 0;

	vector<bool> importIndex;
	string prevID="",prevFullID="";
//	bool isInDel = false;
	inFile.setSiteOnly(true);
	numSamples = header.getNumSamples();


	for (int i = 0; i < numSamples; i++)
	{
		string tempName(header.getSampleName(i));
		individualName.push_back(tempName);
	}

    if(individualName.size()==0)
    {
        cout << "\n No Sample Names found in VCF Input File : "<<filename<<endl;
		cout << " Please check the file properly..\n";
		cout << " Program Aborting ... "<<endl;
		return false;
    }


    string refAllele,altAllele,PrefrefAllele,PrevaltAllele,cno,fixCno,id;
    int bp;
	cout << "\n Reading VCF File to calculate number of records ... \n";
    cout<<endl;


	std::cout << " Number of Samples in VCF File                       : " << numSamples << endl;
	// pre-calculate number of samples read and number of markers read to allocate memory.
	while (inFile.readRecord(record))
	{

        int flag = 1;
		++numReadRecords;

		flag=1-FilterIndicator(record);

        cno=record.getChromStr();
        bp=record.get1BasedPosition();


        // CHECK FOR MULTIPLE CHROMOSOMES, M3VCF TOOLS CAN ONLY ANALYZE SINGLE CHROMOSOMES
        if(FilterCHR=="" && fixCno!=cno && numReadRecords>1)
        {
            cout << "\n Error !!! Reference VCF File contains multiple chromosomes : "<<cno<<", "<<fixCno<<", ... "<<endl;
            cout << " Please use VCF file with single chromosome or specify chromosome using \"--chr\" option !!! "<<endl;
            cout << " Program Aborting ... "<<endl;
            return false;
        }
        fixCno=cno;


        string currID=record.getIDStr();
        stringstream strs1,strs2;
        strs1<<(cno);
        strs2<<(bp);

        // IMPORT FINAL RECORD
		if (flag == 0)
		{
			++numtoBeWrittenRecords;
			markerName.push_back((string)strs1.str()+":"+(string)strs2.str());
            variant thisVariant(currID,cno,bp);
            thisVariant.Qual=record.getQual();
            thisVariant.NoAltAlelles=record.getNumAlts();
//            cout<<thisVariant.NoAltAlelles<<endl;
            thisVariant.Filter=record.getFilter().getString();
//            *thisVariant.FullInfo=record.getInfo();

// VcfRecordInfo Temp=record.getInfo();


            VariantList.push_back(thisVariant);
			importIndex.push_back(true);

		}
		else
		{
			importIndex.push_back(false);
		}

	}


	inFile.close();


    if(FilterCHR=="")
        finChromosome=cno;
    else
        finChromosome=FilterCHR.c_str();


    std::cout << " Number of Markers in VCF File                       : " << numReadRecords << endl;

    if(numReadRecords==0)
    {
        cout << "\n No Variants found in VCF Input File : "<<filename<<endl;
		cout << " Please check the file properly..\n";
		cout << " Program Aborting ... "<<endl;
		return false;
    }


	numMarkers = numtoBeWrittenRecords;
    OrignumHaplotypes=0,OrigNumSamples=numSamples;

    if(!CreateSampleNoHaplotypes(OrignumHaplotypes,filename))
        return false;


    cout<<"\n After variant and sample filtering ..."<<endl<<endl;

    std::cout << " Number of Markers to be read                        : " << numtoBeWrittenRecords<< endl;

    if(numtoBeWrittenRecords<2)
    {
        cout << "\n None/Single marker left after filtering from Input File : "<<filename<<endl;
		cout << " Please check the file or the filtering options properly ...\n";
		cout << " Program Aborting ... "<<endl;
		return false;
    }


    UpdateInputSampleList(filename);


    std::cout <<"\n Starting to load data ...";

    std::cout <<"\n\n Number of Threads to be Used = "<<CPU<<endl;



    int    blockSize = 500;
    int    bufferSize = 100;
//    individualName.resize(numSamples);
    // open file again to start loading the data in the variable haplotype.

    cout<<endl;


//     #pragma omp parallel for

    VcfFileReader inFileBuffer;
    VcfHeader headerBuffer;
    VcfRecord recordBuffer;
    int readIndex = -1;

    inFileBuffer.open(filename, headerBuffer);
    inFileBuffer.setSiteOnly(false);

    vector<vector<String> > Haplotypes(CPU);
    vector<vector<int> > MaxAllele(CPU);



    int currentPiece=0;

    for(currentPiece=0;currentPiece<CPU;currentPiece++)
    {
        Haplotypes[currentPiece].resize(numHaplotypes);
        MaxAllele[currentPiece].resize(bufferSize);
    }
    currentPiece=0;


    vector<int> BufferPosList(CPU);


    int NoMarkersWritten=1;
    int NoMarkersRead=0;

    BufferPosList[0]=1;
    for(int BufferNo=1;BufferNo<CPU;BufferNo++)
    {
        BufferPosList[BufferNo]= BufferPosList[BufferNo-1]+(bufferSize-1);

    }
    ReducedStructureInfoBuffer.clear();
    ReducedStructureInfoBuffer.resize(CPU);

    ReducedStructureInfo.clear();
    int CurrentVariant=0;

//cout<<" LAORA2 "<<endl;
    while (inFileBuffer.readRecord(recordBuffer) && NoMarkersRead<numMarkers)
	{
		// work only with bi-allelic markers
		readIndex++;


		if (importIndex[readIndex])
		{

//cout<<" LAORA "<<endl;
        if (Haplotypes[currentPiece][0].Length()==1)
                {
                    printf("  Loading markers %d - %d  out of %d markers to be loaded... [%.1f%%] ",
                           NoMarkersWritten+currentPiece*(bufferSize-1),
                           min(NoMarkersWritten+(currentPiece+1)*(bufferSize-1),numMarkers),numMarkers,
                           100*(double)(NoMarkersWritten+currentPiece*(bufferSize-1))/numMarkers);
                    cout<<endl;
                }

                NoMarkersRead++;

                string refAllele = recordBuffer.getRefStr();
                string altAllele = recordBuffer.getAltStr();

                VariantList[NoMarkersRead-1].refAlleleString=refAllele;
                VariantList[NoMarkersRead-1].altAlleleString=altAllele;
                int haplotype_index = 0;
                int Inputhaplotype_index = 0;
                bool IsMissing=false;
                MaxAllele[currentPiece][CurrentVariant]=VariantList[NoMarkersRead-1].NoAltAlelles+1;

//                cout<<" WOW = "<<MaxAllele[currentPiece][CurrentVariant]<<"\t"<<VariantList[NoMarkersRead-1].NoAltAlelles+1<<endl;
                for (int i = 0; i<(OrigNumSamples); i++)
                {
                    int NoGt=recordBuffer.getNumGTs(i);

                    if(NoGt==0)
                    {
                        std::cout << "\n Empty Value for Individual : " << individualName[i] << " at Marker : " << VariantList[NoMarkersRead-1].name << endl;
                        std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
                        abort();
                    }
                    else if(NoGt==1 && finChromosome!="X")
                    {
                        std::cout << "\n Single Autosomal Haplotype for Individual : " << individualName[i] << " at Marker : " << VariantList[NoMarkersRead-1].name << endl;
                        std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
                        abort();
                    }
                    for (int j = 0; j<NoGt; j++)
                    {

                        int alleleIndex = recordBuffer.getGT(i, j);
                        if(haplotype_index>=OrignumHaplotypes && finChromosome=="X")
                        {
                            cout << "\n Error in Reference VCF File for Chromosome X !!! "<<endl;
                            cout << "\n Marker : " << VariantList[0].name << " has "<<OrignumHaplotypes <<" haplotypes while ";
                            cout << "Marker : " << VariantList[NoMarkersRead-1].name << " has "<< haplotype_index+1<<" haplotypes."<<endl;
                            cout << " VCF file seems to have both Pseudo-Autosomal region (PAR) and non-PAR of chromosome X. \n";
                            cout << " Please use only either of the two regions ... \n";
                            cout << " See web-page for m3vcftools (Chromosome X) for more details ...\n";
                            cout << " Program Aborting ... "<<endl;
                            abort();
                        }


                        if(SampleInputIndex[i])
                        {
                            if (alleleIndex<0)
                            {

                                cout<<" NO MISSING VALUE ALLOWED"<<endl;
                                abort();
                                if(!IsMissing)
                                {
                                    IsMissing=true;
//                                    MaxAllele[currentPiece][CurrentVariant]++;

                                }

//                                Haplotypes[currentPiece][Inputhaplotype_index]+= (char)('0');
                            }
                            else
                            {
                                Haplotypes[currentPiece][Inputhaplotype_index]+= (char)('0'+alleleIndex);

                            }
                            Inputhaplotype_index++;
                        }

                        haplotype_index++;
                    }


                }

                CurrentVariant++;
                if(haplotype_index<(OrignumHaplotypes-1) && finChromosome=="X")
                {
                    cout << "\n Error in Reference VCF File for Chromosome X !!! "<<endl;
                    cout << "\n Marker : " << VariantList[0].name << " has "<<OrignumHaplotypes <<" haplotypes while ";
                    cout << "Marker : " << VariantList[NoMarkersRead-1].name << " has "<< haplotype_index+1<<" haplotypes."<<endl;
                    cout << " VCF file seems to have both Pseudo-Autosomal region (PAR) and non-PAR of chromosome X. \n";
                    cout << " Please use only either of the two regions ... \n";
                    cout << " See web-page for m3vcftools (Chromosome X) for more details ...\n";
                    cout << " Program Aborting ... "<<endl;
                    abort();
                }



//                cout<< " THIS = "<<MaxAllele[currentPiece][CurrentVariant]<<endl;


                int NewPiece=CPU;
                if (Haplotypes[currentPiece][0].Length() == bufferSize)
                {

                    NewPiece=(currentPiece+1)%CPU;
                    vector<String> tempHaplotypes(numHaplotypes);
                    CurrentVariant=0;
                    if(NewPiece!=0)
                    {
                        for (int i = 0; i < numHaplotypes; i++)
                        {
                            tempHaplotypes[i]=Haplotypes[currentPiece][i][Haplotypes[currentPiece][i].Length()-1];
                            Haplotypes[NewPiece][i]=tempHaplotypes[i];
                        }
                        currentPiece=NewPiece;
                        MaxAllele[currentPiece][CurrentVariant]=VariantList[NoMarkersRead-1].NoAltAlelles+1;

                    }

                }

//abort();

                if(NewPiece==0 || NoMarkersRead==numMarkers)
                {

                    #pragma omp parallel for
                    for(int ThisPiece=0;ThisPiece<=currentPiece;ThisPiece++)
                    {

                        int LastflushPos=BufferPosList[ThisPiece]-1;
                        printf("     Processing Chunk %d for M3VCF ... \n",ThisPiece+1);


                        vector<int> index(numHaplotypes),oldIndex;
                        vector<int> previousDifference(numHaplotypes);
                        vector<int> previousPredecessor(numHaplotypes);
                        vector<int> firstDifference(numHaplotypes-1,0);
                        vector<int> cost(bufferSize+1,0);
                        vector<int> bestSlice(bufferSize+1,0);
                        vector<int> bestComplexity(bufferSize+1,0);
                        vector<vector<int> > bestIndex(bufferSize+1);

                //        vector<ReducedHaplotypeInfo> ReducedStructureInfoTemp;
                        ReducedStructureInfoBuffer[ThisPiece].clear();

                        findUnique RefUnique;

                        RefUnique.updateCoeffs(transFactor,cisFactor);
                        double blockedCost = 0.0;

                        for(int i=0;i<numHaplotypes;i++)
                            index[i]=i;


                        for(int length=1;length<=Haplotypes[ThisPiece][0].Length();length++)
                        {

                            int maxAllele=2; //MaxAllele[ThisPiece][length-1];


                            vector<int> offsets(maxAllele+1,0);
                            for (int i = 0; i < numHaplotypes; i++)
                                offsets[Haplotypes[ThisPiece][i][length - 1] - '0' + 1]++;
                            for (int i = 1; i < maxAllele; i++)
                                offsets[i] += offsets[i - 1];

//abort();
//                            offsets[2]+=offsets[1];


                            oldIndex = index;
                            for (int i = 0; i < numHaplotypes; i++)
                                {
                                    index[offsets[Haplotypes[ThisPiece][oldIndex[i]][length - 1] - '0']++] = oldIndex[i];
                                }

//abort();

                            RefUnique.UpdateDeltaMatrix(Haplotypes[ThisPiece], index, firstDifference, length, blockSize,
                                   oldIndex, previousPredecessor, previousDifference);

                            RefUnique.AnalyzeBlocks(index, firstDifference, length, blockSize,
                               cost, bestSlice, bestComplexity, bestIndex);

                        }


                        if(Haplotypes[ThisPiece][0].Length()>1)
                            blockedCost += RefUnique.FlushBlocks(ReducedStructureInfoBuffer[ThisPiece],
                                                                VariantList,LastflushPos,
                                                                Haplotypes[ThisPiece], MaxAllele[ThisPiece],
                                                                cost,
                                                                bestComplexity, bestSlice, bestIndex);

                    }


                    cout<<endl;
                    NoMarkersWritten+=(CPU*(bufferSize-1));


                    BufferPosList[0]=NoMarkersWritten;
                    for(int BufferNo=1;BufferNo<CPU;BufferNo++)
                    {
                        BufferPosList[BufferNo]= BufferPosList[BufferNo-1]+(bufferSize-1);

                    }

                    vector<String> tempHaplotypes(numHaplotypes);
                    for (int i = 0; i < numHaplotypes; i++)
                    {
                        tempHaplotypes[i]=Haplotypes[currentPiece][i][Haplotypes[currentPiece][i].Length()-1];
                        Haplotypes[0][i]=tempHaplotypes[i];
                    }

                    for(int ThisPiece=0;ThisPiece<CPU;ThisPiece++)
                    {
                        for(int jj=0;jj<(int)ReducedStructureInfoBuffer[ThisPiece].size();jj++)
                            {
                                ReducedStructureInfo.push_back(ReducedStructureInfoBuffer[ThisPiece][jj]);
                            }
                        ReducedStructureInfoBuffer[ThisPiece].clear();

                    }

                    currentPiece=0;


                }

        }


    }



//abort();



    std::cout << "\n Number of Markers Recorded                          : " << markerName.size() << endl;
    std::cout << " Number of Haplotypes Recorded                       : " << numHaplotypes << endl;




    optEndPoints.clear();
    int i;
    for(i=0;i<(int)ReducedStructureInfo.size();i++)
        {
            optEndPoints.push_back(ReducedStructureInfo[i].startIndex);
        }
    optEndPoints.push_back(ReducedStructureInfo[i-1].endIndex);


    if(individualName.size()==0)
    {
        cout << "\n No haplotypes recorded from VCF Input File : "<<filename<<endl;
		cout << " Please check the file properly..\n";
		cout << " Program Aborting ... "<<endl;
		return false;
    }

    numMarkers = markerName.size();
    std::cout << "\n Haplotype Set successfully loaded from VCF File     : " << filename << endl;
	inFile.close();

	return true;

}
void HaplotypeSet::readFile(vector<string> &Vector, String FileName)
{

    Vector.clear();

    if(FileName=="")
        {
            return;
        }

    IFILE fileStream = ifopen(FileName, "r");
    string line;

    if(fileStream)
    {
        while(fileStream->readLine(line)==0)
        {
            Vector.push_back(line.c_str());
            line.clear();

        }
    }
    ifclose(fileStream);
}


bool HaplotypeSet::UpdateInputSampleList(String filename)
{

    SampleInputIndex.clear();
    SampleInputIndex.resize(numSamples,true);


    vector<string> KeepList,RemoveList;

    readFile(KeepList,KeepSample);
    readFile(RemoveList,RemoveSample);
    maleindividualName.clear();
    femaleindividualName.clear();


    for(int index=0;index<numSamples;index++)
    {
        if(KeepList.size()>0)
            if(find(KeepList.begin(), KeepList.end(), individualName[index]) == KeepList.end())
                SampleInputIndex[index]=false;

        if(RemoveList.size()>0)
            if(find(RemoveList.begin(), RemoveList.end(), individualName[index]) != RemoveList.end())
                SampleInputIndex[index]=false;
    }


    vector<string> NewindividualName(0);
    vector<int> NewSampleNoHaplotypes(0);


    if(m3vcfxType)
    {

        HapInputIndex.clear();
        HapInputIndex.resize(numHaplotypes,false);
        int hapIndex=0;
        for(int index=0;index<numSamples;index++)
        {


            if(SampleInputIndex[index])
            {
                HapInputIndex[hapIndex]=true;
                if(SampleNoHaplotypes[index]==2)
                    HapInputIndex[hapIndex+1]=true;
            }
            hapIndex+=SampleNoHaplotypes[index];
        }
    }


    for(int index=0;index<numSamples;index++)
    {
        if(SampleInputIndex[index])
        {
            NewindividualName.push_back(individualName[index]);
            NewSampleNoHaplotypes.push_back(SampleNoHaplotypes[index]);


            if(finChromosome=="X" && !PseudoAutosomal){
                if(SampleNoHaplotypes[index]==1)
                    maleindividualName.push_back(individualName[index]);
                else
                    femaleindividualName.push_back(individualName[index]);
            }
        }
    }


    std::cout << " Number of Samples to be read                        : " << (int)NewindividualName.size() << endl;



    if(NewindividualName.size()==0)
    {
        cout << "\n No haplotypes to be recorded from Input File : "<<filename<<endl;
		cout << " Please check Sample Filtering options properly ..\n";
		cout << " Program Aborting ... "<<endl;
		return false;
    }

    if(finChromosome=="X" && !PseudoAutosomal)
    {
        std::cout << " Number of MALE (FEMALE) Samples to be read          : " << (int)maleindividualName.size();
        std::cout <<" ("<<(int)femaleindividualName.size() <<")"<<endl;
    }


    individualName=NewindividualName;
    SampleNoHaplotypes=NewSampleNoHaplotypes;
    numSamples=(int)SampleNoHaplotypes.size();

    numHaplotypes=0;
    for (int i = 0; i<numSamples; i++)
        numHaplotypes+=SampleNoHaplotypes[i];

    return true;

}









































void HaplotypeSet::PrintDosageForVcfOutputForIDMaleSamples(IFILE vcfdose,
                                                           int MarkerIndex,bool majorIsReference,char refAllele)
{

    bool colonIndex;
    for(int hapId=0;hapId<(int)Dosage.size();hapId++)
        {
            bool a1=ImputedAlleles[hapId][MarkerIndex];
            colonIndex=false;
            ifprintf(vcfdose,"\t");
            if(GT)
            {
                int outAllele1=false;
                if(a1)
                    outAllele1=1;
                ifprintf(vcfdose,"%d",outAllele1);
                colonIndex=true;

            }


            double x=Dosage[hapId][MarkerIndex];

            if(DS)
            {

                if(colonIndex)
                    ifprintf(vcfdose,":");
                colonIndex=false;
                if(majorIsReference)
                    ifprintf(vcfdose,"%.3f",1-x);
                else
                    ifprintf(vcfdose,"%.3f",x);
                colonIndex=true;
            }


            if(GP)
            {

                if(colonIndex)
                    ifprintf(vcfdose,":");
                colonIndex=false;
                double p1,p3;
                if(majorIsReference)
                {
                    p1=x;
                    p3=(1-x);
                }
                else
                {
                    p3=x;
                    p1=(1-x);
                }

                ifprintf(vcfdose,"%.3f,%.3f",p1,p3);

            }


        }


}


void HaplotypeSet::PrintDosageForVcfOutputForID(IFILE vcfdose, int MarkerIndex,bool majorIsReference,char refAllele)
{

    bool colonIndex;
    for(int hapId=0;hapId<(int)Dosage.size()/2;hapId++)
        {

            bool a1=ImputedAlleles[2*hapId][MarkerIndex];
            bool a2=ImputedAlleles[2*hapId+1][MarkerIndex];
            colonIndex=false;
            ifprintf(vcfdose,"\t");
            if(GT)
            {
                int outAllele1=0,outAllele2=0;
                if(a1)
                    outAllele1=1;

                if(a2)
                    outAllele2=1;

                if(!unphasedOutput)
                    ifprintf(vcfdose,"%d|%d",outAllele1,outAllele2);
                else
                {
                    if((a1^a2)==1)
                        ifprintf(vcfdose,"0/1");
                    else if(a1 && a2)
                        ifprintf(vcfdose,"1/1");
                    else
                        ifprintf(vcfdose,"0/0");

                }
                colonIndex=true;

            }


            double x=Dosage[2*hapId][MarkerIndex];
            double y=Dosage[2*hapId+1][MarkerIndex];

            if(DS)
            {

                if(colonIndex)
                    ifprintf(vcfdose,":");
                colonIndex=false;
                if(majorIsReference)
                    ifprintf(vcfdose,"%.3f",1-x+1-y);
                else
                    ifprintf(vcfdose,"%.3f",x+ y);
                colonIndex=true;
            }


            if(GP)
            {

                if(colonIndex)
                    ifprintf(vcfdose,":");
                colonIndex=false;
                double p1,p2,p3;
                if(majorIsReference)
                {
                    p1=x*y;
                    p2=x*(1-y)+y*(1-x);
                    p3=(1-x)*(1-y);
                }
                else
                {
                    p3=x*y;
                    p2=x*(1-y)+y*(1-x);
                    p1=(1-x)*(1-y);
                }

                ifprintf(vcfdose,"%.3f,%.3f,%.3f",p1,p2,p3);

            }


        }


}



void HaplotypeSet::InitializePartialDosageForVcfOutputMaleSamples(int NHaps,
                                                       int NMarkers,vector<bool> &Format)
{
    numHaplotypes = NHaps;
    numMarkers = NMarkers;
    Dosage.resize(numHaplotypes);

    ImputedAlleles.resize(numHaplotypes);
    individualName.resize(numHaplotypes);
    for(int i=0;i<numHaplotypes;i++)
        {
            Dosage[i].resize(NMarkers,3.0);
            ImputedAlleles[i].resize(NMarkers,false);

        }

    GT=Format[0];
    DS=Format[1];
    GP=Format[2];



}

void HaplotypeSet::InitializePartialDosageForVcfOutput(int NHaps,
                                                       int NMarkers,vector<bool> &Format)
{
    numHaplotypes = NHaps;
    numMarkers = NMarkers;
    Dosage.resize(numHaplotypes);
    ImputedAlleles.resize(numHaplotypes);
    individualName.resize(numHaplotypes/2);
    for(int i=0;i<numHaplotypes;i++)
        {
            Dosage[i].resize(NMarkers,3.0);
            ImputedAlleles[i].resize(NMarkers,false);

        }

    GT=Format[0];
    DS=Format[1];
    GP=Format[2];



}


void HaplotypeSet::SaveDosageForVcfOutputSampleWise(int SamID,string &SampleName, vector<float> &dose1,vector<float> &dose2,
                                                    vector<bool> &impAlleles1,vector<bool> &impAlleles2)
{
    individualName[SamID]=SampleName;
    Dosage[2*SamID]=dose1;
    ImputedAlleles[2*SamID]=impAlleles1;
    Dosage[2*SamID+1]=dose2;
    ImputedAlleles[2*SamID+1]=impAlleles2;
}

void HaplotypeSet::SaveDosageForVcfOutputSampleWiseChrX(int SamID,string &SampleName, vector<float> &dose1,
                                                    vector<bool> &impAlleles1)
{
    individualName[SamID]=SampleName;
    Dosage[SamID]=dose1;
    ImputedAlleles[SamID]=impAlleles1;
}


void HaplotypeSet::SaveDosageForVcfOutput(int hapID,vector<float> dose,vector<bool> impAlleles)
{

    Dosage[hapID]=dose;
    ImputedAlleles[hapID]=impAlleles;

}




void HaplotypeSet::reconstructHaplotype(vector<bool> &reHaplotypes,int &index)
{
    if((int)reHaplotypes.size()!=numMarkers)
        {
            cout<<" SIZE MISMATCH "<<endl;
            abort();
        }

    int markerIndex=0,k;
    bool checkAllele=false;

    for(int j=0;j<(int)ReducedStructureInfo.size();j++)
    {

        for(k=0;k<(ReducedStructureInfo[j].endIndex-ReducedStructureInfo[j].startIndex);k++)
        {

            reHaplotypes[markerIndex++]=ReducedStructureInfo[j].uniqueHaps[ReducedStructureInfo[j].uniqueIndexMap[index]][k];
            if(k==0 && j>1)
            {
                if(checkAllele!=(ReducedStructureInfo[j].uniqueHaps[ReducedStructureInfo[j].uniqueIndexMap[index]][k]))
                {
                    cout<<index<<"\t"<<j<<"\t"<<k<<"\t"<<(int)checkAllele<<"\t"<<(int)(ReducedStructureInfo[j].uniqueHaps[ReducedStructureInfo[j].uniqueIndexMap[index]][k])<<endl;
                    cout<<" CHECK ALLELEE  MISMATCH "<<endl;
                    abort();
                }

            }

        }

        if(j==((int)ReducedStructureInfo.size()-1))
        {

            reHaplotypes[markerIndex]=ReducedStructureInfo[j].uniqueHaps[ReducedStructureInfo[j].uniqueIndexMap[index]][k];
        }
        else
        {
            checkAllele=ReducedStructureInfo[j].uniqueHaps[ReducedStructureInfo[j].uniqueIndexMap[index]][k];
        }

    }
}


bool HaplotypeSet::CheckValidChrom(string chr)
{
    bool result=false;

    if(MyChromosome!="" && chr==MyChromosome)
        return true;

    string temp[]={"1","2","3","4","5","6","7","8","9","10","11"
                    ,"12","13","14","15","16","17","18","19","20","21","22","X"};
    std::vector<string> ValidChromList (temp, temp + sizeof(temp) / sizeof(string) );

    for(int counter=0;counter<(int)ValidChromList.size();counter++)
        if(chr==ValidChromList[counter])
            result=true;

    return result;

}


string HaplotypeSet::DetectTargetFileType(String filename)
{
    IFILE fileStream = ifopen(filename, "r");

    string line;
    if(fileStream)
    {
        fileStream->readLine(line);
        if(line.length()<1)
            return "Invalid";
        string tempString;
        tempString=(line.substr(0,17));

        char temp[tempString.length() + 1];
        std::strcpy(temp,tempString.c_str());
        for (char *iter = temp; *iter != '\0'; ++iter)
        {
           *iter = std::tolower(*iter);
        }

        if(((string)temp).compare("##fileformat=vcfv")==0)
            {
                return "vcf";
            }
       return "Invalid";

    }
    else
    {
        return "NA";
    }

    ifclose(fileStream);

    return "NA";

}


bool HaplotypeSet::BasicCheckForTargetHaplotypes(String filename)
{

	VcfFileReader inFile;
	VcfHeader header;
	VcfRecord record;
	std::cout << "\n Performing basic file check on target/GWAS haplotype file : "<<filename << endl;

    std::cout << "\n Checking File ..." << endl;

	inFile.setSiteOnly(true);
    IFILE fileStream = ifopen(filename, "r");

    string line;
    if(!fileStream)
    {
		cout << "\n Program could NOT open file : " << filename << endl;
		return false;
	}
	else
    {
        std::cout << " File Exists ..." << endl;

        std::cout << "\n Checking File Type ..." << endl;

        fileStream->readLine(line);
        if(line.length()<1)
        {
            cout << "\n Target File provided by \"--haps\" must be a VCF file !!! \n";
            cout << " Please check the following file : "<<filename<<endl;
            return false;
        }

        string tempString;
        tempString=(line.substr(0,17));

        char temp[tempString.length() + 1];
        std::strcpy(temp,tempString.c_str());
        for (char *iter = temp; *iter != '\0'; ++iter)
        {
           *iter = std::tolower(*iter);
        }


        if(((string)temp).compare("##fileformat=vcfv")!=0)
        {
            cout << "\n Target File provided by \"--haps\" must be a VCF file !!! \n";
            cout << " Please check the following file : "<<filename<<endl<<endl;
            return false;
        }

    }

    ifclose(fileStream);

    std::cout << " VCF File Type Detected ..." << endl;

    std::cout <<"\n NOTE: Samples will be assumed to be phased irrespective "<<endl;
    cout<<       "       of GT delimiter (| or /) in VCF file !!! " << endl;


    std::cout << "\n Checking variant information ..." << endl;

    int numActualRecords=0,numReadRecords=0;
    int failFilter=0,duplicates=0,notBiallelic=0,inconsistent=0;
    string prevID,currID, refAllele,altAllele,PrefrefAllele,PrevaltAllele,cno,fixCno,id;
    inFile.open(filename, header);



    while (inFile.readRecord(record))
    {

        int flag=0;
        cno=record.getChromStr();
        int bp=record.get1BasedPosition();
        id=record.getIDStr();
        refAllele = record.getRefStr();
        altAllele = record.getAltStr();

        if(numActualRecords==0)
        {
            if(!CheckValidChrom(cno))
            {
                cout << "\n Error !!! Target VCF File contains chromosome : "<<cno<<endl;
                cout << " VCF File can only contain chromosomes 1-22 and X !!! "<<endl;
                cout << " Program Aborting ... "<<endl;
                return false;
            }

        }

        if (record.getNumAlts()>1)
		{
			notBiallelic++;
			flag = 1;
		}
		if (record.getFilter().getString(0).compare("PASS") != 0)
		{
			failFilter++;
			flag = 1;
		}

        stringstream strs3,strs4;
        strs3<<(cno);
        strs4<<(bp);

        currID=(string)strs3.str()+":"+(string)strs4.str();


        if(!CheckValidChrom(cno))
        {
            cout << "\n Error !!! Target VCF File contains chromosome : "<<cno<<endl;
            cout << " VCF File can only contain chromosomes 1-22 and X !!! "<<endl;
            cout << " Program Aborting ... "<<endl;
            return false;
        }
        if(fixCno!=cno && numReadRecords>0)
        {
            cout << "\n Error !!! Target VCF File contains multiple chromosomes : "<<cno<<", "<<fixCno<<", ... "<<endl;
            cout << " Please use VCF file with single chromosome !!! "<<endl;
            cout << " Program Aborting ... "<<endl;
            return false;
        }

        fixCno=cno;
        if (strlen(refAllele.c_str()) == 1 && strlen(altAllele.c_str()) == 1)
        {
            switch (refAllele[0])
            {
                case 'A': case 'a': ; break;
                case 'C': case 'c': ; break;
                case 'G': case 'g': ; break;
                case 'T': case 't': ; break;
                case 'D': case 'd': ; break;
                case 'I': case 'i': ; break;
                case 'R': case 'r': ; break;
                default:
                {
                        flag=1;
                        inconsistent++;
                }
            }
            if(flag==0)
                switch (altAllele[0])
                {
                    case '0':  ; break;
                    case 'A': case 'a': ; break;
                    case 'C': case 'c': ; break;
                    case 'G': case 'g': ; break;
                    case 'T': case 't': ; break;
                    case 'D': case 'd': ; break;
                    case 'I': case 'i': ; break;
                    case 'R': case 'r': ; break;
                    default:
                    {
                        flag=1;
                        inconsistent++;
                    }
                }
        }


        if(prevID==currID)
        {
            if(refAllele==PrefrefAllele && altAllele==PrevaltAllele)
            {
                duplicates++;
                cout << "\n Error !!! Duplicate Variant found chr:"<<cno<<":"<<bp<<" with identical REF = "<<refAllele <<" and ALT = "<<altAllele <<"\n";
                cout << " Program Aborting ... "<<endl;
                return false;
            }

        }

        prevID=currID;
        PrefrefAllele=refAllele;
        PrevaltAllele=altAllele;

        if(flag==0)
        {
            ++numReadRecords;
        }
        numActualRecords++;

    }



	std::cout << " "<<numActualRecords<<" variants in file with "<<numReadRecords<<" variants passing filters ..."<<endl<<endl;
    std::cout << " Checking sample information ..." << endl;

	inFile.close();

    inFile.open(filename, header);
    inFile.setSiteOnly(false);

	inFile.readRecord(record);

    std::cout << " "<<header.getNumSamples()<<" samples found in file ..."<<endl;

    if(header.getNumSamples()==0)
    {
        cout << "\n No haplotypes recorded from VCF Input File : "<<filename<<endl;
		cout << " Please check the file properly..\n";
		cout << " Program Aborting ... "<<endl;
		return false;
    }

    inFile.close();

    finChromosome=cno;

    if(finChromosome=="X")
    {
        cout<<"\n Checking Chromosome X information !!! \n";
        cout << " NOTE: All Samples in target VCF file must be either only MALE or only FEMALE ...  " << endl;

        inFile.open(filename, header);
        inFile.setSiteOnly(false);
        inFile.readRecord(record);
        int tempHapFlag=0;
        int tempHapCount=0;
        for (int i = 0; i<(header.getNumSamples()); i++)
        {
            if(record.getNumGTs(i)==0)
            {
                std::cout << "\n Empty Value for Individual at first variant !!! " << endl;
                std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
                return false;
            }
            else
            {
                if(tempHapFlag!=0 && tempHapFlag!=record.getNumGTs(i))
                {
                    cout << "\n ERROR: Both haploid and diploid samples found in target VCF file "<<filename<<endl;
                    cout << " For chromosome X imputation, Males and Females should be imputed separately.\n";
                    cout << " Program Aborting ... "<<endl;
                    return false;

                }
                tempHapFlag=record.getNumGTs(i);
                tempHapCount+=record.getNumGTs(i);
            }

        }
        inFile.close();
        std::cout << " "<<tempHapCount<<" haplotypes found in file "<<endl;
    }


    std::cout << "\n Initial basic file check on target/GWAS haplotype file successful !!!" << endl<<endl;
	return true;

}


bool HaplotypeSet::LoadTargetHaplotypes(String filename, String targetSnpfile, vector<string> &refSnpList,HaplotypeSet &rHap)
{
	string FileType=DetectTargetFileType(filename);

    if(FileType.compare("vcf")==0)
    {
        cout<<"\n Format = VCF (Variant Call Format) "<<endl;
        return LoadVcfTargetHaplotypes(filename, targetSnpfile, refSnpList,rHap);
    }
    else if(FileType.compare("Invalid")==0)
    {

        cout << "\n Target File provided by \"--haps\" must be a VCF or MaCH file !!! \n";
        cout << " Please check the following file : "<<filename<<endl<<endl;
        return false;
    }
    else if(FileType.compare("NA")==0)
    {
        cout<<"\n Following File File Not Available : "<<filename<<endl<<endl;
        return false;
    }
    return false;
}

bool HaplotypeSet::LoadVcfTargetHaplotypes(String filename, String snpNames, vector<string> &refSnpList, HaplotypeSet &rHap)
{


	bool rsid=false;
	VcfFileReader inFile;
	VcfHeader header;
	VcfRecord record;
	vector<string> tempMarkerNames;
	tempMarkerNames.clear();

	inFile.setSiteOnly(true);
	if (!inFile.open(filename, header))
	{
		cout << "\n Program could NOT open file : " << filename << endl;
		return false;
	}
    vector<int> importIndexList;

    int numReadRecords = 0,numActualRecords=0;
    int failFilter=0,notBiallelic=0,inconsistent=0;
    std::cout << "\n Loading Target Haplotype SNP List from VCF File     : " << filename << endl<<endl;
    string cno;
    while (inFile.readRecord(record))
    {

        int flag=0;
        cno=record.getChromStr();
        int bp=record.get1BasedPosition();
        string id=record.getIDStr();
        string currID;
        string refAllele = record.getRefStr();
        string altAllele = record.getAltStr();
        char Rallele=0,Aallele=0;


        if (record.getNumAlts()>1)
		{
			notBiallelic++;
			flag = 1;
		}
		if (record.getFilter().getString(0).compare("PASS") != 0)
		{
			failFilter++;
			flag = 0;
		}
        stringstream strs3,strs4;
        strs3<<(cno);
        strs4<<(bp);


        if(rsid)
            currID=record.getIDStr();
        else
            currID=(string)strs3.str()+":"+(string)strs4.str();

        if (strlen(refAllele.c_str()) == 1 && strlen(altAllele.c_str()) == 1)
        {
            switch (refAllele[0])
            {
                case 'A': case 'a': Rallele = 1; break;
                case 'C': case 'c': Rallele = 2; break;
                case 'G': case 'g': Rallele = 3; break;
                case 'T': case 't': Rallele = 4; break;
                case 'D': case 'd': Rallele = 5; break;
                case 'I': case 'i': Rallele = 6; break;
                case 'R': case 'r': Rallele = 7; break;
                default:
                {
                        cout << " WARNING !!! Reference allele for SNP for "<<currID<<"is "<<refAllele<<". Will be ignored ..." <<endl;
                        flag=1;
                        inconsistent++;
                }
            }
            if(flag==0)
                switch (altAllele[0])
                {
                    case '0':  Aallele = 0; break;
                    case 'A': case 'a': Aallele = 1; break;
                    case 'C': case 'c': Aallele = 2; break;
                    case 'G': case 'g': Aallele = 3; break;
                    case 'T': case 't': Aallele = 4; break;
                    case 'D': case 'd': Aallele = 5; break;
                    case 'I': case 'i': Aallele = 6; break;
                    case 'R': case 'r': Aallele = 7; break;
                    default:
                    {
                        cout << " WARNING !!! Alternate allele for SNP for "<<currID<<"is "<<altAllele<<". Will be ignored ..." <<endl;
                        flag=1;
                        inconsistent++;
                    }
                }
        }
        else
            Rallele = 7;
        variant thisVariant(currID,cno,bp);
        if(flag==0)
        {
            VariantList.push_back(thisVariant);

            VariantList[numReadRecords].refAllele=Rallele;
            VariantList[numReadRecords].altAllele=Aallele;
            VariantList[numReadRecords].refAlleleString=refAllele;
            VariantList[numReadRecords].altAlleleString=altAllele;
            markerName.push_back(currID);
            ++numReadRecords;
        }
        numActualRecords++;
        importIndexList.push_back(flag);

    }

	std::cout << "\n Number of Markers read from VCF File                : " << numActualRecords << endl;
	std::cout << " Number of Markers with more than One Allele         : " << notBiallelic << endl;
	std::cout << " Number of Markers failing FILTER = PASS             : " << failFilter << endl;
    std::cout << " Number of Markers with inconsistent Ref/Alt Allele  : " << inconsistent << endl;

    std::cout << "\n Number of Markers to be Recorded                    : " << numReadRecords << endl;
    int refMarkerCount=(int)rHap.VariantList.size();

    finChromosome=cno;
	vector<int> knownPosition;
	knownPosition.clear();
	int counter = 0;
	missing.resize(refMarkerCount, true);
	ScaffoldIndex.resize(refMarkerCount, -1);
    UnScaffoldIndex.resize(VariantList.size(), -1);
    AllMaleTarget=false;

    int flag;
	int markerIndex=0;
	vector<string> newMarkerName;
	newMarkerName.clear();
    vector<bool> RefAlleleSwap;

    RefAlleleSwap.clear();

	for (int j = 0; j<(int)VariantList.size(); j++)
	{

		int prevCounter = counter;
		flag=0;
		while(counter<refMarkerCount && flag==0 && rHap.VariantList[counter].bp<=VariantList[j].bp)
        {

            if(rHap.VariantList[counter].chr==VariantList[j].chr
             && rHap.VariantList[counter].bp==VariantList[j].bp)
            {
                prevCounter = counter;

                if(rHap.VariantList[counter].refAlleleString==VariantList[j].refAlleleString
                        && rHap.VariantList[counter].altAlleleString==VariantList[j].altAlleleString)
                    flag=1;
                else if(rHap.VariantList[counter].refAlleleString==VariantList[j].altAlleleString
                        && rHap.VariantList[counter].altAlleleString==VariantList[j].refAlleleString)
                    flag=1;
                else if (VariantList[j].refAlleleString==VariantList[j].altAlleleString
                        && rHap.VariantList[counter].refAlleleString==VariantList[j].refAlleleString)
                    flag=1;
                else if (VariantList[j].refAlleleString==VariantList[j].altAlleleString
                        && rHap.VariantList[counter].altAlleleString==VariantList[j].refAlleleString)
                    flag=1;
                else
                    counter++;
            }
            else
                counter++;
        }
        if(flag==1)
        {
            knownPosition.push_back(counter);
            newMarkerName.push_back(markerName[j]);

            if(rHap.VariantList[counter].refAlleleString==VariantList[j].refAlleleString)
                RefAlleleSwap.push_back(false);
            else
                {

                    RefAlleleSwap.push_back(true);
                    string tempa=VariantList[j].refAlleleString;
                    VariantList[j].refAlleleString=VariantList[j].altAlleleString;
                    VariantList[j].altAlleleString=tempa;
                }

            missing[counter] = false;
            UnScaffoldIndex[markerIndex]=counter;
            ScaffoldIndex[counter]=markerIndex++;

            counter++;
		}
		else
        {
            knownPosition.push_back(-1);
			counter = prevCounter;
        }

	}

	numHaplotypes = 0;
    numMarkers = 0;

	std::cout << " Number of Markers overlapping with Reference List   : " << newMarkerName.size() << endl << endl;

	if (newMarkerName.size() == 0)
	{

		cout << "\n No overlap between Target and Reference markers !!!\n";
		cout << " Please check for consistent marker identifer in reference and target input files..\n";
		cout << " Program Aborting ... \n";
		return false;

	}

	markerName = newMarkerName;
    numHaplotypes = 2 * header.getNumSamples();
    numSamples=header.getNumSamples();
    if(numHaplotypes==0)
    {
        cout << "\n No haplotypes recorded from VCF Input File : "<<filename<<endl;
		cout << " Please check the file properly..\n";
		cout << " Program Aborting ... "<<endl;
		return false;
    }


    for (int i = 0; i < numSamples; i++)
	{
		string tempName(header.getSampleName(i));
		individualName.push_back(tempName);
	}

    if(finChromosome=="X")
    {
        inFile.open(filename, header);
        inFile.setSiteOnly(false);
        inFile.readRecord(record);
        int tempHapCount=0;
        for (int i = 0; i<(numSamples); i++)
        {
            if(record.getNumGTs(i)==0)
            {
                std::cout << "\n Empty Value for Individual : " << individualName[i] << " at Marker : " << VariantList[0].name << endl;
                std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
                return false;
            }
            else
            {
                SampleNoHaplotypes.push_back(record.getNumGTs(i));
                tempHapCount+=record.getNumGTs(i);
            }
            if(i==0)
            {
                if(record.getNumGTs(i)==1)
                    AllMaleTarget=true;
                else
                    AllMaleTarget=false;
            }
        }
        inFile.close();
        numHaplotypes=tempHapCount;
    }



	inFile.close();

	std::cout << " Loading Target Haplotype Set from VCF File          : " << filename << endl<<endl;


	inFile.open(filename, header);
    inFile.setSiteOnly(false);

	inFile.readRecord(record);


	//haplotypes.resize(numHaplotypes);
	haplotypesUnscaffolded.resize(numHaplotypes);
	MissingSampleUnscaffolded.resize(numHaplotypes);




	for (int i = 0; i<numHaplotypes; i++)
	{
		haplotypesUnscaffolded[i].resize(markerName.size(), false);
        MissingSampleUnscaffolded[i].resize(markerName.size(), false);
    }

	// initialize this to 0 again to read each marker from each line of vcf file.
    int importIndex=-1;
    //ScaffoldIndex=knownPosition;
	int readIndex = -1;
	int numtoBeWrittenRecords = 0;
	do
	{
		// work only with bi-allelic markers
		readIndex++;
        if (importIndexList[readIndex] ==0)
		{
		    importIndex++;

            if (knownPosition[importIndex] != -1)
            {
                if (numtoBeWrittenRecords % 10000 == 0)
                       printf("  Loading markers %d out of %d markers to be loaded... [%.1f%%] \n",numtoBeWrittenRecords+1,(int)markerName.size(),100*(double)(numtoBeWrittenRecords + 1)/(int)markerName.size());

                tempMarkerNames.push_back(record.getIDStr());
                string refAllele = record.getRefStr();
                string altAllele = record.getAltStr();
                char Rallele;
                if (strlen(refAllele.c_str()) == 1 && strlen(altAllele.c_str()) == 1)
                {
                    switch (refAllele[0])
                    {
                        case 'A': case 'a': Rallele = 1; break;
                        case 'C': case 'c': Rallele = 2; break;
                        case 'G': case 'g': Rallele = 3; break;
                        case 'T': case 't': Rallele = 4; break;
                        case 'D': case 'd': Rallele = 5; break;
                        case 'I': case 'i': Rallele = 6; break;
                        case 'R': case 'r': Rallele = 7; break;
                        default:
                        {
                                   cout << "\n Data Inconsistency !!! \n";
                                   cout << " Error with reference allele for marker : " << record.getIDStr() << " in VCF File : " << filename;
                                   cout << "\n VCF reference alleles for SNPs can only be A(a), C(c), G(g), or T(t).\n";
                                   cout << " " << record.getIDStr() << " has " << refAllele << endl;
                                   cout << " Program Aborting ... \n\n";
                                   return false;
                        }
                    }
                }
                else
                    Rallele = 7;

                refAlleleList.push_back(Rallele);

                int haplotype_index = 0;
                for (int i = 0; i<(numSamples); i++)
                {

                    if(record.getNumGTs(i)==0)
                    {
                        std::cout << "\n Empty Value for Individual : " << individualName[i] << " at Marker : " <<
                         VariantList[numtoBeWrittenRecords].name << endl;
                        std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
                        return false;
                    }
                    if(record.getNumGTs(i)==1 && finChromosome!="X")
                    {
                        std::cout << "\n Single Autosomal Haplotype for Individual : " << individualName[i] << " at Marker : "
                         << VariantList[numtoBeWrittenRecords].name << endl;
                        std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
                        return false;
                    }

                    if(rHap.PseudoAutosomal)
                    {
                        if(record.getNumGTs(i)!=2)
                        {
                            cout << "\n ERROR: Reference Panel on chromosome X seems to be on Pseudo-Autosomal Region "<<endl;
                            cout <<   "        since all samples have two alleles on Chromosome X. "<<endl;
                            cout <<   "        However, "<<VariantList[numtoBeWrittenRecords].name<<" in target panel has SINGLE haplotypes !!!\n";
                            cout << " Program Aborting ... "<<endl;
                            return false;
                        }
                    }
                    else
                    {
                        if(AllMaleTarget && record.getNumGTs(i)!=1)
                        {
                            cout << "\n ERROR: Reference Panel on chromosome X seems to be on NON-Pseudo-Autosomal Region "<<endl;
                            cout <<   "        since MALE samples had single alleles on Chromosome X. "<<endl;
                            cout <<   "        However, "<<VariantList[numtoBeWrittenRecords].name<<" in target panel (all MALE samples) \n";
                            cout <<   "        has DOUBLE haplotypes !!!\n";
                            cout << " Program Aborting ... "<<endl;
                            return false;
                        }
                        if(!AllMaleTarget && record.getNumGTs(i)!=2)
                        {
                            cout << "\n ERROR: Reference Panel on chromosome X seems to be on NON-Pseudo-Autosomal Region "<<endl;
                            cout <<   "        since MALE samples had single alleles on Chromosome X. "<<endl;
                            cout <<   "        However, "<<VariantList[numtoBeWrittenRecords].name<<" in target panel (all FEMALE samples) \n";
                            cout <<   "        has SINGLE haplotypes !!!\n";
                            cout << " Program Aborting ... "<<endl;
                            return false;
                        }



                    }

                    for (int j = 0; j<record.getNumGTs(i); j++)
                    {

                        int alleleIndex = record.getGT(i, j);
                        if (alleleIndex<0)
                        {
                            MissingSampleUnscaffolded[haplotype_index][numtoBeWrittenRecords] = true;
                        }
                        else
                        {

                            if(!RefAlleleSwap[numtoBeWrittenRecords])
                            {
                                if(alleleIndex==1)
                                {
                                    haplotypesUnscaffolded[haplotype_index][numtoBeWrittenRecords] = true;
                                }
                            }
                            else
                            {
                                if(alleleIndex==0)
                                {
                                    haplotypesUnscaffolded[haplotype_index][numtoBeWrittenRecords] = true;
                                }
                            }


                        }

                        haplotype_index++;



                    }
                }
            ++numtoBeWrittenRecords;
            }
		}
	} while (inFile.readRecord(record));

	std::cout << "\n\n Number of Markers Recorded                          : " << markerName.size() << endl;
	std::cout << " Number of Haplotypes Recorded                       : " << (haplotypesUnscaffolded.size()) << endl;
    numMarkers = markerName.size();

    std::cout << "\n Haplotype Set successfully loaded from VCF File     : " << filename << endl;
	inFile.close();



	return true;

}


bool HaplotypeSet::getScaffoldedHaplotype(int sample,int marker)
{

   if(missing[marker]==true)
        {
            return false;
        }
   else
        return haplotypesUnscaffolded[sample][ScaffoldIndex[marker]];
}

void HaplotypeSet::Create(vector<bool> &tempHaplotype)
{

    numHaplotypes=1;
    numMarkers=(int)tempHaplotype.size();
    missing.resize(numMarkers,false);
    ScaffoldIndex.resize(numMarkers);
    vector<bool> tempMissin(numMarkers,false);


    haplotypesUnscaffolded.clear();
    MissingSampleUnscaffolded.clear();

    haplotypesUnscaffolded.push_back(tempHaplotype);
    MissingSampleUnscaffolded.push_back(tempMissin);


    for(int i=0;i<numMarkers;i++)
        ScaffoldIndex[i]=i;


}


bool HaplotypeSet::getMissingScaffoldedHaplotype(int sample,int marker)
{

   if(missing[marker]==true)
        {
            return true;
        }
   else
        return MissingSampleUnscaffolded[sample][ScaffoldIndex[marker]];
}
//
//
//bool HaplotypeSet::FastLoadHaplotypes(String filename, int maxIndiv, int maxMarker,String CHR,
//                                      int START,int END,int WINDOW,bool rsid,bool compressOnly,bool filter)
//{
//
//
//    string FileType=DetectReferenceFileType(filename);
//
//    if(FileType.compare("m3vcf")==0)
//    {
//        cout<<"\n Format = M3VCF (Minimac3 VCF File) "<<endl;
////        if(compressOnly)
////        {
////            cout << "\n Reference File provided by \"--refHaps\" is an M3VCF file !!! \n";
////            cout << " M3VCF files cannot be processed further !!! "<<endl;
////            return false;
////        }
//
//
//        return readm3vcfFile(filename,CHR,START,END,WINDOW);
//    }
//    else if(FileType.compare("NA")==0)
//    {
//        cout<<"\n Following File File Not Available : "<<filename<<endl;
//        return false;
//    }
//    else if(FileType.compare("Invalid")==0)
//    {
//
//        cout << "\n Reference File provided by \"--refHaps\" must be a VCF or M3VCF file !!! \n";
//        cout << " Please check the following file : "<<filename<<endl;
//        return false;
//    }
//    cout<<"\n Format = VCF (Variant Call Format) "<<endl;
//    vcfType=true;
//
//    optEndPoints.clear();
//	VcfFileReader inFile;
//	VcfHeader header;
//	VcfRecord record;
//
//	if (!inFile.open(filename, header))
//	{
//		cout << "\n Program could NOT open file : " << filename << endl;
//		return false;
//	}
//
//    std::cout << "\n Loading Reference Haplotype Set from VCF File       : " << filename << endl;
////    if(END==0)
////        END=300000000;
//    int OrigStartPos=START;
//    int OrigEndPos=END;
//    if (WINDOW > 0)
//    {
//        if (START-WINDOW < 0)
//            START = 0;
//        else
//            START -= WINDOW;
//
//        END += WINDOW;
//    }
//    stringstream strs;
//    strs<<(END);
//
//    if(CHR!="")
//    {
//        std::cout << "\n Region specified by user (including window = "<<WINDOW <<" bp) : chr" << CHR
//        <<":"<<START <<"-"<< (END > 0 ? (string)(strs.str()) :"END") << endl;
//    }
//
//
//
//    PrintStartIndex=0;
//    PrintEndIndex=0;
//    int     numReadRecords = 0;
//	int numtoBeWrittenRecords = 0;
//	vector<bool> importIndex;
//    int haplotype_index = 0;
//	int numSamplesRead = 0;
//	int notBiallelic = 0;
//	int failFilter = 0;
//	int duplicates = 0;
//	int insertions = 0;
//	int deletions = 0;
//	int inconsistent=0;
//	string prevID="";
////	bool isInDel = false;
//	inFile.setSiteOnly(true);
//	numSamplesRead = header.getNumSamples();
//
//	for (int i = 0; i < numSamplesRead; i++)
//	{
//		string tempName(header.getSampleName(i));
//		individualName.push_back(tempName);
//	}
//
//    string refAllele,altAllele,PrefrefAllele,PrevaltAllele,cno,fixCno,id;
//    int bp;
//	cout << "\n Reading VCF File to calculate number of records ... \n";
//    cout<<endl;
//	// pre-calculate number of samples read and number of markers read to allocate memory.
//	while (inFile.readRecord(record))
//	{
//
//        int flag = 0;
//        if (maxMarker != 0 && numtoBeWrittenRecords>=maxMarker)
//            break;
//
//		++numReadRecords;
//
//		if (record.getNumAlts()>1)
//		{
//			notBiallelic++;
//			flag = 1;
//		}
//		if (record.getFilter().getString(0).compare("PASS") != 0)
//		{
//			failFilter++;
//			if(filter)
//                flag = 1;
//		}
////		isInDel=false;
//
//		refAllele = record.getRefStr();
//		cno=record.getChromStr();
//        bp=record.get1BasedPosition();
//        id=record.getIDStr();
//        altAllele = record.getAltStr();
//
//
//        if(numReadRecords==1 && CHR=="")
//        {
//            if(!CheckValidChrom(cno))
//            {
//                cout << "\n Error !!! Reference VCF File contains chromosome : "<<cno<<endl;
//                cout << " VCF File can only contain chromosomes 1-22 and X !!! "<<endl;
//                cout << " Program Aborting ... "<<endl;
//                return false;
//            }
//
//        }
//        if(CHR=="" && fixCno!=cno && numReadRecords>1)
//        {
//            cout << "\n Error !!! Reference VCF File contains multiple chromosomes : "<<cno<<", "<<fixCno<<", ... "<<endl;
//            cout << " Please use VCF file with single chromosome or specify chromosome using \"--chr\" option !!! "<<endl;
//            cout << " Program Aborting ... "<<endl;
//            return false;
//        }
//
//
//        fixCno=cno;
//
//        string currID;
//
//        stringstream strs1,strs2;
//        strs1<<(cno);
//        strs2<<(bp);
//
//
//        if(rsid)
//            currID=record.getIDStr();
//        else
//            currID=(string)strs1.str()+":"+(string)strs2.str();
//
//        //cout<<prevID<<endl;
//        if (strlen(refAllele.c_str()) == 1 && strlen(altAllele.c_str()) == 1)
//        {
//            switch (refAllele[0])
//            {
//                case 'A': case 'a': break;
//                case 'C': case 'c': break;
//                case 'G': case 'g': break;
//                case 'T': case 't': break;
//                case 'D': case 'd': break;
//                case 'I': case 'i': break;
//                case 'R': case 'r': break;
//                default:
//                {
//                   cout << "\n WARNING !!! ";
//                   cout << " Reference allele for marker : " << currID << " is : " <<refAllele<<endl;
////                   cout<<" In VCF File : " << filename;
////                   cout << "\n VCF reference alleles for SNPs (not INDELs/SVs) can only be A(a), C(c), G(g), or T(t).";
//                   cout << " Variant will be ignored... \n";
//                   flag=1;
//                   inconsistent++;
//                }
//            }
//
//            if(flag==0)
//                switch (altAllele[0])
//                    {
//                    case 'A': case 'a': break;
//                    case 'C': case 'c': break;
//                    case 'G': case 'g': break;
//                    case 'T': case 't': break;
//                    case 'D': case 'd': break;
//                    case 'I': case 'i': break;
//                    case 'R': case 'r': break;
//                    default:
//                    {
//                        cout << "\n WARNING !!! ";
//                        cout << " Alternate allele for marker : " <<currID << " is : " <<altAllele<<endl;
////                        cout<<" In VCF File : " << filename;
////                        cout << "\n VCF alternate alleles for SNPs (not INDELs/SVs) can only be A(a), C(c), G(g), or T(t).";
//                        cout << " Variant will be ignored... \n";
//                        flag=1;
//                    }
//                }
//        }
//        else if(strlen(refAllele.c_str())<strlen(altAllele.c_str()))
//        {
////			isInDel = true;
//			insertions++;
//		}
//		else
//        {
////			isInDel = true;
//			deletions++;
//        }
//
//        stringstream strs3,strs4;
//        strs3<<(cno);
//        strs4<<(bp);
//
//        if(prevID==currID)
//        {
//            if(refAllele==PrefrefAllele && altAllele==PrevaltAllele)
//            {
//
//                cout << " WARNING !!! Duplicate Variant found chr:"
//                <<(string)strs3.str()+":"+
//                (string)strs4.str()<<" with identical REF = "
//                <<refAllele
//                 <<" and ALT = "<<altAllele <<"\n";
//                duplicates++;
//            }
//
//        }
//
//
////        prevSNP=1-isInDel;
//        prevID=currID;
//        PrefrefAllele=refAllele;
//        PrevaltAllele=altAllele;
//
//		if(CHR!="")
//        {
//            if(cno.compare(CHR.c_str())!=0)
//                flag=1;
//            else
//            {
//                if(END>0)
//                {
//                    if(bp>END || bp<START)
//                        flag=1;
//                }
//                else
//                    if(bp<START)
//                        flag=1;
//            }
//
//        }
//
////cout<<fixCno<<"\t"<<cno<<endl;
//
//
//		if (flag == 0)
//		{
//			if(bp<OrigStartPos)
//                PrintStartIndex++;
//
//            if(CHR=="" || bp<=OrigEndPos)
//                PrintEndIndex=numtoBeWrittenRecords;
//			++numtoBeWrittenRecords;
//			markerName.push_back(currID);
//            variant thisVariant(currID,cno,bp);
//            VariantList.push_back(thisVariant);
//			importIndex.push_back(true);
//
//		}
//		else
//		{
//			importIndex.push_back(false);
//		}
//
//	}
//
//
//	inFile.close();
//
//
//    if(CHR=="")
//        finChromosome=cno;
//    else
//        finChromosome=CHR.c_str();
//
//
//
//
//	std::cout << "\n Number of Markers read from VCF File                : " << numReadRecords << endl;
//	std::cout << " Number of Markers with more than Two Alleles        : " << notBiallelic << endl;
//	std::cout << " Number of Markers failing FILTER = PASS             : " << failFilter << endl;
//    std::cout << " Number of Markers with inconsistent Ref/Alt Allele  : " << inconsistent << endl;
//    std::cout << " Number of Markers with duplicate ID/Position        : " << duplicates << endl;
//    std::cout << " Number of Insertions                                : " << insertions << endl;
//	std::cout << " Number of Deletions                                 : " << deletions << endl;
//
//	if (maxIndiv == 0)
//		maxIndiv = numSamplesRead;
//	numMarkers = numtoBeWrittenRecords;
//	numHaplotypes = (maxIndiv<numSamplesRead) ? (2 * maxIndiv) : (2 * numSamplesRead);
//    numSamples=numHaplotypes/2;
//    if(finChromosome!="X")
//    {
//        std::cout << "\n Number of Markers to be Recorded                    : " << numtoBeWrittenRecords << endl;
//        std::cout << " Number of Haplotypes to be Recorded                 : " << (numHaplotypes) << endl;
//    }
//    else
//    {
//        cout<<"\n Chromosome X Detected !!! \n";
//        std::cout << "\n Number of Markers to be Recorded                    : " << numtoBeWrittenRecords << endl;
//        PseudoAutosomal=false;
//        inFile.open(filename, header);
//        inFile.setSiteOnly(false);
//        int numWrittenRecords = 0;
////        int readIndex = -1;
//        inFile.readRecord(record);
//        int tempHapCount=0,MaleCount=0,FemaleCount=0;
//        for (int i = 0; i<(numSamples); i++)
//        {
//            if(record.getNumGTs(i)==0)
//            {
//                std::cout << "\n Empty Value for Individual : " << individualName[i] << " at Marker : " << VariantList[numWrittenRecords].name << endl;
//                std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
//                return false;
//            }
//            else
//            {
//                if(record.getNumGTs(i)==1)
//                    MaleCount++;
//                else
//                    FemaleCount++;
//                SampleNoHaplotypes.push_back(record.getNumGTs(i));
//
//                tempHapCount+=record.getNumGTs(i);
//            }
//
//        }
//        inFile.close();
//        std::cout << " Number of Samples (Haplotypes) to be Recorded       : " << numSamples << " ("<<tempHapCount <<") "<<endl;
//        numHaplotypes=tempHapCount;
//
//        if(MaleCount>0)
//        {
//        std::cout << " Number of MALE Samples (Haplotypes)                 : " << MaleCount << " ("<<MaleCount <<") "<<endl;
//        std::cout << " Number of FEMALE Samples (Haplotypes)               : " << FemaleCount<< " ("<<FemaleCount*2 <<") "<<endl;
//        }
//        else
//        {
//            std::cout << "\n All " << numSamples<<" samples have two alleles on Chromosome X !!! "<< endl;
//            cout<<" Cannot determine number of male/female samples from Pseudo-Autosomal Region ..."<<endl;
//            PseudoAutosomal=true;
//        }
//
//
//    }
//
//
//   if(numtoBeWrittenRecords<2)
//    {
//        cout << "\n None/Single marker left after filtering from Input File : "<<filename<<endl;
//		cout << " Please check the file or the filtering options properly ...\n";
//		cout << " Program Aborting ... "<<endl;
//		return false;
//    }
//
//    std::cout <<"\n Starting to load data ...";
//
//
//
//    int    blockSize = 500;
//    int    bufferSize = 1000;
//    int LastflushPos=0;
//    vector<int> index(numHaplotypes),oldIndex;
//    vector<int> previousDifference(numHaplotypes);
//    vector<int> previousPredecessor(numHaplotypes);
//    vector<int> firstDifference(numHaplotypes-1,0);
//    vector<int> cost(bufferSize+1,0);
//    vector<int> bestSlice(bufferSize+1,0);
//    vector<int> bestComplexity(bufferSize+1,0);
//    vector<vector<int> > bestIndex(bufferSize+1);
//
//    findUnique RefUnique;
//
//    RefUnique.updateCoeffs(transFactor,cisFactor);
//    vector<String> Haplotypes(numHaplotypes);
//    double blockedCost = 0.0;
//    for(int i=0;i<numHaplotypes;i++)
//        index[i]=i;
//    individualName.resize(numSamples);
//    // open file again to start loading the data in the variable haplotype.
//
//    cout<<endl;
//	inFile.open(filename, header);
//	inFile.setSiteOnly(false);
//	int numWrittenRecords = 0;
//	int readIndex = -1;
//
//
//
//
//
//
//
//	while (inFile.readRecord(record) && numWrittenRecords<numMarkers)
//	{
//		// work only with bi-allelic markers
//		readIndex++;
//
//		if (importIndex[readIndex])
//		{
//            if (Haplotypes[0].Length()==1)
//                {
////                    cout<<numWrittenRecords<<"\t"<<readIndex<<"\t"<<Haplotypes[0].Length()<<endl;
//                    printf("  Loading markers %d - %d  out of %d markers to be loaded... [%.1f%%] "
//                    , LastflushPos + 1,
//                    min(LastflushPos+bufferSize,numMarkers),numMarkers,100*(double)(LastflushPos + 1)/numMarkers);
//                    cout<<endl;
//                }
//
//
//
//            string refAllele = record.getRefStr();
//			string altAllele = record.getAltStr();
//			char Rallele,Aallele;
//
//			if (strlen(refAllele.c_str()) == 1 && strlen(altAllele.c_str()) == 1)
//				{
//				    switch (refAllele[0])
//                    {
//                        case 'A': case 'a': Rallele = 1; break;
//                        case 'C': case 'c': Rallele = 2; break;
//                        case 'G': case 'g': Rallele = 3; break;
//                        case 'T': case 't': Rallele = 4; break;
//                        case 'D': case 'd': Rallele = 5; break;
//                        case 'I': case 'i': Rallele = 6; break;
//                        case 'R': case 'r': Rallele = 7; break;
//                        default:
//                        {
//                                   cout << "\n\n Data Inconsistency !!! \n";
//                                   cout << " Error with reference allele for marker : " << record.getIDStr() << " in VCF File : " << filename;
//                                   cout << "\n VCF reference alleles for SNPs can only be A(a), C(c), G(g), or T(t).\n";
//                                   cout << " " << record.getIDStr() << " has " << refAllele << endl;
//                                   cout << "\n Program Aborting ... \n\n";
//                                   return false;
//                        }
//                    }
//
//                    switch (altAllele[0])
//                    {
//                        case 'A': case 'a': Aallele = 1; break;
//                        case 'C': case 'c': Aallele = 2; break;
//                        case 'G': case 'g': Aallele = 3; break;
//                        case 'T': case 't': Aallele = 4; break;
//                        case 'D': case 'd': Aallele = 5; break;
//                        case 'I': case 'i': Aallele = 6; break;
//                        case 'R': case 'r': Aallele = 7; break;
//                        default:
//                        {
//                                   cout << "\n\n Data Inconsistency !!! \n";
//                                   cout << " Error with alternate allele for marker : " << record.getIDStr() << " in VCF File : " << filename;
//                                   cout << "\n VCF alternate alleles for SNPs can only be A(a), C(c), G(g), or T(t).\n";
//                                   cout << " " << record.getIDStr() << " has " << altAllele << endl;
//                                   cout << "\n Program Aborting ... \n\n";
//                                   return false;
//                        }
//                    }
//
//
//				}
//
//			else
//				{
//				    Rallele = 7;
//				    if(strlen(refAllele.c_str())<strlen(altAllele.c_str()))
//                        Aallele=6;
//                    else
//                        Aallele=5;
//				}
//
//            refAlleleList.push_back(Rallele);
//            VariantList[numWrittenRecords].refAllele=Rallele;
//            VariantList[numWrittenRecords].altAllele=Aallele;
//            VariantList[numWrittenRecords].refAlleleString=refAllele;
//            VariantList[numWrittenRecords].altAlleleString=altAllele;
//            haplotype_index = 0;
//            for (int i = 0; i<(numSamples); i++)
//            {
//                if(record.getNumGTs(i)==0)
//                {
//                    std::cout << "\n Empty Value for Individual : " << individualName[i] << " at Marker : " << VariantList[numWrittenRecords].name << endl;
//                    std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
//                    return false;
//                }
//                if(record.getNumGTs(i)==1 && finChromosome!="X")
//                {
//                    std::cout << "\n Single Autosomal Haplotype for Individual : " << individualName[i] << " at Marker : " << VariantList[numWrittenRecords].name << endl;
//                    std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
//                    return false;
//                }
//                for (int j = 0; j<record.getNumGTs(i); j++)
//                {
//
//                    int alleleIndex = record.getGT(i, j);
//                    if (alleleIndex<0)
//                    {
//                        std::cout << "\n Missing Value for Individual : " << individualName[i] << " at Marker : " << VariantList[numWrittenRecords].name << endl;
//                        return false;
//                    }
//
//                    else
//                    {
//                        if(haplotype_index>=numHaplotypes && finChromosome=="X")
//                        {
//                            cout << "\n Error in Reference VCF File for Chromosome X !!! "<<endl;
//                            cout << "\n Marker : " << VariantList[0].name << " has "<<numHaplotypes <<" haplotypes while ";
//                            cout << "Marker : " << VariantList[numWrittenRecords].name << " has "<< haplotype_index+1<<" haplotypes."<<endl;
//                            cout << " VCF file seems to have both Pseudo-Autosomal region (PAR) and non-PAR of chromosome X. \n";
//                            cout << " Please use only either of the two regions ... \n";
//                            cout << " See web-page for Minimac3 (Chromosome X Imputation) for more details ...\n";
//                            cout << " Program Aborting ... "<<endl;
//                            return false;
//                        }
//                        Haplotypes[haplotype_index]+= (char)('0'+record.getGT(i, j));
//                        haplotype_index++;
//                    }
//
//                }
//
//
//            }
//
//            if(haplotype_index<(numHaplotypes-1) && finChromosome=="X")
//            {
//                cout << "\n Error in Reference VCF File for Chromosome X !!! "<<endl;
//                cout << "\n Marker : " << VariantList[0].name << " has "<<numHaplotypes <<" haplotypes while ";
//                cout << "Marker : " << VariantList[numWrittenRecords].name << " has "<< haplotype_index+1<<" haplotypes."<<endl;
//                cout << " VCF file seems to have both Pseudo-Autosomal region (PAR) and non-PAR of chromosome X. \n";
//                cout << " Please use only either of the two regions ... \n";
//                cout << " See web-page for Minimac3 (Chromosome X Imputation) for more details ...\n";
//                cout << " Program Aborting ... "<<endl;
//                return false;
//            }
//            numWrittenRecords++;
//            int length = Haplotypes[0].Length();
//            vector<int> offsets(3,0);
//			for (int i = 0; i < numHaplotypes; i++)
//                offsets[Haplotypes[i][length - 1] - '0' + 1]++;
//            offsets[2]+=offsets[1];
////            cout<<VariantList[numWrittenRecords].name<<endl;
//
////            int bal=*max_element(index.begin(),index.end());
////            cout<<bal<<endl;
//
//
//            oldIndex = index;
//            for (int i = 0; i < numHaplotypes; i++)
//                {
//                    index[offsets[Haplotypes[oldIndex[i]][length - 1] - '0']++] = oldIndex[i];
////                    cout<<" INDEX [ "<<i<<
//                }
////            for (int i = 0; i < numHaplotypes; i++)
////                {
//////                    index[offsets[Haplotypes[oldIndex[i]][length - 1] - '0']++] = oldIndex[i];
////                    cout<<" INDEX [ "<<i<<" = "<<index[i]<<endl;
////                }
////
////                bal=*max_element(oldIndex.begin(),oldIndex.end());
//
//            RefUnique.UpdateDeltaMatrix(Haplotypes, index, firstDifference, length, blockSize,
//                           oldIndex, previousPredecessor, previousDifference);
//            RefUnique.AnalyzeBlocks(index, firstDifference, length, blockSize,
//                       cost, bestSlice, bestComplexity, bestIndex);
////             bal=*max_element(oldIndex.begin(),oldIndex.end());
////            cout<<bal<<endl;
//
//            if (Haplotypes[0].Length() == bufferSize)
//            {
//                blockedCost += RefUnique.FlushBlocks(optEndPoints,ReducedStructureInfo,VariantList,LastflushPos, Haplotypes, cost,
//                                            bestComplexity, bestSlice, bestIndex);
//                LastflushPos+=bufferSize;
//                LastflushPos--;
//                vector<String> tempHaplotypes(numHaplotypes);
//                for (int i = 0; i < numHaplotypes; i++)
//                {
//                    tempHaplotypes[i]=Haplotypes[i][Haplotypes[i].Length()-1];
//                    Haplotypes[i]=tempHaplotypes[i];
//                }
//
//
//                vector<int> tempoffsets(3,0);
//                for (int i = 0; i < numHaplotypes; i++)
//                    tempoffsets[Haplotypes[i][0] - '0' + 1]++;
//                tempoffsets[2]+=tempoffsets[1];
//
//                for (int i = 0; i < numHaplotypes; i++)
//                    index[tempoffsets[Haplotypes[i][0] - '0']++] = i;
//
//
//
//
//            }
//        }
//	}
//
//    if (Haplotypes[0].Length() > 1)
//      {
//      blockedCost +=  RefUnique.FlushBlocks(optEndPoints, ReducedStructureInfo,VariantList, LastflushPos, Haplotypes, cost,
//                                  bestComplexity, bestSlice, bestIndex);
//
////      double originalCost = 1.0 * (double) markerName.size() * (double)Haplotypes.size();
//
////      printf("\n After %d markers ...\n", (int)markerName.size());
////      printf("   Original cost = %.0f\n", originalCost);
////      printf("        New cost = %.0f\n", blockedCost);
////      printf("        Speed-up = %.1f\n", originalCost / (blockedCost + 1e-6));
//      }
//
//    std::cout << "\n Number of Markers Recorded                          : " << markerName.size() << endl;
//    std::cout << " Number of Haplotypes Recorded                       : " << numHaplotypes << endl;
//
//    optEndPoints.clear();
//    int i;
//    for(i=0;i<(int)ReducedStructureInfo.size();i++)
//        {
////            cout<<i<<"\t"<<ReducedStructureInfo[i].startIndex<<"\t"<<ReducedStructureInfo[i].endIndex<<endl;
//            optEndPoints.push_back(ReducedStructureInfo[i].startIndex);
//        }
//    optEndPoints.push_back(ReducedStructureInfo[i-1].endIndex);
//
//    //cout<<" NOW WHAT = "<<optEndPoints.back()<<endl;
//
//    if(individualName.size()==0)
//    {
//        cout << "\n No haplotypes recorded from VCF Input File : "<<filename<<endl;
//		cout << " Please check the file properly..\n";
//		cout << " Program Aborting ... "<<endl;
//		return false;
//    }
//
//    numMarkers = markerName.size();
//    std::cout << "\n Haplotype Set successfully loaded from VCF File     : " << filename << endl;
//	inFile.close();
//
//	return true;
//
//}
//
