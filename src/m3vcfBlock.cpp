#include "m3vcfBlock.h"

m3vcfBlock::m3vcfBlock()
{
    reset();
}


m3vcfBlock::~m3vcfBlock()
{
}


bool m3vcfBlock::read(IFILE filePtr, m3vcfHeader &ThisHeader,
                      bool InfoOnly,
                      bool siteOnly,
                     VcfRecordDiscardRules* discardRules,
                     VcfSubsetSamples* sampleSubset)
{
    // Clear out any previously set values.
    reset();
    const char *equalSep="=",*semicolSep=";";
    std::string tempString;

    if(filePtr == NULL)
    {
        myStatus.setStatus(StatGenStatus::FAIL_ORDER,
                           "Error reading VCF record before opening the file.");
        return(false);
    }

    if(ifeof(filePtr))
    {
        // End of file, just return false.
        return(false);
    }

    // Read the chromosome.
    if(!readTilTab(filePtr, myChrom))
    {
        if(myChrom.empty())
        {
            // EOF.
            return(false);
        }
        // Not an empty line.
        myStatus.setStatus(StatGenStatus::FAIL_PARSE,
                           "Error reading M3VCF Record CHROM.");
        return(false);
    }


    // Read the Block Boundaries
    if(!readTilTab(filePtr, tempString))
    {
        myStatus.setStatus(StatGenStatus::FAIL_PARSE,
                           "Error reading M3VCF Record BLOCK BOUNDARIES.");
        return(false);
    }
    else
    {
        // Read the positions, so convert them to an integer.
        char *end_str;
//        cout<<" STRING = "<<tempString<<endl;
        char  *pch  = strtok_r ((char*)tempString.c_str(),"-", &end_str);
        startPosition = atoi(pch);
        pch  = strtok_r (NULL,"-", &end_str);
//        cout<<" WAIT "<<pch<<endl;
        if(!pch)
        {
            myStatus.setStatus(StatGenStatus::FAIL_PARSE,
                                "Error parsing M3VCF Block BOUNDARIES.");
            return(false);
        }
        endPosition = atoi(pch);
    }

    readTilTab(filePtr, tempString);
    readTilTab(filePtr, tempString);
    readTilTab(filePtr, tempString);
    readTilTab(filePtr, tempString);
    readTilTab(filePtr, tempString);


    // Read the Number of Markers, Number of UniqueReps from Block ID.
    tempString.clear();
    if(!readTilTab(filePtr, tempString))
    {
        myStatus.setStatus(StatGenStatus::FAIL_PARSE,
                           "Error reading M3VCF Block INFO.");
        return(false);
    }
    else
    {
        // Read Number of Markers
        string InfoString = FindTokenWithPrefix(tempString.c_str(),semicolSep, "VARIANTS=");
        if(InfoString!="") { numMarkers = atoi(MyTokenize(InfoString.c_str(), equalSep, 2).c_str());}
        else
        {
            myStatus.setStatus(StatGenStatus::FAIL_PARSE,
                               "Error reading M3VCF Block VARIANTS.");
            return(false);
        }

        // Read Number of Unique Representatives
        InfoString = FindTokenWithPrefix(tempString.c_str(),semicolSep, "REPS=");
        if(InfoString!="") { numUniqueReps = atoi(MyTokenize(InfoString.c_str(), equalSep, 2).c_str());}
        else
        {
            myStatus.setStatus(StatGenStatus::FAIL_PARSE,
                               "Error reading M3VCF Block REPS.");
            return(false);
        }
    }


    readTilTab(filePtr, tempString);
    if(InfoOnly)
    {
        // Do not store genotypes, so just consume the rest of the line.
        filePtr->readTilChar("\n");
        return true;
    }




    int index = 0;
    UniqueIndexMap.resize(ThisHeader.getNumHaplotypes());
    while(index<ThisHeader.getNumHaplotypes())
    {
//        cout<<tempString<<endl;
        tempString.clear();
        if(!readTilTab(filePtr, tempString) && index<ThisHeader.getNumHaplotypes()-1)
        {
            myStatus.setStatus(StatGenStatus::FAIL_PARSE,
                               "Error reading M3VCF UNIQUE INDEX MAP.");
            return(false);
        }
        else
        {
            UniqueIndexMap[index++]=atoi(tempString.c_str());
        }
    }
    return(true);
}



void m3vcfBlock::reset()
{
    NoMarkersRead=0;
    FinishedReadingBlock=false;
    startPosition = -2;
    endPosition = -1;
    numMarkers = 0;
    numUniqueReps = 0;
    myChrom.clear();
    myStatus = StatGenStatus::SUCCESS;

}


// Return the error after a failed call.
const StatGenStatus& m3vcfBlock::getStatus()
{
    return(myStatus);
}

bool m3vcfBlock::write(IFILE filePtr, bool siteOnly)
{
    if(filePtr == NULL)
    {
        myStatus.setStatus(StatGenStatus::FAIL_ORDER,
                           "Error writing M3VCF record before opening the file.");
        return(false);
    }

    if(myChrom.length() == 0)
    {
        ifprintf(filePtr, ".\t");
    }
    else
    {
        ifprintf(filePtr, "%s\t", myChrom.c_str());
    }

    if(startPosition == -2)
    {
        ifprintf(filePtr, ".-.\t");
    }
    else
    {
        ifprintf(filePtr, "%d-%d\t", startPosition, endPosition);
    }

    ifprintf(filePtr, "<BLOCK>\t.\t.\t.\t.\tVARIANTS=%d,REPS=%d\t.", numMarkers, numUniqueReps);
    for(int i=0;i<(int)UniqueIndexMap.size();i++)
    {
        ifprintf(filePtr, "\t%d",UniqueIndexMap[i]);
    }
    ifprintf(filePtr, "\n");

   return true;
}



//
//
//const char* m3vcfBlock::getAlleles(unsigned int index)
//{
//    if(index == 0)
//    {
//        return(myRef.c_str());
//    }
//    if(index > getNumAlts())
//    {
//        // Index out of range.
//        // Throw an exception.
//        throw(std::runtime_error("m3vcfBlock::getAlleles called with an index that is greater than the number of alternates."));
//        return(NULL);
//    }
//    // Alternate allele, so return the alternate.
//    return(myAltArray.get(index-1).c_str());
//}

//
//int m3vcfBlock::getIntAllele(unsigned int index)
//{
//    const char* alleles = getAlleles(index);
//    switch(alleles[0])
//    {
//        case 'A':
//            return(1);
//            break;
//        case 'C':
//            return(2);
//            break;
//        case 'G':
//            return(3);
//            break;
//        case 'T':
//            return(4);
//            break;
//        default:
//            std::cerr << "m3vcfBlock::getIntAllele, unknown allele, "
//                      << alleles[0] << std::endl;
//    }
//    return(0);
//}
//
//
//unsigned int m3vcfBlock::getNumAlts()
//{
//    int numAlts = myAltArray.size();
//    if(numAlts != 0)
//    {
//        // Already parsed so just return the number of alternates.
//        return(numAlts);
//    }
//    // Check if it is just '.'.
//    if((myAlt.length() == 1) && (myAlt == "."))
//    {
//        // No alternates.
//        return(0);
//    }
//
//    // Parse the alternates by looping looking for commas.
//    std::string* altStr = &(myAltArray.getNextEmpty());
//    for(std::string::iterator iter = myAlt.begin(); iter != myAlt.end(); iter++)
//    {
//        if(*iter == ',')
//        {
//            altStr = &(myAltArray.getNextEmpty());
//        }
//        else
//        {
//            altStr->push_back(*iter);
//        }
//    }
//    return(myAltArray.size());
//}
//
//
//int m3vcfBlock::getAlleleCount(unsigned int index,
//                              VcfSubsetSamples* sampleSubset)
//{
//    unsigned int numAlts = getNumAlts();
//    if(index > numAlts)
//    {
//        // Index out of range.
//        // Throw an exception.
//        throw(std::runtime_error("m3vcfBlock::getAlleles called with an index that is greater than the number of alternates."));
//        return(-1);
//    }
//
//    if(myAlleleCount.size() == 0)
//    {
//        unsigned int gt = 0;
//        myAlleleCount.resize(numAlts+1, 0);
//
//        // Loop through the samples, counting the number of each allele.
//        for(int sampleNum = 0; sampleNum < getNumSamples();
//            sampleNum++)
//        {
//            if((sampleSubset != NULL) &&
//               !(sampleSubset->keep(sampleNum)))
//            {
//                // Skip this sample.
//                continue;
//            }
//            for(int gtNum = 0; gtNum < getNumGTs(sampleNum); gtNum++)
//            {
//                gt = getGT(sampleNum, gtNum);
//                if((gt < 0) || (gt > numAlts))
//                {
//                    // Out of range GT, so continue to the next gt
//                    continue;
//                }
//                // Increment the minor allele count
//                ++myAlleleCount[gt];
//            }
//        }
//    }
//
//    // Alternate allele, so return the alternate.
//    return(myAlleleCount[index]);
//}
//
//


bool m3vcfBlock::readTilTab(IFILE filePtr, std::string& stringRef)
{
    int charRead = 0;
    while(1)
    {
        charRead = ifgetc(filePtr);

        if((charRead == '\n') || (charRead == EOF))
        {
            // Didn't find a tab, found a '\n' or eof
            // It still populated the string with values up
            // until the tab.
            return(false);
        }
        if(charRead == '\t')
        {
            // hit the tab character, so exit the loop.
            break;
        }
        stringRef += charRead;
    }
    return(true);
}
