#include "m3vcfRecord.h"

void m3vcfRecord::copyStartInfotoBlock(m3vcfBlock &thisBlock)
{
    thisBlock.setChrom(myChrom);
    thisBlock.setStartBasePosition(BasePositionVal);
}


void m3vcfRecord::copyEndInfotoBlock(m3vcfBlock &thisBlock)
{
    thisBlock.setEndBasePosition(BasePositionVal);
}


void m3vcfRecord::copyFromVcfRecord(VcfRecord &thisRecord)
{
    myChrom = thisRecord.getChromStr();
    BasePositionVal = thisRecord.get1BasedPosition();
    varID = thisRecord.getIDStr();
    refAlleleString = thisRecord.getRefStr();
    allAltAlleleString = thisRecord.getAltStr();
    infoString.clear();
    for (int i = 0; i < thisRecord.getInfo().getNumInfoFields(); ++i)
    {
        std::pair<std::string, std::string> p = thisRecord.getInfo().getInfoPair(i);
        if (i != 0)
            infoString += ";";
        infoString += p.first;
        if (p.second.size())
            infoString += "=" + p.second;
    }

}


m3vcfRecord::m3vcfRecord()
{
    reset();
}


m3vcfRecord::~m3vcfRecord()
{
}


bool m3vcfRecord::read  (IFILE filePtr,
                          m3vcfBlock &ThisBlock,
                        bool siteOnly,
                        VcfRecordDiscardRules* discardRules,
                        VcfSubsetSamples* sampleSubset)
{
    // Clear out any previously set values.
    reset();
    std::string tempString;

    if(filePtr == NULL)
    {
        myStatus.setStatus(StatGenStatus::FAIL_ORDER,
                           "Error reading M3VCF record before opening the file.");
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


    // Read the Base Position
    if(!readTilTab(filePtr, BasePosition))
    {
        myStatus.setStatus(StatGenStatus::FAIL_PARSE,
                           "Error reading M3VCF Record VARIANT POSITION.");
        return(false);
    }
    else
    {
        // Read the positions, so convert them to an integer.
        BasePositionVal = atoi(BasePosition.c_str());
    }

    // Read the Variant ID.
    if(!readTilTab(filePtr, varID))
    {
        myStatus.setStatus(StatGenStatus::FAIL_PARSE,
                           "Error reading M3VCF Record ID.");
        return(false);
    }

    // Read the Ref.
    if(!readTilTab(filePtr, refAlleleString))
    {
        myStatus.setStatus(StatGenStatus::FAIL_PARSE,
                           "Error reading M3VCF Record REF.");
        return(false);
    }

    // Read the Alt.
    altAlleleStringArray.clear();
    if(!readTilTab(filePtr, allAltAlleleString))
    {
        myStatus.setStatus(StatGenStatus::FAIL_PARSE,
                           "Error reading M3VCF Record ALT.");
        return(false);
    }

    readTilTab(filePtr, tempString);
    readTilTab(filePtr, tempString);

    // Read the Info String.
    if(!readTilTab(filePtr, infoString))
    {
        myStatus.setStatus(StatGenStatus::FAIL_PARSE,
                           "Error reading M3VCF Record INFO.");
        return(false);
    }

    ThisBlock.AnotherMarkerRead();

//abort();
    if(siteOnly)
    {
        // Do not store genotypes, so just consume the rest of the line.
        filePtr->readTilChar("\n");
        return true;
    }

    // Get Unique Representative Alleles
    numUniqueReps=ThisBlock.getNumUniqueReps();
    tempString.clear();
    if(readTilTab(filePtr, tempString))
    {
        cout<<" WHOA = "<<tempString<<endl;
        myStatus.setStatus(StatGenStatus::FAIL_PARSE,
                           "Error reading M3VCF UNIQUE REP ALLELES");
        return(false);
    }
    int index = 0;
    UniqueRepAllele.resize(numUniqueReps,'0');
    while(index<numUniqueReps)
    {
        UniqueRepAllele[index]=tempString[index];
        index++;
    }
    return(true);
}



void m3vcfRecord::reset()
{
    myChrom.clear();
    BasePosition.clear();
    varID.clear();
    refAlleleString.clear();
    allAltAlleleString.clear();
    altAlleleStringArray.clear();
    infoString.clear();
    UniqueRepAllele.clear();

    BasePositionVal = -2;
    numUniqueReps = 0;

    myStatus = StatGenStatus::SUCCESS;
}


// Return the error after a failed call.
const StatGenStatus& m3vcfRecord::getStatus()
{
    return(myStatus);
}

bool m3vcfRecord::write(IFILE filePtr, bool siteOnly)
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

    if(BasePositionVal == 0)
    {
        ifprintf(filePtr, ".\t");
    }
    else
    {
        ifprintf(filePtr, "%d\t", BasePositionVal);
    }

    if(varID.length() == 0)
    {
        ifprintf(filePtr, ".\t");
    }
    else
    {
        ifprintf(filePtr, "%s\t", varID.c_str());
    }

    if(refAlleleString.length() == 0)
    {
        ifprintf(filePtr, ".\t");
    }
    else
    {
        ifprintf(filePtr, "%s\t", refAlleleString.c_str());
    }

    if(allAltAlleleString.length() == 0)
    {
        ifprintf(filePtr, ".\t");
    }
    else
    {
        ifprintf(filePtr, "%s\t", allAltAlleleString.c_str());
    }

    ifprintf(filePtr, ".\t");
    ifprintf(filePtr, ".\t");


    if(infoString.length() == 0)
    {
        ifprintf(filePtr, ".\t");
    }
    else
    {
        ifprintf(filePtr, "%s\t", infoString.c_str());
    }


    for(int i=0; i<numUniqueReps; i++)
    {
        ifprintf(filePtr, "%c", UniqueRepAllele[i]);
    }
    ifprintf(filePtr, "\n");

   return true;
}




//
//const char* m3vcfRecord::getAlleles(unsigned int index)
//{
//    if(index == 0)
//    {
//        return(myRef.c_str());
//    }
//    if(index > getNumAlts())
//    {
//        // Index out of range.
//        // Throw an exception.
//        throw(std::runtime_error("m3vcfRecord::getAlleles called with an index that is greater than the number of alternates."));
//        return(NULL);
//    }
//    // Alternate allele, so return the alternate.
//    return(myAltArray.get(index-1).c_str());
//}

//
//int m3vcfRecord::getIntAllele(unsigned int index)
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
//            std::cerr << "m3vcfRecord::getIntAllele, unknown allele, "
//                      << alleles[0] << std::endl;
//    }
//    return(0);
//}

//
//unsigned int m3vcfRecord::getNumAlts()
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
//int m3vcfRecord::getAlleleCount(unsigned int index,
//                              VcfSubsetSamples* sampleSubset)
//{
//    unsigned int numAlts = getNumAlts();
//    if(index > numAlts)
//    {
//        // Index out of range.
//        // Throw an exception.
//        throw(std::runtime_error("m3vcfRecord::getAlleles called with an index that is greater than the number of alternates."));
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


bool m3vcfRecord::readTilTab(IFILE filePtr, std::string& stringRef)
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

