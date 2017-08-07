#ifndef M3VCFTOOLS_M3VCFBLOCK_H
#define M3VCFTOOLS_M3VCFBLOCK_H

#include "StatGenStatus.h"
#include "m3vcfBlockHeader.h"
#include "m3vcfRecord.h"

/// This header file provides interface to read/write M3VCF whole blocks with the header and number of records
///
class m3vcfBlock
{

private:
    m3vcfBlock(const m3vcfBlock& thisRecord);
    m3vcfBlock& operator=(const m3vcfBlock& thisRecord);

    m3vcfBlockHeader myBlockHeader;
    vector<m3vcfRecord> myRecords;



    // The status of the last failed command.
    StatGenStatus myStatus;
public:
    /// Default Constructor, initializes the variables, but does not open
    /// any files.
    m3vcfBlock();
    /// Destructor
    virtual ~m3vcfBlock();

    /// Reset function
    void reset()
    {
        myBlockHeader.reset();
        myRecords.clear();
        myStatus = StatGenStatus::SUCCESS;
    }


    bool readHeader(IFILE filePtr,  m3vcfHeader &ThisHeader,
              bool InfoOnly=false,
              bool siteOnly=false)
    {
        if(myBlockHeader.read(filePtr, ThisHeader, InfoOnly, siteOnly))
        {
            myRecords.resize(myBlockHeader.getNumMarkers());
            return true;
        }
        return false;
    }

    void CopyBlockHeader(m3vcfBlockHeader &ThisHeader) { myBlockHeader=ThisHeader; myRecords.resize(ThisHeader.getNumMarkers()); }
    bool isBlockFinished() { return myBlockHeader.isBlockFinished(); }
    bool readRecord(IFILE filePtr) { return myRecords[myBlockHeader.LastReadVariant()].read(filePtr,myBlockHeader); }
    int getEndBasePosition() { return myBlockHeader.getEndBasePosition(); }
    int getStartBasePosition() { return myBlockHeader.getStartBasePosition(); }
    int getNumMarkers(){return myBlockHeader.getNumMarkers();};
    void CopyToBlockHeader(m3vcfBlockHeader &ThisHeader) { ThisHeader = myBlockHeader;}
    m3vcfRecord *getM3vcfRecord(int index){if(index<(int)myRecords.size())return &myRecords[index];else return NULL;}


};







#endif //M3VCFTOOLS_M3VCFBLOCK_H
