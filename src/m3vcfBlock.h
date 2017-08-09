#ifndef M3VCFTOOLS_M3VCFBLOCK_H
#define M3VCFTOOLS_M3VCFBLOCK_H

#include "StatGenStatus.h"
#include "m3vcfBlockHeader.h"
#include "m3vcfRecord.h"
#include "m3vcfFileWriter.h"

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
    
    void DeleteVariantsFrom(int Number)
    {
        vector<int> tempIndex(Number);
        for(int i=0; i<Number; i++)tempIndex[i]=i;
        DeleteVariantsList(tempIndex);
    }
    void DeleteVariantsTill(int Number)
    {
        vector<int> tempIndex(myBlockHeader.getNumMarkers() - Number - 1);
        for(int i=0; i<(int)tempIndex.size(); i++)
            tempIndex[i]=Number+i;
        DeleteVariantsList(tempIndex);
    }
    
    void write (m3vcfFileWriter<m3vcfHeader> &filePtr)
    {
        filePtr.writeBlock(myBlockHeader);
        for(int i=0; i<myBlockHeader.getNumMarkers(); i++)
            filePtr.writeRecord(myRecords[i]);
    }
    
    void DeleteVariantsList(vector<int> &ImportIndex)
    {
        myBlockHeader.NoMarkersRead=0;
        myBlockHeader.FinishedReadingBlock=false;
        myBlockHeader.startPosition=myRecords[ImportIndex[0]].getBasePosition();
        myBlockHeader.endPosition=myRecords[ImportIndex.back()].getBasePosition();
        myBlockHeader.numMarkers=(int)ImportIndex.size();
        
        for(int i=0; i<(int)ImportIndex.size(); i++)
        {
            myRecords[i]=myRecords[ImportIndex[i]];            
        }
        myRecords.resize(myBlockHeader.numMarkers);
    }
    
    void writeVcfRecordGenotypes(IFILE filePtr, int index)
    {
        myRecords[index].write(filePtr, true);
        myRecords[index].writeVcfRecordGenotypes(filePtr, myBlockHeader);
    }
    
           bool read(IFILE filePtr,  m3vcfHeader &ThisHeader,
              bool siteOnly=false)
       {
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

            if(!myBlockHeader.read(filePtr, ThisHeader))
            {
                // End of file, just return false.
                return(false);
            }
            
            myRecords.resize(myBlockHeader.getNumMarkers());
            int index=0;
            while(!myBlockHeader.isBlockFinished())
            {
                myRecords[index++].read(filePtr,myBlockHeader);
            }
            return true;         
       }
           
        

    bool copyRecord(int index, VcfRecord &thisRecord)
    {
        if(index>=myBlockHeader.getNumMarkers())
            return false;
        else return myRecords[index].copyRecord(thisRecord);
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
