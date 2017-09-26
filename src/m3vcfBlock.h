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
    vector<int> SampleToHaplotypeMapper;


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
        SampleToHaplotypeMapper.clear();
        myBlockHeader.reset();
        myRecords.clear();
        myStatus = StatGenStatus::SUCCESS;
    }

    void SetHeader(m3vcfBlockHeader &ThisHeader)
    {
        myBlockHeader = ThisHeader;
        myRecords.resize(ThisHeader.getNumMarkers());
    }
    
    void SetMapper(vector<int> &ThisMapper)
    {
        SampleToHaplotypeMapper=ThisMapper;
    }
   
    bool readHeader(IFILE filePtr,  m3vcfHeader &ThisHeader,
              bool InfoOnly=false,
              bool siteOnly=false)
    {
        SampleToHaplotypeMapper.clear();
        if(myBlockHeader.read(filePtr, ThisHeader, InfoOnly, siteOnly))
        {
            myRecords.resize(myBlockHeader.getNumMarkers());
            SampleToHaplotypeMapper.resize(ThisHeader.getNumSamples(),0);
            for(int i=1; i<ThisHeader.getNumSamples(); i++)
            {
                SampleToHaplotypeMapper[i]=SampleToHaplotypeMapper[i-1] + myBlockHeader.getSamplePloidy(i-1);
            }
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
            tempIndex[i]=Number+i+1;
        DeleteVariantsList(tempIndex);
    }
    
    void write (m3vcfFileWriter<m3vcfHeader> &filePtr)
    {
        filePtr.writeBlock(myBlockHeader);
        for(int i=0; i<myBlockHeader.getNumMarkers(); i++)
            filePtr.writeRecord(myRecords[i]);
    }
    
    void writeHeader (m3vcfFileWriter<m3vcfHeader> &filePtr)
    {
        filePtr.writeBlock(myBlockHeader);
    }
    void writeRecord (m3vcfFileWriter<m3vcfHeader> &filePtr, int index)
    {
       filePtr.writeRecord(myRecords[index]);
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
    
    void writeVcfRecordGenotypes(IFILE filePtr, int index, bool siteOnly)
    {
        myRecords[index].write(filePtr, true);
        if(siteOnly==false)
            myRecords[index].writeVcfRecordGenotypes(filePtr, myBlockHeader);
        else
            ifprintf(filePtr, "\n");
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
                myRecords[index++].read(filePtr,myBlockHeader,siteOnly);
            }
            return true;         
       }



    
//    bool copyRecord(int index, VcfRecord &thisRecord)
//    {
//        if(index>=myBlockHeader.getNumMarkers())
//            return false;
//        else return myRecords[index].copyRecord(thisRecord);
//    }

   


    void CopyBlockHeader(m3vcfBlockHeader &ThisBlockHeader, m3vcfHeader &ThisHeader) { 
        myBlockHeader=ThisBlockHeader; 
        myRecords.resize(ThisBlockHeader.getNumMarkers());
        SampleToHaplotypeMapper.clear();
        SampleToHaplotypeMapper.resize(ThisHeader.getNumSamples(),0);
        for(int i=1; i<ThisHeader.getNumSamples(); i++)
        {
            SampleToHaplotypeMapper[i]=SampleToHaplotypeMapper[i-1] + myBlockHeader.getSamplePloidy(i-1);
        }
        
    }
    bool isBlockFinished() { return myBlockHeader.isBlockFinished(); }
    bool readRecord(IFILE filePtr) { return myRecords[myBlockHeader.LastReadVariant()].read(filePtr,myBlockHeader); }
    int getEndBasePosition() { return myBlockHeader.getEndBasePosition(); }
    int getStartBasePosition() { return myBlockHeader.getStartBasePosition(); }
    int getNumMarkers(){return myBlockHeader.getNumMarkers();};
    int getNumHaplotypes(){return myBlockHeader.getNumHaplotypes();};
    int getSamplePloidy(int index) {return(myBlockHeader.SampleNoHaplotypes[index]);}
    
    AlleleType getAllele(int RecordIndex, int haploIndex)
    {
        return myRecords[RecordIndex].getAllele(myBlockHeader, haploIndex);
    }
    void CopyToBlockHeader(m3vcfBlockHeader &ThisHeader) { ThisHeader = myBlockHeader;}
    void CopyToBlockHeader(m3vcfBlock &ThisBlock) 
    { 
        ThisBlock.SetHeader(myBlockHeader);  
        ThisBlock.SetMapper(SampleToHaplotypeMapper);
    }
    m3vcfRecord *getM3vcfRecord(int index){if(index<(int)myRecords.size())return &myRecords[index];else return NULL;}
    
    
    void swapPhase(vector<int> SampleIndices)
    {
        int HapMappingIndex = 0;
        for(int i=0;i<(int)SampleIndices.size(); i++)
        {   
            if(myBlockHeader.SampleNoHaplotypes[i]==2)
            {
                int index = SampleToHaplotypeMapper[SampleIndices[i]];
                int temp = myBlockHeader.UniqueIndexMap[index];
                myBlockHeader.UniqueIndexMap[index] = myBlockHeader.UniqueIndexMap[index+1];
                myBlockHeader.UniqueIndexMap[index+1] = temp;
            }              
        }
    }

};







#endif //M3VCFTOOLS_M3VCFBLOCK_H
