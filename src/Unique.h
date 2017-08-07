#ifndef UNIQUE_H_INCLUDED
#define UNIQUE_H_INCLUDED
#include<cmath>
#include<fstream>
#include "StringBasics.h"
#include<vector>
#include "assert.h"
#include "m3vcfRecord.h"
#include "m3vcfBlockHeader.h"
using namespace std;

template <class HaploType> class Compressor
{

    private:

        bool FIXED_LENGTH;
        int transFactor,cisFactor;
        int numHaplotypes, bufferSize;
        int maxAllele;
        vector<int> index,oldIndex;
        vector<int> previousDifference;
        vector<int> previousPredecessor;
        vector<int> firstDifference;
        vector<int> cost;
        vector<int> bestSlice;
        vector<int> bestComplexity;
        vector<vector<int> > bestIndex;
        vector<int> blockBoundaries;
        int BlockStartPointer;
        int blockEnd;
        vector<int> examplars;
        int flushMarkerIndex;

    public:
    /// Default Constructor, initializes the variables, but does not open
    /// any files.
    Compressor()
    {
    }
    /// Destructor
    virtual ~Compressor()
    {
    }
 void       CompressChunk       (HaploType & Haplotypes, bool fixedLength);
 void       Initialize          ();
 void       FlushBlocks         (HaploType & haplotypes);
 void       UpdateDeltaMatrix   (HaploType & haplotypes,int length);
    void       AnalyzeBlocks       (int length);
    void       SimpleAnalyzeBlocks       (int length);
 bool       NewBlockReady       ()
 {
    if(flushMarkerIndex==0) return true;
    return (BlockStartPointer<flushMarkerIndex);
 }
    int        getBlockHeaderEndPosition() {return blockEnd;};
    void       GetM3vcfRecord               (m3vcfRecord &thisM3vcfRecord, HaploType & Haplotypes);
    void       GetBlockHeader               (HaploType & Haplotypes, m3vcfBlockHeader &thisBlock);
    void       CreateBoundaries             (HaploType & Haplotypes);
    void       CreateSimpleBoundaries       (HaploType & Haplotypes);


};


template <class HaploType> void Compressor<HaploType>::CompressChunk
                            (HaploType & Haplotypes, bool fixedLength)
{
    numHaplotypes=Haplotypes.size();
    bufferSize=Haplotypes[0].length();
    FIXED_LENGTH = fixedLength;
    Initialize();

    for(int length=1;length<=Haplotypes[0].length();length++)
    {

        vector<int> offsets(maxAllele+1,0);
        for (int i = 0; i < numHaplotypes; i++)
            offsets[Haplotypes[i][length - 1] - '0' + 1]++;
        for (int i = 1; i < maxAllele; i++)
            offsets[i] += offsets[i - 1];

        oldIndex = index;
        for (int i = 0; i < numHaplotypes; i++)
        {
            index[offsets[Haplotypes[oldIndex[i]][length - 1] - '0']++] = oldIndex[i];
        }

        UpdateDeltaMatrix(Haplotypes, length);
        if(!FIXED_LENGTH)
            AnalyzeBlocks(length);
//        else
//            SimpleAnalyzeBlocks(length);

    }

    if(Haplotypes[0].length()>1)
    {
        if(FIXED_LENGTH)
            CreateSimpleBoundaries(Haplotypes);
        else
            CreateBoundaries(Haplotypes);
    }
    else
    {
        abort();
    }

}



template <class HaploType> void Compressor<HaploType>::CreateBoundaries
                            (HaploType & Haplotypes)
{
    int where   = Haplotypes[0].length()-1;
    blockBoundaries.clear();
    while (where != 0)
    {
        blockBoundaries.push_back(where);
        where = where - bestSlice[where+1]+1;
    };
    BlockStartPointer=0;
    flushMarkerIndex=0;
}


template <class HaploType> void Compressor<HaploType>::CreateSimpleBoundaries
        (HaploType & Haplotypes)
{

    int where   = Haplotypes[0].length()-1;

    bestIndex[where] = index;
    blockBoundaries.clear();
    blockBoundaries.push_back(where);
    BlockStartPointer=0;
    flushMarkerIndex=0;
}

template <class HaploType> void Compressor<HaploType>::GetBlockHeader
                            (HaploType & Haplotypes,
                             m3vcfBlockHeader &thisBlock)
{
    blockEnd = blockBoundaries.back();
    blockBoundaries.pop_back();
    examplars.clear();
    examplars.push_back(bestIndex[blockEnd][0]);
    flushMarkerIndex=BlockStartPointer;

    int zero=0;



    thisBlock.initUniqueIndexMap(numHaplotypes);
    thisBlock.assignUniqueIndexMap(bestIndex[blockEnd][0], zero);
    int countCard=0;

    for (int i = 1; i < numHaplotypes; i++)
    {
        int previous = bestIndex[blockEnd][i-1];
        int current  = bestIndex[blockEnd][i];
        countCard++;
        for (int j = BlockStartPointer; j <= blockEnd; j++)
            if (Haplotypes[previous][j] != Haplotypes[current][j])
            {
                examplars.push_back(current);
                break;
            }
        thisBlock.assignUniqueIndexMap(current, examplars.size() - 1);

     }

    thisBlock.setNumMarkers(blockEnd-BlockStartPointer+1);
    thisBlock.setNumUniqueReps((int)examplars.size());


    BlockStartPointer = blockEnd;
    return;
}


template <class HaploType> void Compressor<HaploType>::GetM3vcfRecord
                            (m3vcfRecord &thisM3vcfRecord,
                             HaploType & Haplotypes)
{
    thisM3vcfRecord.setNumUniqueReps(examplars.size());
    thisM3vcfRecord.initUniqueRepAllele(examplars.size());
    thisM3vcfRecord.initAltHapIndex();
    for (int j = 0; j < (int)examplars.size(); j++)
    {
        thisM3vcfRecord.assignUniqueRepAllele(j,Haplotypes[examplars[j]][flushMarkerIndex]);

        if(Haplotypes[examplars[j]][flushMarkerIndex]=='1')
        {

            thisM3vcfRecord.PushThisIndex(j);
        }
    }
    flushMarkerIndex++;
    return;
}


template <class HaploType> void Compressor<HaploType>::AnalyzeBlocks
                                (int length)
{
    // Try to figure out optimal block split
    int blockSize=length;
    for (int i = 1; i <= blockSize && i <= length; i++)
    {
        int distinctHaplos = 1;

        for (int j = 0; j < numHaplotypes - 1; j++)
            if (i > firstDifference[j])
                distinctHaplos++;

      int currentCost=1;

      if(i>1)
        currentCost= transFactor * numHaplotypes + (i) * distinctHaplos  * cisFactor + cost[length - i+1];

      if (i==2)
         {
             cost[length] = currentCost;
             bestSlice[length] = 2;
//             bestComplexity[length] = distinctHaplos;
         }
      else if (cost[length] > currentCost)
         {
         cost[length] = currentCost;
         bestSlice[length] = i;
//         bestComplexity[length] = distinctHaplos;
         }
      else if (cost[length] + transFactor * numHaplotypes < currentCost)
         break;
      }
    bestIndex[length-1] = index;

}

template <class HaploType> void Compressor<HaploType>::SimpleAnalyzeBlocks
        (int length)
{
    bestIndex[length-1] = index;
}



template <class HaploType> void Compressor<HaploType>::Initialize()
{
    index.resize(numHaplotypes);
    oldIndex.resize(numHaplotypes);
    previousDifference.resize(numHaplotypes);
    previousPredecessor.resize(numHaplotypes);
    firstDifference.resize(numHaplotypes-1,0);
    cost.resize(bufferSize+1,0);
    bestSlice.resize(bufferSize+1,0);
    bestComplexity.resize(bufferSize+1,0);
    bestIndex.resize(bufferSize+1);
    transFactor = 3;
    cisFactor = 2;

    for(int i=0;i<numHaplotypes;i++)
        index[i]=i;

    maxAllele=2; //MaxAllele[ThisPiece][length-1];

}


template <class HaploType> void Compressor<HaploType>::UpdateDeltaMatrix
                            (HaploType & haplotypes,
                             int length)
{
    int blockSize=length;

    previousPredecessor[oldIndex[0]] = -1;
    for (int i = 1; i < numHaplotypes; i++)
    {
        previousPredecessor[oldIndex[i]] = oldIndex[i-1];
        previousDifference[oldIndex[i]]  = firstDifference[i-1];
    }

    for (int i = 1; i < numHaplotypes; i++)
    {

        if (index[i-1] == previousPredecessor[index[i]])
        {
            firstDifference[i-1] =
            haplotypes[index[i]][length - 1] ==
            haplotypes[index[i-1]][length - 1] ?
            previousDifference[index[i]] + 1 : 0;
            continue;
        }

        int diff = 0;
        while (diff < length && diff < blockSize &&
               haplotypes[index[i]][length - diff - 1] ==
               haplotypes[index[i-1]][length - diff - 1])
        diff++;

        firstDifference[i - 1] = diff;
    }
}





#endif // UNIQUE_H_INCLUDED
