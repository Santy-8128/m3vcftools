#ifndef __M3VCF_BLOCK_H__
#define __M3VCF_BLOCK_H__

#include <vector>
#include <stdlib.h>
#include "Tokenize.h"
#include "m3vcfHeader.h"
#include "VcfRecordFilter.h"
#include "VcfRecordInfo.h"
#include "VcfRecordGenotype.h"
#include "StatGenStatus.h"
#include "VcfRecordDiscardRules.h"
#define LENGTH_OF_BLOCK_ID 7
using namespace std;

/// This header file provides interface to read/write VCF files.
class m3vcfBlock
{

private:
    m3vcfBlock(const m3vcfBlock& thisRecord);
    m3vcfBlock& operator=(const m3vcfBlock& thisRecord);

    static const char ALT_DELIM = ',';

    int NoMarkersRead;
    bool FinishedReadingBlock;
    std::string myChrom;
    int startPosition;
    int endPosition;
    int numMarkers;
    int numUniqueReps;
    vector<int> UniqueIndexMap;

    // The status of the last failed command.
    StatGenStatus myStatus;

    std::string myDummyString;
    // Set Pointers to each piece.


public:
    /// Default Constructor, initializes the variables, but does not open
    /// any files.
    m3vcfBlock();
    /// Destructor
    virtual ~m3vcfBlock();

    /// Read the next Vcf data line from the file.
    /// \param filePtr IFILE to read from.
    /// \param siteOnly only store the first 8 columns
    /// \param sampleSubset pointer to sample subset information,
    ///        but NULL if sample subsetting is not to be done.
    /// \param discardRules pointer to record discard information,
    ///        but NULL if record discarding is not to be done.
    /// \return true if a line was successfully read from the specified filePtr,
    /// false if not.
    bool read(IFILE filePtr,  m3vcfHeader &ThisHeader,
              bool InfoOnly=false,
              bool siteOnly=false,
              VcfRecordDiscardRules* discardRules = NULL,
              VcfSubsetSamples* sampleSubset = NULL);


    /// Reset this header, preparing for a new one.
    void reset();

    /// Return if block is finished reading
    bool isBlockFinished()
    {
        return FinishedReadingBlock;
    }

    /// Increments NoMarkersRead by 1
    void AnotherMarkerRead()
    {
        NoMarkersRead++;
//        cout<< " CHANGED TEXT  = "<<NoMarkersRead<<"\t"<<numMarkers<<endl;
        if(NoMarkersRead>=numMarkers)
        {
            FinishedReadingBlock=true;
        }
    }

    /// Returns the status associated with the last method that sets the status.
    /// \return StatGenStatus of the last command that sets status.
    const StatGenStatus& getStatus();







    /// Write this data line to the file (including the newline).
    /// \param filePtr IFILE to write to.
    /// \param siteOnly only write the first 8 columns
    /// \return true if a line was successfully written to the specified filePtr,
    /// false if not.
    bool write(IFILE filePtr, bool siteOnly);

//

////     bool isValid();
//
//    ///////////////////////
//    /// @name  Get Vcf Fields
//    /// Get methods for record fields (do not set status).
//    //@{
        const char* getChromStr() {return(myChrom.c_str());}
        int getStartBasePosition() {return(startPosition);}
        int getEndBasePosition() {return(endPosition);}
        int getNumMarkers() {return(numMarkers);}
        int getNumUniqueReps() {return(numUniqueReps);}

        void initUniqueIndexMap(int n){UniqueIndexMap.resize(n);return;}
        void assignUniqueIndexMap(int &pos, int val){UniqueIndexMap[pos]=val;}

        void setChrom(string &chrom) {myChrom =chrom;};//myChrom = chrom;}
        void setStartBasePosition(int &pos) {startPosition = pos;}
        void setEndBasePosition(int &pos) {endPosition = pos;}
        void setNumMarkers(int pos) {numMarkers = pos;}
        void setNumUniqueReps(int pos) {numUniqueReps = pos;}

//            myChrom = thisRecord.getChromStr();

//    string infoString;
//    vector<AlleleType> UniqueRepAllele;
//    int numUniqueReps;







 //    const char* getIDStr() {return(myID.c_str());}
//    const char* getRefStr() {return(myRef.c_str());}
//    int getNumRefBases() {return(myRef.size());}
//    const char* getAltStr() {return(myAlt.c_str());}
//    /// Return a pointer to the alleles at the specified index with index 0
//    /// being the reference string for this position and index 1 starting
//    /// the alternate alleles, throwing an exception if the index is out of
//    /// range.
//    /// \param index allele index (0 for reference, 1 for first alt,
//    ///  2 for second, etc)
//    /// \return string of the alleles at the specified index
//    const char* getAlleles(unsigned int index);
//    /// Return the int value of the first allele in the string at the
//    /// specified index with index 0 being the reference string for
//    /// this position and index 1 starting the alternate alleles,
//    /// throwing an exception if the index is out of range.
//    /// \param index allele index (0 for reference, 1 for first alt,
//    ///  2 for second, etc)
//    /// \return int allele at the specified index, 1=A, 2=C, 3=G, 4=T
//    int getIntAllele(unsigned int index);
//    /// Return the number of alternates listed in the Alts string.
//    unsigned int getNumAlts();
//
//    float getQual() {return(myQualNum);}
//    const char* getQualStr() {return(myQual.c_str());}
//
//    /// Return a reference to the filter information.
//    VcfRecordFilter& getFilter(){return(myFilter);}
//    /// Return whether or not all filters were passed.
//    int passedAllFilters() { return(myFilter.passedAllFilters()); }
//
//    /// Get a reference to the information field.
//    VcfRecordInfo& getInfo() {return myInfo;}
//
//    /// Get a reference to the genotype fields.
//    VcfRecordGenotype& getGenotypeInfo() {return myGenotype;}
//
//    inline int getNumSamples() { return(myGenotype.getNumSamples()); }
//
//    inline int getNumGTs(int index) { return(myGenotype.getNumGTs(index)); }
//
//    inline int getGT(int sampleNum, unsigned int gtIndex)
//    { return(myGenotype.getGT(sampleNum, gtIndex)); }
//
//    inline void setGT(int sampleNum, unsigned int gtIndex, int newGt)
//    { myGenotype.setGT(sampleNum, gtIndex, newGt); }
//
//    /// Return true if all of the samples are phased and none are unphased,
//    /// false if any are unphased or not phased.
//    inline bool allPhased() { return(myGenotype.allPhased()); }
//
//    /// Return true if all of the samples are unphased and none are phased,
//    /// false if any are phased or not unphased.
//    inline bool allUnphased() { return(myGenotype.allUnphased()); }
//
//    /// Return true if all samples of all records have all the genotype alleles
//    /// specified, false if not or if any GT field is missing.
//    bool hasAllGenotypeAlleles() { return(myGenotype.hasAllGenotypeAlleles()); }
//
//    /// Return the number of occurances of the specified allele index in the
//    /// genotypes for this record.  Index 0 for the reference.  The alternate
//    /// alleles start with index 1.  An exception is thrown if the index is
//    /// out of range.  Optionally, the specified subset of samples can be
//    /// skipped when determining allele counts.  (If the record is read with
//    /// just a subset of samples, those are automatically excluded here
//    /// regardless of the passed in sampleSubset.
//    /// \param index allele index (0 for reference, 1 for first alt,
//    ///  2 for second, etc)
//    /// \param sampleSubset pointer to sample subset information,
//    ///        but NULL if additional sample subsetting is not to be done.
//    /// \return int allele count for the specified ref/alt.
//    int getAlleleCount(unsigned int index,
//                       VcfSubsetSamples* sampleSubset = NULL);
//
//    //@}
//
//
//    ///////////////////////
//    /// @name  Set Vcf Fields
//    /// Set methods for record fields (do not set status).
//    //@{
//    void setChrom(const char* chrom) {myChrom = chrom;}
//    void set1BasedPosition(int pos) {my1BasedPosNum = pos;}
//    void setID(const char* id) {myID = id;}
//    void setRef(const char* ref) {myRef = ref;}
//    void setAlt(const char* alt) {myAlt = alt; myAltArray.clear();}
//    //    void setQual(float qual) {myQualNum = qual; }
//    void setQual(const char* qual)
//    {
//        myQual = qual;
//        if(myQual != ".")
//        {
//            myQualNum = atof(qual);
//        }
//        else
//        {
//            myQualNum = -1;
//        }
//    }

protected:

    /// Read the specified file until a tab, '\n', or EOF is found
    /// and appending the read characters to stringRef (except for the
    /// stopping character).
    /// \param filePtr open file to be read.
    /// \param stringRef reference to a string that should be appended to
    /// with the characters read from the file until a '\t', '\n', or EOF
    /// is found.
    /// \return true if a '\t' stopped the reading, false if '\n' or
    /// EOF stopped the reading.
    bool readTilTab(IFILE filePtr, std::string& stringRef);

};





#endif
