
#ifndef __M3VCF_HEADER_H__
#define __M3VCF_HEADER_H__

#include <vector>
#include <VcfHeader.h>
#include "StringArray.h"
#include "StatGenStatus.h"
#define NUM_M3VCF_NON_SAMPLE_HEADER_COLS 9
#define M3VCF_VERSION "2.0"
using namespace std;

/// This header file provides interface for dealing with M3VCF Meta/Header lines.
class m3vcfHeader
{
    
    
private:
    m3vcfHeader(const m3vcfHeader& thisHeader);
    m3vcfHeader& operator=(const m3vcfHeader& thisHeader);

    // Make sure the last header line is synced with the parsed header line.
    // This is used when samples are removed.
    void syncHeaderLine();

    static const int NUM_NON_SAMPLE_HEADER_COLS = 9;

    // Is set to true once the header line has been set, false until then.
    bool myHasHeaderLine;

    int numSamples;
    std::vector<string> myHeaderLines;
    std::vector<string> SampleNames;
    std::vector<int> SamplePloidy;

    StringArray myParsedHeaderLine;

    // The status of the last failed command.
    StatGenStatus myStatus;
    
    
public:
    /// Default Constructor, initializes the variables.
    m3vcfHeader();
    /// Destructor
    virtual ~m3vcfHeader();

    /// Read the header from the specified file replacing any previous header
    /// contents.
    /// \param filePtr IFILE to read from.
    /// \return true if an entire meta/header was successfully read from
    /// the specified filePtr, false if not.
    bool read(IFILE filePtr);

    /// This function copies header information from VcfHeader
    /// and updates to m3vcf header
    void copyHeader(VcfHeader &thisHeader);
    
    /// Write the header to the specified file.
    /// \param filePtr IFILE to write to.
    /// \return true if an entire meta/header was successfully written to
    /// the specified filePtr, false if not.
    bool write(IFILE filePtr);

    /// Reset this header, preparing for a new one.
    void reset();

     /// Sparse the Header Line and get names of samples and their ploidy
    bool SparseSampleNames();

    
    /// Returns the status associated with the last method that sets the status.
    /// \return StatGenStatus of the last command that sets status.
    const StatGenStatus& getStatus();

    /// Return the number of meta-lines (lines starting with ##)
    int getNumMetaLines();

    /// Return the specified meta-line (index starting at 0)
    /// or NULL if out of range.
    /// Will return the headerline if the header line's index is specified.
    const char* getMetaLine(unsigned int index);

    /// Return the header line, the line containing #chrom...
    const char* getHeaderLine();

    /// Returns the number of samples in the header line or 0 if the header
    /// line has not yet been read.
    int getNumSamples()  {return numSamples;} ;

    /// Returns the name of the specified sample or NULL if the sample number
    /// is out of range (first sample is index 0).
    const char* getSampleName(unsigned int index) const;

    /// Returns the index of the specified sample or -1 if the sample name
    /// is not found (first sample is index 0).
    int getSampleIndex(const char* sampleName) const;

    /// Remove the sample at the specified index.
    void removeSample(unsigned int index);

    /////////////////
    /// Add Lines

    /// Add a Meta Line to the end of the currently specified meta lines.
    /// Return false if the meta line is invalid (does not start with ##)
    /// A return of false means the line was not added.
    bool appendMetaLine(const char* metaLine);

    /// Replace the header line if one exists or add it if one does not.
    /// Return false if the header line is invalid (does not start with #).
    /// A return of false means the line was not added.
    bool addHeaderLine(const char* headerLine);

    /// Check if metaLine already exists in Header
    bool ifExistsMetaLine(const char* metaLine);

    /// Check Merge if possible before merging
    bool checkMergeHeader(m3vcfHeader &Header);

    /// Merge Header Information
    void mergeHeader(m3vcfHeader &Header);

protected:

};


#endif
