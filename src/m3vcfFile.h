#ifndef __M3VCF_FILE_H__
#define __M3VCF_FILE_H__

#include "InputFile.h"
#include "m3vcfHeader.h"

/// This header file provides interface to read/write VCF files.
class m3vcfFile {
public:
    /// Default Constructor, initializes the variables, but does not open
    /// any files.
    m3vcfFile();
    /// Destructor
    virtual ~m3vcfFile();

    /// Close the file if it is open.
    void close();

    /// When set to true, read only the first 8 columns, skipping the format
    /// and genotype fields, so when reading do not store them, and when
    /// writing do not write them.  Defaults to read/write all columns.
    /// This setting is maintained even when the file is reset/closed.
    /// \param siteOnly process only the first 8 columns
    void setSiteOnly(bool siteOnly) {mySiteOnly = siteOnly;}

    /// Get the number of m3vcf records that have been processed (read/written)
    /// so far including any filtered records.
    int getNumRecords() {return(myNumRecords);}

    // Get the Status of the last call that sets status.
    //    inline StatGenStatus::Status getStatus()
    //    {
    //        return(myStatus.getStatus());
    //    }

protected:
    // Open the m3vcf file with the specified filename
    // with the specified mode.
    // \param  filename the m3vcf file to open.
    // \param  mode how to open (r/w).
    // \return true = success; false = failure.
    bool open(const char* filename, const char* mode,
              InputFile::ifileCompression compressionMode = InputFile::DEFAULT);

    void reset();
//    virtual void resetFile() = 0;

    IFILE  myFilePtr;

    StatGenStatus myStatus;

    bool mySiteOnly;

    // Number of records read/written so far.  Child classes need to set this.
    int myNumRecords;

private:
    m3vcfFile(const m3vcfFile& vcfFile);
    m3vcfFile& operator=(const m3vcfFile& vcfFile);
};

#endif
