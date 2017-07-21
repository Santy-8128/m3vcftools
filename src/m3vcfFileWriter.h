
#ifndef __M3VCF_FILE_WRITER_H__
#define __M3VCF_FILE_WRITER_H__

#include "m3vcfFile.h"
#include "m3vcfBlock.h"
#include "m3vcfRecord.h"

/// This header file provides interface to write M3VCF files.

template <class HeaderType> class m3vcfFileWriter : public m3vcfFile
{
public:
    /// Default Constructor, initializes the variables, but does not open
    /// any files.
    m3vcfFileWriter();
    /// Destructor
    virtual ~m3vcfFileWriter();

    /// Open the m3vcf file with the specified filename for writing.
    /// \param filename the m3vcf file to open for writing.
    /// \param header to be written the file
    /// \param compressionMode type of compression to use for writing
    /// \return true = success; false = failure.
    bool open(const char* filename, HeaderType& header,
              InputFile::ifileCompression compressionMode);

    /// Open the m3vcf file with the specified filename for writing using the
    /// default compression (BGZF).
    /// \param filename the m3vcf file to open for writing.
    /// \param header to be written the file
    /// \return true = success; false = failure.
    bool open(const char* filename, HeaderType& header);

    /// Write the M3VCF data line to the file.
    /// \param record record to write to the file.
    /// \return true if successfully wrote, false if not.
    bool writeRecord(m3vcfRecord& record);

    /// Write the M3VCF data line to the file.
    /// \param record record to write to the file.
    /// \return true if successfully wrote, false if not.
    bool writeBlock(m3vcfBlock& block);

protected:
    virtual void resetFile() {}

private:
    m3vcfFileWriter(const m3vcfFileWriter& ThisFileWriter);
    m3vcfFileWriter& operator=(const m3vcfFileWriter& ThisFileWriter);
};


template <class HeaderType> m3vcfFileWriter<HeaderType> ::m3vcfFileWriter()
    : m3vcfFile()
{
}


template <class HeaderType> m3vcfFileWriter<HeaderType>::~m3vcfFileWriter()
{
}


template <class HeaderType> bool m3vcfFileWriter<HeaderType>::open(const char* filename, HeaderType& header,
                         InputFile::ifileCompression compressionMode)
{
//    myStatus = StatGenStatus::SUCCESS;
    if(m3vcfFile::open(filename, "w", compressionMode))
    {
        // Successfully opened, so write the header.
        if(!header.write(myFilePtr))
        {
            // Failed, so copy the status.
            myStatus = header.getStatus();
            return(false);
        }
    }
    else
    {
        // Failed, status set by m3vcfFile::open.
        return(false);
    }

    // Successfully opened and read the header.
    return(true);
}


template <class HeaderType> bool m3vcfFileWriter<HeaderType>::open(const char* filename, HeaderType& header)
{
    return(open(filename, header, InputFile::UNCOMPRESSED));
}

template <class HeaderType> bool m3vcfFileWriter<HeaderType>::writeBlock(m3vcfBlock& block)
{
    if(!block.write(myFilePtr, mySiteOnly))
    {
        myStatus = block.getStatus();
        return(false);
    }
    ++myNumRecords;
    return(true);
}
template <class HeaderType> bool m3vcfFileWriter<HeaderType>::writeRecord(m3vcfRecord& record)
{
    if(!record.write(myFilePtr, mySiteOnly))
    {
        myStatus = record.getStatus();
        return(false);
    }
    ++myNumRecords;
    return(true);
}




#endif
