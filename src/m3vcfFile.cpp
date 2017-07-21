#include "m3vcfFile.h"

m3vcfFile::m3vcfFile()
{
  myFilePtr = NULL;
  mySiteOnly = false;
  myNumRecords = 0;
}

m3vcfFile::~m3vcfFile()
{
  // Close the file.
  if (myFilePtr != NULL)
  {
      // If we already have an open file, close it.
      ifclose(myFilePtr);
      myFilePtr = NULL;
  }
}


bool m3vcfFile::open(const char* filename, const char* mode,
                   InputFile::ifileCompression compressionMode)
{
    // Reset for any previously operated on files.
    reset();

    myFilePtr = ifopen(filename, mode, compressionMode);

    if(myFilePtr == NULL)
    {
        std::string errorMessage = "Failed to Open ";
        errorMessage += filename;
        errorMessage += " for ";
        errorMessage += mode;
        myStatus.setStatus(StatGenStatus::FAIL_IO, errorMessage.c_str());
        return(false);
    }

    return(true);
}


void m3vcfFile::close()
{
    reset();
}


void m3vcfFile::reset()
{
    // Reset the child class.
//    resetFile();

    // Close the file.
    if (myFilePtr != NULL)
    {
        // If we already have an open file, close it.
        ifclose(myFilePtr);
        myFilePtr = NULL;
    }
    myNumRecords = 0;
}
