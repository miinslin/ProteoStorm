////////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <stdio.h>

#include "StatusFile.h"
#include "Logger.h"
////////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace specnets;
////////////////////////////////////////////////////////////////////////////////
bool writeStatusFile(string fileName, string status)
{

  FILE * f = fopen(fileName.c_str(), "wb");

  if (f == NULL)
  {
    perror("The following error occurred");
    ERROR_MSG("Problem encountered creating status file");
    exit(-1);
  }

  int val = fwrite(status.c_str(), sizeof(char), status.length(), f);
  if (val != status.length()) {
    ERROR_MSG("Problem encountered writing status file");
    ERROR_MSG(val);
    ERROR_MSG(status.length());
    exit(-1);
  }
  fclose(f);
  DEBUG_MSG("Updated program status to " + status);
}
////////////////////////////////////////////////////////////////////////////////
