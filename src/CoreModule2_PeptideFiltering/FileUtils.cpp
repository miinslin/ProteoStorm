// Header Include
#include "FileUtils.h"
#include "Logger.h"

// System Includes
#include <fstream>
#include <sys/stat.h>

using namespace std;

namespace specnets
{
  bool fileExists(const std::string & filename)
  {
    struct stat fileInfo;
    // stat() returns 0 for success
    if (stat(filename.c_str(), &fileInfo) != 0)
    {
      return false;
    }

    return true;
#if 0
    ifstream is(filename.c_str(), ios::in | ios::binary);
    if (!is)
    {
      return false;
    }
    is.close();
    return true;
#endif
  }

  string getPath(const string & basePath,
                 const string & relPath,
                 bool relative)
  {
    //DEBUG_VAR(basePath);
    //DEBUG_VAR(relPath);

    if (relPath.empty()) {
      //DEBUG_VAR(basePath);
      return basePath;
    }

    if (relative) {
      //DEBUG_VAR(relPath);
      return relPath;
    }

    string returnPath(basePath);
    string addPath(relPath);

    bool up = false;
    if (addPath[0] == '.') {
      addPath = addPath.substr(1, addPath.length() - 1);
      //DEBUG_VAR(addPath);
      if (!addPath.empty() && addPath[0] == '.') {
        addPath = addPath.substr(1, addPath.length() - 1);
        //DEBUG_VAR(addPath);
        up = true;
      }
    }
    //DEBUG_VAR(addPath);
    //DEBUG_VAR(up);

    returnPath = basePath;
    if (returnPath[returnPath.length() - 1] == '/') {
      returnPath = returnPath.substr(1, returnPath.length() - 1);
    }
    //DEBUG_VAR(returnPath);

    if (up) {
      size_t lastSlashPos = returnPath.find_last_of('/');
      if (lastSlashPos != string::npos) {
        returnPath = returnPath.substr(0, lastSlashPos);
      }
    }
    //DEBUG_VAR(returnPath);

    if (addPath[0] != '/') {
      returnPath += '/';
    }
    returnPath += addPath;

    //DEBUG_VAR(returnPath);

    return returnPath;
  }

} //namespace specnets
