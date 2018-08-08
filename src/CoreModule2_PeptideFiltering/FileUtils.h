#ifndef _FileUtils_H_
#define _FileUtils_H_

// System Includes
#include <string>

namespace specnets
{
  bool fileExists(const std::string & filename);
  
  std::string getPath(const std::string & basePath, const std::string & addPath, bool relative);

} //namespace specnets

#endif // _FileUtils_H_
