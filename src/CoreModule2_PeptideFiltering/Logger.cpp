// Header Include
#include "DelimitedTextReader.h"
#include "Logger.h"

// System Includes
#include <iostream>
#include <map>

using namespace std;

namespace specnets
{
  static map<string,Logger *> m_loggerMap;
  static Logger * m_defaultLogger = 0;

  //--------------------------------------------------------------------------
  Logger & Logger::makeLogger(const string & logName)
  {    
    if (m_loggerMap.find(logName) == m_loggerMap.end())
    {
      if (logName.empty())
      {
        m_loggerMap[logName] = new Logger();
      } 
      else
      {
        m_loggerMap[logName] = new Logger(logName);
      } 
    }
    return *m_loggerMap[logName];
  }

  //--------------------------------------------------------------------------
  void Logger::setDefaultLogger(Logger & defaultLogger)
  {
    m_defaultLogger = &defaultLogger;
  }

  //--------------------------------------------------------------------------
  Logger & Logger::getDefaultLogger(void)
  {
    if (m_defaultLogger != 0)
    {
      return *m_defaultLogger;
    }
    m_defaultLogger = &makeLogger("");
  }

  //--------------------------------------------------------------------------
  Logger & Logger::getLogger()
  {
    Logger & logger = makeLogger("");
    return logger;
  }

  //--------------------------------------------------------------------------
  Logger & Logger::getLogger(unsigned int reportingLevel)
  {
    Logger & logger = makeLogger("");
    logger.setReportingLevel(reportingLevel);
    return logger;
  }

  //--------------------------------------------------------------------------
  Logger & Logger::getLogger(const std::string & fileName)
  {
    Logger & logger = makeLogger(fileName);
    return logger;
  }

  //--------------------------------------------------------------------------
  Logger & Logger::getLogger(const std::string & fileName, unsigned int reportingLevel)
  {
    Logger & logger = makeLogger(fileName);
    logger.setReportingLevel(reportingLevel);
    return logger;
  }

  //--------------------------------------------------------------------------
  void Logger::closeAllLoggers(void)
  {
    map<string,Logger *>::iterator itr = m_loggerMap.begin();
    map<string,Logger *>::iterator itr_end = m_loggerMap.end();
    for ( ; itr != itr_end; itr++)
    {
      delete itr->second;
    }
    m_defaultLogger = 0x0;
  }

  //--------------------------------------------------------------------------
  Logger::Logger(unsigned int reportingLevel)
    : logStream(cerr)
    , m_reportingLevel(reportingLevel)
  {
    // EMPTY
  }

  //--------------------------------------------------------------------------
  Logger::Logger(const std::string & fileName, unsigned int reportingLevel)
    : fileStream(fileName.c_str())
    , logStream(fileStream)
    , m_reportingLevel(reportingLevel)
  {
    // EMPTY
  }

  //--------------------------------------------------------------------------
  Logger::~Logger(void)
  {
    if (fileStream.is_open())
    {
      fileStream.close();
    }
  }

  //--------------------------------------------------------------------------
  unsigned int Logger::getReportingLevel(void)
  {
    return m_reportingLevel;
  }

  //--------------------------------------------------------------------------
  void Logger::output(unsigned int level, const std::string & message)
  {
    if (level >= m_reportingLevel)
    {
      logStream << message << endl;
      logStream.flush();
#if 0      
      if(level >= LOG_CERR_LEVEL) {
        cerr << message << endl;
        cerr.flush();
      } else
      if(level >= LOG_COUT_LEVEL) {
        cout << message << endl;
        cout.flush();
      }
#endif
    }
  }

  //--------------------------------------------------------------------------
  void Logger::outputMemoryUsage(unsigned int level, const std::string & message)
  {
    if (level < m_reportingLevel) {
      return;
    }

    ifstream ifsMem("/proc/meminfo");
    if (!ifsMem.is_open() || !ifsMem.good()) {
      logStream << message << "Unable to open /proc/meminfo for memory info" << endl;
      logStream.flush();
    } else {
      while (!ifsMem.eof()) {
        string lineBuffer;
        DelimitedTextReader::getlineChomp(ifsMem, lineBuffer);
        if (lineBuffer.find("MemTotal") != string::npos ||
            lineBuffer.find("MemFree") != string::npos) {
          logStream << message << lineBuffer << endl;
        }
      }
      logStream.flush();
    }
    ifsMem.close();
    return;
  }

  //--------------------------------------------------------------------------
  void Logger::setReportingLevel(unsigned int level)
  {
    if (level > 10)
    {
      return;
    }   
    m_reportingLevel = level;
  }

} //namespace specnets


