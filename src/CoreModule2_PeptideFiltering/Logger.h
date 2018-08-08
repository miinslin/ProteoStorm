#ifndef _Logger_H_
#define _Logger_H_

// System Includes
#include <fstream>
#include <sstream>
#include <string>


#define LOG_CERR_LEVEL  10
#define LOG_COUT_LEVEL  5



#define LOG_MESSAGE0(logger, level, type)         \
{                                                \
  std::stringstream ss;                          \
  ss << type << __FILE__ << " : " << __LINE__ ;  \
  logger.output(level,ss.str());                 \
}
#define LOG_MESSAGE(logger, level, type, msg)    \
{                                                \
  std::stringstream ss;                          \
  ss << type << __FILE__                         \
     << " : " << __LINE__  <<  ": " << msg;      \
  logger.output(level,ss.str());                 \
}
#define MEM_MESSAGE(logger, level, type)         \
{                                                \
  std::stringstream ss;                          \
  ss << type << __FILE__                         \
     << " : " << __LINE__  <<  ": ";             \
  logger.outputMemoryUsage(level,ss.str());      \
}

#define LOG_TRACE(logger)            LOG_MESSAGE0(logger, 1, "DEBUG: ")
#define LOG_DEBUG_MSG(logger, msg)   LOG_MESSAGE(logger, 1, "DEBUG: ", msg)
#define LOG_DEBUG_VAR(logger, var)   LOG_MESSAGE(logger, 1, "DEBUG: ", #var << " = " << var)
#define LOG_WARN_MSG(logger, msg)    LOG_MESSAGE(logger, 5, "WARNING: ", msg)
#define LOG_WARN_VAR(logger, var)    LOG_MESSAGE(logger, 5, "WARNING: ", #var << " = " << var)
#define LOG_ERROR_MSG(logger, msg)   LOG_MESSAGE(logger, 10, "ERROR: ", msg)
#define LOG_ERROR_VAR(logger, var)   LOG_MESSAGE(logger, 10, "ERROR: ", #var << " = " << var)

#define DEBUG_TRACE              LOG_MESSAGE0(Logger::getDefaultLogger(), 1, "DEBUG: ")
#define DEBUG_MSG(msg)           LOG_MESSAGE(Logger::getDefaultLogger(), 1, "DEBUG: ", msg)
#define DEBUG_VAR(var)           LOG_MESSAGE(Logger::getDefaultLogger(), 1, "DEBUG: ", #var << " = " << var)
#define WARN_MSG(msg)            LOG_MESSAGE(Logger::getDefaultLogger(), 5, "WARNING: ", msg)
#define WARN_VAR(var)            LOG_MESSAGE(Logger::getDefaultLogger(), 5, "WARNING: ", #var << " = " << var)
#define ERROR_MSG(msg)           LOG_MESSAGE(Logger::getDefaultLogger(), 10, "ERROR: ", msg)
#define ERROR_VAR(var)           LOG_MESSAGE(Logger::getDefaultLogger(), 10, "ERROR: ", #var << " = " << var)

#define DEBUG_MEM                MEM_MESSAGE(Logger::getDefaultLogger(), 1, "MEMORY: ")


namespace specnets
{
  class Logger
  {
  public:

    static Logger & getLogger(void);
    static Logger & getLogger(unsigned int reportingLevel);
    static Logger & getLogger(const std::string & fileName);
    static Logger & getLogger(const std::string & fileName, unsigned int reportingLevel);
    static void closeAllLoggers(void);
  
    static Logger & getDefaultLogger();
    static void setDefaultLogger(Logger & defaultLogger);
    
    unsigned int getReportingLevel();
    void output(unsigned int level, const std::string & message);
    void outputMemoryUsage(unsigned int level, const std::string & message);

  private:
    static Logger & makeLogger(const std::string & logName);

    Logger(unsigned int reportingLevel = 5);
    Logger(const std::string & fileName, unsigned int reportingLevel = 5);
    ~Logger(void);
    
    void setReportingLevel(unsigned int level);
    
    std::ostream & logStream;
    std::ofstream fileStream;
    
    unsigned int m_reportingLevel;
    
    friend class LoggerCleaner;
  };

  // Simply object to automatically clean up loggers upon exit from scope
  class LoggerCleaner
  {
  public:
    LoggerCleaner(void) {};
    ~LoggerCleaner(void) {Logger::closeAllLoggers();}
  };
  
} //namespace specnets

#endif // _Logger_H_
