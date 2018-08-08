////////////////////////////////////////////////////////////////////////////////
#ifndef __SPECIFIC_H__
#define __SPECIFIC_H__
////////////////////////////////////////////////////////////////////////////////
#include <string>
#include <dirent.h>

using namespace std;

#if !defined(ACCESSPERMS)
#define ACCESSPERMS (S_IRWXU|S_IRWXG|S_IRWXO)
#endif

////////////////////////////////////////////////////////////////////////////////
// This method diverts the segfault to a method that outputs
// a message to the logfile
int addSegFaultDivert(void);

////////////////////////////////////////////////////////////////////////////////
// set environment variables.
void mysetenv(const char *name, const char *value);

////////////////////////////////////////////////////////////////////////////////
// wait for a process to finish
void waitProcess(void);

////////////////////////////////////////////////////////////////////////////////
// is fork() available?
bool forkable(void);

////////////////////////////////////////////////////////////////////////////////
// fork. Empty function if fork() is not available
int myfork(void);

////////////////////////////////////////////////////////////////////////////////
// processPath
string processPath(const char *path);
string processPath2(const char *path);

////////////////////////////////////////////////////////////////////////////////
// system() replacement
int spsSystem(const char *cmd, const char *args = NULL);

/**
 * Adds an environment variable at the beggining of a command. Used on a system call when in need to previously define the environment variable. The original value is not destroyd; the needed value is added to the existing one.
 *@param: Variable containing the command
 *@param: Variable to be added
 *@param: Variable value
 *@retur : None
 */

void addEnvironmentVariable(string &command,
                            const char *variable,
                            const string &variableValue,
                            bool atBegin = true);

////////////////////////////////////////////////////////////////////////////////
// Copy directory recursively
int FileCopy(string &source, string &dest);
int FileCopy(const char *source, const char *dest);
int CopyDirectory(string source, string dest);
////////////////////////////////////////////////////////////////////////////////
// sleep is a system depedent method
#if !defined(__linux__) && !defined(__CYGWIN__)
#include <windows.h>
#define sleep(n) Sleep(1000 * n)
#endif

///////////////////////////////////////////////////////////////////////////////
#if defined(__linux__) || defined(__CYGWIN__)
#define MKDIR(__A) mkdir(__A, ACCESSPERMS)
#else
#define MKDIR(__A) mkdir(__A)
#endif
////////////////////////////////////////////////////////////////////////////////
//#if defined(__MINGW32__)
//#define abort() return -9
//#endif
////////////////////////////////////////////////////////////////////////////////
#if defined(__MINGW32__)
#include <direct.h>
#include <process.h>
#define GetCurrentDir _getcwd
#define spsGetPID _getpid
#else
#include <unistd.h>
#include <sys/types.h>
#include <iostream>
#define GetCurrentDir getcwd
#define spsGetPID getpid
#endif // defined(__MINGW32__) || defined(__CYGWIN__)
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////

