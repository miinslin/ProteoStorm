////////////////////////////////////////////////////////////////////////////////
// This file contains all functions that are, for some reason, system dependent
////////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <dirent.h>
#include <iostream>
#include <fstream>

#include "Specific.h"
#include "Logger.h"
#include "StatusFile.h"

////////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace specnets;
////////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// Catch signals, like segfault
//-----------------------------------------------------------------------------
#if defined(__linux__) || defined(__CYGWIN__)
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void segfault_sigaction(int signal, siginfo_t *si, void *arg)
{
  string statusFileName("status.txt");
  ERROR_MSG("Caught segfault at address " << si->si_addr);
  writeStatusFile(statusFileName, "Error");
  exit(0);
}
#endif
//-----------------------------------------------------------------------------
// The segfault divert function: to where the segfaults are diverted to
//-----------------------------------------------------------------------------
int addSegFaultDivert(void)
{
#if defined(__linux__) || defined(__CYGWIN__)
  // catch signals (redirect)
  struct sigaction sa;

  memset(&sa, 0, sizeof(struct sigaction));
  sigemptyset(&sa.sa_mask);
  sa.sa_sigaction = segfault_sigaction;
  sa.sa_flags   = SA_SIGINFO;

  sigaction(SIGSEGV, &sa, NULL);
  return 1;
  // end catch signals (redirect)
#endif
  return 0;
}
////////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// cwait - waits for a forked process to finish
//-----------------------------------------------------------------------------
#if defined(__linux__) || defined(__CYGWIN__)
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif
//-----------------------------------------------------------------------------
#if defined(__linux__) || defined(__CYGWIN__)
struct cwait_t {

  int operator () ()
  {
    pid_t nw;
    int nstatus;

    do {
      nw = waitpid(-1, & nstatus, WUNTRACED | WCONTINUED);

      if (nw == -1)
        return nstatus;
      if (! WIFEXITED(nstatus))
        return nstatus;
    } while (! WIFEXITED(nstatus) && ! WIFSIGNALED(nstatus));

    return 0;
  }
} cwait;
#endif
////////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// Used to test if fork() methods are present.
//-----------------------------------------------------------------------------
bool forkable(void)
{
#if defined(__linux__) || defined(__CYGWIN__)
  return true;
#else
  return false;
#endif
}
//-----------------------------------------------------------------------------
// fork()
//-----------------------------------------------------------------------------
int myfork(void)
{
#if defined(__linux__) || defined(__CYGWIN__)
  return fork();
#endif
}
//-----------------------------------------------------------------------------
// wait for a process
//-----------------------------------------------------------------------------
void waitProcess(void)
{
#if defined(__linux__) || defined(__CYGWIN__)
  cwait();
#endif
}
////////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// processPath
//-----------------------------------------------------------------------------
string processPath(const char *path)
{
  string aux = path, aux2;
#if defined(__MINGW32__)
  for(int i = 0 ; i < aux.length() ; i++)
    if(aux[i] == '/')
      aux2 += "\\";
    else
      aux2 += aux[i];
#else
  aux2 = aux;
#endif
  return aux2;
}

string processPath2(const char *path)
{
  string aux = path, aux2;
#if defined(__MINGW32__)
  for(int i = 0 ; i < aux.length() ; i++)
    if(aux[i] == '/')
      aux2 += "//";
    else
      aux2 += aux[i];
#else
  aux2 = aux;
#endif
  return aux2;
}
////////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// spsSystem - replace the system() call
//-----------------------------------------------------------------------------
int spsSystem(const char *cmd, const char *args)
{
  string aux;

#if defined(__MINGW32__)
  aux += '"';
#endif

  aux += processPath(cmd);
  if(args) {
    //string aa = processPath2(args);
    aux += ' '; aux += args;
  }

#if defined(__MINGW32__)
  aux += '"';
#endif

  //cout << aux << endl;
  return system(aux.c_str());
}

////////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// addEnvironmentVariable()
//-----------------------------------------------------------------------------
void addEnvironmentVariable(string &command,
                            const char *variable,
                            const string &variableValue,
                            bool atBegin)
{
#if defined(__linux__) || defined(__CYGWIN__)

  // get variable from system
  char *curEnvValue = getenv(variable);

  // Set the final variable value as the default variable value
  string finalEnvValue = variableValue;

  // check if the environment variable is already defined.
  if (curEnvValue) {
    // if it is, check if it contains the value we want.
    finalEnvValue = curEnvValue;
    // if the variable does not contained the value we need, add it.
    if (finalEnvValue.find(variableValue) == string::npos) {
      finalEnvValue += ':';
      finalEnvValue += variableValue;
    }
  }

  // Add to the command
  string finalCommand = variable;
  finalCommand += '=';

  if (atBegin) {
    finalCommand += finalEnvValue;
    finalCommand += ":$";
    finalCommand += variable;
  } else {
    finalCommand += '$';
    finalCommand += variable;
    finalCommand += ':';
    finalCommand += finalEnvValue;
  }
  //finalCommand += " ; export ";
  //finalCommand += variable;
  //finalCommand += " ; ";
  finalCommand += ' ';
  finalCommand += command;

  command = finalCommand;

#endif
}
////////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// setenv
//-----------------------------------------------------------------------------
void mysetenv(const char *name, const char *value)
{
#ifdef HAVE_SETENV
  setenv(name, value, 1);
#else
  int len = strlen(value)+1+strlen(value)+1;
  char *str = (char *)malloc(len + 10);
  if(str) {
    sprintf(str, "%s=%s", name, value);
    if(putenv(str) != EXIT_SUCCESS)
      ;
    //free(str); // memcheck says this is unneccessary
  }
#endif
}
////////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// Directory structure copy
//-----------------------------------------------------------------------------
#include <sys/stat.h>
#define BUFFER_SIZE 32768

int FileCopy(string &source, string &dest)
{
  return FileCopy(source.c_str(), dest.c_str());
}
//-----------------------------------------------------------------------------
int FileCopy(const char *source, const char *dest)
{
  char buf[BUFFER_SIZE];
  size_t size;

  FILE* src = fopen(source, "rb");
  FILE* dst = fopen(dest, "wb");

  //if(!dst)
  //  cout << "Failed to open dst file" << endl;

  //if(!src)
  //  cout << "Failed to open src file" << endl;

  if(!src || !dst)
    return 0;
  // clean and more secure
  // feof(FILE* stream) returns non-zero if the end of file indicator for stream is set

  while(size = fread(buf, 1, BUFFER_SIZE, src))
    fwrite(buf, 1, size, dst);

  fclose(src);
  fclose(dst);

  //ifstream in(source, ios_base::in | ios_base::binary);
  //ofstream out(dest, ios_base::out | ios_base::binary);
  //out << in.rdbuf();
  return 1;
}
//-----------------------------------------------------------------------------
int CopyDirectory(string source, string dest)
{
  DIR *dirp;
  unsigned char isFile =0x8;

  if((dirp = opendir(source.c_str())) == NULL) {
    cout << "Opendir( " << source << " returned NULL" << endl;
    return 0;
  }

  struct dirent *dptr;
  struct stat dir_stat;

  while(dptr = readdir(dirp)) {

    //cout << "Found " << dptr->d_name << endl;
    /* skip the "." and ".." entries, to avoid loops. */
    if (strcmp(dptr->d_name, ".") == 0)
      continue;
    if (strcmp(dptr->d_name, "..") == 0)
      continue;
    // compose full name
    string s1(source), d1(dest);
    s1 += '/'; s1 += dptr->d_name;
    d1 += '/'; d1 += dptr->d_name;
    // check if it is a directory
    if (stat(s1.c_str(), &dir_stat) == -1) {
      //cout << "stat failed" << endl;
      continue;
    }
    if (S_ISDIR(dir_stat.st_mode)) {
      //cout << "Copy directory: " << s1 << " -> " << d1 << endl;
      MKDIR(d1.c_str());
      CopyDirectory(s1, d1);
    } else {
      //cout << "Copy file: " << s1 << " -> " << d1 << endl;
      FileCopy(s1, d1);
    }
  }

  closedir(dirp);
  return 1;
}
////////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// set the stack size
//-----------------------------------------------------------------------------
#if defined(__linux__) || defined(__CYGWIN__)
#include <sys/resource.h>
#include <errno.h>
#endif
int setStackSize(int newSize)
{

#if defined(__linux__) || defined(__CYGWIN__)
  const rlim_t kStackSize = newSize * 1024 * 1024;   // min stack size in MB
  struct rlimit rl;
  int result;

  result = getrlimit(RLIMIT_STACK, &rl);
  DEBUG_MSG("Stack soft limit: " << (int)(rl.rlim_cur) );
  DEBUG_MSG("Stack hard limit: " << (int)(rl.rlim_max) );

  if (result == 0) {
    if (rl.rlim_cur < kStackSize)  {
      rl.rlim_cur = kStackSize;
      if(rl.rlim_cur > rl.rlim_max)
        rl.rlim_cur = rl.rlim_max;
      DEBUG_MSG("Changing stack soft limit to : " << rl.rlim_cur );
      result = setrlimit(RLIMIT_STACK, &rl);
      if (result != 0)  {
        int aa = errno;
        WARN_MSG("setrlimit returned result = " << result);
        WARN_MSG("errno = " << aa);
        string aux = strerror(aa);
        WARN_MSG("error = " << aux);
        return 0;
      }
      DEBUG_MSG("Stack soft limit changed to : " << kStackSize );
    }
  }

  // ...
#endif

  return 1;
}
////////////////////////////////////////////////////////////////////////////////
