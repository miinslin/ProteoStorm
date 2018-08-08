#ifndef UTILS_H
#define UTILS_H

#include "twovalues.h"
#include "Logger.h"

#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <set>
#include <map>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>
#include <time.h>

using namespace std;

//delimiter for .csv files
extern const char* CSV_SEP;

extern const float POS_INF;


////////////////////////////////////////////////////////////////////////////////
// Class to hold a file name, and split its components
////////////////////////////////////////////////////////////////////////////////
 /*! \brief FilenameManager class

   Provides funtionalities assotiated to filename management.

   */
class FilenameManager {

 public:

  /*! \brief full file name
       The initial filename string
   */
  string filenameFull;

  /*! \brief path
  The path section, after filename splitting.
   */
  string path;

  /*! \brief the file name (excluding extension)
    The filename section, after filename splitting, excluding the extension.
   */
  string filename;

  /*! \brief the filename extension
    The filename extension
   */
  string extension;


  //! \name CONSTRUCTORS
  //@{
  /*! \brief The default constructor
   */
  FilenameManager() {};

  /*! \brief The constructor
   @param the filename
   */
  FilenameManager(const char *str) {initialize(str);};

  /*! \brief The constructor
   @param the filename
   */
  FilenameManager(const string &str) {initialize(str.c_str());};
  //@}

  /*! \brief Initialization method
    The initialization method splits the filename into it's components and coverts them to lower case.
   @param the filename
   */
  void initialize(const char *str) {filenameFull=str;splitFilename();};

  // split the full filename into its components
  /*! \brief splitFilename
   Splits the filename into it's components.
   */
  void splitFilename(void);

  /*! \brief joinFilename
   builds the full filename from its components.
   */
  void joinFilename(void);

  /*! \brief lowerCaseExtension
   Convert the extension to lower case.
   */
  void lowerCaseExtension(void);

  /*! \brief checkExtension
    checks if extension is known, by comparing agains a list of hardcoded known extensions.
   @param includeExternal Also checks agains the extensions used by external tools (e.x. proteowizard)
   @return Returns true if the extension is known.
   */
  bool checkExtension(bool includeExternal = true);

  // get filename with extension
  inline string getFilenameWithExtension(void)
  {
    string str(filename);
    str += ".";
    str += extension;
    return str;
  }

};

int removeFolderFiles(const string &dir, const vector<string> &files);

int removeFolder(const string &dir);

////////////////////////////////////////////////////////////////////////////////
  /*! \brief directoryContents
    Gets the files names in a given directory that have the given preffix and suffix. Files list may be sorted.
   @param dir Directory to list
   @param prefix file prefix
   @param sufix file suffix
   @param outputFiles list of files found
   @param sort if true, the file list is sorted
   @return Returns true if there were no errors
   */
bool directoryContents(const string & dir,
                                 const string & prefix,
                                 const string &sufix,
                                 vector<string> &outputFiles,
                                 bool sort);

////////////////////////////////////////////////////////////////////////////////
  /*! \brief directoryContents
    Gets the files names in a given directory that have the given preffix and suffix. Files list may be sorted.
   @param dir Directory to list
   @param prefix file prefix
   @param sufix file suffix
   @param sort if true, the file list is sorted
   @return outputFiles list of files found
   */
vector<string> directoryContents(const string & dir,
                                 const string & prefix,
                                 const string &sufix,
                                 bool sort);



/**
 * Compares two float numbers.
 *@Param: first string
 *@Param: second string
 *@Param: sensitivity
 *@return: Boolean result
 */
bool floatCompare(float A, float B, int maxUlps);

/**
 * Compares two strings by size, then alphabetically.
 *@Param: first string
 *@Param: second string
 *@return: Boolean result
 */
bool compare_nocase(string first, string second);

string intToString(int i);
/**
 * string right trim. Aplies to string objects.
 */
#define rtrim(_STR) _STR.erase(_STR.find_last_not_of(string(" \t\f\v\n\r"))+1)
#define ltrim(_STR) _STR.erase(0, _STR.find_first_not_of(string(" \t\f\v\n\r"))+1)

/**
 * Writes a file index into a file
 *@param: file name
 *@param: strings list (file index)
 *@retur : return status: true = ok ; false = error
 */
bool writeFileIndex(const string &fileName, vector<string> &index);

/**
 * Reads a text file content and returns it's contents in lines
 Reads a file index stored in a text file and return it as a vector of strings.
 *@param: file name
 *@param: returned set of sub strings
 *@return: return status: true = ok ; false = error
 */
bool readFilesFromFile(const string& str,
                       vector<string>& subStrings,
                       bool fromCurrentDir = false);
/**
 * Reads a file into a buffer
 *@param: file name
 *@param: file size
 *@return: File buffer; null if error
 */
char *readFile(const string &iFileName, unsigned &length);
/**
 * Splits a string into a vector of strings, not discarding empty string
 *@param: string to split
 *@param: returned set of sub strings
 *@param: delimiters; defaultings to space
 */
void stringSplit2(string str, vector<string> &results, const string &delim =
    " ");

void stringSplit2(string str, list<string> &results, const string &delim =
    " ");

void stringSplit2(string str, vector<char *> &results, const string &delim =
    " ");

    
/**
 * Splits a string into a vector of strings, discarding empty strings
 *@param: string to split
 *@param: returned set of sub strings
 *@param: delimiters; defaultings to space
 */
void stringSplit(const string& str,
                 vector<string>& subStrings,
                 const string& delimiters = " ");
/**
 * Joins a vector of strings into a single string
 *@param: returned string
 *@param:  set of sub strings
 *@param: delimiter; defaults to space
 */
void stringJoin(string &outputStr,
                const vector<string> &subStrings,
                const string & delimiters = " ");

/**
 * Generates the list of prefixes of a given string in order of increasing length
 * @param str
 * @param outputPrefixes
 */
void getPrefixes(const string& str, list<string>& outputPrefixes);

/**
 * Generates the list of suffixes of a given string in order of decreasing length
 * @param str
 * @param outputSuffixes
 */
void getSuffixes(const string& str, list<string>& outputSuffixes);

/**
 * Returns true if shortStr is a case-sensitive prefix of longStr, false is not
 * @param shortStr
 * @param longStr
 * @return true
 */
bool isPrefix(const string& shortStr, const string& longStr);

/**
 * Reverses a string
 * @param str
 * @return reversed string
 */
string reverseString(const string& str);

/**
 * Writes a vector of strings to a file in binary format
 * @param fp file ptr
 * @param strings vector of strings
 * @return true if written successfully, false otherwise
 */
bool writeStringsToBinaryStream(FILE* fp, vector<string>& strings);

template<class T> bool writeStringMapToBinaryStream(FILE* fp,
                                                    map<string, T>& stringMap)
{

  vector<string> strings(stringMap.size());
  T* values = (unsigned short*)malloc(sizeof(T) * stringMap.size());

  unsigned int idxUse = 0;
  for (typename map<string, T>::iterator mapIt = stringMap.begin(); mapIt
      != stringMap.end(); mapIt++)
  {
    strings[idxUse] = mapIt->first;
    values[idxUse] = mapIt->second;
    idxUse++;
  }

  if (!writeStringsToBinaryStream(fp, strings))
  {
    free(values);
    return false;
  }

  if (strings.size() == 0)
  {
    return true;
  }

  unsigned int count = fwrite(values, sizeof(T), idxUse, fp);

  free(values);
  return (bool)count;

  //if (count == 0)
  //{
  //  free(values);
  //  return false;
  //}
  //free(values);
  //return true;
}

/**
 * Reads a vector of strings from a file in binary format
 * @param fp file ptr
 * @param strings vector of strings
 * @return true if read successfully, false otherwise
 */
bool readStringsFromBinaryStream(FILE* fp, vector<string>& strings);

template<class T> bool readStringMapFromBinaryStream(FILE* fp,
                                                     map<string, T>& stringMap)
{

  vector<string> strings;
  stringMap.clear();

  if (!readStringsFromBinaryStream(fp, strings))
  {
    return false;
  }

  if (strings.size() == 0)
  {
    return true;
  }

  T* values = (unsigned short*)malloc(sizeof(T) * strings.size());

  unsigned int count = fread(values, sizeof(T), strings.size(), fp);
  if (count == 0)
  {
    free(values);
    return false;
  }

  for (unsigned int i = 0; i < strings.size(); i++)
  {
    stringMap[strings[i]] = values[i];
  }
  free(values);
  return true;
}

/**
 * Gets the filename from a full path
 *@param outputFilename - Input path string
 *@param includeExt - include the file extension
 *@return: string containing file name
 */
void extractFileName(const string & path,
                     string &outputFilename,
                     bool includeExt = true);

/**
 * Copies a file to another location
 * @param src Source filename
 * @param dest destination filename
 * @return Returns true, if the file is copied successfully.  Returns
 * false if the src file does not exists
 */
bool copyFile(string src, string dest);

/**
 * Concatenates 2 files and writes them to a new file
 * @param f1 Filename 1
 * @param f2 Filename 2
 * @param dest The file to create from the concatenation of f1 and f2
 * @return Returns true if concatenation succeeds, false otherwise
 */
bool concatenateFiles(string f1, string f2, string dest);

/**
 * Gets file extension from full path
 * @param filename - Input filename
 * @param fileExt - output file extension
 *
 */
void extractFileExt(const string & filename, string & fileExt);

/**
 * Replace occurences in a string
 *@param string context (the string to search. This string WILL be changed)
 *@param string from (what to change)
 *@param string to (replace pattern)
 *@return string: a copy of context string
 */
string& replaceAll(string& context, const string& from, const string& to);

string makeBracketsMods(string& peptide);

/**
 * Makes a directory/folder if it can be made (path to folder is valid and it does  not exist)
 *@param path folder path to create
 *@return true is folder was created, false if not
 */
bool mkdir_if_not_exist(const char* path);

/**
 * Recursively generate directory/folder and all subdirectories if they can be made.
 *@param path folder path to create
 *@return true is folder was created or already exists, false if not
 */
bool mkdirRecurse(const string &path);

/**
 * Maps command-line arguments to their appropriate flags
 *@param argv array of command-line arguments
 *@param argc number of arguments in argv
 *@param flags argument flags to look for
 *@param parsedArgs where to map found flags to their value
 *@return
 */
void parseArguments(char *argv[], int argc, set<string>& flags, map<string,
    string>& parsedArgs);

void
getHistogramInfo(vector<float>& data, vector<TwoValues<float> >& binFreq);

float getResolution(float peakTol);

/**
 * Parses a c string to a float
 *@param str char*
 *@return float value
 */
float getFloat(const char* str);

/**
 * Rounds a decimal number for a # of decimal places
 *@param double value to round
 *@param int number of decimal places
 *@return double rounded value
 */
double roundDouble(double x, int places);

/**
 * Rounds and concatenates a float to an int
 *@param float x
 *@return int value
 */
int floatToInt(float x);

/**
 * Rounds and concatenates a double to an int
 *@param double x
 *@return int value
 */
int doubleToInt(double x);

/**
 * Parses a c string to an int
 *@param str char*
 *@return int value
 */
int getInt(const char* str);

/**
 * Depreciated. See MZRange::EqualWithinRange
 */
bool isEqual(float f1, float f2, float range);

/**
 * Parses an int to a string
 *@param x integer
 *@param equalizeLength pad zeros to beginning of string
 *  if its length is less than this
 *@return string
 */
string parseInt(int x, int equalizeLength = 0);

/**
 * Parses a float, which is rounded, to a string
 *@param x float
 *@param prec number of decimal points to round to
 *@return string
 */
string parseFloat(float x, int prec);

/**
 * Parses a double, which is rounded, to a string
 *@param x double
 *@param prec number of decimal points to round to
 *@return string
 */
string parseDouble(double x, int prec);

/**
 * Parses a double to a string while displaying exponent
 *@param x double
 *@param prec number of decimal points to round to
 *@return string
 */
string parseDoubleSci(double x, int prec);

/**
 * Splits text by delim to a vector of string
 *@param text char*
 *@param vec output string vector
 *@param delim char*
 */
bool splitText(const char* text, list<string>& vec, const char* delim);

/**
 * Splits string by delim to a vector of strings
 *@param text string
 *@param vec output string vector
 *@param delim string
 */
bool splitText(const char* text, vector<string>& vec, const char* delim);

/**
 * Converts text containing numbers delimeted by delim to a vector of floats
 *@param text char*
 *@param vec output string list
 *@param delim char*
 */
//bool readTextToVector(const char* text, list<float>& vec, const char* delim);

/**
 * From a string of positive integers and integer ranges (ex. "-3,5,7-9"), retrieves all
 *   numbers specified (ex. a set of 0,1,2,3,5,7,8,9).
 *@param text input string of integers and/or integer ranges separated by commas
 *@param nums output set of all numbers specifed by text
 *@return true if parsing suceeded, false if not
 */
bool getRanges(const char* text, set<short>& nums);

bool getdir(string dir, list<string> &files);

class ProgressDisplay
{
public:
  long end;
  int len;
  ProgressDisplay(ostream& output, long e, const char* msg = "Progress: % ");
  void showProgress(ostream& output, long pos);
  void clear(ostream& output);
};

/**
 * TODO: add description
 */
class Utils
{
public:

  /**
   * TODO: add description
   */
// const double FLOAT_ERR = 0.0001;
  static double FLOAT_ERR;
  /**
   * TODO: add description
   */
  static vector<float> &unique(vector<float> &v, float resolution, vector<
      unsigned int> *idx = 0);

  /**
   * Intersects two sorted sets of masses.
   *
   *@param v1
   *@param v2
   *@param putHere
   *@param tolerance
   *@return
   */
  static int intersect(vector<float> &v1,
                       vector<float> &v2,
                       vector<float> &putHere,
                       float tolerance);

  /**
   * Saves a set of pairs of floats to binary file.
   * Format: 1 int with number of pairs followed by numEntries pairs of floats.
   *
   *@param filename
   *@param data
   *@return
   */
  static int save_tcb(char *filename, vector<TwoValues<float> > data);

  /**
   * TODO: add description
   *
   *@param filename
   *@param data
   *@return
   */
  static int save_tcb_list(char *filename, list<TwoValues<float> > data);

  /**
   * Converts a set of values to an indicator boolean sps::vector with true
   * whenever a value is present in values.
   *
   *@param values
   *@param output
   *@param tolerance
   *@param multFactor
   *@return
   */
  static unsigned int list_to_vector_bool(vector<float> &values,
                                          vector<bool> &output,
                                          float tolerance,
                                          float multFactor);

  /**
   * Computes a table of binomial probabilities for n=1..maxN, k=1..n;
   * probs[i][j]=p(k=j|n=j+1).
   *
   *@param maxN
   *@param p
   *@param probs
   */
  static void binomial(unsigned short maxN,
                       float p,
                       vector<vector<float> > &probs);

  /**
   * Computes the Gaussian P(X<=x).
   *
   *@param x
   *@param mean
   *@param stddev
   *@return
   */
  static double gaussiancdf(double x, double mean = 0.0, double stddev = 1.0);

  /**
   * Computes scores[i][j]=log(probsSignal[i][j]/probsNoise[i][j]) for all
   * probsNoise[i][j]>0.
   *
   *@param probsSignal
   *@param probsNoise
   *@param scores
   */
  static void logscores(vector<vector<float> > &probsSignal, vector<vector<
      float> > &probsNoise, vector<vector<float> > &scores);

};

/**
 * Filter values from the input argument values - keep only those within
 * tolerance of some entry in validValues. NOTE: Assumes that both values and
 * validValues are sorted in increasing order.
 *
 *@param values
 *@param validValues
 *@param tolerance
 *@return
 */
template<class T> void filterValues(vector<T> &values,
                                    vector<T> &validValues,
                                    T &tolerance,
                                    vector<T> &finalValues)
{
  finalValues.resize(values.size());
  unsigned idxVal = 0, idxValid, idxF = 0;
  while (idxVal < values.size() and idxValid < validValues.size())
  {
    if (fabs(values[idxVal] - validValues[idxValid]) <= tolerance + 0.0001)
    {
      finalValues[idxF++] = values[idxVal++];
      continue;
    }
    if (values[idxVal] < validValues[idxValid])
      idxVal++;
    else
      idxValid++;
  }
  finalValues.resize(idxF);
}

/**
 * TODO: add description
 *
 *@param meansIndexFN
 *@param varsIndexFN
 *@param params
 *@return
 */
bool LoadGaussianParams(const char *meansIndexFN,
                        const char *varsIndexFN,
                        vector<TwoValues<float> > &params);

bool loadCsvToTxt(const char* infilecsv,
                  const char* outfiletxt,
                  float threshold = 0.1);
/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T> int Save_binArray(const char *filename,
                                    vector<vector<T> > &data)
{
  FILE *fp;
  unsigned int numLines = data.size(), numCols;

  if (numLines == 0)
    return -1;
  else
    numCols = data[0].size();

  T *streamData = (T *)malloc(sizeof(T) * numLines * numCols);
  unsigned int i, j;

  fp = fopen(filename, "wb");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numLines, sizeof(int), 1, fp); // Number of lines in data
  fwrite(&numCols, sizeof(int), 1, fp); // Number of cols in data

  for (i = 0; i < numLines; i++)
  {
    if (data[i].size() != numCols)
    {
      cerr << "ERROR in Save_binArray: Not enough data at position " << i
          << "!\n";
      exit(-1);
    }
    for (j = 0; j < numCols; j++)
      streamData[i * numCols + j] = data[i][j];
  }

  fwrite(streamData, sizeof(T), numLines * numCols, fp);
  free(streamData);

  fclose(fp);
  return 1;
}

/**
 * Used when needing to map each input spectrum to 2 final spectrum (when GenoMS
 * is being merged with CSPS).
 *
 *@param filename
 *@param data
 *@return
 */
template<class T> int Save_binArrayDouble(const char *filename,
                                    vector<vector<T> > &data)
{
  FILE *fp;
  unsigned int numLines = data.size()*2, numCols;

  if (numLines == 0)
    return -1;
  else
    numCols = data[0].size();

  T *streamData = (T *)malloc(sizeof(T) * numLines * numCols);
  unsigned int i, j;

  fp = fopen(filename, "wb");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numLines, sizeof(int), 1, fp); // Number of lines in data
  fwrite(&numCols, sizeof(int), 1, fp); // Number of cols in data

  for(int copy = 0; copy < 2; ++copy) {
    for (i = 0; i < numLines/2; i++)
      {
	if (data[i].size() != numCols)
	  {
	    cerr << "ERROR in Save_binArray: Not enough data at position " << i
		 << "!\n";
	    exit(-1);
	  }
	for (j = 0; j < numCols; j++)
	  streamData[(copy * numLines/2) + i * numCols + j] = data[i][j];
      }
  }
  cout << "Wrote " << numLines << " from " << data.size() << endl;
  fwrite(streamData, sizeof(T), numLines * numCols, fp);
  free(streamData);

  fclose(fp);
  return 1;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T> int Save_binArray(const char *filename, vector<T> &data)
{
  FILE *fp;
  unsigned int numLines = data.size(), numCols = 1;
  T *streamData = (T *)malloc(sizeof(T) * numLines);
  unsigned int i;

  fp = fopen(filename, "wb");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numLines, sizeof(int), 1, fp); // Number of lines in data
  fwrite(&numCols, sizeof(int), 1, fp); // Number of cols in data

  for (i = 0; i < numLines; i++)
    streamData[i] = data[i];

  fwrite(streamData, sizeof(T), numLines, fp);
  free(streamData);

  fclose(fp);
  return 1;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T> int Save_binArray(const char *filename,
                                    vector<TwoValues<T> > &data)
{
  FILE *fp;
  unsigned int numLines = data.size(), numCols = 2;
  T *streamData = (T *)malloc(sizeof(T) * numCols * numLines);
  unsigned int i;
  fp = fopen(filename, "wb");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numLines, sizeof(int), 1, fp); // Number of lines in data
  fwrite(&numCols, sizeof(int), 1, fp); // Number of cols in data

  for (i = 0; i < numLines; i++)
  {
    streamData[i << 1] = data[i][0];
    streamData[(i << 1) + 1] = data[i][1];
  }

  fwrite(streamData, sizeof(T), numLines * numCols, fp);
  free(streamData);

  fclose(fp);
  return 1;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T> int Save_binArray(const char *filename,
                                    list<TwoValues<T> > &data)
{
  FILE *fp;
  unsigned int numLines = data.size(), numCols = 2;
  T *streamData = (T *)malloc(sizeof(T) * numCols * numLines);

  fp = fopen(filename, "wb");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numLines, sizeof(int), 1, fp); // Number of lines in data
  fwrite(&numCols, sizeof(int), 1, fp); // Number of cols in data

  typename list<TwoValues<T> >::iterator iter;
  unsigned int i = 0;
  for (iter = data.begin(); iter != data.end(); iter++)
  {
    streamData[i++] = (*iter)[0];
    streamData[i++] = (*iter)[1];
  }

  fwrite(streamData, sizeof(T), numLines * numCols, fp);
  free(streamData);

  fclose(fp);
  return 1;
}

/**
 * Writes two binary arrays as if they were 1.  So if matrix 1 is aXb and matrix 2 is a'xb, then the
 * final array will be (a+a')xb
 */
template<class T> int Save_doubleBinArray(const char * filename, vector<vector<
    T> > & data1, vector<vector<T> > & data2)
{
  FILE * fp;

  unsigned int numLines = data1.size() + data2.size();
  if (numLines == 0)
    return -1;

  typename vector<vector<T> >::iterator iter1 = data1.begin();
  typename vector<vector<T> >::iterator iter2 = data2.begin();

  //num cols must be the same!!
  if (iter1->size() != iter2->size())
  {
    return -1;
  }

  unsigned int numCols = iter1->size();
  T *streamData = (T *)malloc(sizeof(T) * numLines * numCols);
  unsigned int i, j;

  //cout << "We have 2 arrays (" << data1.size() << "x" << numCols << ") and (" << data2.size() << "x" << numCols << ")" << endl;

  fp = fopen(filename, "wb");

  fwrite(&numLines, sizeof(int), 1, fp);
  fwrite(&numCols, sizeof(int), 1, fp);
  for (i = 0; iter1 != data1.end(); iter1++, i++)
    for (j = 0; j < numCols; j++)
      streamData[i * numCols + j] = (*iter1)[j];
  //cout << "i = " << i << endl;
  for (; iter2 != data2.end(); iter2++, i++)
    for (j = 0; j < numCols; j++)
      streamData[i * numCols + j] = (*iter2)[j];

  fwrite(streamData, sizeof(T), numLines * numCols, fp);
  free(streamData);

  fclose(fp);
  return 1;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T> int Save_binArray(const char *filename,
                                    list<vector<T> > &data)
{
  FILE *fp;
  unsigned int numLines = data.size(), numCols;
  typename list<vector<T> >::iterator iter = data.begin();
  if (numLines == 0)
    return -1;
  else
    numCols = iter->size();
  T *streamData = (T *)malloc(sizeof(T) * numLines * numCols);
  unsigned int i, j;

  fp = fopen(filename, "wb");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numLines, sizeof(int), 1, fp); // Number of lines in data
  fwrite(&numCols, sizeof(int), 1, fp); // Number of cols in data

  for (i = 0; iter != data.end(); iter++, i++)
    for (j = 0; j < numCols; j++)
      streamData[i * numCols + j] = (*iter)[j];

  fwrite(streamData, sizeof(T), numLines * numCols, fp);
  free(streamData);

  fclose(fp);
  return 1;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */

template<class T, template < typename U, typename = std::allocator<U> > class C> int
Load_binArray(const char *filename, C< C<T> > &data)
{
  FILE *fp; int numLines=0, numCols=0;
  unsigned int i,j; int n;

  fp = fopen(filename,"rb");
  if (fp==0)
  { std::cerr << "ERROR: cannot open " << filename << "\n"; return -1;}

  n=fread(&numLines,sizeof(int),1,fp); if(numLines<0 | n!=1)
  { std::cerr << "Invalid number of lines ("<<numLines<<") in "<<filename<<"\n";return -1;}
  n=fread(&numCols,sizeof(int),1,fp); if(numCols<0 | n!=1)
  { std::cerr << "Invalid number of columns ("<<numCols<<") in "<<filename<<"\n";return -1;}
  if(numLines==0 or numCols==0)
  { data.resize(0); return 0;}

  T *streamData = (T *)malloc(sizeof(T)*numLines*numCols);
  if(!streamData)
  { cerr<<"(Load_binArray) Error allocating memory to read "<<numLines<<"/"<<numCols<<" lines/cols of "<<sizeof(T)<<"-byte elements!\n"; return -1;}
  n=fread(streamData,sizeof(T),numLines*numCols,fp);
  fclose(fp);
  if(n!=(numLines*numCols))
  { std::cerr << "Error reading binary array: not enough elements in the file.\n"; free(streamData); return -1;}

  data.resize(numLines);
  if(data.size()!=numLines)
  { cerr<<"(Load_binArray) Error allocating memory to store "<<numLines<<"/"<<numCols<<" lines/cols of "<<sizeof(T)<<"-byte elements!\n"; free(streamData); return(-1);}
  unsigned int streamIdx=0;
  for(i=0;i<(unsigned int)numLines;i++)
  { data[i].resize(numCols);
    for(j=0;j<(unsigned int)numCols;j++) data[i][j] = streamData[streamIdx++];
  }

  free(streamData);
  return 1;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T, template < typename U, typename = std::allocator<U> > class C> int
Load_binArray(const char *filename, C<T> &data, bool outputErrors = true)
{
  FILE *fp; int numLines=0, numCols=0;
  unsigned int i,j; int n;

  fp = fopen(filename,"rb");
  if (fp==0)
  { if(outputErrors) std::cerr << "ERROR: cannot open " << filename << "\n"; return -1;}

  n=fread(&numLines,sizeof(int),1,fp); if(numLines<0 | n!=1)
  { if(outputErrors) std::cerr << "Invalid number of lines ("<<numLines<<") in "<<filename<<"\n";return -1;}
  n=fread(&numCols,sizeof(int),1,fp); if(numCols!=1 | n!=1)
  { if(outputErrors) cerr << "Invalid number of columns ("<<numCols<<") in "<<filename<<"\n";return -1;}
  if(numLines==0)
  { data.resize(0); return 0;}

  T *streamData = (T *)malloc(sizeof(T)*numLines);
  if(!streamData)
  { if(outputErrors) cerr<<"(Load_binArray) Error allocating memory to read "<<numLines<<" lines of "<<sizeof(T)<<"-byte elements!\n"; return -1;}
  n=fread(streamData,sizeof(T),numLines,fp);
  fclose(fp);
  if(n!=(numLines))
  { if(outputErrors) std::cerr << "Error reading binary array: not enough elements in the file.\n"; free(streamData); return -1;}

  data.resize(numLines);
  if(data.size()!=numLines)
  { if(outputErrors) cerr<<"(Load_binArray) Error allocating memory to store "<<numLines<<" lines of "<<sizeof(T)<<"-byte elements!\n"; free(streamData); return -1;}
  for(i=0;i<(unsigned int)numLines;i++) data[i]=streamData[i];

  free(streamData);
  return 1;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T> int Load_binArray(const char *filename,
                                    vector<TwoValues<T> > &data)
{

  FILE *fp;
  int numLines = 0, numCols = 0, n = 0;

  fp = fopen(filename, "rb");
  if (fp == 0)
  {
    std::cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  n = fread(&numLines, sizeof(int), 1, fp);
  if (numLines < 1)
  {
    std::cerr << "Invalid number of lines (" << numLines << ") in " << filename
        << "\n";
    return -1;
  }
  n = fread(&numCols, sizeof(int), 1, fp);
  if (numCols != 2)
  {
    cerr << "Invalid number of columns (" << numCols << ") in " << filename
        << "\n";
    return -1;
  }
  if (numLines == 0)
  {
    data.resize(0);
    return 0;
  }

  data.resize(numLines);
  if (data.size() != numLines)
  {
    cerr << "(Load_binArray) Error allocating memory to store " << numLines
        << " lines of " << sizeof(T) << "-byte elements!\n";
    return -1;
  }

  char * streamData = (char*)malloc(numLines * numCols * sizeof(T));
  if (streamData == 0x0)
  {
    cerr << "(Load_binArray) Error allocating memory to store " << numLines
        << " lines of " << sizeof(T) << "-byte elements!\n";
    return -1;
  }
  n = fread(streamData, sizeof(T), numLines * numCols, fp);
  fclose(fp);

  unsigned int streamIdx = 0;
  for (unsigned int i = 0; i < numLines; i++)
  {
    data[i][0] = *((T*)&streamData[streamIdx]);
    streamIdx += sizeof(T);
    data[i][1] = *((T*)&streamData[streamIdx]);
    streamIdx += sizeof(T);
  }
  free(streamData);
  return 1;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param numLinesCols
 *@return
 */
inline int Load_binArraySize(const char *filename,
                             TwoValues<unsigned int> &numLinesCols)
{
  FILE *fp;
  numLinesCols.set(0, 0);
  fp = fopen(filename, "rb");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }
  TwoValues<int> sz(0, 0);

  unsigned int n;
  n = fread(&sz[0], sizeof(int), 1, fp);
  if ((sz[0] < 0) | (n != 1))
  {
    cerr << "Invalid number of lines (" << sz[0] << ") in " << filename << "\n";
    return -1;
  }
  n = fread(&sz[1], sizeof(int), 1, fp);
  if ((sz[1] < 0) | (n != 1))
  {
    cerr << "Invalid number of columns (" << sz[1] << ") in " << filename
        << "\n";
    return -1;
  }
  numLinesCols.set((unsigned int)sz[0], (unsigned int)sz[1]);

  fclose(fp);
  return numLinesCols[0];
}
//template<class T1,class T2,class T3> int Save_binListArray(char *filename, vector<T2> &data);

// T1 = element type, T2 = vector<T1>/list<T1>, T3 = T2::iterator

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T1, class T2, class T3> int Load_binListArray(const char *filename,
                                                             vector<T2> &data)
{
  FILE *fp;
  unsigned int numLists = 0, numElems = 0;
  unsigned int *index;
  unsigned int i, streamIdx;
  int n;

  fp = fopen(filename, "rb");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  // Read number of lists
  n = fread(&numLists, sizeof(unsigned int), 1, fp);
  if (numLists <= 0 | n != 1)
  {
    cerr << "Invalid number of lists (" << numLists << ")\n";
    return -1;
  }

  // Read index
  index = (unsigned int*)malloc(sizeof(unsigned int) * numLists);
  n = fread(index, sizeof(unsigned int), numLists, fp);
  if (n != numLists)
  {
    cerr << "Invalid index!\n";
    free(index);
    return -1;
  }
  numElems = 0;
  for (i = 0; i < numLists; i++)
    numElems += index[i];

  // Read data
  T1 *streamData = (T1 *)malloc(sizeof(T1) * numElems);
  n = fread(streamData, sizeof(T1), numElems, fp);
  fclose(fp);
  if (n != (int)numElems)
  {
    cerr
        << "Error reading binary list array: not enough elements in the file.\n";
    free(streamData);
    free(index);
    return -1;
  }

  // Convert back to C++ structure
  data.resize(numLists);
  T3 iter;
  streamIdx = 0;
  for (i = 0; i < numLists; i++)
  {
    data[i].resize(index[i]);
    for (iter = data[i].begin(); iter != data[i].end(); iter++)
      (*iter) = streamData[streamIdx++];
  }

  free(index);
  free(streamData);
  return numLists;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T1, class T2, class T3> int Load_binListArray(const char *filename,
                                                             list<T2> &data)
{
  FILE *fp;
  unsigned int numLists = 0, numElems = 0;
  unsigned int *index;
  int n;

  fp = fopen(filename, "rb");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  // Read number of lists
  n = fread(&numLists, sizeof(unsigned int), 1, fp);
  if (numLists <= 0 | n != 1)
  {
    cerr << "Invalid number of lists (" << numLists << ")\n";
    return -1;
  }

  // Read index
  index = (unsigned int*)malloc(sizeof(unsigned int) * numLists);
  n = fread(index, sizeof(unsigned int), numLists, fp);
  if (n != numLists)
  {
    cerr << "Invalid index!\n";
    free(index);
    return -1;
  }
  numElems = 0;
  for (unsigned int i = 0; i < numLists; i++)
    numElems += index[i];

  // Read data
  T1 *streamData = (T1 *)malloc(sizeof(T1) * numElems);
  n = fread(streamData, sizeof(T1), numElems, fp);
  fclose(fp);
  if (n != (int)numElems)
  {
    cerr
        << "Error reading binary list array: not enough elements in the file.\n";
    free(streamData);
    free(index);
    return -1;
  }

  // Convert back to C++ structure
  data.resize(numLists);
  T3 iter;
  unsigned int i = 0;
  for (iter = data.begin(); iter != data.end(); iter++)
  {
    iter->resize(index[i]);
    for (unsigned int j = 0; j < index[i]; j++)
    {
      (*iter)[j] = streamData[j];
    }
    i++;
  }

  free(index);
  free(streamData);
  return numLists;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T1, class T1b, class T2, class T3> int Load_binListArray(const char *filename,
                                                                        vector<
                                                                            T2> &data)
{
  FILE *fp;
  unsigned int numLists = 0, numElems = 0;
  unsigned int *index;
  unsigned int i, streamIdx;
  int n;

  fp = fopen(filename, "rb");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  // Read number of lists
  n = fread(&numLists, sizeof(unsigned int), 1, fp);
  if (numLists <= 0 | n != 1)
  {
    cerr << "Invalid number of lists (" << numLists << ")\n";
    return -1;
  }

  // Read index
  index = (unsigned int*)malloc(sizeof(unsigned int) * numLists);
  n = fread(index, sizeof(unsigned int), numLists, fp);
  if (n != numLists)
  {
    cerr << "Invalid index!\n";
    free(index);
    return -1;
  }
  numElems = 0;
  for (i = 0; i < numLists; i++)
  {
    numElems += index[i];
  }

  // Read data
  char *streamData = (char *)malloc((sizeof(T1) + sizeof(T1b)) * numElems);
  n = fread(streamData, sizeof(T1) + sizeof(T1b), numElems, fp);
  fclose(fp);
  if (n != (int)numElems)
  {
    cerr
        << "Error reading binary list array: not enough elements in the file.\n";
    free(streamData);
    free(index);
    return -1;
  }

  // Convert back to C++ structure
  data.resize(numLists);
  T3 iter;
  streamIdx = 0;
  for (i = 0; i < numLists; i++)
  {
    data[i].resize(index[i]);
    for (iter = data[i].begin(); iter != data[i].end(); iter++)
    {
      iter->first = *((T1*)&streamData[streamIdx]);
      streamIdx += sizeof(T1);
      iter->second = *((T1b*)&streamData[streamIdx]);
      streamIdx += sizeof(T1b);
    }
  }

  free(index);
  free(streamData);
  return numLists;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T1, class T2, class T3> int Save_binListArray(const char *filename,
                                                             vector<T2> &data)
{
  unsigned int numEntries1 = data.size(), numEntries2 = 0;
  unsigned int *index = (unsigned int*)malloc(sizeof(unsigned int)
      * numEntries1);
  for (unsigned int i = 0; i < numEntries1; i++)
  {
    index[i] = data[i].size();
    numEntries2 += index[i];
  }
  T1 *streamData = (T1*)malloc(sizeof(T1) * numEntries2);

  int dataArrayIdx = 0;
  T3 iter;
  for (unsigned int i = 0; i < numEntries1; i++)
    for (iter = data[i].begin(); iter != data[i].end(); iter++)
      streamData[dataArrayIdx++] = *iter;

  FILE *fp = fopen(filename, "wb");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numEntries1, sizeof(int), 1, fp); // Number of lists in data
  fwrite(index, sizeof(int), numEntries1, fp); // List sizes
  fwrite(streamData, sizeof(T1), numEntries2, fp); // List elements

  fclose(fp);
  free(index);
  free(streamData);
  return 1;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T1, class T2, class T3> int Save_binListArray(const char *filename,
                                                             list<T2> &data)
{
  unsigned int numEntries1 = data.size(), numEntries2 = 0;
  unsigned int *index = (unsigned int*)malloc(sizeof(unsigned int)
      * numEntries1);
  typename list<T2>::iterator dataIter;
  //	for(unsigned int i=0;i<numEntries1;i++) { index[i]=data[i].size(); numEntries2+=index[i]; }
  unsigned int indexIdx = 0;
  for (dataIter = data.begin(); dataIter != data.end(); dataIter++, indexIdx++)
  {
    index[indexIdx] = dataIter->size();
    numEntries2 += index[indexIdx];
  }
  T1 *streamData = (T1*)malloc(sizeof(T1) * numEntries2);

  int dataArrayIdx = 0;
  T3 iter;
  for (dataIter = data.begin(); dataIter != data.end(); dataIter++)
    for (iter = dataIter->begin(); iter != dataIter->end(); iter++)
      streamData[dataArrayIdx++] = *iter;

  FILE *fp = fopen(filename, "wb");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numEntries1, sizeof(int), 1, fp); // Number of lists in data
  fwrite(index, sizeof(int), numEntries1, fp); // List sizes
  fwrite(streamData, sizeof(T1), numEntries2, fp); // List elements

  fclose(fp);
  free(index);
  free(streamData);
  return 1;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T1, class T1b, class T2, class T3> int Save_binListArray(const char *filename,
                                                                        vector<
                                                                            list<
                                                                                pair<
                                                                                    T1,
                                                                                    T1b> > > &data)
{
  unsigned int numEntries1 = data.size(), numEntries2 = 0;
  unsigned int *index = (unsigned int*)malloc(sizeof(unsigned int)
      * numEntries1);
  for (unsigned int i = 0; i < numEntries1; i++)
  {
    index[i] = data[i].size();
    numEntries2 += index[i];
  }
  char * streamData = (char*)malloc((sizeof(T1) + sizeof(T1b)) * numEntries2);

  int dataArrayIdx = 0;
  T3 iter;
  for (unsigned int i = 0; i < numEntries1; i++)
  {
    for (iter = data[i].begin(); iter != data[i].end(); iter++)
    {
      strncpy(&streamData[dataArrayIdx], (char*)&iter->first, 4);
      dataArrayIdx += sizeof(T1);
      strncpy(&streamData[dataArrayIdx], (char*)&iter->second, 4);
      dataArrayIdx += sizeof(T1b);
    }
  }

  FILE *fp = fopen(filename, "wb");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numEntries1, sizeof(int), 1, fp); // Number of lists in data
  fwrite(index, sizeof(int), numEntries1, fp); // List sizes
  fwrite(streamData, sizeof(T1) + sizeof(T1b), numEntries2, fp); // List elements

  fclose(fp);
  free(index);
  free(streamData);
  return 1;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T1, class T1b, class T2, class T3> int Save_binListArray(const char *filename,
                                                                        list<
                                                                            vector<
                                                                                pair<
                                                                                    T1,
                                                                                    T1b> > > &data)
{
  unsigned int numEntries1 = data.size(), numEntries2 = 0;
  unsigned int *index = (unsigned int*)malloc(sizeof(unsigned int)
      * numEntries1);
  typename list<T2>::iterator dataIter;
  unsigned int indexIdx = 0;
  for (dataIter = data.begin(); dataIter != data.end(); dataIter++, indexIdx++)
  {
    index[indexIdx] = dataIter->size();
    numEntries2 += index[indexIdx];
  }

  char * streamData = (char*)malloc((sizeof(T1) + sizeof(T1b)) * numEntries2);

  int dataArrayIdx = 0;
  T3 iter;
  for (dataIter = data.begin(); dataIter != data.end(); dataIter++)
  {
    for (iter = dataIter->begin(); iter != dataIter->end(); iter++)
    {
      strncpy(&streamData[dataArrayIdx], (char*)&iter.first, 4);
      dataArrayIdx += sizeof(T1);
      strncpy(&streamData[dataArrayIdx], (char*)&iter.second, 4);
      dataArrayIdx += sizeof(T1b);
    }
  }

  FILE *fp = fopen(filename, "wb");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numEntries1, sizeof(int), 1, fp); // Number of lists in data
  fwrite(index, sizeof(int), numEntries1, fp); // List sizes
  fwrite(streamData, sizeof(T1) + sizeof(T1b), numEntries2, fp); // List elements

  fclose(fp);
  free(index);
  free(streamData);
  return 1;
}

/**
 * TODO: add description
 */
class BufferedLineReader
{
  char *lines;
  vector<unsigned int> linesIdx;
public:
  BufferedLineReader()
  {
    lines = NULL;
  }
  ~BufferedLineReader()
  {
    reset();
  }

  /**
   * TODO: add description
   */
  void reset()
  {
    if (lines)
      free(lines);
    lines = NULL;
    linesIdx.resize(0);
  }

  /**
   * TODO: add description
   *
   *@return the number of lines currently in the reader.
   */
  unsigned int size()
  {
    return linesIdx.size();
  }

  /**
   * TODO: add description
   *
   *@param lineNum the line number to get.
   *@return
   */
  char *getline(unsigned int lineNum)
  {
    if (lineNum >= linesIdx.size())
      return NULL;
    else
      return &lines[linesIdx[lineNum]];
  }

  /**
   * TODO: add description
   *
   *@param filename
   *@return
   */
  short Load(const char *filename);
};

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T> bool Load_binArray_multiple(const char *filename, vector<
    vector<T> > &data)
{
  BufferedLineReader blrNames;
  unsigned int fIdx, numFiles;
  data.resize(0);
  if (blrNames.Load(filename) <= 0)
  {
    cerr << "ERROR reading " << filename << "!\n";
    return false;
  }
  numFiles = blrNames.size();

  // Count total number of entries
  TwoValues<unsigned int> totLinesCols(0, 0), curLinesCols(0, 0);
  for (fIdx = 0; fIdx < numFiles; fIdx++)
  {
    if (strlen(blrNames.getline(fIdx)) == 0)
      continue;
    const char* fName = blrNames.getline(fIdx);
    if (Load_binArraySize(fName, curLinesCols) < 0)
      return false;
    if (fIdx == 0)
      totLinesCols[1] = curLinesCols[1];
    else if (totLinesCols[1] != curLinesCols[1])
    {
      cerr << "ERROR: Inconsistent number of columns in "
          << blrNames.getline(fIdx) << " (Load_binArray_multiple)\n";
      return false;
    }
    totLinesCols[0] += curLinesCols[0];
  }

  // Load data
  FILE *fp;
  vector<vector<T> > curData;
  data.resize(totLinesCols[0]);
  unsigned int pivot, idxData = 0;
  for (fIdx = 0; fIdx < numFiles; fIdx++)
  {
    if (strlen(blrNames.getline(fIdx)) == 0)
      continue;
    const char* fName = blrNames.getline(fIdx);
    if (Load_binArray(fName, curData) < 0)
      return false;
    if (idxData + curData.size() > data.size())
    {
      cerr << "ERROR: Too many lines in " << blrNames.getline(fIdx)
          << " (Load_binArray_multiple)\n";
      return false;
    }
    for (pivot = 0; pivot < curData.size(); pivot++)
      data[idxData++] = curData[pivot];
  }

  return true;
}


/**
 * compareBinArray: compares two bin arrays and returns a bin array comparison structure
 *
 *@param binArray1
 *@param binArray2
 *@return number of different elements found
 */
template<class T> int Compare_binArray(vector<vector<T> > &a1, vector<vector<T> > &a2, vector<pair<int,int> > &res)
{
  int numRows1 = a1.size();
  int numRows2 = a2.size();

  // check bin array length
  if(numRows1 != numRows2) {
    res.push_back(pair<int,int>(numRows1, numRows2));
    return -1;
  }

  // check for empty arrays.
  if(numRows1 == 0)
    return -2;

  int numCols1 = a1[0].size();
  int numCols2 = a2[0].size();

  // check bin array width
  if(numCols1 != numCols2) {
    res.push_back(pair<int,int>(numCols1, numCols2));
    return -3;
  }

  // check for empty rows.
  if(numCols1 == 0)
    return -4;

  // cycle thru the rows
  for(int i = 0 ; i < numRows1 ; i++)

    // cycle thru the columns
    for(int j = 0 ; j < numCols1 ; j++)

      // compare elements
      if(a1[i][j] != a2[i][j])

        // set result element if array elements are different
        res.push_back(pair<int,int>(i,j));

  return res.size();
}


std::string get_time_string();

/**
 * trim: trims string of spaces and tabs from both ends
 *
 *@param str
 *@return trimmed string
 */
std::string trim(std::string str);


string string_replace_all(std::string s,
                      std::string toReplace,
                      std::string replaceWith);

#endif
