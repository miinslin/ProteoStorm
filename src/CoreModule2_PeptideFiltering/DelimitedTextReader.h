/*
 * delimitedtextreader.h
 *
 *  Created on: Jan 2, 2011
 *      Author: jsnedecor
 */

#ifndef DELIMITEDTEXTREADER_H_
#define DELIMITEDTEXTREADER_H_

//System includes
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

//External includes
#include "Logger.h"
#include "utils.h"

namespace specnets
{


  void insertValue(vector<string> & fields, vector<vector<string> > &lines);
  void insertValue(vector<string> & fields, vector<vector<float> > &lines);
  void insertValue(vector<string> & fields, vector<vector<int> > &lines);
  void insertValue(vector<string> & fields, vector<vector<double> > &lines);

  void insertValue(list<string> & fields, list<list<string> > &lines);
  void insertValue(list<string> & fields, list<list<float> > &lines);
  void insertValue(list<string> & fields, list<list<int> > &lines);
  void insertValue(list<string> & fields, list<list<double> > &lines);


  class DelimitedTextReader
  {
  public:
    /*! \brief Load delimited text file with no header.
     *
     * @param filename - Path to delimited text file.
     * @param delim - String delimiter
     * @param delim - Comment delimiter. Assume this is at the beginning of the string
     * @param lines - vector of vectors to contain delimited text.
     */


    static bool loadDelimitedFileNoHeader(const char * filename,
                                   const string &delim,
                                   const string &commentDelim,
                                   vector<vector<string> > &lines);

    static bool loadDelimitedFileNoHeader(const char * filename,
                                   const string &delim,
                                   const string &commentDelim,
                                   vector<vector<float> > &lines);

    static bool loadDelimitedFileNoHeader(const char * filename,
                                   const string &delim,
                                   const string &commentDelim,
                                   vector<vector<double> > &lines);

    static bool loadDelimitedFileNoHeader(const char * filename,
                                   const string &delim,
                                   const string &commentDelim,
                                   vector<vector<int> > &lines);

    static bool loadDelimitedFileNoHeader(ifstream &filename,
                                   const string &delim,
                                   const string &commentDelim,
                                   vector<vector<int> > &lines);

    static bool loadDelimitedFileNoHeader(ifstream &filename,
                                   const string &delim,
                                   const string &commentDelim,
                                   vector<vector<string> > &lines);

    static bool loadDelimitedFileNoHeader(ifstream &filename,
                                   const string &delim,
                                   const string &commentDelim,
                                   vector<vector<double> > &lines);

    static bool loadDelimitedFileNoHeader(ifstream &filename,
                                   const string &delim,
                                   const string &commentDelim,
                                   vector<vector<float> > &lines);



    /*! \brief Load delimited text file.
     *
     * @param filename - Path to delimited text file.
     * @param delim - String delimiter
     * @param delim - Comment delimiter. Assume this is at the beginning of the string
     * @param lines - vector of vectors to contain delimited text.
     * @param header - vector to contain header information
     * @param requiredHeader - vector containing required header names
     * @param requiredIndices - indices to the required header.
     */

     static bool loadDelimitedFile(const char * filename,
                           const string &delim,
                           const string &commentDelim,
                           map<string, unsigned int> &header,
                           vector<vector<string> > &lines,
                           vector<string> &requiredHeader,
                           vector<int> &requiredHeaderIndex);


    static bool loadDelimitedFile(const char * filename,
                           const string &delim,
                           const string &commentDelim,
                           map<string, unsigned int> &header,
                           vector<vector<float> > &lines,
                           vector<string> &requiredHeader,
                           vector<int> &requiredHeaderIndex);

    static bool loadDelimitedFile(const char * filename,
                           const string &delim,
                           const string &commentDelim,
                           map<string, unsigned int> &header,
                           vector<vector<double> > &lines,
                           vector<string> &requiredHeader,
                           vector<int> &requiredHeaderIndex);

    static bool loadDelimitedFile(const char * filename,
                           const string &delim,
                           const string &commentDelim,
                           map<string, unsigned int> &header,
                           vector<vector<int> > &lines,
                           vector<string> &requiredHeader,
                           vector<int> &requiredHeaderIndex);

    static bool loadDelimitedFile(const char * filename,
                           const string &delim,
                           const string &commentDelim,
                           vector<string> &header,
                           vector<vector<int> > &lines,
                           vector<string> &requiredHeader,
                           vector<int> &requiredHeaderIndex);

    static bool loadDelimitedFile(const char * filename,
                           const string &delim,
                           const string &commentDelim,
                           vector<string> &header,
                           vector<vector<float> > &lines,
                           vector<string> &requiredHeader,
                           vector<int> &requiredHeaderIndex);

    static bool loadDelimitedFile(const char * filename,
                           const string &delim,
                           const string &commentDelim,
                           vector<string> &header,
                           vector<vector<double> > &lines,
                           vector<string> &requiredHeader,
                           vector<int> &requiredHeaderIndex);

    static bool loadDelimitedFile(const char * filename,
                           const string &delim,
                           const string &commentDelim,
                           vector<string> &header,
                           vector<vector<string> > &lines,
                           vector<string> &requiredHeader,
                           vector<int> &requiredHeaderIndex);
    /*! \brief Use getline and clear off carriage returns.
     *
     * Note: This situation only occurs if DOS files are being read in on Linux,
     * ordinarily getline correctly chomps off the terminating character
     * @param fileHandle - ifstream handle to file
     * @param lineBuffer - string to read into
     */
    static void getlineChomp(ifstream &fileHandle, string &lineBuffer);

    /*! \brief Check required headers
     *
     * @param header - map of header to check
     * @param requiredHeader - vector of values to check
     * @param missingIndex - integer containing idex of header we're missing
     */
    static bool checkHeader(map<string, unsigned int> &header,
                     vector<string> &requiredHeader,
                     vector<int> &requiredHeaderIndex,
                     int &missingIndex);
    /*! \brief Check required headers
     *
     * @param header - vector of header to check
     * @param requiredHeader - vector of values to check
     * @param missingIndex - integer containing idex of header we're missing
     */
    static bool checkHeader(vector<string> &header,
                     vector<string> &requiredHeader,
                     vector<int> &requiredHeaderIndex,
                     int &missingIndex);
		     
  //========================================================================================================
  
    static bool loadDelimitedFileNoHeader(const char * filename,
                                   const string &delim,
                                   const string &commentDelim,
                                   list<list<string> > &lines);

    static bool loadDelimitedFileNoHeader(const char * filename,
                                   const string &delim,
                                   const string &commentDelim,
                                   list<list<float> > &lines);

    static bool loadDelimitedFileNoHeader(const char * filename,
                                   const string &delim,
                                   const string &commentDelim,
                                   list<list<double> > &lines);

    static bool loadDelimitedFileNoHeader(const char * filename,
                                   const string &delim,
                                   const string &commentDelim,
                                   list<list<int> > &lines);

    static bool loadDelimitedFileNoHeader(ifstream &filename,
                                   const string &delim,
                                   const string &commentDelim,
                                   list<list<int> > &lines);

    static bool loadDelimitedFileNoHeader(ifstream &filename,
                                   const string &delim,
                                   const string &commentDelim,
                                   list<list<string> > &lines);

    static bool loadDelimitedFileNoHeader(ifstream &filename,
                                   const string &delim,
                                   const string &commentDelim,
                                   list<list<double> > &lines);

    static bool loadDelimitedFileNoHeader(ifstream &filename,
                                   const string &delim,
                                   const string &commentDelim,
                                   list<list<float> > &lines);

     static bool loadDelimitedFile(const char * filename,
                           const string &delim,
                           const string &commentDelim,
                           map<string, unsigned int> &header,
                           list<list<string> > &lines,
                           vector<string> &requiredHeader,
                           vector<int> &requiredHeaderIndex);


    static bool loadDelimitedFile(const char * filename,
                           const string &delim,
                           const string &commentDelim,
                           map<string, unsigned int> &header,
                           list<list<float> > &lines,
                           vector<string> &requiredHeader,
                           vector<int> &requiredHeaderIndex);

    static bool loadDelimitedFile(const char * filename,
                           const string &delim,
                           const string &commentDelim,
                           map<string, unsigned int> &header,
                           list<list<double> > &lines,
                           vector<string> &requiredHeader,
                           vector<int> &requiredHeaderIndex);

    static bool loadDelimitedFile(const char * filename,
                           const string &delim,
                           const string &commentDelim,
                           map<string, unsigned int> &header,
                           list<list<int> > &lines,
                           vector<string> &requiredHeader,
                           vector<int> &requiredHeaderIndex);

    static bool loadDelimitedFile(const char * filename,
                           const string &delim,
                           const string &commentDelim,
                           vector<string> &header,
                           list<list<int> > &lines,
                           vector<string> &requiredHeader,
                           vector<int> &requiredHeaderIndex);

    static bool loadDelimitedFile(const char * filename,
                           const string &delim,
                           const string &commentDelim,
                           vector<string> &header,
                           list<list<float> > &lines,
                           vector<string> &requiredHeader,
                           vector<int> &requiredHeaderIndex);

    static bool loadDelimitedFile(const char * filename,
                           const string &delim,
                           const string &commentDelim,
                           vector<string> &header,
                           list<list<double> > &lines,
                           vector<string> &requiredHeader,
                           vector<int> &requiredHeaderIndex);

    static bool loadDelimitedFile(const char * filename,
                           const string &delim,
                           const string &commentDelim,
                           vector<string> &header,
                           list<list<string> > &lines,
                           vector<string> &requiredHeader,
                           vector<int> &requiredHeaderIndex);
  
    static bool loadHeader(const char * filename,
			    const string &delim,
			    const string &commentDelim,
			    vector<string> &header,
			    vector<string> &requiredHeader,
			    vector<int> &requiredHeaderIndex,
			    ifstream & delimitedFileHandle);
		     
    static bool getNextLine(ifstream & delimitedFileHandle,
			    const string &delim,
			    const string &commentDelim,
			    vector<string> & fields);
			    
    static bool getNextLine(ifstream & delimitedFileHandle,
			    const string &delim,
			    const string &commentDelim,
			    vector<char *> & fields);
  };
}
#endif /* DELIMITEDTEXTREADER_H_ */
