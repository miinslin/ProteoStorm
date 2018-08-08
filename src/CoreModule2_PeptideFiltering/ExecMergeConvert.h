/*
 * ExecMergeConvert.h
 *
 *  Created on: Oct 5, 2011
 *      Author: aguthals
 */

#ifndef EXECMERGECONVERT_H_
#define EXECMERGECONVERT_H_

#include "ExecBase.h"

#include "utils.h"
#include "spectrum.h"
#include "SpecSet.h"
#include "Specific.h"

// System Includes
#include <string>
#include <vector>

using namespace std;

namespace specnets
{
  class ExecMergeConvert : public ExecBase
  {
  public:

    static bool
    loadSpecsetMultiple(const string &exeDir,
                        const string &specFilePaths,
                        SpecSet* spectra,
                        vector<string>* fileNames = 0,
                        vector<pair<int, int> >* fileSpecIdxs = 0,
                        const char* separator = ";");

    static bool convertSpecset(const string& exeDir,
                               const string& inFile,
                               const string& outFile);

    static bool loadSingleSpecset(const string &exeDir,
                                  const string &specFilePath,
                                  SpecSet* spectra,
                                  FilenameManager &fm);

    static bool loadSpecset(const string& exeDir,
                            const string &specFilePath,
                            SpecSet* spectra,
                            vector<string>* fileNames = 0,
                            vector<pair<int, int> >* fileSpecIdxs = 0);

    static bool saveSpecset(const string& specFilePath, 
                            SpecSet* spectra, 
                            bool outputScans = true);

    static bool saveSpecsetMultiple(const string& dirName,
                                    const string& specFilePaths,
                                    SpecSet* spectra,
                                    vector<string>* fileNames = 0,
                                    const char* separator = ";");
                                    
    static bool saveSpecsetMultiple(const string& specFilePaths,
                                    SpecSet* spectra,
                                    vector<string>* fileNames = 0,
                                    const char* separator = ";");

    static pair<bool, list<string> > getFileList(const string& listFilePath);

    stringstream m_recordedProcessing;
    string m_separator;

    ExecMergeConvert(void);

    ExecMergeConvert(const ParameterList & inputParams);

    ExecMergeConvert(const ParameterList & inputParams,
                     vector<pair<int, int> >* inputSpecsetIdx,
                     SpecSet * inputSpecset);

        ExecMergeConvert(const ParameterList & inputParams,
                         SpecSet * outputSpecset);

    ExecMergeConvert(const ParameterList & inputParams,
                     vector<pair<int, int> >* inputSpecsetIdx,
                     SpecSet * inputSpecset,
                     SpecSet * outputSpecset);

    virtual ~ExecMergeConvert(void);

    virtual ExecBase * clone(const ParameterList & input_params) const;

    virtual bool invoke(void);

    virtual bool loadInputData(void);

    virtual bool saveOutputData(void);

    virtual bool saveInputData(std::vector<std::string> & filenames);

    virtual bool loadOutputData(void);

    virtual std::vector<ExecBase *> const & split(int numSplit);

    virtual bool merge(void);

    virtual bool validateParams(std::string & error);

    virtual bool saveActivity(const string& infoFilePath);

  private:
    vector<pair<int, int> >* m_inputSpecsetIdx; // Set of input specsets
    bool ownInput; //! Does this object "own" the input data structures (and hence have to free them)

    SpecSet * m_outputSpecset; //! Output specset
    SpecSet * m_inputSpecset; //! Input specset
    bool ownOutput; //! Does this object "own" the output data pointers (and hence have to free them)

    vector<pair<string, int> >* m_recordedInput;
    vector<pair<string, int> >* m_recordedOutput;

    bool initialize();
  };
}

#endif /* EXECMERGECONVERT_H_ */
