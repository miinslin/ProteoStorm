/*
 * PeptideSpectrumMatchSet.h
 *
 *  Created on: Apr 27, 2011
 *      Author: jsnedecor
 */

#ifndef PEPTIDESPECTRUMMATCHSET_H_
#define PEPTIDESPECTRUMMATCHSET_H_

//Module includes
#include "DelimitedTextReader.h"
#include "PeptideSpectrumMatch.h"
#include "SpecSet.h"

//TR1 includes. GCC 4.0 and above only!
#ifdef __GLIBCXX__
#  include <tr1/memory>
#  include <tr1/unordered_map>
#  include <tr1/unordered_set>
#else
#  ifdef __IBMCPP__
#    define __IBMCPP_TR1__
#  endif
#  include <memory>
#  include <unordered_map>
#  include <unordered_set>
#endif

//System includes
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cstring>
#include <map>

namespace specnets
{

  /**
   * @see PeptideSpectrumMatch.h
   */
  class PeptideSpectrumMatch;

  /**
   * @see SpecSet
   */
  class SpecSet;

  typedef std::tr1::shared_ptr<PeptideSpectrumMatch> psmPtr;

  class PeptideSpectrumMatchSet
  {
  public:

    /*! \brief Return associated  PeptideSpectrumMatch for
     that vector position.
     */
    psmPtr & operator[](unsigned int i);

    /*! \brief Return associated  PeptideSpectrumMatch
     for that vector position.
     */
    const psmPtr & operator[](unsigned int i) const;

    /*! \brief Set value of one PeptideSpectrumMatchSet to another

     */
    PeptideSpectrumMatchSet & operator=(PeptideSpectrumMatchSet &other);

    /*! \brief Adds to psmSet vector
     *
     */
    unsigned int push_back(const psmPtr &other);

    /*! \brief Returns size of m_psmSet vector

     */
    inline unsigned int size() const
    {
      return (unsigned int)m_psmSet.size();
    }

    /*! \brief Resizes m_psmSet vector

     @param newSize the new size of the parameter vector
     */
    unsigned int resize(unsigned int newSize);

    /** Counts modification frequencies
     *
     */
    void getModCounts(map<string, map<float, float> > & mapModCount,
                      bool useOrig = false,
                      vector<pair<string,float> > * pexclusionList = 0x0);
    void getModCounts(map<string, map<float, float> > & mapModCount,
                      map<float, vector<string> > mapPrefered,
                      int nPref,
                      bool useOrig,
                      vector<pair<string,float> > * pexclusionList = 0x0);
    void getRecursiveModCounts(map<string, map<float, float> > & mapModCount,
                               PeptideSpectrumMatchSet & returnSet,
                               bool useOrig,
                               vector<pair<string,float> > * pexclusionList = 0x0);

    void recomputeGapsByPreferredPtm(map<float, vector<string> > & mapPrefered,
                                     int nPref,
                                     bool useOrig);

    /** Saves the modification frequency table to a file
     *
     */
    bool saveModMatrix(const char * filename,
                       bool useOrig = false,
                       bool recurse = false,
                       vector<pair<string,float> > * pexclusionList = 0x0);

    /** Loads the PSM set from a standard file
     *
     */
    bool loadFromFile(const char * filename, int firstScan = -1, int lastScan =
                          -1);

    /** Loads only the numbers of PSMs from a standard file
     *
     */
    bool loadSizesFromFile(const char * filename,
                           map<int, int> & mapPsmCounts,
                           int & totalPsms,
                           int firstScan = -1,
                           int lastScan = -1);

    /** Load multiple results files.
     *
     *@param resultsFileList = tab delimited file which contains columns Path.
     */
    bool loadFromFiles(const char * resultsFileList);

    bool saveToFile(const char * filename, 
                    bool includeHeader = true,
                    bool variantCol = false);

    /** Inspect result file parsing
     *
     *@param resultsFile = input inspect results
     *@param zeroIndexed = is file zero indexed or not (yes for MGF, no for mzXML)
     */
    bool loadInspectResultsFile(const char * resultsFile, bool zeroIndexed = 1);
    /** Load multiple inspect results files.
     *
     *@param resultsFileList = tab delimited file which contains columns Path and (optionally) isZeroIndexed.
     */
    bool loadInspectResultsFiles(const char * resultsFileList);

    bool loadMSGFPlusResultsFile(const char * resultsFile, bool zeroIndexed =
                                     false);

    /** MSGFDB result file parsing
     *
     *@param resultsFile = input MSGFDB results
     *@param zeroIndexed = is file zero indexed or not (should always be false for MSGFDB)
     */
    bool loadMSGFDBResultsFile(const char * resultsFile, bool zeroIndexed =
                                   false);

    /** Load multiple msgf results files.
     *
     *@param resultsFileList = tab delimited file which contains columns Path and (optionally) isZeroIndexed.
     */

    bool loadMSGFDBResultsFiles(const char * resultsFileList);

    /** MODa result file parsing
     *
     *@param resultsFile = input MODa results
     *@param zeroIndexed = is file zero indexed or not (should always be false for MODa)
     */
    bool loadModaResultsFile(const char * resultsFile, 
                             bool zeroIndexed = false);

    /** MSplit result file parsing
     *
     *@param resultsFile = input MSplit results
     *@param zeroIndexed = is file zero indexed or not (should always be false for MSplit?)
     */
    bool loadMsplitResultsFile(const char * resultsFile, 
                               bool zeroIndexed = false);

    /** MSplit result file parsing
     *
     *@param resultsFile = input MSplit results
     *@param zeroIndexed = is file zero indexed or not (should always be false for MSplit?)
     */
    bool loadMsplitResultsFile2(const char * resultsFile,
                                bool zeroIndexed = false);

    /** Multipass result file parsing
     *
     *@param resultsFile = input Multipass results
     *@param zeroIndexed = is file zero indexed or not (should always be false for MODa)
     */
    bool loadMultipassResultsFile(const char * resultsFile, 
                                  bool zeroIndexed = false);

    /** SpecNets report file parsing
     *
     *@param resultsFile = input SpecNets report results
     *@param zeroIndexed = is file zero indexed or not (should always be false for MODa)
     */
    bool loadSpecnetsReportFile(const char * resultsFile);

    /** Specnets result file parsing
     *
     *@param resultsFile = input inspect results
     *@param zeroIndexed = is file zero indexed or not (yes for MGF, no for mzXML)
     */
    bool
    loadSpecnetsResultsFile(const char * resultsFile, bool zeroIndexed = 0);

    /**
     * Loads the appropriate results file depending on the format
     * @param resultsFile
     * @param fileType
     *           "" --> PSMSet
     *           "msgfdb" --> MSGFDB
     *           "inspect" --> Inspect
     *           "moda" --> ModA
     *           "msgfplus" --> MSGF+
     * @return true if file was successfully loaded, false if not
     */
    bool Load(const string& resultsFile, const string& fileType);

    /**
     * Match up annotations with spectra
     * @param spectra spectra to add to psms
     *  @param addToSpectra = choose whether to update psmList on
     * spectra
     */
    unsigned int addSpectra(SpecSet * spectra,
                            bool addToSpectra = true,
                            bool suppressWarnings = false);

    /**
     * Match up annotations with spectra
     *
     * @param spectra spectra to add to psms
     * @param filename = original file name
     * @param addToSpectra = choose whether to update psmList on
     * spectra
     */
    unsigned int addSpectra(SpecSet * spectra,
                            string filename,
                            bool addToSpectra = true,
                            bool suppressWarnings = false);

    /**
     * Match up annotations with spectra based on filenames defined in spectra.
     *
     * @param spectra spectra to add to psms
     * @param filename = original file name
     * @param addToSpectra = choose whether to update psmList on
     * spectra
     */
    unsigned int addSpectraByFilename(SpecSet * spectra, bool addToSpectra =
                                          true,
                                      bool suppressWarnings = false);

    /**
     * Adds the most spectra, either by scan # only or scan # + filename. Tie goes to scan # only
     * @param spectra spectra to add to psms
     * @param addToSpectra = choose whether to update psmList on
     * @return number of annotated spectra
     */
    unsigned int addMostSpectra(SpecSet * spectra,
                                bool addToSpectra = true,
                                bool suppressWarnings = false);

    /**
     * DEPRECATED. THIS ASSUMES THERE WAS ONE INPUT FILE, please use the other cluster() method below.
     *
     * Clusters PSMs. Scan #s in each mapped list of clusterInfo must match to scan numbers loaded into
     *   this PSM set. Each cluster of PSMs is replaced a single consensus PSM
     *
     * @param clusterInfo indexed by scan number of clustered spectrum. Each value is
     *   a list containing all scan numbers that were merged into the corresponding cluster
     * @param mergeType details how to choose consensus PSM for each cluster
     *        0 -> PSM with most abundant peptide match
     *        (add more here)
     * @return number of PSMs that were clustered
     */
    int cluster(map<int, list<int> >& clusterInfo, short mergeType = 0);

    /**
     * Clusters PSMs. Pairs of scan #s and filenames in each mapped list of clusterInfo must match to scan
     *   #s/filenames loaded into this PSM set. Each cluster of PSMs is replaced a single consensus PSM
     *
     * @param clusterInfo indexed by scan number of clustered spectrum. Each value is
     *   a list containing all pairs of scan #s/filenames that were merged into the corresponding cluster
     * @param mergeType details how to choose consensus PSM for each cluster
     *        0 -> PSM with most abundant peptide match
     *        (add more here)
     * @return number of PSMs that were clustered
     */
    int cluster(map<int, list<pair<int, string> > >& clusterInfo,
                short mergeType = 0);

    /** Return PSM for a particular scan number
     *
     */
    bool getResultsByScan(int scan,
                          vector<PeptideSpectrumMatchSet>& output_results);

    /** Build PSMSet from input SpecSet
     *  Note: this will overwrite any existing PSMs.
     */
    void getPSMSet(SpecSet * spectra);

    /** Make PSMSet from input SpecSet
     *     This is similar to getPSMSet except it does not make new PSMs but uses shared pointers
     *  Note: this will overwrite any existing PSMs.
     */
    void makePSMSetShared(SpecSet * spectra);

    void removePsmSetItem(psmPtr);

    void maximumParsimony(void);

    void sortBySpecIndex();//actually by m_scanNum

    vector<psmPtr> m_psmSet;

  protected:
    int m_doubleLoad;
    void createMapIndex(vector<string> & header, vector<int> & mapIndex);

    bool load(const char * fileName,
              vector<string> & requiredHeader,
              string fieldDelim,
              string commentDelim,
              bool zeroIndexed,
              bool isInspect,
              int firstScan = -1,
              int lastScan = -1);

    bool parseLine(psmPtr currMatch,
                   vector<string> & line,
                   vector<int> & mapIndex,
                   bool zeroIndexed,
                   bool isInspect);

    bool parseLine(psmPtr currMatch,
                   vector<char *> & line,
                   vector<int> & mapIndex,
                   bool zeroIndexed,
                   bool isInspect);

    bool readHeader(const char * filename,
                    ifstream & ifs,
                    string fieldDelim,
                    string commentDelim,
                    vector<int> & index);

    bool readNextPsm(ifstream & ifs,
                     vector<int> & index,
                     string fieldDelim,
                     string commentDelim,
                     bool zeroIndexed,
                     bool isInspect,
                     psmPtr nextPsm);
  };

  class PeptideSpectrumMatchSetSpectralLibraryLoader : public PeptideSpectrumMatchSet
  {
  public:
    bool
    loadSpecnetsResultsFile(const char * resultsFile, bool zeroIndexed = 0);
  };
}

#endif /* PEPTIDESPECTRUMMATCHSET_H_ */
