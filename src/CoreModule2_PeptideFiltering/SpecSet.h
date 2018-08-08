/*
 * SpecSet.h
 *
 *  Created on: Aug 30, 2011
 *      Author: jsnedecor
 */

#ifndef SPECSET_H_
#define SPECSET_H_

#include "spectrum.h"
#include "PeptideSpectrumMatchSet.h"

//TR1 includes. GCC 4.0 and above only!
#ifdef __GLIBCXX__
#  include <tr1/unordered_map>
#else
#  ifdef __IBMCPP__
#    define __IBMCPP_TR1__
#  endif
#  include <unordered_map>
#endif

namespace specnets
{

  typedef vector<pair<int, vector<int> > > SpecCompareData;

  typedef tr1::unordered_map<string, unsigned int> SpectrumMap;

  class Spectrum;

  class SpecSet
  {
    friend class SpecPlot;
    friend class ContPlot;
    friend class SpsPlot;
    friend class ReportTableGenerator;

  public:

    /**
     * Standard constructor
     * @param sz initial number of spectra to make room for
     * @param enableIndexing if true, index spectra by their getUniqueID() value as they are added.
     *   Index becomes invalid if Spectrum scan numbers or filenames are modified outside of this class.
     */
    SpecSet(unsigned int sz = 0, bool enableIndexing = false);

    /**
     * TODO: add description
     *
     *@param filename
     */
    SpecSet(char *filename);

    /**
     * Deconstructor. Deletes index if necessary
     */
    ~SpecSet();

    /**
     * TODO: add description
     *
     *@param i
     *@return
     */
    inline Spectrum &operator[](unsigned int i)
    {
      return specs[i];
    }

    inline const Spectrum & operator[](unsigned int i) const
    {
      return specs[i];
    }


    /**
     * Copy contents from another SpecSet
     *
     *@param other
     *@return
     */
    SpecSet &operator=(const SpecSet &other);

    /**
     * Add Spectrum to back of SpecSet
     *
     *@param other
     *@return
     */
    void push_back(const Spectrum & x);

    /**
     * Insert Spectrum into SpecSet
     */
    void insert(vector<Spectrum>::iterator position,
                vector<Spectrum>::iterator start,
                vector<Spectrum>::iterator end);

    vector<Spectrum>::iterator begin(void);

    vector<Spectrum>::iterator end(void);

    /**
     * Get scan at scan number
     * @param scan_num - scan number (indexed from 1!)
     */
    Spectrum* getScan(int scan_num);

    /**
     * Gets spectrum index with specified ID
     * @param specID Spectrum ID (see Spectrum::getUniqueID())
     * @return -1 if spectrum is not found, otherwise its index
     */
    int getIndexFromID(const string& specID);

    /**
     * Gets spectrum with specified ID
     * @param specID Spectrum ID (see Spectrum::getUniqueID())
     * @return 0 if spectrum is not found, otherwise its pointer
     */
    Spectrum* getIndex(const string& specID);

    /**
     * Gets spectrum with specified ID
     * @param specID Spectrum ID (see Spectrum::getUniqueID())
     * @return 0 if spectrum is not found, otherwise its pointer
     */
    const Spectrum* getIndex(const string& specID) const;

    /**
     * Indexes the SpecSet by each Spectrum's ID (see Spectrum::getUniqueID()).
     *   If the SpecSet was constructed with enableIndexing=false, the index map is initialized.
     *   The index becomes invalid if Spectrum scan numbers or filenames are modified outside of this class.
     */
    void index();

    /**
     * TODO: add description
     *
     *@return
     */
    inline unsigned int size() const
    {
      return specs.size();
    }

    /**
     * Swaps the contents of this SpecSet with other w/o making an in-memory copy. Runs in constant time.
     */
    void swap(SpecSet& other);

    /**
     * Clears the set of spectra
     */
    void clear();

    /**
     * TODO: add description
     */
    unsigned int resize(unsigned int newSize);

    /**
     * Sets the peak tolerance in classic or ppm form
     */
    void setPeakTolerance(float tolerance, bool applyPPM = false);

    /**
     * Sets the parent mass tolerance in classic or ppm form
     */
    void setParentMassTolerance(float tolerance, bool applyPPM = false);

    /**
     * Removes or resize spectra with less than minimum peak count
     */
    void removeSpectraBelowMinPeaks(int min_peaks, bool remove = false);

    /**
     * Appends the set of spectra from another SpecSet to this one by swapping in the other SpecSet's memory.
     * @param other set of Spectra to swap. Since it is being swapped with empty spectra,
     *   it will have size 0 when this function returns.
     * @param updateScans if true, increment all scan #s of the appended spectra by the current size of this SpecSet
     */
    void swapAppendSpecSet(SpecSet& other, bool updateScans = true);

    /**
     * Appends the set of spectra from another SpecSet to this one.
     * @param other set of spectra to append to this one.
     * @param updateScans if true, increment all scan #s of the appended spectra by the current size of this SpecSet
     */
    void appendSpecSet(const SpecSet& other, bool updateScans = true);

    /**
     * Clear all associated psmList values in SpecSet
     *
     */
    void clearPsms(void);

    /**
     * TODO: add description
     *
     *@param results
     *@param resolution
     *@return
     */
    vector<list<int> > &getMassesHistogram(vector<list<int> > &results,
                                           float resolution = 1.0) const;

    /**
     * TODO: add description
     *
     *@param newResolution
     *@param enforceRounding
     */
    void setResolution(float newResolution, bool enforceRounding);

    /**
     * Set fileName parameter on each spectrum in specset.
     *
     *@param filename Filename to set for all spectra
     */
    void setFilename(string &filename);

    /**
     * TODO: add description
     *
     *@param tolerance
     *@param ionOffset
     *@param includeY0k
     */
    void addZPMpeaks(float tolerance, float ionOffset, bool includeY0k);

    /**
     * Compute the average intensity of all peaks in all spectrum
     */
    float averageIntensity(void);

    /**
     * TODO: add description
     *
     */
    bool maximumParsimony(void);

    /**
     * Copies spectra from this set to outputSpecs that match a fragmentation method
     * @param outputSpecs
     * @param fragType
     * @return # of spectra matching fragmentation method
     */
    unsigned int extractSpectra(SpecSet& outputSpecs,
                                Spectrum::FragType fragType);

    /**
     * Swaps out spectra from this set to outputSpecs that match a fragmentation method.
     *   Any spectra in this set matching fragType will be empty when this function returns.
     * @param outputSpecs
     * @param fragType
     * @return # of spectra matching fragmentation method
     */
    unsigned int swapExtractSpectra(SpecSet& outputSpecs,
                                    Spectrum::FragType fragType);

    /**
     * Attempts to extract scan #s of spectra acquired from the same precursor ion.
     *   Paired (or tripled, quadrupled, etc) spectra must be found in consecutive order and
     *   have the same parent mass within 1 ppm -> the accuracy of float.
     * @param outputPairs sequence of clusters, each holding the scans of all spectra belonging
     *   to the same precursor
     * @return # of clusters of size 2 or more
     */
    unsigned int extractPairedScans(vector<vector<unsigned int> >& outputPairs,
                                    int numConsec = -1,
                                    int maxConsec = -1);

    /**
     * TODO: add description
     *
     *@param newTotalInt
     */
    void normalize(float newTotalInt = 1);

    /**
     * TODO: add description
     *
     *@param features
     *@param featureValue
     *@param output
     *@return
     */
    template<class T> unsigned int extract(vector<T> &features,
                                           T featureValue,
                                           SpecSet &output);

    /**
     * TODO: add description
     *
     *@param filename
     *@return
     */
    unsigned int Load(const char *filename, const char *ext = NULL);
    /**
     * TODO: add description
     *
     *@param filename
     *@return
     */
    bool LoadSpecSet_pkl_mic(const char* filename);
    /**
     * TODO: add description
     *
     *@param filename
     *@return
     */
    unsigned int LoadSpecSet_pkl(const char *filename);

    /**
     *  LoadSpecSet_ms2: Loads a set of spectra from a .ms2 file. Spectra must be separated
     *  by at least one empty line.
     *
     *   Note: If a spectrum has more than 1 header then this function only uses the last header.
     *
     *@param filename
     *@return
     */
    unsigned int LoadSpecSet_ms2(const char *filename);

    /**
     *  LoadSpectrum_mgf: Loads a sing spectrum, if use_index = 1, then it attempts to use index, if index doesnt exist, then we will create it
     *
     *@param filename
     *@return
     */
    unsigned int LoadSpectrum_mgf(const char *filename,
                                  unsigned int scan,
                                  int index,
                                  int use_index);

    /**
     * LoadSpecSet_mgf: Loads a set of spectra from a .mgf file. Recognizes and processes
     * these MGF tags: BEGIN IONS, END IONS, CHARGE, PEPMASS, PEPTIDE
     *
     * Note: If a spectrum has more than 1 header then this function only uses the last header.
     *
     *@param filename
     *@return
     */
    unsigned int LoadSpecSet_mgf(const char *filename);

    /**
     * LoadSpecSet_mgf: Loads a set of spectra from a .mgf file. Recognizes and processes
     * these MGF tags: BEGIN IONS, END IONS, CHARGE, PEPMASS, PEPTIDE
     *
     * Note: If a spectrum has more than 1 header then this function only uses the last header.
     *
     *@param filename
     *@param scan
     *@param index
     *@return
     */
    unsigned int LoadSpecSet_mgf(const char *filename,
                                 unsigned int scan,
                                 int index);

    /**
     * Loads a specset with the accompanying peptide spectrum matches with input being pklbin
     *
     *@param spectra_filename
     *@param annotation_filename
     *@return
     */
    unsigned int
    LoadSpecSet_pklbin_with_annotation(const char * spectra_filename,
                                       const char * annotation_filename);

    /**
     * Loads a specset with the accompanying peptide spectrum matches with input being mgf
     *
     *@param spectra_filename
     *@param annotation_filename
     *@return
     */
    unsigned int
    LoadSpecSet_mgf_with_annotation(const char * spectra_filename,
                                    const char * annotation_filename);

    /**
     * LoadSpecSet_prms: Loads a set of spectra from a .prms file, as output but the
     * current version of pepnovo_prms (2007/04/09). Text file with multiple spectra,
     * each delimited in the following format:
     *
     * - ">> <original msms file index> <scan/index # in msms file> <single-spectrum file name>"
     *    -- scan/index # in msms file is stored in the Spectrum.scan field.
     *    -- The whole info line is stored in the Spectrum.info field.
     *  - Set of peaks in "<peak mass> <peak intensity/score>" format
     *  - Spectrum terminated by empty line
     *  - Last pair of mass/intensity corresponds to the corrected sum-of-amino-acid-masses parent mass
     *
     *  Note: If a spectrum has more than 1 header then this function only uses the last header.
     *
     *@param filename
     *@return
     */
    unsigned int LoadSpecSet_prms(const char *filename);

    /**
     * Loads spectra in format saved by PepNovo
     *
     *@param filename
     *@param prmOrigins Filled by 3rd column of each spectrum if specified. In the latest version of PepNovo, each PRM's origin MS/MS peak is saved in this 3rd column.
     *@return
     */
    unsigned int LoadSpecSet_prmsv3(const char *filename,
                                    vector<vector<string> >* prmOrigins = 0);

    /** DEPRECATED
     * LoadSpecSet_pklbin - loads a set of spectra in binary format. File format
     * 1 int - number of spectra in the file
     * numSpecs shorts - number of peaks per spectrum in the file.
     * arrays of numPeaks+1 pairs of floats - [parentMass charge] + array of [mass intensity]
     *
     *@param filename Name of the pklbin file to load
     *@param countSpectraOnly If then only counts number of spectra in the file without loading the spectra (defaults to false)
     *@return
     */
    //unsigned int LoadSpecSet_pklbin(const char *filename,
    //                                bool countSpectraOnly = false);
    /**
     * VERSION 1: Loads all data fields of each spectrum from an open pklbin file
     * @param fp Open file pointer where the number of spectra and version of
     *    already been read from. This will be closed upon return
     * @param numSpecs Number of spectra to read. This has already been read from the file
     * @param oldVersion true if reading in the original pklbin format (pre version 1).
     * @param subversion Sub-version of this format. Currently not used (at version 1, subversion 0)
     * @return true if the SpecSet was successfully loaded, false if not. In either case the file will
     *    be closed
     */
    bool loadPklBin_1(FILE* fp, int numSpecs, bool oldVersion, char subversion);

    /**
     * VERSION 2: Loads all data fields of each spectrum from an open pklbin file
     * @param fp Open file pointer where the number of spectra and version of
     *    already been read from. This will be closed upon return
     * @param numSpecs Number of spectra to read. This has already been read from the file
     * @param subversion Sub-version of this format. Currently not used (at version 2, subversion 0)
     * @return true if the SpecSet was successfully loaded, false if not. In either case the file will
     *    be closed
     */
    bool loadPklBin_2(FILE* fp, int numSpecs, char subversion);

    /**
     * VERSION 3: Loads all data fields of each spectrum from an open pklbin file
     * @param fp Open file pointer where the number of spectra and version of
     *    already been read from. This will be closed upon return
     * @param numSpecs Number of spectra to read. This has already been read from the file
     * @param subversion Sub-version of this format. Currently not used (at version 2, subversion 0)
     * @return true if the SpecSet was successfully loaded, false if not. In either case the file will
     *    be closed
     */
    bool loadPklBin_3(FILE* fp, int numSpecs, char subversion);

    /**
     * Loads all data fields of each spectrum in a pklbin file. Backwards compatible with all old versions.
     *
     *@param filename
     *@param psmFileame
     *@param peaksFilename
     *@return
     */
    int loadPklBin(const char * filename,
                   const char * psmFileame = 0x0,
                   const char * peaksFilename = 0x0);

    /**
     * TODO: add description
     *
     *@param filename
     *@return
     */
    short SaveSpecSet(const char *filename, bool outputScans = true);

    /**
     * TODO: add description
     *
     *@param filename
     *@return
     */
    short SaveSpecSet_mgf(const char *filename, bool outputScans = true);

    /**
     *@param base_filename
     *@return
     */
    bool SaveSpecSet_dta(const char *filename);

    /**
     * TODO: add description
     *
     *@param filename
     *@return
     */
    short SaveSpecSet_pkl(const char *filename);

    /**
     * Saves the annotation_peptide field of each spectrum to the indexed line in output file
     *@param filename output filename, if a spectrum has an annotation_peptide it will be written to the same line # as its index
     *@return true if file was written successfully, false otherwise
     */
    bool SaveAnnotations(const char* filename);

    /**
     * TODO: add description
     *
     *@param filename
     *@return
     */
    short SaveSpecSet_ms2(const char* filename);

    /**
     * saves the .bin component of a set of spectra in .bin format containing (per row)
     *
     * 1 int  - scan #
     * 1 int - ms level
     *
     *@param filename
     *@return
     */
    short SaveSpecSet_bin(const char *filename);

    /** DEPRECATED
     * TODO: add description
     *
     *@param filename
     *@param filename
     *@return
     */
    //short SaveSpecSet_pklbin(const char *filename, const char *binFilename =
    //    NULL);
    /**
     * TODO: add description
     *
     *@param filename
     *@return
     */
    void SaveScanNums(const char *filename);

    /** DEPRECATED
     * SaveSpecSet_pklbin - saves a set of spectra in binary format. File format
     * 1 int - number of spectra in the file
     * numSpecs shorts - number of peaks per spectrum in the file
     * arrays of numPeaks+1 pairs of floats - [parentMass charge] + array of [mass intensity]
     *
     *@param filename
     *@param specs
     *@return
     */
    // short SaveSpecSet_pklbin(char *filename, vector<Spectrum *> specs);
    /**
     * TODO: add description
     *
     *@param filename
     *@param specs
     *@return
     */
    int saveMatchedProts(const char *filename);

    /**
     * TODO: Current pklbin format. Includes scan # and msLevel.
     *
     *@param filename
     *@return
     */
    int savePklBin(const char *filename,
                   const char * psmFileame = 0x0,
                   const char * peaksFilename = 0x0);

    /**
     * compare - Compares two specsets
     *
     *@param toSet specset to compare to
     *@param SpecCompareData structure to hold compare results
     *@return number of different spectra
     */
    int compare(SpecSet &toSet, SpecCompareData &cd);

    void saveTags(const std::string & filename, int tagLen, int gap, int tagNum, float fragTol, bool isPRMSpec);

  protected:
    vector<Spectrum> specs;

    SpectrumMap* m_index;

    void computeIndex();

    unsigned FindSpectrumInMgfByScan(BufferedLineReader &blr,
                                     unsigned &lineIdx,
                                     unsigned scan);

    unsigned FindSpectrumInMgfByScanOrIndex(const char * filename,
                                            unsigned scan,
                                            int index,
                                            unsigned int &filestartposition);

    unsigned int Calculate_MGF_index(const char *filename,
                                     vector<vector<unsigned int> > &file_index);

    bool FindSpectrumInMgfByIndex(BufferedLineReader &blr,
                                  unsigned &lineIdx,
                                  unsigned index);

  };
}
#endif /* SPECSET_H_ */
