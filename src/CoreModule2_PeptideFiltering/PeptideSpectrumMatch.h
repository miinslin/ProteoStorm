/*
 * PeptideSpectrumMatch.h
 *
 *  Created on: Apr 5, 2011
 *      Author: jsnedecor
 */

#ifndef PEPTIDESPECTRUMMATCH_H_
#define PEPTIDESPECTRUMMATCH_H_

//Module includes
#include "utils.h"
#include "spectrum.h"
#include "spectrum_scoring.h"
#include "aminoacid.h"
#include "MathUtils.h"
#include "IonMass.h"

//System includes
#include <iostream>
#include <fstream>
#include <vector>
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
  /**
   * @see spectrum.h
   */
  class Spectrum;

  struct ftIonFragment;

  class MS2ScoringModel;

  class MS2ScoringModelSet;

  class IonMass;

  /*! \brief Peptide spectrum match class. A single spectrum is matched to a single annotation
   *
   * It is possible to have more than one annotation per spectrum, however, they will show up
   * as different PSMs.
   */
  class PeptideSpectrumMatch
  {
  public:
    //! \name CONSTRUCTORS
    //@{
    /*! \brief constructor for PSM class
     *
     */
    PeptideSpectrumMatch();
    //@}

    //! \name DESTRUCTOR
    //@{
    virtual ~PeptideSpectrumMatch();
    //@}

    /*! \brief set PSM to another PSM, copies peak annotations, spectrum pointers, etc.
     *
     */
    PeptideSpectrumMatch(const PeptideSpectrumMatch &other);

    /*! \brief set PSM to another PSM, copies peak annotations, spectrum pointers, etc.
     *
     */
    virtual PeptideSpectrumMatch &operator=(const PeptideSpectrumMatch &other);

    /** helper function for annotate. Runs through match vector and
     *sets annotation vector
     *@param matches - vector of ion indices (from 0) to set matches for
     *@param annotation - vector of fragments to set
     *@param ionIdx - index of ion 0..N-1 where N is length of peptide. i.e. for b1, ionIdx = 1
     *@param currIonFrag - pointer to current ftIonFragment we're currently annotating
     */
    virtual void
    setAnnotationToMatches(vector<int> &matches,
                           vector<pair<const ftIonFragment*, short> > &annotation,
                           int ionIdx,
                           const ftIonFragment* currIonFrag);

    /** helper function for annotate. Runs through match vector and
     *sets annotation vector, taking only the higher intensity match in case of duplicates
     *@param matches - vector of ion indices (from 0) to set matches for
     *@param annotation - vector of fragments to set
     *@param ionIdx - index of ion 0..N-1 where N is length of peptide. i.e. for b1, ionIdx = 1
     *@param currIonFrag - pointer to current ftIonFragment we're currently annotating
     */
    virtual void
    setAnnotationToMatchesNoDups(vector<int> &matches,
                                 vector<pair<const ftIonFragment*, short> > &annotation,
                                 int ionIdx,
                                 const ftIonFragment* currIonFrag);

    /**
     * Removes the +8 or +10 at the C-terminus (if it exists)
     * @param modifySpectrum modify the spectrum's parent mass accordingly. WARNING: the spectrum
     *   pointer must be properly set or this will throw an error
     * @return
     */
    virtual float removeSILACMod(const bool modifySpectrum = false);

    void setDbMatch(const string & protein,
                    const int dbIndex,
                    const float startMass);

    void getDbMatch(string & protein, int & dbIndex, float & startMass) const;

    /**
     * add annotations for matched peaks by looking up the Spectrum's fragmentation type in MS2ScoringModelSet (using per-peak tolerances)
     * @param peptide amino acid sequence used to determine peak annotations
     * @param ionNamesInclude  comma delimited string indicating which ions with which names should be annotated. If string is simply "all"
     * then just include all ions in MS2Model. ex. "y,b,y++,b++"
     * @param inputIonTypes definition of ion types: mass offsets, ion probabilities, prefix/suffix.
     * is copied into ionTypes
     * @param prmOffset sets applied to theoretical ion masses to locate Prefix Residue Masses
     * @param srmOffset sets applied to theoretical ion masses to locate Suffix Residue Masses
     * @param removeDuplicates - if true, no duplicate annotations are allowed for the same peak
     * @param jumps - AAJumps indicating amino acid masses
     */
    virtual bool annotate(const string &peptide,
                          const string &ionNamesInclude,
                          const MS2ScoringModelSet &inputIonTypes,
                          const float prmOffset,
                          const float srmOffset,
                          const AAJumps &jumps,
                          const bool removeDuplicates = true,
                          const bool ignoreParentCharge = false,
                          const bool retainOldAnnotations = false);

    /**
     * add annotations for matched peaks from MS2ScoringModel (using per-peak tolerances)
     * @param peptide amino acid sequence used to determine peak annotations
     * @param ionNamesInclude  comma delimited string indicating which ions with which names should be annotated. If string is simply "all"
     * then just include all ions in MS2Model. ex. "y,b,y++,b++"
     * @param inputIonTypes definition of ion types: mass offsets, ion probabilities, prefix/suffix.
     * is copied into ionTypes
     * @param prmOffset sets applied to theoretical ion masses to locate Prefix Residue Masses
     * @param srmOffset sets applied to theoretical ion masses to locate Suffix Residue Masses
     * @param removeDuplicates - if true, no duplicate annotations are allowed for the same peak
     * @param jumps - AAJumps indicating amino acid masses
     */
    virtual bool annotate(const string &peptide,
                          const string &ionNamesInclude,
                          const MS2ScoringModel &inputIonTypes,
                          const float prmOffset,
                          const float srmOffset,
                          const AAJumps &jumps,
                          const bool removeDuplicates = true,
                          const bool ignoreParentCharge = false,
                          const bool retainOldAnnotations = false);

    /** DEPRECATED - spectrum per-peak tolerances should already be set
     * add annotations for matched peaks from MS2ScoringModel
     * @param peptide amino acid sequence used to determine peak annotations
     * @param ionNamesInclude  comma delimited string indicating which ions with which names should be annotated. If string is simply "all"
     * then just include all ions in MS2Model. ex. "y,b,y++,b++"
     * @param inputIonTypes definition of ion types: mass offsets, ion probabilities, prefix/suffix.
     * is copied into ionTypes
     * @param prmOffset sets applied to theoretical ion masses to locate Prefix Residue Masses
     * @param srmOffset sets applied to theoretical ion masses to locate Suffix Residue Masses
     * @param peakTol - mass tolerance (in Da) when matching ion masses to peak masses
     * @param removeDuplicates - if true, no duplicate annotations are allowed for the same peak
     * @param jumps - AAJumps indicating amino acid masses
     */
    virtual bool annotate(const string &peptide,
                          const string &ionNamesInclude,
                          const MS2ScoringModel &inputIonTypes,
                          const float prmOffset,
                          const float srmOffset,
                          const float peakTol,
                          const AAJumps &jumps,
                          const bool removeDuplicates = true,
                          const bool ignoreParentCharge = false,
                          const bool retainOldAnnotations = false);

    /** DEPRECATED - spectrum per-peak tolerances should already be set
     * add annotations for matched peaks from MS2ScoringModel
     * @param peptide amino acid sequence used to determine peak annotations
     * @param ionNamesInclude  comma delimited string indicating which ions with which names should be annotated. If string is simply "all"
     * then just include all ions in MS2Model. ex. "y,b,y++,b++"
     * @param inputIonTypes definition of ion types: mass offsets, ion probabilities, prefix/suffix.
     * is copied into ionTypes
     * @param prmOffset sets applied to theoretical ion masses to locate Prefix Residue Masses
     * @param srmOffset sets applied to theoretical ion masses to locate Suffix Residue Masses
     * @param peakTol - mass tolerance (in Da) when matching ion masses to peak masses
     * @param removeDuplicates - if true, no duplicate annotations are allowed for the same peak
     */
    virtual bool annotate(const string &peptide,
                          const string &ionNamesInclude,
                          const MS2ScoringModel &inputIonTypes,
                          const float prmOffset,
                          const float srmOffset,
                          const float peakTol,
                          const bool removeDuplicates = true,
                          const bool ignoreParentCharge = false,
                          const bool retainOldAnnotations = false);

    /** DEPRECATED - spectrum per-peak tolerances should already be set
     * add annotations for matched peaks from ionInputMasses
     * @param ionInputMasses ion masses
     * @param ionInputFragments  ion fragments
     * @param peakTol - mass tolerance (in Da) when matching ion masses to peak masses
     * @param jumps - AAJumps
     * @param removeDuplicates - if true, no duplicate annotations are allowed for the same peak
     */
    virtual bool annotate(const vector<IonMass> &ionInputMasses,
                          const MS2ScoringModel &model,
                          const float peakTol,
                          const AAJumps &jumps,
                          const bool removeDuplicates = true,
                          const bool ignoreParentCharge = false,
                          const bool retainOldAnnotations = false);

    /** DEPRECATED - spectrum per-peak tolerances should already be set
     * add annotations for matched peaks from ionInputMasses
     * @param ionInputMasses ion masses
     * @param ionInputFragments  ion fragments
     * @param peakTol - mass tolerance (in Da) when matching ion masses to peak masses
     * @param removeDuplicates - if true, no duplicate annotations are allowed for the same peak
     */
    virtual bool annotate(const vector<IonMass> &ionInputMasses,
                          const MS2ScoringModel &model,
                          const float peakTol,
                          const bool removeDuplicates = true,
                          const bool ignoreParentCharge = false,
                          const bool retainOldAnnotations = false);

    /**
     * Count the number of annotated peaks and the total explained intensity (must call annotate() to annotate peaks first)
     * @param ionNamesInclude comma delimited string indicating which ions with which names should be annotated. If string is simply "all"
     *   then just include all ions in MS2Model. ex. "y,b,y++,b++".
     * @param outputIonCounts optional output map pairing each ion type's name w/ annotated peaks to a pair of: pair.first = number of peaks
     *   annotated to ion type, pair.second = total intensity of matched peaks
     * @return pair.first = total number of matching peaks, pair.second = total matched intensity

     virtual pair<int, float> countAnnotatedPeaks(string &ionNamesInclude,
     map<string, pair<int, float> >* outputIonCounts = 0);*/

    /** Sets m_charge on psm based on m_annotation and parentMZ from the associated spectrum
     *
     */
    bool setChargeByAnnotation(void);
    /** Sets m_charge on psm based on m_annotation and parentMZ from the associated spectrum
     *
     * @param jumps AAJumps object
     */
    bool setChargeByAnnotation(const AAJumps &jumps);

    /** Returns an annotation of the spectrum in the PSM using the matched peaks and the DB sequence/spectrum passed in
     *
     * @param Spectrum dbSpec database spectrum
     * @param jumps AAJumps object
     * @param bool useOriginal use m_origAnnotation (instead of m_annotation)
     */
    void getMatchedPeaksFromAnnotation(Spectrum & dbSpec,
                                       AAJumps & aaJumps,
                                       bool useOriginal = false);

    /** Returns an annotation of the spectrum in the PSM using the matched peaks and the DB sequence/spectrum passed in
     *
     */
    void getAnnotationFromMatchedPeaks(Spectrum & dbSpec,
                                       const string & dbSeqStr,
                                       string & annotation);

    /** Replaces all gap annotations with single AA annotations (on first AA)
     *
     */
    void changeGapAnnosToSingle(void);

    /** Adds fixed cysteine mods to all cysteines in annotation
     *     Call changeGapAnnosToSingle() first to pull all mods on single AA
     *
     */
    void addFixedCysteineMods(void);

    /** Maps peakList indices to ion names.
     *
     *
     * @param outputMatches output for matched indices
     * @param ionMap mapping between ion name + ion index keyed to
     * index of outputMatches
     */
    virtual void mapIons(vector<vector<int> > &outputMatches,
                         std::tr1::unordered_map<string, int> &ionMap) const;

    /** Strips off the leading and following AA specifiers
     * For example: Changes M.ACDE.R to ACDE
     * @param original = original string
     * @param stripped = output string without leading and trailing AA specification
     */
    void stripPrecedingAndFollowing(const string &original, string &specnets);

    /** Inspect annotation parsing
     * Changes Inspect annotations to SpecNets (should add converse function)
     * @param inspect = inspect string
     * @param specnets = output string for specnets
     */
    static void inspectToSpecNets(const string &inspect, string &specnets);

    /** Return whether an annotation contains a modification
     *
     *  This simply looks for the first "(" character in the string
     *  if it's there, we assume the peptide is modified
     */
    bool isModified(void) const;

    /** Return whether an annotation is tryptic
     *
     *  Checks for proper tryiptyc cleavage
     */
    bool isTryptic(bool useOrig = false);

    /** Returns all associated modifications for annotation
     *
     * This might be worth replacing with something that parses modifications
     * when m_annotation is changed.
     * @param modifications = output vector for modifications
     * @return = returns false when modification is unable to be parsed
     */
    bool getModifications(vector<float> & modifications, bool useOriginal =
                              false) const;

    /** Returns all associated modifications for annotation
     *
     * This might be worth replacing with something that parses modifications
     * when m_annotation is changed.
     * @param modifications = output vector for modifications
     * @return = returns false when modification is unable to be parsed
     */
    bool getModificationsAndPositions(vector<float> &modifications,
                                      vector<unsigned int> &positions,
                                      bool useOriginal = false) const;

    /** Returns all associated modifications for annotation
     *
     * This might be worth replacing with something that parses modifications
     * when m_annotation is changed.
     * @param modifications = output vector for modifications
     * @return = returns false when modification is unable to be parsed
     */
    bool getModificationsAndPositions(vector<float> & modifications,
                                      vector<unsigned int> & positions,
                                      vector<unsigned int> & lengths,
                                      bool useOriginal = false) const;

    /** Creates an annotation from the modification data vectors
     *
     */
    static void makeAnnotationFromData(string & cleanAnnotation,
                                       vector<float> & modifications,
                                       vector<unsigned int> & positions,
                                       vector<unsigned int> & lengths,
                                       string & completeAnnotation);

    /**
     * Returns unmodified version of the peptide
     *
     * @param inputPeptide = modified input peptide sequence
     * @param outputPeptide = output unmodified annotation
     */
    static void getUnmodifiedPeptide(const string &inputPeptide,
                                     string &outputPeptide);

    /**
     * Counts number of gaps in annotation
     * @param psm - PeptideSpectrumMatch with pointer to spectrum (@see PeptideSpectrumMatch)
     */
    int countGaps();

    /**
     * Returns unmodified version of the peptide.
     *
     * Note that this does not remove gaps, only modifications. (denoted by ())
     * @param outputPeptide = output unmodified annotation.
     */
    void getUnmodifiedPeptide(string &outputPeptide) const;

    /**
     * Inserts modifications into peptide
     *
     * @param modifications = array of modifications to put into peptide
     * @param positions = locations to insert modifications into
     * @param outputPeptide = output annotation
     * @return = returns false when modifications aren't able to be inserted.
     */
    bool insertModifications(const vector<float> &modifications,
                             const vector<unsigned int> &positions,
                             string &outputPeptide) const;
    /**
     * Inserts modifications into peptide
     *
     * @param unmodifiedPeptide = input string
     * @param modifications = array of modifications to put into peptide
     * @param positions = locations to insert modifications into
     * @param outputPeptide = output annotation
     * @return = returns false when modifications aren't able to be inserted.
     */
    static bool insertModifications(const string &unmodifiedPeptide,
                                    const vector<float> &modifications,
                                    const vector<unsigned int> &positions,
                                    string &outputPeptide);

    /**
     * Generates all theoretical masses from model and peptide sequence
     *
     * @param peptide amino acid sequence used to determine peak annotations
     * @param ionNamesInclude  comma delimited string indicating which ions with which names should be annotated. If string is simply "all"
     * then just include all ions in MS2Model. ex. "y,b,y++,b++"
     * @param inputIonTypes definition of ion types: mass offsets, ion probabilities, prefix/suffix.
     * is copied into ionTypes
     * @param prmOffset sets applied to theoretical ion masses to locate Prefix Residue Masses
     * @param srmOffset sets applied to theoretical ion masses to locate Suffix Residue Masses
     * @param peakTol - mass tolerance (in Da) when matching ion masses to peak masses
     * @param jumps - AAJumps indicating amino acid masses
     * @param ionNames - output ion names, i.e. b2, y1 etc.
     * @param theoreticalMasses - output theoretical masses for ionNames
     */
    static bool generateTheoreticalMasses(const string &peptide,
                                          const string &ionNamesInclude,
                                          const MS2ScoringModel &inputIonTypes,
                                          const float prmOffset,
                                          const float srmOffset,
                                          const AAJumps &jumps,
                                          vector<string> &ionNames,
                                          vector<float> &theoreticalMasses);

    static bool generateTheoreticalMasses(const string &peptide,
                                          const string &ionNamesInclude,
                                          const MS2ScoringModel &inputIonTypes,
                                          const float prmOffset,
                                          const float srmOffset,
                                          vector<string> &ionNames,
                                          vector<float> &theoreticalMasses);

    /**
     * Generates all theoretical internal fragment masses from model and peptide sequence including all
     * possible modified and unmodified fragments. Note that any modifications on peptide are ignored.
     *
     * @param unmodPeptide amino acid sequence used to determine peak annotations
     * @param ionNamesInclude  comma delimited string indicating which ions with which names should be annotated. If string is simply "all"
     * then just include all ions in MS2Model. ex. "y,b,y++,b++"
     * @param inputIonTypes definition of ion types: mass offsets, ion probabilities, prefix/suffix.
     * is copied into ionTypes
     * @param prmOffset sets applied to theoretical ion masses to locate Prefix Residue Masses
     * @param srmOffset sets applied to theoretical ion masses to locate Suffix Residue Masses
     * @param peakTol - mass tolerance (in Da) when matching ion masses to peak masses
     * @param jumps - AAJumps indicating amino acid masses
     * @param ionOutputMasses - output ion mass including fragment type
     * @param ionOutputFragments - ion fragments used in ionMass.
     */
    static bool
    generateTheoreticalMassesWithMods(const string &unmodPeptide,
                                      const string &ionNamesInclude,
                                      const MS2ScoringModel &inputIonTypes,
                                      const float prmOffset,
                                      const float srmOffset,
                                      const AAJumps &jumps,
                                      const vector<vector<float> > &massShifts,
                                      const vector<float> &fixedMods,
                                      const vector<unsigned int> &fixedPositions,
                                      vector<IonMass> &ionOutputMasses);

    bool loadFromFile(std::ifstream & ifs);
    static bool saveHeaderToFile(std::ofstream & ofs, bool variantCol = false);
    bool saveToFile(std::ofstream & ofs, bool variantCol = false);
    int getPeptideLength();
    
    string m_spectrumFile; //!original spectrum file from which annotation was loaded. Optional
    int m_scanNum; //! Scan number of spectrum.
    string m_annotation; //! Peptide annotation. Uses SpecNets format regardless of original input format.
    string m_origAnnotation; //! Original annotation, if original format has been changed.
    char   m_precedingAA;
    char   m_followingAA;
    /**
     * Annotation for each peak in peakList.
     * annotation[i] - contains the annotation for i-th spectrum peak
     * annotation[i].first is a pointer to a structure defining the type of ion (ftIonFragment *)
     * annotation[i].second is the index of the ion in that corresponding ion series (e.g., first b-ion has index 1)
     */
    vector<pair<const ftIonFragment*, short> > m_peakAnnotations;

    vector<TwoValues<int> > m_matchedPeaks;

    /**
     *  vector loaded in from MS2Model for annotation;
     */
    vector<ftIonFragment> m_ionTypes;
    string m_protein; //! Protien annotation.
    int m_dbIndex; //! Protein index in DB (if needed)
    int m_numMods; //! Number of modifications
    int m_matchOrientation; //! Forward (0) or reverse (1) tag match orientation
    float m_startMass; //! Starting mass of protein match
    float m_fdr;
    int m_charge; //! Peptide charge
    float m_score; //! associated score
    float m_compScore;//complementary score
    double m_pValue; //! associated p-value.
    bool m_isDecoy; //! indicates whether hit is decoy or not
    float m_strict_envelope_score; //! strict envelope score of MS1
    float m_unstrict_envelope_score; //! unstrict envelope score of MS1
    Spectrum * m_spectrum; //! associated spectrum match
    float m_mz; //! Theoretical MZ of the ID

    int m_variantGroup;           //! The variant group (unmod peptide + mass) that it belongs too
    string m_peptideRegionGroup;     //! The peptide variant region group number
    string m_peptideRegion;       //! The peptide variant region from the DB
    int m_peptideRegionLength; //! Length of the peptide variant region
    int m_startIndex; //! The starting index of the annotation in the proetin
    int m_endIndex; //! The ending index of the annotation in the proetin

    /**
     *  Below is data for Spectral Library Searches on CCMS to display the appropriate Library File
     **/
    string m_library_name;

    /**
     *  Below is meta data regarding the spectrum used for loading spectra
     *  into CCMS for creating spectral libraries
     */
    string m_submission_metadata;
    string m_organism;
    string m_compound_name;
    string m_smiles;
    string m_InChI;
    string m_InChI_Aux;

    int m_useYendPts; //! Whether or not the spectrum needs to be reversed (0/1) to match the maximum number of PRMs
    float m_pepFdr;

    string m_notes;
    string m_ionmode;
    float m_mz_error_ppm;
    float m_exactmass;
    float m_parentmass_difference;
    int m_shared_peaks;

    int m_abundance; //1 is low, 2 is high

    /**
     * 1 = Gold
     * 2 = Silver
     * 3 = Bronze
     **/
    int m_library_quality;

    string m_spectrumID;

    /**
     *  Below are fields related to Spectrum Spectrum Matches
     *  used for spectral library searches
     */
    vector<float> m_ion_extraction;

    /**
     *  SLGF
     */
    vector<pair<float, float> > SLGF_distribution;

  private:
    void internalCopy(const PeptideSpectrumMatch &other);
  };
}
#endif /* PEPTIDESPECTRUMMATCH_H_ */
