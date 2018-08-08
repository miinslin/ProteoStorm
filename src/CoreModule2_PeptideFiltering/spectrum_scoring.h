#ifndef SPECTRUM_SCORING_H
#define SPECTRUM_SCORING_H

#include <string>
#include <vector>
#include <list>
#include <cstring>
#include "spectrum.h"

namespace specnets
{
  using namespace std;

  /**
   * @see spectrum.h
   */
  class Spectrum;

  struct ftIonFragment
  {
    bool isNTerm, isIF; // true if ion is N-terminal (e.g. b-ion, b-H2O, etc), is internal fragment vs. modified total peptide.
    float massOffset, // mass offset (in Da) from the corresponding (prefix/suffix) sum of amino acid masses. Note that this includes the mass of charge!
        prob; // probability of finding this type of fragment in a spectrum
    string name; // name of the ion (e.g., "b", "y-H2O", "a-NH3-H2O", "y++", "b++", etc)
    int charge; //charge of ion type

    bool isMainIon; //
    bool isBreakIon; // if it is a break ion. Drawn in spectra images, on top, if so.
    bool hasLabel; // if the label is drwan
    string label; // label used to draw the peaks. If empty, uses 'name'


    ftIonFragment()
    {
      isNTerm = false;
      massOffset = 0;
      prob = 0;
      charge = 1;
      isIF = true;
      isMainIon = false;
      isBreakIon = false;
      hasLabel = false;
    }
  };

  struct ftCorrelatedIntensity
  {
    unsigned short fragIndex1, // Index of the first correlated fragment in a companion vector of ftIonFragment
        fragIndex2; // Index of the second correlated fragment in a companion vector of ftIonFragment
    float prob; // probability of seeing intensity(fragIndex1) >= intensity(fragIndex2)
    ftCorrelatedIntensity()
    {
      fragIndex1 = 0;
      fragIndex2 = 0;
      prob = 0;
    }
  };

  // NOTE: This structure is used differently when scoring MS3 spectra - fragIndex1 will correspond
  //       to an MS3 peak but fragIndex2 will correspond to a peak in the parent MS2 spectrum.
  struct ftCorrelatedOccurrence
  {
    unsigned short fragIndex1, // Index of the first correlated fragment in a companion vector of ftIonFragment
        fragIndex2; // Index of the second correlated fragment in a companion vector of ftIonFragment
    float probIfPresent, // probability of seeing fragIndex1 when fragIndex2 is present
        probIfAbsent; // probability of seeing fragIndex1 when fragIndex2 is not present
    ftCorrelatedOccurrence()
    {
      fragIndex1 = 0;
      fragIndex2 = 0;
      probIfPresent = 0;
      probIfAbsent = 0;
    }
  };

  class MS2ScoringModel
  {
  protected:
    virtual bool ParseParameter(char *parameterLine); // Parses scoring model parameters from a text string
    list<ftIonFragment> tmp_probs; // Temporary variable used by LoadModel()

  public:
    float peakTolPPM, peakTolDA; // Peak mass tolerance in PPMs or Daltons (PPM is used if >0)
    float noiseProb;
    vector<ftIonFragment> probs;
    list<ftCorrelatedIntensity> probs_ci;
    list<ftCorrelatedOccurrence> probs_co;

    MS2ScoringModel();
    virtual ~MS2ScoringModel()
    {
    }

    bool LoadModel(const char *filename);

    void getParentMassIonNames(vector<string> &parentIons) const;

    void getBreakIonsNames(vector<string> &bIons, vector<string> &yIons) const;

  };

  class MS2ScoringModelSet
  {
  protected:

    map<string, MS2ScoringModel> m_modelMap;

  public:
    MS2ScoringModelSet()
    {
      m_modelMap.clear();
    }

    inline void addModel(const string &fragType, MS2ScoringModel &model)
    {
      m_modelMap[fragType] = model;
    }

    inline const MS2ScoringModel& getModel(const string &fragType) const
    {
      return m_modelMap.at(fragType);
    }

    inline void removeModel(const string &fragType)
    {
      m_modelMap.erase(fragType);
    }

    inline bool containsModel(const string &fragType) const
    {
      return m_modelMap.count(fragType) > 0;
    }

    inline void clearModels()
    {
      m_modelMap.clear();
    }
  };

  class MS3ScoringModel : public MS2ScoringModel
  {
  protected:
    virtual bool ParseParameter(char *parameterLine);

  public:
    list<ftCorrelatedOccurrence> probs_co23; // Correlated peak occurrences between the MS2/MS3

  };

  /*! \brief Gets names of parent mass ions for further processing.

   @param parentIons vector to hold parent ions.
   */

  void ScoreSpectrum(Spectrum &in,
                     MS2ScoringModel &model,
                     bool removeNegativeScores = false,
                     Spectrum *out = NULL);
  void ScoreSpectrumMS3(Spectrum &ms3in,
                        Spectrum &ms2in,
                        MS3ScoringModel &model,
                        float ms3shift,
                        bool removeNegativeScores = false,
                        Spectrum *ms3out = NULL);
}
#endif
