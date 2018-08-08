#include "db_fasta.h"
#include "spectrum.h"
#include "label.h"
#include "Logger.h"
#include "IsoEnvelope.h"
#include "utils.h"
#include "alignment_scoring.h"
#include "AlignmentUtils.h"

#include <errno.h>
#include <string>
#include <map>

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <list>
#include <limits>
#include <stdio.h>
#include <stdlib.h>

const float MAX_SPEC_PROB_PEAK = 100.0;
#define DEBUG_SPEC_PROB  0
#define DEBUG_SPRINKLE 0

using namespace std;

namespace specnets
{

  const unsigned short Spectrum::BIN_VERSION = 1;
  const unsigned short Spectrum::BIN_SUBVERSION = 4;

  const string Spectrum::BIN_VERSION_ID = "Spectrum_binVersion";
  const string Spectrum::BIN_SUBVERSION_ID = "Spectrum_binSubVersion";

  Spectrum::FragType Spectrum::parseActivation(const string& activString)
  {
    string myCopy = activString;
    std::transform(myCopy.begin(), myCopy.end(), myCopy.begin(), ::toupper);

    if (myCopy == "CID")
    {
      return FragType_CID;
    }
    else if (myCopy == "HCD")
    {
      return FragType_HCD;
    }
    else if (myCopy == "PRM")
    {
      return FragType_PRM;
    }
    else if (myCopy == "PRM_ETD")
    {
      return FragType_PRM_ETD;
    }
    else if (myCopy.compare(0, 3, "ETD") == 0)
    {
      return FragType_ETD;
    }

    else
    {
      WARN_MSG("Found unknown fragmentation ID \'" << activString << "\', assuming CID");
      return FragType_CID;
    }
  }

  string Spectrum::activationToString(Spectrum::FragType& activ)
  {
    if (activ == FragType_CID)
    {
      return "CID";
    }
    else if (activ == FragType_HCD)
    {
      return "HCD";
    }
    else if (activ == FragType_ETD)
    {
      return "ETD";
    }
    else if (activ == FragType_PRM)
    {
      return "PRM";
    }
    else if (activ == FragType_PRM_ETD)
    {
      return "PRM_ETD";
    }
    else
    {
      ERROR_MSG("Found unknown fragmentation type \'" << activ << "\', returning empty string");
      return "";
    }
  }

  Spectrum::MassAnalyzerType Spectrum::parseMassAnalyzer(const string& analyzerString)
  {
    string myCopy = analyzerString;
    std::transform(myCopy.begin(), myCopy.end(), myCopy.begin(), ::tolower);

    if (myCopy.find("qtof") != string::npos
        || myCopy.find("quadrupole") != string::npos)
    {
      return MassAnalyzer_QTOF;
    }
    else if (myCopy.find("orbitrap") != string::npos
        || myCopy.find("ftms") != string::npos)
    {
      return MassAnalyzer_ORBI_TRAP;
    }
    else if (myCopy.find("ion trap") != string::npos)
    {
      return MassAnalyzer_ION_TRAP;
    }
    else
    {
      WARN_MSG("Found unknown mass analyzer type \'" << analyzerString << "\', assuming ion trap");
      return MassAnalyzer_ION_TRAP;
    }
  }

  string Spectrum::massAnalyzerToString(Spectrum::MassAnalyzerType& massAnalyzer)
  {
    if (massAnalyzer == MassAnalyzer_QTOF)
    {
      return "qtof";
    }
    else if (massAnalyzer == MassAnalyzer_ORBI_TRAP)
    {
      return "orbitrap";
    }
    else if (massAnalyzer == MassAnalyzer_ION_TRAP)
    {
      return "ion trap";
    }
    else
    {
      ERROR_MSG("Found unknown mass analyzer type \'" << massAnalyzer << "\', returning empty string");
      return "";
    }
  }

  // -------------------------------------------------------------------------
  Spectrum::Spectrum()
  {
    initialize();
  }
  // -------------------------------------------------------------------------
  Spectrum::~Spectrum()
  {
    //DEBUG_MSG("destructor");
    //	std::cerr << "destructor" << endl;

    list<psmPtr>::iterator it;
    //set Spectrum pointer to null for associated PSMs
    for (it = psmList.begin(); it != psmList.end(); it++)
    {
      (*it)->m_spectrum = (Spectrum *)NULL;
    }
  }
  // -------------------------------------------------------------------------
  Spectrum::Spectrum(const Spectrum &other)
  {
    initialize();
    internalCopy(other);
  }

  void Spectrum::initialize()
  {
    // DEBUG_MSG("constructor");
    //	std::cerr << "constructor" << endl;
    parentMass = 0;
    parentMassTol = 0;
    oldParentMassTol = 0;
    parentCharge = 0;
    parentMZ = 0;
    scan = 0;
    msLevel = 0;
    msFragType = FragType_CID;
    msMassAnalyzerType = MassAnalyzer_ION_TRAP;
    peakList.resize(0);
    psmList.resize(0);
    peakTols.resize(0);
    oldPeakTols.resize(0);
    m_gfProbs.resize(0);
    resolution = 1;
    idDist = 1;
    instrument_name = "";
    ITOL = -1.0;
    ITOLU = "";
    TOL = -1.0;
    TOLU = "";
    spectrum_quality = -1;
    fileName = "";
    fileIndex = -1;
    retention_time = 0.0;
    precursor_intensity = 0.0;
    collision_energy = 0.0;
    m_reversed = false;
    precursor_kl = 0.0;
    precursor_scan = 0;
  }

  void Spectrum::swap(Spectrum& other)
  {
    Spectrum temp;
    temp.copyNP(other);
    other.copyNP(*this);
    this->copyNP(temp);
    peakList.swap(other.peakList);
    peakTols.swap(other.peakTols);
    oldPeakTols.swap(other.oldPeakTols);
    psmList.swap(other.psmList);
    m_gfProbs.swap(other.m_gfProbs);
  }

  // -------------------------------------------------------------------------
  Spectrum &Spectrum::operator=(const Spectrum &other)
  {
    //	std::cerr << "assignment" << endl;

    if (this == &other)
      return *this;
    internalCopy(other);
    return (*this);
  }





  // -------------------------------------------------------------------------
  void Spectrum::push_back(const TwoValues<float> & x)
  {
    peakList.push_back(x);
    resize(peakList.size());
  }
  // -------------------------------------------------------------------------
  Spectrum &Spectrum::copyNP(const Spectrum &other)
  {
    //	std::cerr << "copynp" << endl;
    parentMass = other.parentMass;
    parentMassTol = other.parentMassTol;
    parentCharge = other.parentCharge;
    parentMZ = other.parentMZ;
    scan = other.scan;
    collision_energy = other.collision_energy;
    msMassAnalyzerType = other.msMassAnalyzerType;
    msLevel = other.msLevel;
    msFragType = other.msFragType;
    resolution = other.resolution;
    idDist = other.idDist;
    instrument_name = other.instrument_name;
    ITOL = other.ITOL;
    ITOLU = other.ITOLU;
    TOL = other.TOL;
    TOLU = other.TOLU;
    spectrum_quality = other.spectrum_quality;
    fileName = other.fileName;
    fileIndex = other.fileIndex;
    m_reversed = other.m_reversed;
    precursor_intensity = other.precursor_intensity;
    retention_time = other.retention_time;
    precursor_kl = other.precursor_kl;
    precursor_scan = other.precursor_scan;

    return (*this);
  }

  // -------------------------------------------------------------------------
  void Spectrum::copyWithShallowPsms(const Spectrum &other)
  {
    copyNP(other);

    peakList.resize(other.peakList.size());
    for (unsigned int i = 0; i < other.peakList.size(); i++)
      peakList[i] = other.peakList[i];
    peakTols.resize(other.peakTols.size());
    for (unsigned int i = 0; i < other.peakTols.size(); i++)
    {
      peakTols[i] = other.peakTols[i];
    }
    oldPeakTols = other.oldPeakTols;

    //copy psms
    psmList.clear();
    for (list<psmPtr>::const_iterator it = other.psmList.begin();
        it != other.psmList.end(); it++)
    {
      psmList.push_back(*it);
    }

    //copy spectral probabilities
    m_gfProbs = other.m_gfProbs;

    return;
  }

  // -------------------------------------------------------------------------
  void Spectrum::internalCopy(const Spectrum &other)
  {
    copyNP(other);

    peakList.resize(other.peakList.size());
    for (unsigned int i = 0; i < other.peakList.size(); i++)
      peakList[i] = other.peakList[i];
    peakTols.resize(other.peakTols.size());
    for (unsigned int i = 0; i < other.peakTols.size(); i++)
    {
      peakTols[i] = other.peakTols[i];
    }
    oldPeakTols = other.oldPeakTols;

    //copy psms
    psmList.clear();
    for (list<psmPtr>::const_iterator it = other.psmList.begin();
        it != other.psmList.end(); it++)
    {
      psmPtr newPsm(new PeptideSpectrumMatch);
      *newPsm = **it;
      psmList.push_back(newPsm);
    }

    //copy spectral probabilities
    m_gfProbs = other.m_gfProbs;

    return;
  }

  // -------------------------------------------------------------------------
  inline void Spectrum::weightedAverage(int minIdx,
                                        int maxIdx,
                                        float &avgMass,
                                        float &totScore)
  {
    int i = 0;
    totScore = 0;
    avgMass = 0;
    for (i = minIdx; i <= maxIdx; i++)
      totScore += peakList[i][1];
    if (totScore <= 0)
      for (i = minIdx; i <= maxIdx; i++)
        avgMass += peakList[i][0] / (maxIdx - minIdx + 1); // regular average
    else
      for (i = minIdx; i <= maxIdx; i++)
        avgMass += peakList[i][0] * peakList[i][1] / totScore; // weighted average
  }

  void Spectrum::setResolution(float newResolution, bool enforceRounding)
  {
    resolution = newResolution;
    idDist = 1 / resolution;
    parentMass =
        enforceRounding ? round(parentMass / resolution) :
            parentMass / resolution;
    if (!enforceRounding)
      for (unsigned int i = 0; i < peakList.size(); i++)
        peakList[i][0] = peakList[i][0] / resolution;
    else
    {
      unsigned int j = 0;
      for (unsigned int i = 0; i < peakList.size(); i++)
      {
        peakList[i][0] = round(peakList[i][0] / resolution);
        if (j < i)
        {
          if (abs(peakList[j][0] - peakList[i][0]) < 0.0001)
            peakList[j][1] += peakList[i][1];
          else
          {
            j++;
            peakList[j] = peakList[i];
            peakTols[j] = peakTols[i];
          }
        }
      }
      resize(j + 1);
    }
  }
  // -------------------------------------------------------------------------
  /**
   * Sets the tolerance of each peak in Da
   */
  void Spectrum::setPeakTolerance(float peakTol)
  {
    peakTols.resize(size());
    for (int i = 0; i < size(); i++)
    {
      peakTols[i] = peakTol;
    }
  }
  // -------------------------------------------------------------------------
  /**
   * Sets the ppm tolerance of each peak
   */
  void Spectrum::setPeakTolerancePPM(float peakTolPPM)
  {
    for (int i = 0; i < size(); i++)
    {
      peakTols[i] = peakList[i][0] * peakTolPPM * PPM_FACTOR;
    }
  }
  // -------------------------------------------------------------------------
  TwoValues<float>* Spectrum::front()
  {
    return &peakList.front();
  }
  // -------------------------------------------------------------------------
  TwoValues<float>* Spectrum::back()
  {
    return &peakList.back();
  }
  // -------------------------------------------------------------------------
  /**
   * Sets the peak tolerance of the parent mass
   */
  void Spectrum::setParentMassTol(float pmTol)
  {
    parentMassTol = pmTol;
  }
  // -------------------------------------------------------------------------
  /**
   * Sets the ppm tolerance of the parent mass
   */
  void Spectrum::setParentMassTolPPM(float pmTolPPM)
  {
    parentMassTol = (parentMZ * pmTolPPM * PPM_FACTOR) * ((float)parentCharge);
  }
  // -------------------------------------------------------------------------
  void Spectrum::setCharge(short charge)
  {
    float newChargeF = (float)charge;
    float prevTol = parentMassTol / newChargeF;
    parentCharge = charge;
    parentMZ = (parentMass + ((newChargeF - 1.0) * AAJumps::massHion))
        / newChargeF;
    parentMassTol = prevTol * newChargeF;
  }

  void Spectrum::setFilename(const string &specsFilename)
  {
    FilenameManager fm(specsFilename);
    fileName = fm.getFilenameWithExtension();
  }

  /**
   * Sets the parentMass and parentMZ given the parent charge and input parent mass
   */
  void Spectrum::setParentMass(float pm)
  {
    parentMass = pm;
    float chargeF = (float)parentCharge;
    parentMZ = (parentMass + ((chargeF - 1.0) * AAJumps::massHion)) / chargeF;
    if (parentCharge == 0)
    {
      parentMZ = pm;
    }
  }

  void Spectrum::setParentMZ(float mz)
  {
    parentMZ = mz;
    float chargeF = (float)parentCharge;
    parentMass = (parentMZ * chargeF) - (AAJumps::massHion * (chargeF - 1.0));
    //parentMass = (parentMZ * chargeF) -  chargeF ;
    if (parentCharge == 0)
    {
      parentMass = mz;
    }
  }


  // -------------------------------------------------------------------------
  /**
   * Copies all per-peak tolerances into an internal vector
   */
  void Spectrum::rememberTolerances()
  {
    oldPeakTols = peakTols;
    oldParentMassTol = parentMassTol;
  }
  // -------------------------------------------------------------------------
  /**
   * Restores peak tolerances to as they were before last call to rememberTolerances()
   */
  void Spectrum::revertTolerances()
  {
    if (oldPeakTols.size() == peakTols.size())
    {
      peakTols = oldPeakTols;
      parentMassTol = oldParentMassTol;
    }
    else
    {
      WARN_MSG("Could not revert peak tolerances because peak list size has changed!");
      return;
    }
  }

  void Spectrum::roundPeaks(float roundFactor, bool roundPM, bool roundToInts)
  {
    for (unsigned int i = 0; i < peakList.size(); i++)
    {
      peakList[i][0] *= roundFactor;
      if (roundToInts)
      {
        peakList[i][0] = round(peakList[i][0]);
      }
      peakTols[i] *= roundFactor;
    }
    if (roundPM)
    {
      parentMass *= roundFactor;
      if (roundToInts)
      {
        parentMass = round(parentMass);
      }
      parentMassTol *= roundFactor;
    }
  }
  // -------------------------------------------------------------------------
  float Spectrum::getTolerance(int index) const
  {
    return peakTols[index];
  }
  // -------------------------------------------------------------------------
  void Spectrum::setTolerance(int index, float tolerance)
  {
    peakTols[index] = tolerance;
  }
  // -------------------------------------------------------------------------
  void Spectrum::sortPeaks()
  {
    map<float, float> peakTolMap;

    for (int i = 0; i < peakList.size(); i++)
    {
      peakTolMap[peakList[i][0]] = peakTols[i];
    }

    sort(peakList.begin(), peakList.end());

    for (int i = 0; i < peakList.size(); i++)
    {
      peakTols[i] = peakTolMap[peakList[i][0]];
    }
  }
  // -------------------------------------------------------------------------
  int Spectrum::insertPeak(float mass,
                           float intensity,
                           float tolerance,
                           int atIndex)
  {
    vector<TwoValues<float> >::iterator peakIt = peakList.begin();
    vector<float>::iterator tolIt = peakTols.begin();
    int insertIdx = 0;

    if (atIndex >= 0)
    {
      if (atIndex >= peakList.size() || atIndex >= peakTols.size())
      {
        atIndex = min(peakList.size(), peakTols.size()) - 1;
      }
      peakIt += atIndex;
      tolIt += atIndex;
      insertIdx = atIndex;
    }
    else
    {
      if (size() > 0)
      {
        insertIdx = findClosest(mass);
        if (peakList[insertIdx][0] < mass)
        {
          insertIdx++;
        }
      }
      else
      {
        insertIdx = 0;
      }
      peakIt += insertIdx;
      tolIt += insertIdx;
    }

    peakList.insert(peakIt, TwoValues<float>(mass, intensity));
    peakTols.insert(tolIt, tolerance);

    return insertIdx;
  }
  // -------------------------------------------------------------------------
  int Spectrum::insertPeak(MZRange* peak, int atIndex)
  {
    return insertPeak(peak->getMass(),
                      peak->getIntensity(),
                      peak->getTolerance(),
                      atIndex);
  }
  // -------------------------------------------------------------------------
  int Spectrum::insertPeaks(list<MZRange>& newPeaks)
  {

    int idxUse = size();
    resize(size() + newPeaks.size());
    for (list<MZRange>::iterator peakIt = newPeaks.begin();
        peakIt != newPeaks.end(); peakIt++)
    {
      peakList[idxUse][0] = peakIt->getMass();
      peakList[idxUse][1] = peakIt->getIntensity();
      peakTols[idxUse] = peakIt->getTolerance();
      ++idxUse;
    }
    sortPeaks();
  }
  // -------------------------------------------------------------------------
  int Spectrum::findPeaks(float mass,
                          float tolerance,
                          list<int>* matches,
                          int* startIdx) const
  {
    if (size() == 0)
    {
      return -1;
    }
    int idx;

    if (startIdx == (int*)0)
    {
      idx = findClosest(mass);
    }
    else
    {
      idx = *startIdx;
      if (idx < 0)
      {
        idx = 0;
      }
      while (idx < size() && peakList[idx][0] < mass)
      {
        ++idx;
      }
      if (idx > 0 && peakList[idx][0] > mass && peakList[idx - 1][0] < mass
          && peakList[idx][0] - mass > mass - peakList[idx - 1][0])
      {
        --idx;
      }
    }

    if (idx >= size() || idx < 0)
    {
      return -1;
    }

    if (MZRange::EqualWithinRange(peakList[idx][0],
                                  mass,
                                  peakTols[idx] + tolerance))
    {
      if (matches != 0)
      {
        matches->clear();
        matches->push_back(idx);
        float lowerBound = mass - tolerance;
        for (int i = idx - 1;
            i >= 0 && peakList[i][0] + peakTols[i] >= lowerBound; i--)
        {
          matches->push_back(i);
        }
        float upperBound = mass + tolerance;
        for (int i = idx + 1;
            i < size() && peakList[i][0] - peakTols[i] <= upperBound; i++)
        {
          matches->push_back(i);
        }
        matches->sort();
      }
      return idx;
    }
    else
    {
      if (matches != 0)
      {
        matches->clear();
      }
      return -1;
    }
  }
  // -------------------------------------------------------------------------
  int Spectrum::findPeaks(const MZRange& range,
                          list<int>* matches,
                          int* startIdx) const
  {
    return findPeaks(range.getMass(), range.getTolerance(), matches, startIdx);
  }
  // -------------------------------------------------------------------------
  bool Spectrum::removePeak(int peakIndex)
  {
    if (peakIndex < 0 || peakIndex >= peakList.size()
        || peakIndex >= peakTols.size())
    {
      WARN_MSG("Invalid peak index " << peakIndex);
      return false;
    }
    vector<TwoValues<float> >::iterator peakIt = peakList.begin();
    vector<float>::iterator tolIt = peakTols.begin();

    peakList.erase(peakIt + peakIndex);
    peakTols.erase(tolIt + peakIndex);

    return true;
  }
  // -------------------------------------------------------------------------
  bool Spectrum::removePeaks(list<int>& peakIndices)
  {

    if (peakIndices.size() == 0)
    {
      return true;
    }
    if (size() < peakIndices.size())
    {
      WARN_MSG("Remove peak list size [" << peakIndices.size() <<
               "] is greater than number of peaks [" << size() << "]");
    }
    vector<bool> removed(size(), false);

    for (list<int>::iterator peakIt = peakIndices.begin();
        peakIt != peakIndices.end(); peakIt++)
    {
      if (*peakIt < 0)
      {
        ERROR_MSG("Peak index below zero [" << *peakIt << "] in removePeaks");
        return false;
      } else  if (*peakIt >= size()) {
        WARN_MSG("Peak index out of bounds in removePeaks " <<
                 *peakIt << "[ size is ]" << size() << "]");
      } else {
        removed[*peakIt] = true;
      }
    }

    vector<TwoValues<float> > newPeakList(size());
    vector<float> newTols(size());
    int idxUse = 0;
    for (int i = 0; i < removed.size(); i++)
    {
      if (!removed[i])
      {
        newPeakList[idxUse] = peakList[i];
        newTols[idxUse] = peakTols[i];
        ++idxUse;
      }
    }
    newPeakList.resize(idxUse);
    newTols.resize(idxUse);
    resize(idxUse);

    peakList.assign(newPeakList.begin(), newPeakList.end());
    peakTols.assign(newTols.begin(), newTols.end());

    return true;
  }
  // -------------------------------------------------------------------------
  float Spectrum::getPeakDensity()
  {
    return ((float)size()) / parentMass;
  }
  // -------------------------------------------------------------------------
  void Spectrum::addZPMpeaks(float tolerance,
                             float ionOffset,
                             bool includeY0k,
                             bool ctermH2O,
                             SpectrumPeakLabels *labels)
  {
    float massB0 = ionOffset * idDist, massBk = parentMass
        - ((ctermH2O ? AAJumps::massMH : AAJumps::massHion) - ionOffset) * idDist,
        massY0 = ((ctermH2O ? AAJumps::massH2O : 0) + ionOffset) * idDist,
        massYk = parentMass - (AAJumps::massHion - ionOffset) * idDist;

    float endptIntensity = 0.1;

    bool modTols = false;

    if (tolerance >= 0)
    {
      modTols = true;
      rememberTolerances();
      setPeakTolerance(tolerance);
    }

    if (findPeaks(massB0) < 0)
    {
      int idx = insertPeak(massB0, endptIntensity, 0);
      if (modTols)
      {
        vector<float>::iterator peakIt = oldPeakTols.begin();
        peakIt += idx;
        oldPeakTols.insert(peakIt, 0);
      }
      if (labels != 0)
      {
        vector<PeakLabel>::iterator labelIt = labels->peakLabels.begin();
        labelIt += idx;
        PeakLabel newLabel(0, 0, 0, 0, 1);
        labels->peakLabels.insert(labelIt, newLabel);
      }
    }

    if (findPeaks(massBk) < 0)
    {
      int idx = insertPeak(massBk, endptIntensity, 0);
      if (modTols)
      {
        vector<float>::iterator peakIt = oldPeakTols.begin();
        peakIt += idx;
        oldPeakTols.insert(peakIt, 0);
      }
      if (labels != 0)
      {
        vector<PeakLabel>::iterator labelIt = labels->peakLabels.begin();
        labelIt += idx;
        PeakLabel newLabel(0, 0, 0, 0, 1);
        labels->peakLabels.insert(labelIt, newLabel);
      }
    }

    if (includeY0k && findPeaks(massY0) < 0)
    {
      int idx = insertPeak(massY0, endptIntensity, 0);
      if (modTols)
      {
        vector<float>::iterator peakIt = oldPeakTols.begin();
        peakIt += idx;
        oldPeakTols.insert(peakIt, 0);
      }
      if (labels != 0)
      {
        vector<PeakLabel>::iterator labelIt = labels->peakLabels.begin();
        labelIt += idx;
        PeakLabel newLabel(0, 0, 0, 0, 1);
        labels->peakLabels.insert(labelIt, newLabel);
      }
    }

    if (includeY0k && findPeaks(massYk) < 0)
    {
      int idx = insertPeak(massYk, endptIntensity, 0);
      if (modTols)
      {
        vector<float>::iterator peakIt = oldPeakTols.begin();
        peakIt += idx;
        oldPeakTols.insert(peakIt, 0);
      }
      if (labels != 0)
      {
        vector<PeakLabel>::iterator labelIt = labels->peakLabels.begin();
        labelIt += idx;
        PeakLabel newLabel(0, 0, 0, 0, 1);
        labels->peakLabels.insert(labelIt, newLabel);
      }
    }

    if (modTols)
    {
      revertTolerances();
    }
  }

  // -------------------------------------------------------------------------
  void Spectrum::addZPMpeaks22(float tolerance,
                               float ionOffset,
                               bool includeY0k,
                               float endptIntensity,
                               bool ctermH2O,
                               SpectrumPeakLabels *labels)
  {
	float mH2O = ctermH2O? AAJumps::massH2O : 0;
    float massB0 = ionOffset * idDist,
    	  massBk = parentMass + (ionOffset-mH2O-AAJumps::massHion) * idDist,
          massY0 = (ionOffset + mH2O) * idDist,
          massYk = parentMass + (ionOffset-AAJumps::massHion) * idDist;

    bool modTols = false;

    if (tolerance >= 0)
    {
      modTols = true;
      rememberTolerances();
      setPeakTolerance(tolerance);
    }

    int findedIndex = findPeaks(massB0);
    if (findedIndex < 0)
    {
      int idx = insertPeak(massB0, endptIntensity, 0);
      if (modTols)
      {
        vector<float>::iterator peakIt = oldPeakTols.begin();
        peakIt += idx;
        oldPeakTols.insert(peakIt, 0);
      }
      if (labels != 0)
      {
        vector<PeakLabel>::iterator labelIt = labels->peakLabels.begin();
        labelIt += idx;
        PeakLabel newLabel(0, 0, 0, 0, 1);
        labels->peakLabels.insert(labelIt, newLabel);
      }
    }
    else {
    	peakList[findedIndex][1] = endptIntensity;
    }

    findedIndex = findPeaks(massBk);
    if (findedIndex < 0)
    {
      int idx = insertPeak(massBk, endptIntensity, 0);
      if (modTols)
      {
        vector<float>::iterator peakIt = oldPeakTols.begin();
        peakIt += idx;
        oldPeakTols.insert(peakIt, 0);
      }
      if (labels != 0)
      {
        vector<PeakLabel>::iterator labelIt = labels->peakLabels.begin();
        labelIt += idx;
        PeakLabel newLabel(0, 0, 0, 0, 1);
        labels->peakLabels.insert(labelIt, newLabel);
      }
    }
    else {
    	peakList[findedIndex][1] = endptIntensity;
    }

    if (includeY0k){
		findedIndex = findPeaks(massY0);
		if (findedIndex < 0)
		{
		  int idx = insertPeak(massY0, endptIntensity, 0);
		  if (modTols)
		  {
			vector<float>::iterator peakIt = oldPeakTols.begin();
			peakIt += idx;
			oldPeakTols.insert(peakIt, 0);
		  }
		  if (labels != 0)
		  {
			vector<PeakLabel>::iterator labelIt = labels->peakLabels.begin();
			labelIt += idx;
			PeakLabel newLabel(0, 0, 0, 0, 1);
			labels->peakLabels.insert(labelIt, newLabel);
		  }
		}
		else {
			peakList[findedIndex][1] = endptIntensity;
		}

		findedIndex = findPeaks(massYk);
		if (findedIndex < 0)
		{
		  int idx = insertPeak(massYk, endptIntensity, 0);
		  if (modTols)
		  {
			vector<float>::iterator peakIt = oldPeakTols.begin();
			peakIt += idx;
			oldPeakTols.insert(peakIt, 0);
		  }
		  if (labels != 0)
		  {
			vector<PeakLabel>::iterator labelIt = labels->peakLabels.begin();
			labelIt += idx;
			PeakLabel newLabel(0, 0, 0, 0, 1);
			labels->peakLabels.insert(labelIt, newLabel);
		  }
		}
		else {
			peakList[findedIndex][1] = endptIntensity;
		}
    }

    if (modTols)
    {
      revertTolerances();
    }
  }


  void Spectrum::removeZPMpeaks()
  {
    list<int> peaksToRemove;
    for (int i = 0; i < peakList.size(); i++)
    {
      if (peakList[i][0] + 1.0 < AAJumps::minAAmass
          || peakList[i][0] > parentMass - 37.0)
      {
        peaksToRemove.push_back(i);
      }
    }

    removePeaks(peaksToRemove);
  }

  // -------------------------------------------------------------------------
  void Spectrum::maximizeZPMpeaks(float tolerance,
                                  float ionOffset,
                                  bool includeY0k)
  {

    bool modTols = false;

    if (tolerance >= 0)
    {
      modTols = true;
      rememberTolerances();
      setPeakTolerance(tolerance);
    }

    float massB0 = ionOffset * idDist, massBk = parentMass
        - (AAJumps::massMH - ionOffset) * idDist, massY0 = (AAJumps::massH2O
        + ionOffset) * idDist, massYk = parentMass
        - (AAJumps::massHion - ionOffset) * idDist;

    unsigned int pivot;
    float maxScore = 0;
    for (pivot = 0; pivot < peakList.size(); pivot++)
      if (peakList[pivot][1] > maxScore)
        maxScore = peakList[pivot][1];
    for (pivot = 0;
        pivot < peakList.size()
            and peakList[pivot][0] <= massB0 + getTolerance(pivot); pivot++)
      peakList[pivot][1] = maxScore;
    if (includeY0k)
    {
      for (;
          pivot < peakList.size()
              and peakList[pivot][0] < massY0 - getTolerance(pivot); pivot++)
        ;
      for (;
          pivot < peakList.size()
              and peakList[pivot][0] <= massY0 + getTolerance(pivot); pivot++)
        peakList[pivot][1] = maxScore;
    }
    for (; pivot < peakList.size() and peakList[pivot][0] < massBk; pivot++)
      ;
    for (; pivot < peakList.size() and peakList[pivot][0] <= massBk; pivot++)
      peakList[pivot][1] = maxScore;
    if (includeY0k)
    {
      for (;
          pivot < peakList.size()
              and peakList[pivot][0] < massYk - getTolerance(pivot); pivot++)
        ;
      for (;
          pivot < peakList.size()
              and peakList[pivot][0] <= massYk + getTolerance(pivot); pivot++)
        peakList[pivot][1] = maxScore;
    }

    if (modTols)
    {
      revertTolerances();
    }
  }
  // -------------------------------------------------------------------------
  void Spectrum::normalize(float newTotalInt, bool removeNegatives)
  { // Normalizes total intensity to 100
    float totalIntensity = 0;
    for (unsigned int i = 0; i < peakList.size(); i++)
      totalIntensity += peakList[i][1];
    for (unsigned int i = 0; i < peakList.size(); i++)
      peakList[i][1] = newTotalInt * peakList[i][1] / totalIntensity;
  }
  // -------------------------------------------------------------------------
  void Spectrum::normalize2(float newNorm2)
  { // Normalizes to Euclidian norm newNorm2
    float factor = 0;
    for (unsigned int i = 0; i < peakList.size(); i++)
      factor += peakList[i][1] * peakList[i][1];
    factor = sqrt(factor) / newNorm2;
    for (unsigned int i = 0; i < peakList.size(); i++)
      peakList[i][1] = peakList[i][1] / factor;
  }

  // -------------------------------------------------------------------------
  float Spectrum::cosine(Spectrum &withSpec,
                         vector<pair<unsigned int, unsigned int> > *matchedPeaksIdx)
  {
    // Greedy calculation of highest-possible cosine (maximum bipartite heuristic)
    /*				FindMatchPeaksAll2(*this,
     withSpec,
     0,
     peakTol,
     idxMatched1_zero,
     idxMatched2_zero);

     FindMatchPeaksAll2(*this,
     withSpec,
     pmDiff,
     peakTol,
     idxMatched1_other,
     idxMatched2_other);

     //cerr<<"Got "<<idxMatched1_zero.size()<<"/"<<idxMatched2_zero.size()<<" (at 0) + "<<idxMatched1_other.size()<<"/"<<idxMatched2_other.size()<<" (at "<<pmDiff<<") matched peaks\n";

     if( idxMatched1_zero.size() + idxMatched1_other.size() < minNumMatchedPeaks )
     continue;
     peakMatches.resize( idxMatched1_zero.size() + idxMatched1_other.size() );
     for(unsigned int i=0; i < this->size(); i++) peakUsed1[i]=0;
     for(unsigned int i=0; i < withSpec.size(); i++) peakUsed2[i]=0;
     for(unsigned int i=0; i < idxMatched1_zero.size(); i++) {
     peakMatches[i].i = idxMatched1_zero[i];
     peakMatches[i].j = idxMatched2_zero[i];
     peakMatches[i].score = (*this)[idxMatched1_zero[i]][1] * withSpec[idxMatched2_zero[i]][1];
     }
     for(unsigned int i=0, j=idxMatched1_zero.size(); i < idxMatched1_other.size(); i++, j++) {
     peakMatches[j].i = idxMatched1_other[i];
     peakMatches[j].j = idxMatched2_other[i];
     peakMatches[j].score = (*this)[idxMatched1_other[i]][1] * withSpec[idxMatched2_other[i]][1];
     }
     std::sort(peakMatches.begin(), peakMatches.end(), GPCAux_cmp);

     cosine = 0;
     score1 = 0;
     score2 = 0;
     unsigned int numMatchedPeaks = 0;
     for(int i=(int)peakMatches.size()-1; i >= 0; i--) {
     if( peakUsed1[ peakMatches[i].i ]==0 and peakUsed2[ peakMatches[i].j ]==0 ) {
     cosine += peakMatches[i].score;
     peakUsed1[ peakMatches[i].i ] = 1;
     peakUsed2[ peakMatches[i].j ] = 1;
     score1 += (*this)[peakMatches[i].i][1];
     score2 += withSpec[peakMatches[i].j][1];
     numMatchedPeaks++;
     //cerr<<" --- used ("<<peakMatches[i].i<<","<<specSet[spec1][peakMatches[i].i][0]<<") / ("<<peakMatches[i].j<<","<<specSet[spec2][peakMatches[i].j][0]<<")\n";
     } else {
     //cerr<<" --- skipped ("<<peakMatches[i].i<<","<<specSet[spec1][peakMatches[i].i][0]<<") / ("<<peakMatches[i].j<<","<<specSet[spec2][peakMatches[i].j][0]<<")\n";
     }
     }
     //cerr<<"Got "<<numMatchedPeaks<<" matched cosine peaks with cosine "<<cosine<<" and matched intensities "<<score1<<"/"<<score2<<"\n";

     return cosine; */
  }

  // -------------------------------------------------------------------------
  void Spectrum::guessPrecursorZPM(Spectrum &parentMS1,
                                   float peakTol,
                                   short maxZ,
                                   IsoEnvelope &isoEnvs,
                                   bool strictMode)
  {

    WARN_MSG("guessPrecursorZPM has not been updated for per-peak tolerances");

    // Look for the best isotopic envelope match over all possible charges and monoisotopic masses
    float curMonoMass, chargeIncrement, bestPM = 0, curMatch, bestMatch =
        1000000; // Match score of zero is optimal
    short bestZ = 0;
    vector<float> massEnvelope;
    //cerr<<"Processing scan "<<scan<<endl;
    for (short charge = 1; charge <= maxZ; charge++)
    {
      chargeIncrement = 1.0 / ((float)charge);
      for (short massOffset = -2; massOffset <= 0; massOffset++)
      {
        curMonoMass = parentMass + ((float)massOffset) * chargeIncrement;
        isoEnvs.ExtractEnvelope(curMonoMass,
                                charge,
                                parentMS1,
                                peakTol,
                                massEnvelope,
                                strictMode);
        //cerr<<"massEnvelope ["<<charge<<","<<massOffset<<"]: "; for(unsigned int i=0; i<massEnvelope.size(); i++) cerr<<massEnvelope[i]<<" "; cerr<<"\n"; cerr.flush();
        isoEnvs.normalize(massEnvelope);
        //cerr<<"normalized massEnvelope: "; for(unsigned int i=0; i<massEnvelope.size(); i++) cerr<<massEnvelope[i]<<" "; cerr<<"\n"; cerr.flush();
        curMatch = isoEnvs.ScoreEnvelope(curMonoMass, massEnvelope, strictMode);
        //cerr<<"score = "<<curMatch<<endl;
        if (curMatch < bestMatch)
        {
          bestMatch = curMatch;
          bestPM = curMonoMass * charge - charge + 1;
          bestZ = charge;
        }
      }
    }
    //cerr<<" --> Chose parentCharge = "<<bestZ<<", parentMass = "<<bestPM<<" (was "<<parentMass<<")\n";
    parentCharge = bestZ;
    parentMass = bestPM;
  }
  // -------------------------------------------------------------------------
  short Spectrum::findMatches(float baseMass,
                              float peakTol,
                              vector<int> &matchesIdx,
                              int startIdx)
  {
    vector<float> tempTols;
    bool modTols = false;

    if (peakTol >= 0)
    {
      modTols = true;
      rememberTolerances();
      setPeakTolerance(peakTol);
    }
    int* startPtr = &startIdx;
    if (startIdx < 0)
    {
      startPtr = (int*)0;
    }

    list<int> matchesList;
    int idx = findPeaks(baseMass, 0, &matchesList, startPtr);
    matchesIdx.resize(matchesList.size());
    int idxUse = 0;
    for (list<int>::iterator mIt = matchesList.begin();
        mIt != matchesList.end(); mIt++)
    {
      if (*mIt >= startIdx)
      {
        matchesIdx[idxUse] = *mIt;
        idxUse++;
      }
    }
    matchesIdx.resize(idxUse);

    if (modTols)
    {
      revertTolerances();
    }

    return (short)idxUse;
    /*
     unsigned int baseIdx = (startIdx >= 0 and startIdx < peakList.size())
     ? (unsigned int)startIdx : 0;
     unsigned int numPeaks = peakList.size();
     unsigned int matchCount = 0;
     matchesIdx.resize(numPeaks);
     for (; baseIdx > 0 and peakList[baseIdx][0] >= baseMass - peakTol; baseIdx--)
     ;
     for (; baseIdx < numPeaks and peakList[baseIdx][0] < baseMass - peakTol; baseIdx++)
     ; // Get to first peak above the lower end of the tolerance window
     for (; baseIdx < numPeaks and peakList[baseIdx][0] <= baseMass + peakTol; baseIdx++)
     matchesIdx[matchCount++] = baseIdx;
     matchesIdx.resize(matchCount);
     return matchCount;
     */
  }
  // -------------------------------------------------------------------------
  int Spectrum::findClosest(float mass) const
  {
    if (peakList.size() == 0)
      return -1;

    int lowerLim = 0, upperLim = (int)peakList.size() - 1, middle;
    while (upperLim > lowerLim + 1)
    { // Set lowerLim to the index of the largest monoisotopic mass < mass
      middle = (lowerLim + upperLim) / 2;
      if (peakList[middle][0] > mass)
        upperLim = middle;
      else
        lowerLim = middle;
    }
    if (fabs(peakList[lowerLim][0] - mass) < fabs(peakList[upperLim][0] - mass))
      return lowerLim;
    else
      return upperLim;
  }
  // -------------------------------------------------------------------------
  void Spectrum::mergePeakList(Spectrum &spec, Spectrum *putHere)
  {
    int idxUpdated, idxOld, idxNew;
    vector<TwoValues<float> > &updatedPeakList =
        (putHere == 0) ? peakList : putHere->peakList;
    vector<float> &updatedPeakTols =
        (putHere == 0) ? peakTols : putHere->peakTols;
    vector<TwoValues<float> > *oldPeakList;
    vector<float> *oldPeakTols;

    if (putHere == 0)
    {
      oldPeakList = new vector<TwoValues<float> >;
      oldPeakList->resize(peakList.size());
      for (int i = 0; i < peakList.size(); i++)
        (*oldPeakList)[i] = peakList[i];

      oldPeakTols = new vector<float>;
      oldPeakTols->resize(peakTols.size());
      for (int i = 0; i < peakTols.size(); i++)
        (*oldPeakTols)[i] = peakTols[i];
    }
    else
    {
      oldPeakList = &peakList; // No copying necessary if results goes into another object (putHere)
      oldPeakTols = &peakTols;
      putHere->copyNP(*this);
    }

    // Sorted merge
    updatedPeakList.resize(oldPeakList->size() + spec.size());
    updatedPeakTols.resize(oldPeakList->size() + spec.size());
    idxUpdated = 0;
    idxOld = 0;
    idxNew = 0;
    while (idxOld < oldPeakList->size() or idxNew < spec.size())
    {
      if (idxOld == oldPeakList->size())
      {
        updatedPeakList[idxUpdated] = spec[idxNew];
        updatedPeakTols[idxUpdated++] = spec.peakTols[idxNew++];
        continue;
      }
      if (idxNew == spec.size())
      {
        updatedPeakList[idxUpdated] = (*oldPeakList)[idxOld];
        updatedPeakTols[idxUpdated++] = (*oldPeakTols)[idxOld++];
        continue;
      }

      if (fabs((*oldPeakList)[idxOld][0] - spec[idxNew][0])
          < min((*oldPeakTols)[idxOld], spec.peakTols[idxNew]) + 0.0001)
      {
        updatedPeakTols[idxUpdated] = min((*oldPeakTols)[idxOld],
                                          spec.peakTols[idxNew]);
        updatedPeakList[idxUpdated][0] = (*oldPeakList)[idxOld][0];
        updatedPeakList[idxUpdated++][1] = (*oldPeakList)[idxOld++][1]
            + spec[idxNew++][1];
      }
      else if ((*oldPeakList)[idxOld][0] < spec[idxNew][0])
      {
        updatedPeakTols[idxUpdated] = (*oldPeakTols)[idxOld];
        updatedPeakList[idxUpdated++] = (*oldPeakList)[idxOld++];
      }
      else
      {
        updatedPeakTols[idxUpdated] = spec.peakTols[idxNew];
        updatedPeakList[idxUpdated++] = spec[idxNew++];
      }
    }
    updatedPeakList.resize(idxUpdated); // Just in case some peaks get merged along the way
    updatedPeakTols.resize(idxUpdated);

    if (putHere == 0)
    {
      delete oldPeakList;
      delete oldPeakTols;
    }
  }
  // -------------------------------------------------------------------------
  void Spectrum::mergePeakList(vector<TwoValues<float> > &newPeaks,
                               Spectrum *putHere)
  {
    WARN_MSG("mergePeakList has not been updated for per-peak tolerances");

    int idxUpdated, idxOld, idxNew;
    vector<TwoValues<float> > &updatedPeakList =
        (putHere == 0) ? peakList : putHere->peakList;
    vector<TwoValues<float> > *oldPeakList;

    if (putHere == 0)
    {
      oldPeakList = new vector<TwoValues<float> >;
      oldPeakList->resize(peakList.size());
      for (int i = 0; i < peakList.size(); i++)
        (*oldPeakList)[i] = peakList[i];
    }
    else
    {
      oldPeakList = &peakList; // No copying necessary if results goes into another object (putHere)
      putHere->parentMass = parentMass;
      putHere->parentCharge = parentCharge;
    }

    // Sorted merge
    updatedPeakList.resize(oldPeakList->size() + newPeaks.size());
    idxUpdated = 0;
    idxOld = 0;
    idxNew = 0;
    while (idxOld < oldPeakList->size() or idxNew < newPeaks.size())
    {
      if (idxOld == oldPeakList->size())
      {
        updatedPeakList[idxUpdated++] = newPeaks[idxNew++];
        continue;
      }
      if (idxNew == newPeaks.size())
      {
        updatedPeakList[idxUpdated++] = (*oldPeakList)[idxOld++];
        continue;
      }

      if (fabs((*oldPeakList)[idxOld][0] - newPeaks[idxNew][0]) < 0.0001)
      {
        updatedPeakList[idxUpdated][0] = (*oldPeakList)[idxOld][0];
        updatedPeakList[idxUpdated++][1] = (*oldPeakList)[idxOld++][1]
            + newPeaks[idxNew++][1];
      }
      else if ((*oldPeakList)[idxOld][0] < newPeaks[idxNew][0])
        updatedPeakList[idxUpdated++] = (*oldPeakList)[idxOld++];
      else
        updatedPeakList[idxUpdated++] = newPeaks[idxNew++];
    }
    updatedPeakList.resize(idxUpdated); // Just in case some peaks get merged along the way

    if (putHere == 0)
    {
      delete oldPeakList;
    }
  }
  // -------------------------------------------------------------------------
  void Spectrum::mergeClosestPeaks(Spectrum& newSpec, short mergeType)
  {
    if (size() == 0)
    {
      peakList = newSpec.peakList;
      peakTols = newSpec.peakTols;
      return;
    }

    list<int> peaksAppend;
    map<int, list<MZRange> > peaksMerge;
    float mass1, intensity1, tolerance1, mass2, intensity2, tolerance2;
    int closestIdx;
    MZRange nextRange1, nextRange2;
    list<MZRange> mergePeaks;

    for (int i = 0; i < newSpec.size(); i++)
    {
      mass1 = newSpec[i][0];
      intensity1 = newSpec[i][1];
      tolerance1 = newSpec.peakTols[i];
      nextRange1.set(mass1, intensity1, tolerance1);

      closestIdx = findClosest(mass1);
      mass2 = peakList[closestIdx][0];
      intensity2 = peakList[closestIdx][1];
      tolerance2 = peakTols[closestIdx];
      nextRange2.set(mass2, intensity2, tolerance2);

      if (MZRange::EqualWithinRange(mass1, mass2, tolerance1 + tolerance2))
      {
        if (peaksMerge.count(closestIdx) == 0)
        {
          mergePeaks.clear();
          mergePeaks.push_back(nextRange2);
          mergePeaks.push_back(nextRange1);
          peaksMerge[closestIdx] = mergePeaks;
        }
        else
        {
          peaksMerge[closestIdx].push_back(nextRange1);
        }
      }
      else
      {
        peaksAppend.push_back(i);
      }
    }

    for (map<int, list<MZRange> >::iterator mergeIt = peaksMerge.begin();
        mergeIt != peaksMerge.end(); mergeIt++)
    {
      mass1 = peakList[mergeIt->first][0];
      intensity1 = peakList[mergeIt->first][1];
      tolerance1 = peakTols[mergeIt->first];
      nextRange1.set(mass1, intensity1, tolerance1);

      switch (mergeType)
      {
      case 0:
        nextRange1.MergeMZRanges(&mergeIt->second, 1);
      case 1:
        nextRange1.MergeMZRanges(&mergeIt->second, 2);
      case 2:
        nextRange1.MergeMZRanges(&mergeIt->second, 3);
      case 3:
        nextRange1.MergeMZRanges(&mergeIt->second, 0);
      }

      peakList[mergeIt->first][0] = nextRange1.getMass();
      peakList[mergeIt->first][1] = nextRange1.getIntensity();
      peakTols[mergeIt->first] = nextRange1.getTolerance();
    }

    int idx = peakList.size();
    peakList.resize(peakList.size() + peaksAppend.size());
    peakTols.resize(peakList.size());
    for (list<int>::iterator peakIt = peaksAppend.begin();
        peakIt != peaksAppend.end(); peakIt++)
    {
      peakList[idx] = newSpec[*peakIt];
      peakTols[idx] = newSpec.peakTols[*peakIt];
      idx++;
    }

    sortPeaks();
  }

  // -------------------------------------------------------------------------
  void Spectrum::addPeaksPartial(Spectrum & addSpectrum,
                                 Spectrum & outSpec,
                                 float multiplier)
  {
    if (DEBUG_SPRINKLE)
      DEBUG_VAR(scan);
    int scanIndex = scan - 1;
    outSpec = *this;

    float avgIntensityStar = 0.0;
    for (int iPeak = 0; iPeak < outSpec.size(); iPeak++)
    {
      avgIntensityStar += outSpec[iPeak][1];
      if (DEBUG_SPRINKLE)
        DEBUG_MSG(outSpec[iPeak][0] << "  " << outSpec[iPeak][1]);
    }
    avgIntensityStar /= (float)outSpec.size();
    if (DEBUG_SPRINKLE)
      DEBUG_VAR(avgIntensityStar);

    float avgIntensityPrm = 0.0;
    if (DEBUG_SPRINKLE)
      DEBUG_VAR(addSpectrum.parentMass);
    addSpectrum.filterLowMassPeaks(57.0);
    addSpectrum.filterHighMassPeaks(addSpectrum.parentMass - 57.0);
    for (int iPeak = 0; iPeak < addSpectrum.size(); iPeak++)
    {
      if (DEBUG_SPRINKLE)
        DEBUG_MSG(addSpectrum[iPeak][0] << "  " << addSpectrum[iPeak][1]);
      avgIntensityPrm += addSpectrum[iPeak][1];
    }
    avgIntensityPrm /= (float)addSpectrum.size();
    if (DEBUG_SPRINKLE)
      DEBUG_VAR(avgIntensityPrm);

    // Remove any peaks that already exist in this spectrum
    list<int> peaksToRemove;
    for (int iPeak = 0; iPeak < outSpec.size(); iPeak++)
    {
      int iIndex = addSpectrum.findPeaks(outSpec[iPeak][0],
                                         outSpec.getTolerance(iPeak));
      if (iIndex != -1)
      {
        peaksToRemove.push_back(iIndex);
      }
    }
    if (DEBUG_SPRINKLE)
      DEBUG_VAR(peaksToRemove.size());
    list<int>::iterator itrList = peaksToRemove.begin();
    list<int>::iterator itrListEnd = peaksToRemove.end();
    for (; itrList != itrListEnd; itrList++)
    {
      if (DEBUG_SPRINKLE)
        DEBUG_VAR(*itrList);
    }

    if (peaksToRemove.size() != 0)
    {
      if (!addSpectrum.removePeaks(peaksToRemove)) {
        WARN_MSG("Remove peaks failed");
      }
    }

    // Smash down the intensities of the add spectra by desired factor
    for (int iPeak = 0; iPeak < addSpectrum.size(); iPeak++)
    {
      addSpectrum[iPeak][1] *= avgIntensityStar / avgIntensityPrm * multiplier;
      if (DEBUG_SPRINKLE)
        DEBUG_MSG(addSpectrum[iPeak][0] << "  " << addSpectrum[iPeak][1]);
    }

    if (DEBUG_SPRINKLE)
      DEBUG_VAR(outSpec.size());
    outSpec.mergeClosestPeaks(addSpectrum, 1);

    if (DEBUG_SPRINKLE)
    {
      DEBUG_VAR(outSpec.size());
      for (int iPeak = 0; iPeak < outSpec.size(); iPeak++)
      {
        DEBUG_MSG(outSpec[iPeak][0] << "  " << outSpec[iPeak][1]);
      }
    }
    return;
  }

  // -------------------------------------------------------------------------
  void Spectrum::massesToIndices(Spectrum &masses,
                                 vector<int> &indices,
                                 float peakTol)
  {
    vector<float> tempTols;
    bool modTols = false;

    if (peakTol >= 0)
    {
      modTols = true;
      rememberTolerances();
      setPeakTolerance(peakTol);
    }

    int idxClosest, idxPeaks = 0, idxMasses;
    float distClosest; // = peakTol + 1; // Supremum for distance to closest peak
    indices.resize(masses.size());
    for (idxMasses = 0; idxMasses < (int)masses.size(); idxMasses++)
    {
      while (idxPeaks > 0
          and peakList[idxPeaks][0]
              > masses[idxMasses][0] - getTolerance(idxPeaks))
        idxPeaks--;
      while (idxPeaks < (int)peakList.size()
          and peakList[idxPeaks][0]
              < masses[idxMasses][0] - getTolerance(idxPeaks))
        idxPeaks++;
      idxClosest = -1;
      distClosest = getTolerance(idxPeaks) + 1.0;
      while (idxPeaks < (int)peakList.size()
          and abs(peakList[idxPeaks][0] - masses[idxMasses][0])
              <= getTolerance(idxPeaks) + 0.0001)
      {
        if (abs(peakList[idxPeaks][0] - masses[idxMasses][0]) < distClosest)
        {
          idxClosest = idxPeaks;
          distClosest = abs(peakList[idxPeaks][0] - masses[idxMasses][0]);
        }
        idxPeaks++;
      }
      if (idxClosest >= 0)
        indices[idxMasses] = idxClosest;
      else
        indices[idxMasses] = -1;
      idxPeaks = min(idxPeaks, (int)peakList.size() - 1);
    }

    if (modTols)
    {
      revertTolerances();
    }
  }
  // -------------------------------------------------------------------------
  void Spectrum::selectIndices(vector<int> &idx)
  {
    WARN_MSG("selectIndices has not been updated for per-peak tolerances");

    unsigned int i;
    for (i = 0; i < idx.size(); i++)
      if (idx[i] < peakList.size())
        peakList[i] = peakList[idx[i]];
      else
        break;
    resize(i);
  }
  // -------------------------------------------------------------------------
  unsigned int Spectrum::filterLowIntensity(float minIntensity)
  {
    list<int> peaksToRemove;

    for (int i = 0; i < size(); i++)
    {
      if (peakList[i][1] < minIntensity)
      {
        peaksToRemove.push_back(i);
      }
    }

    removePeaks(peaksToRemove);

    return size();
  }

  // -------------------------------------------------------------------------
  unsigned int Spectrum::filterLowMassPeaks(float minMass)
  {
    unsigned int minIdxToKeep;
    for (minIdxToKeep = 0; minIdxToKeep < peakList.size(); minIdxToKeep++)
      if (peakList[minIdxToKeep][0] >= minMass)
        break;
    if (minIdxToKeep == 0)
      return 0;

    for (unsigned int peakIdx = 0; peakIdx < peakList.size() - minIdxToKeep;
        peakIdx++)
    {
      peakList[peakIdx] = peakList[minIdxToKeep + peakIdx];
      peakTols[peakIdx] = peakTols[minIdxToKeep + peakIdx];
    }
    peakList.resize(peakList.size() - minIdxToKeep);
    peakTols.resize(peakList.size());
    return (minIdxToKeep);
  }
  // -------------------------------------------------------------------------
  unsigned int Spectrum::filterHighMassPeaks(float maxMass)
  {
    int maxIdxToKeep;
    for (maxIdxToKeep = peakList.size() - 1; maxIdxToKeep > 0; maxIdxToKeep--)
      if (peakList[maxIdxToKeep][0] <= maxMass)
        break;
    if (maxIdxToKeep == peakList.size() - 1)
      return 0;

    unsigned int prevSize = peakList.size();
    peakList.resize(maxIdxToKeep + 1);
    peakTols.resize(maxIdxToKeep + 1);
    return (prevSize - peakList.size());
  }
  // -------------------------------------------------------------------------
  /**
   * Removes peaks that have an intensity ranked worse than maxRank
   *   compared to all neighboring peaks +/- windowRadius
   * @param maxRank maximum allowable rank of each peak
   * @param windowRadius radius of peak comparison in the spectrum
   * @return
   */
  void Spectrum::rankFilterPeaks(int maxRank, float windowRadius)
  {
    //If spectrum has too many peaks (>1000, filter to that first)
    //Filter those so this doesnt take forever
    if (this->size() > 1000)
    {
      DEBUG_MSG("Spectum Contains more the 1000 peaks, globalling filtering before window filtering");
      this->selectTopK(1000);
    }

    list<int> peaksToRemove;

    list<TwoValues<float> > rankedPeaks;
    int j, k, curRank, peakIdx;
    float mass, intensity, minRange, maxRange, prevInten;
    bool foundPeak;
    int numPeaks = size();
    for (int i = 0; i < numPeaks; i++)
    {
      rankedPeaks.clear();
      mass = peakList[i][0];
      intensity = peakList[i][1];
      minRange = mass - windowRadius;
      maxRange = mass + windowRadius;
      rankedPeaks.push_back(TwoValues<float>(intensity, (float)i));
      j = i - 1;
      while (j >= 0 && peakList[j][0] > minRange)
      {
        rankedPeaks.push_back(TwoValues<float>(peakList[j][1], (float)j));
        --j;
      }
      j = i + 1;
      while (j < numPeaks && peakList[j][0] < maxRange)
      {
        rankedPeaks.push_back(TwoValues<float>(peakList[j][1], (float)j));
        ++j;
      }

      rankedPeaks.sort();
      rankedPeaks.reverse();
      curRank = 0;
      foundPeak = false;
      prevInten = -1.0;
      //cout << "ranked: ";
      for (list<TwoValues<float> >::iterator rankIt = rankedPeaks.begin();
          rankIt != rankedPeaks.end(); rankIt++)
      {
        //cout << "(" << (*rankIt)[0] << "," << (*rankIt)[1] << "), ";
        peakIdx = floatToInt((*rankIt)[1]);
        if ((*rankIt)[0] != prevInten)
        {
          curRank++;
          prevInten = (*rankIt)[0];
        }
        if (curRank > maxRank)
        {
          break;
        }
        if (peakIdx == i)
        {
          foundPeak = true;
          break;
        }
      }
      //cout << "\n";
      if (!foundPeak)
      {
        peaksToRemove.push_back(i);
      }
    }

    removePeaks(peaksToRemove);
  }

  void Spectrum::consolidatePeaks(const bool addIntensities)
  {
    list<int> peaksToRemove;
    vector<bool> checked(size(), false);

    for (int i = 0; i < size(); i++)
    {
      if (checked[i])
      {
        continue;
      }
      MZRange basePeak = getPeak(i);
      list<int> peakGroup;
      peakGroup.push_back(i);
      int maxIntenPeak = i;
      float maxIntensity = basePeak.getIntensity();
      float totIntensity = maxIntensity;

      for (int j = i + 1; j < size(); j++)
      {
        if (checked[j])
        {
          continue;
        }

        MZRange nextPeak = getPeak(j);

        if (nextPeak != basePeak)
        {
          continue;
        }
        peakGroup.push_back(j);
        totIntensity += nextPeak.getIntensity();
        checked[j] = true;

        if (nextPeak.getIntensity() > maxIntensity)
        {
          maxIntensity = nextPeak.getIntensity();
          maxIntenPeak = j;
        }
      }

      if (peakGroup.size() < 2)
      {
        continue;
      }

      if (addIntensities)
      {
        peakList[maxIntenPeak][1] = totIntensity;
      }

      for (list<int>::iterator pIt = peakGroup.begin(); pIt != peakGroup.end();
          pIt++)
      {
        if (*pIt != maxIntenPeak)
        {
          peaksToRemove.push_back(*pIt);
        }
      }
    }

    removePeaks(peaksToRemove);
  }

  int Spectrum::removeChargeReducedPrecursors(const float reducedPrecursorTol,
                                              const bool removeNeutralLoss,
                                              const float minKMedianIntensity)
  {
    if (size() == 0)
    {
      return 0;
    }
    Spectrum tempSpec;
    tempSpec.resize(this->size());
    for (int i = 0; i < size(); i++)
    {
      tempSpec[i][0] = peakList[i][1];
      tempSpec[i][1] = (float)i;
    }
    tempSpec.sortPeaks();

    int medianIdx = ((int)tempSpec.size()) / 2;
    float minIntensityRemove = tempSpec[medianIdx][0] * minKMedianIntensity;

    //DEBUG_VAR(minIntensityRemove);
    //DEBUG_VAR(tempSpec[medianIdx][0]);

    set<int> peaksToRemoveSet;
    const float neutralPM = this->parentMass - AAJumps::massHion;
    for (short charge = 1; charge <= this->parentCharge; charge++)
    {
      const float fCharge = (float)charge;
      const float mz = (neutralPM + (AAJumps::massHion * fCharge)) / fCharge;

      list<int> curMatches;
      this->findPeaks(mz, reducedPrecursorTol, &curMatches);

      for (list<int>::const_iterator pIt = curMatches.begin();
          pIt != curMatches.end(); pIt++)
      {
        if (peakList[*pIt][1] >= minIntensityRemove)
        {
          peaksToRemoveSet.insert(*pIt);
        }
      }

      if (!removeNeutralLoss)
      {
        continue;
      }

      float minNeutralLoss = mz;
      if (charge <= 2)
      {
        minNeutralLoss -= 60.0 / fCharge;
      }
      else
      {
        minNeutralLoss -= 18.0 / fCharge;
      }

      for (int i = 0; i < size() && peakList[i][0] <= mz; i++)
      {
        if (peakList[i][0] < minNeutralLoss)
        {
          continue;
        }

        if (peakList[i][1] >= minIntensityRemove)
        {
          peaksToRemoveSet.insert(i);
        }
      }
    }

    list<int> peaksToRemoveList(peaksToRemoveSet.begin(),
                                peaksToRemoveSet.end());
    //DEBUG_VAR(peaksToRemoveList.size());
    this->removePeaks(peaksToRemoveList);

    return peaksToRemoveList.size();
  }

  // -------------------------------------------------------------------------
  void Spectrum::outputDebug()
  {
    DEBUG_MSG(parentMass << " " << parentCharge);
    for (unsigned int i = 0; i < peakList.size(); i++)
    {
      DEBUG_MSG(peakList[i][0] << " " << peakList[i][1] << " " << peakTols[i]);
    }
  }

  // -------------------------------------------------------------------------
  void Spectrum::output(ostream &output)
  {
    output << parentMass << " " << parentCharge << endl;
    for (unsigned int i = 0; i < peakList.size(); i++)
      output << peakList[i][0] << " " << peakList[i][1] << " " << peakTols[i]
          << endl;
  }

  // -------------------------------------------------------------------------
  void Spectrum::output_ms2(ostream &output)
  {
    //output << ":0.0.0\n";
    output << parentMass << " " << parentCharge << endl;
    for (unsigned int i = 0; i < peakList.size(); i++)
      output << peakList[i][0] << " " << peakList[i][1] << endl;
  }
  // -------------------------------------------------------------------------
  void Spectrum::reverse(float pmOffset, Spectrum *putHere)
  {
    // m_reversed = !m_reversed;   LARS: Pretty sure this is wrong (it is done below)
    vector<TwoValues<float> > &updatedPeakList =
        (putHere == 0) ? peakList : putHere->peakList;
    vector<float> &updatedPeakTols =
        (putHere == 0) ? peakTols : putHere->peakTols;

    vector<TwoValues<float> > *oldPeakList;
    vector<float> *oldPeakTols;

    if (putHere == 0)
    {
      oldPeakList = new vector<TwoValues<float> >;
      oldPeakTols = new vector<float>;
      oldPeakList->resize(peakList.size());
      oldPeakTols->resize(peakList.size());
      for (unsigned int i = 0; i < peakList.size(); i++)
      {
        (*oldPeakList)[i] = peakList[i];
        (*oldPeakTols)[i] = peakTols[i];
      }
    }
    else
    {
      oldPeakList = &peakList; // No copying necessary if results goes into another object (putHere)
      oldPeakTols = &peakTols;
      (*putHere) = *this;
    }

    float totMass = parentMass - AAJumps::massHion + pmOffset;
    unsigned int numPeaks = oldPeakList->size();
    for (unsigned int i = 0; i < numPeaks; i++)
    {
      updatedPeakList[numPeaks - i - 1].set(totMass - (*oldPeakList)[i][0],
                                            (*oldPeakList)[i][1]);
      updatedPeakTols[numPeaks - i - 1] = (*oldPeakTols)[i];
    }

    if (putHere == 0)
    {
      m_reversed = !m_reversed;
      delete oldPeakList;
      delete oldPeakTols;
    }
    else
    {
      putHere->m_reversed = !m_reversed;
    }
  }
  // -------------------------------------------------------------------------
  void Spectrum::rotate(float offset, float cterm)
  {
    WARN_MSG("rotate has not been updated for per-peak tolerances");

    vector<TwoValues<float> > tmpPeaks;
    int peakIdx, v;

    if (offset < 0)
    {
      for (peakIdx = 0; peakIdx < (int)peakList.size(); peakIdx++)
        if (peakList[peakIdx][0] > -offset)
          break;
      int firstAbove = peakIdx;

      tmpPeaks.resize(firstAbove);
      for (peakIdx = 0; peakIdx < firstAbove; peakIdx++) // Compute peaks with masses lower than offset
        tmpPeaks[peakIdx].set(parentMass - cterm + peakList[peakIdx][0]
                                  + offset,
                              peakList[peakIdx][1]);
      for (peakIdx = 0; peakIdx < (int)peakList.size() - firstAbove; peakIdx++) // Set peaks with masses higher than offset
      {
        v = peakIdx + firstAbove;
        peakList[peakIdx].set(peakList[v][0] + offset, peakList[v][1]);
      }
      v = peakList.size() - firstAbove;
      for (peakIdx = 0; peakIdx < firstAbove; peakIdx++) // Set peaks with masses lower than offset
        peakList[v + peakIdx] = tmpPeaks[peakIdx];
    }
    else
    {
      for (peakIdx = 0; peakIdx < (int)peakList.size(); peakIdx++)
        if (peakList[peakIdx][0] > parentMass - cterm - offset)
          break;
      int firstWrap = peakIdx, numWrap = peakList.size() - firstWrap;

      tmpPeaks.resize(numWrap);
      for (peakIdx = firstWrap; peakIdx < (int)peakList.size(); peakIdx++) // Compute peaks with masses higher than offset
        tmpPeaks[peakIdx - firstWrap].set(peakList[peakIdx][0] + offset
                                              - (parentMass - cterm),
                                          peakList[peakIdx][1]);
      for (peakIdx = (int)peakList.size() - 1; peakIdx >= numWrap; peakIdx--) // Set peaks with masses higher than offset
      {
        v = peakIdx - numWrap;
        peakList[peakIdx].set(peakList[v][0] + offset, peakList[v][1]);
      }
      for (peakIdx = 0; peakIdx < numWrap; peakIdx++) // Set peaks with masses lower than offset
        peakList[peakIdx] = tmpPeaks[peakIdx];
    }
  }
  // -------------------------------------------------------------------------
  void Spectrum::selectTopK(unsigned int topK, Spectrum *putHere)
  {
    WARN_MSG("selectTopK has not been updated for per-peak tolerances and should be replaced by rankFilterPeaks");

    vector<TwoValues<float> > &updatedPeakList =
        (putHere == 0) ? peakList : putHere->peakList;
    vector<TwoValues<float> > *oldPeakList;
    unsigned int peakIdx;

    if (putHere == 0)
    {
      oldPeakList = new vector<TwoValues<float> >;
      oldPeakList->resize(peakList.size());
      for (unsigned int i = 0; i < peakList.size(); i++)
        (*oldPeakList)[i] = peakList[i];
    }
    else
    {
      oldPeakList = &peakList; // No copying necessary if results goes into another object (putHere)
      (*putHere) = *this;
    }
    vector<TwoValues<float> > sortedInts(oldPeakList->size());

    if (topK < oldPeakList->size())
    {
      for (peakIdx = 0; peakIdx < oldPeakList->size(); peakIdx++)
        sortedInts[peakIdx].set((*oldPeakList)[peakIdx][1], peakIdx);
      sort(sortedInts.begin(), sortedInts.end());

      updatedPeakList.resize(topK);
      for (peakIdx = 1; peakIdx <= topK; peakIdx++)
        updatedPeakList[peakIdx - 1] =
            (*oldPeakList)[(unsigned int)sortedInts[sortedInts.size() - peakIdx][1]];

      sort(updatedPeakList.begin(), updatedPeakList.end());
    }

    if (putHere == 0)
      delete oldPeakList;
  }
  // -------------------------------------------------------------------------
  void Spectrum::selectTopK(unsigned int topK,
                            TwoValues<float> w,
                            Spectrum *putHere)
  {
    WARN_MSG("selectTopK has not been updated for per-peak tolerances and should be replaced by rankFilterPeaks");

    vector<TwoValues<float> > &updatedPeakList =
        (putHere == 0) ? peakList : putHere->peakList;
    vector<TwoValues<float> > *oldPeakList;
    unsigned int peakIdx, idxLow = 0, idxHigh = 0, putIdx = 0;

    if (putHere == 0)
    {
      oldPeakList = new vector<TwoValues<float> >;
      oldPeakList->resize(peakList.size());
      for (unsigned int i = 0; i < peakList.size(); i++)
        (*oldPeakList)[i] = peakList[i];
    }
    else
    {
      oldPeakList = &peakList; // No copying necessary if results goes into another object (putHere)
      (*putHere) = *this;
    }
    vector<TwoValues<float> > sortedInts;

    if (topK < oldPeakList->size())
    {
      for (peakIdx = 0; peakIdx < oldPeakList->size(); peakIdx++)
      {
        while (idxLow < peakIdx
            and (*oldPeakList)[idxLow][0] < (*oldPeakList)[peakIdx][0] + w[0])
          idxLow++;
        while (idxHigh < oldPeakList->size()
            and (*oldPeakList)[idxHigh][0] < (*oldPeakList)[peakIdx][0] + w[1])
          idxHigh++;

        if (topK < idxHigh - idxLow)
        {
          sortedInts.resize(idxHigh - idxLow);
          for (unsigned int neighIdx = idxLow; neighIdx < idxHigh; neighIdx++)
            sortedInts[neighIdx - idxLow].set((*oldPeakList)[neighIdx][1],
                                              neighIdx);
          sort(sortedInts.begin(), sortedInts.end());
          bool keep = false;
          for (unsigned int rank = 1; rank <= topK; rank++)
            if (sortedInts[sortedInts.size() - rank][1] == peakIdx)
            {
              keep = true;
              break;
            }
          if (not keep)
            continue;
        }

        updatedPeakList[putIdx++] = (*oldPeakList)[peakIdx];
      }
      updatedPeakList.resize(putIdx);
    }

    if (putHere == 0)
      delete oldPeakList;
  }

  // -------------------------------------------------------------------------
  void Spectrum::getSymmetricPeakPairs(float pmOffset,
                                       float tolerance,
                                       vector<vector<float> > &pairs,
                                       vector<vector<int> > &pairsIdx)
  {

    vector<float> tempTols;
    bool modTols = false;

    if (tolerance >= 0)
    {
      modTols = true;
      rememberTolerances();
      setPeakTolerance(tolerance);
    }

    float aaMass = parentMass + (pmOffset - idDist) * AAJumps::massHion;
    int curPair, // Index of the current pair being generated
        i, j, k; // iterators over the prefix/suffix peaks

    vector<float> peakScores(peakList.size()); // Peak scores including scores of neighboring peaks within tolerance
    vector<bool> processed(peakList.size()); // Keeps track of which peaks have already been used in some pair

    for (k = 0; k < (int)peakList.size(); k++)
    {
      processed[k] = false;
      peakScores[k] = peakList[k][1];
      // Removed 06/08/30 - It's reasonable but this is not the right place to do this
      //		for(i=k-1; i>=0 && peakList[i][0]>=peakList[k][0]-tolerance; i--) peakScores[k]+=peakList[i][1];
      //		for(i=k+1; i<peakList.size() && peakList[i][0]<=peakList[k][0]+tolerance; i++) peakScores[k]+=peakList[i][1];
    }

    //	for(k=0;k<peakList.size();k++) cerr<<k<<": "<<peakList[k][0]<<"\t"<<peakList[k][1]<<", "<<peakScores[k]<<"\n";

    unsigned int maxNumPairs = (unsigned int)round(((double)peakList.size()
        + 1.0) * (double)peakList.size() / 2);
    pairs.resize(maxNumPairs);
    pairsIdx.resize(maxNumPairs);
    //	pairs.resize(peakList.size()*(1+2*(int)ceil(tolerance/(0.1*idDist))));
    //	pairsIdx.resize(peakList.size()*(1+2*(int)ceil(tolerance/(0.1*idDist))));
    i = 0; // Current peak
    j = peakList.size() - 1; // Candidate symmetric peak for current peak
    curPair = 0;
    while (i <= j)
    {
      if (peakList[i][0] <= aaMass - peakList[j][0])
      {
        for (k = j;
            k > i
                && peakList[k][0]
                    >= aaMass - peakList[i][0] - getTolerance(k)
                        + getTolerance(i); k--)
        {
          pairs[curPair].resize(3);
          pairsIdx[curPair].resize(2);
          pairs[curPair][2] = peakScores[i] + peakScores[k];
          pairs[curPair][0] = peakList[i][0];
          pairs[curPair][1] = peakList[k][0];
          pairsIdx[curPair][0] = i;
          pairsIdx[curPair][1] = k;
          curPair++;
          processed[k] = true;
        }
        if (k == j && !processed[i])
        {
          pairs[curPair].resize(3);
          pairsIdx[curPair].resize(2);
          pairs[curPair][2] = peakScores[i];
          pairs[curPair][0] = peakList[i][0];
          pairs[curPair][1] = aaMass - peakList[i][0];
          pairsIdx[curPair][0] = i;
          pairsIdx[curPair][1] = -1;
          curPair++;
        }
        processed[i] = true;
        i++;
      }
      else
      {
        for (k = i;
            k < j
                && peakList[k][0]
                    <= aaMass - peakList[j][0] + getTolerance(k)
                        + getTolerance(j); k++)
        {
          pairs[curPair].resize(3);
          pairsIdx[curPair].resize(2);
          pairs[curPair][2] = peakScores[j] + peakScores[k];
          pairs[curPair][0] = peakList[k][0];
          pairs[curPair][1] = peakList[j][0];
          pairsIdx[curPair][0] = k;
          pairsIdx[curPair][1] = j;
          curPair++;
          processed[k] = true;
        }
        if (k == i && !processed[j])
        {
          pairs[curPair].resize(3);
          pairsIdx[curPair].resize(2);
          pairs[curPair][2] = peakScores[j];
          pairs[curPair][0] = aaMass - peakList[j][0];
          pairs[curPair][1] = peakList[j][0];
          pairsIdx[curPair][0] = -1;
          pairsIdx[curPair][1] = j;
          curPair++;
        }
        processed[j] = true;
        j--;
      }
    }
    pairs.resize(curPair);
    pairsIdx.resize(curPair);

    if (modTols)
    {
      revertTolerances();
    }
  }
  // -------------------------------------------------------------------------
  void Spectrum::makeSymmetric(float pmOffset,
                               float tolerance,
                               vector<int>* indicies)
  {
    bool modTols = false;

    if (tolerance >= 0)
    {
      modTols = true;
      rememberTolerances();
      setPeakTolerance(tolerance);
    }

    // Find all symmetric pairs in the spectrum
    vector<vector<float> > pairs;
    vector<vector<int> > pairsIdx;
    getSymmetricPeakPairs(pmOffset, -1.0, pairs, pairsIdx);

    map<float, int> indexMap;

    // Check which peaks already have a symmetric peak in the spectrum
    vector<bool> needsPair(size());
    unsigned int peakIdx;
    for (peakIdx = 0; peakIdx < needsPair.size(); peakIdx++)
    {
      needsPair[peakIdx] = true;
      if (indicies)
      {
        indexMap[peakList[peakIdx][0]] = peakIdx;
      }
    }
    for (unsigned int pairIdx = 0; pairIdx < pairsIdx.size(); pairIdx++)
    {
      if (pairsIdx[pairIdx][0] >= 0 and pairsIdx[pairIdx][1] >= 0)
      {
        needsPair[pairsIdx[pairIdx][0]] = false;
        needsPair[pairsIdx[pairIdx][1]] = false;
      }
    }

    // Count how many new peaks are necessary
    unsigned int curCount = peakList.size(), newCount = 0;
    for (peakIdx = 0; peakIdx < needsPair.size(); peakIdx++)
    {
      newCount += (unsigned int)needsPair[peakIdx];
    }
    resize(curCount + newCount);
    if (modTols)
    {
      oldPeakTols.resize(curCount + newCount);
    }

    // Insert new symmetric peaks
    unsigned int newIdx = curCount;
    for (peakIdx = 0; peakIdx < needsPair.size(); peakIdx++)
      if (needsPair[peakIdx])
      {
        peakList[newIdx].set(parentMass + (pmOffset - 1) * AAJumps::massHion
                                 - peakList[peakIdx][0],
                             peakList[peakIdx][1]);
        peakTols[newIdx] = peakTols[peakIdx];
        if (indicies)
        {
          indexMap[peakList[newIdx][0]] = -1;
        }
        if (modTols)
        {
          oldPeakTols[newIdx] = oldPeakTols[peakIdx];
        }
        if (peakList[newIdx][0] > 0)
        {
          newIdx++;
        }
      }
    resize(newIdx);

    if (modTols)
    {
      oldPeakTols.resize(newIdx);
      peakTols = oldPeakTols;
    }

    sortPeaks();

    if (indicies)
    {
      indicies->resize(size());
      for (int i = 0; i < size(); i++)
      {
        (*indicies)[i] = indexMap[peakList[i][0]];
      }
    }
  }
  // -------------------------------------------------------------------------
  float Spectrum::getGFProbability(float score,
                                   unsigned int atMass,
                                   bool roundMasses,
                                   AAJumps *aaMasses,
                                   unsigned int maxNumAAMods,
                                   vector<unsigned int> *modMasses,
                                   vector<pair<unsigned int, bool> > *ntermModMasses,
                                   bool outputGFTables,
                                   bool useExistingRawScore)
  {
    if (DEBUG_SPEC_PROB)
      DEBUG_VAR(useExistingRawScore);
    if (parentMass < AAJumps::massMH)
    {
      WARN_MSG("Warning: bad spectrum with parent mass " << parentMass<< ", skipping GF calculation.");
      return 1.0;
    }
    unsigned int maxMass =
        (unsigned int)roundMasses ?
            round(0.9995 * (parentMass - AAJumps::massMH)) :
            round(parentMass - AAJumps::massMH);
    if (maxMass != atMass)
      maxMass = atMass;

    // Calculate spectral probabilities
    if (aaMasses)
    {

      vector<int> prmScores(maxMass + 1);
      for (unsigned int massIdx = 0; massIdx <= maxMass; massIdx++)
        prmScores[massIdx] = 0;

      // Calculate maximum score, assume minimum score is zero
      int maxScore = 0;
      for (unsigned int peakIdx = 0; peakIdx < peakList.size(); peakIdx++)
      {
        int peakMass =
            (int)roundMasses ? round(0.9995 * peakList[peakIdx][0]) :
                peakList[peakIdx][0];
        int peakScore = (int)round(peakList[peakIdx][1]);
        if (peakScore < 0)
        {
          ERROR_MSG("Spectrum::getGFProbability():Negative PRM scores not supported().");
          return 1.0;
        }
        if (peakMass >= 0 and peakMass <= maxMass)
        {
          //cerr<<" --- kept "<<peakIdx<<": "<<peakMass<<", "<<peakScore<<endl;
          prmScores[peakMass] += peakScore;
          maxScore += peakScore;
        }
      }
      //cerr<<"maxScore = "<<maxScore<<", input score = "<<score<<endl;

      if (DEBUG_SPEC_PROB)
        DEBUG_VAR(aaMasses->size());
      // Initialize structures for spectral probabilities
      float aaProb = 1.0 / ((float)aaMasses->size());
      if (DEBUG_SPEC_PROB)
        DEBUG_VAR(aaProb);
      vector<unsigned int> intAAMasses(aaMasses->size());

      for (unsigned int aaIdx = 0; aaIdx < intAAMasses.size(); aaIdx++)
        intAAMasses[aaIdx] = (unsigned int)round((*aaMasses)[aaIdx]);

      vector<vector<float> > *currentProbs = new vector<vector<float> >,
          *prevProbs = new vector<vector<float> >, *tmp; // Used when reassigning currentProbs to prevProbs
      if (DEBUG_SPEC_PROB)
        DEBUG_VAR(maxMass);
      currentProbs->resize(maxMass + 1);
      prevProbs->resize(maxMass + 1);
      if (DEBUG_SPEC_PROB)
        DEBUG_VAR(maxScore);
      for (unsigned int massIdx = 0; massIdx <= maxMass; massIdx++)
      {
        (*currentProbs)[massIdx].resize(maxScore + 1);
        (*prevProbs)[massIdx].resize(maxScore + 1);
        for (unsigned int scoreIdx = 0; scoreIdx <= maxScore; scoreIdx++)
        {
          (*currentProbs)[massIdx][scoreIdx] = 0;
          (*prevProbs)[massIdx][scoreIdx] = 0;
        }
      }
      if (DEBUG_SPEC_PROB)
        DEBUG_TRACE;

      unsigned int scoreAtZero = 0;
      vector<int> matchesIdx;
      findMatches(0, 0.5, matchesIdx); // 0.5 is because if GF using 1 Da bins, not related to peak mass tolerance
      for (unsigned int i = 0; i < matchesIdx.size(); i++)
        if (peakList[matchesIdx[i]][0] >= 0 and peakList[matchesIdx[i]][1] > 0)
          scoreAtZero += ((unsigned int)round(peakList[matchesIdx[i]][1]));
      (*currentProbs)[0][scoreAtZero] = 1.0;

      // Process N-terminal mods
      if (ntermModMasses)
      {
        unsigned int numGenericMods = 0;
        for (unsigned int deltaIdx = 0; deltaIdx < ntermModMasses->size();
            deltaIdx++)
          if (not (*ntermModMasses)[deltaIdx].second)
            numGenericMods++;

        float ntermProb = 1.0 / (1.0 + ((float)numGenericMods));
        for (unsigned int deltaIdx = 0; deltaIdx < ntermModMasses->size();
            deltaIdx++)
        {
          // Allow starting at the offset corresponding to n-terminal modification
          if ((*ntermModMasses)[deltaIdx].first <= maxMass)
          {
            if ((*ntermModMasses)[deltaIdx].second) // Mod mass includes one amino acid mass (e.g., PyroGlu)
              (*currentProbs)[(*ntermModMasses)[deltaIdx].first][scoreAtZero] =
                  aaProb;
            else
              (*currentProbs)[(*ntermModMasses)[deltaIdx].first][scoreAtZero] =
                  ntermProb;
          }
        }
        (*currentProbs)[0][scoreAtZero] = ntermProb;

      } // if(ntermModMasses) {

      // Initialize all probs to zero
      m_gfProbs.resize(maxScore + 1);
      for (unsigned int scoreIdx = 0; scoreIdx <= maxScore; scoreIdx++)
        m_gfProbs[scoreIdx] = 0.0;

      // Calculate spectral probabilities
      for (unsigned int numMods = 0; numMods <= maxNumAAMods; numMods++)
      {
        //cerr<<"numMods = "<<numMods<<"\n";
        if (numMods > 0)
        {
          tmp = prevProbs;
          prevProbs = currentProbs;
          currentProbs = tmp;
          for (unsigned int massIdx = 0; massIdx <= maxMass; massIdx++)
            for (unsigned int scoreIdx = 0; scoreIdx <= maxScore; scoreIdx++)
              (*currentProbs)[massIdx][scoreIdx] = 0;
        }

        for (unsigned int massIdx = 0; massIdx <= maxMass; massIdx++)
        {
          unsigned int curPRMScore = prmScores[massIdx];

          // Iterate over predecessors at amino acid mass distance
          for (unsigned int deltaIdx = 0; deltaIdx < intAAMasses.size();
              deltaIdx++)
          {
            if (intAAMasses[deltaIdx] <= massIdx)
            {
              unsigned int prevMass = massIdx - intAAMasses[deltaIdx];
              for (unsigned int scoreIdx = curPRMScore; scoreIdx <= maxScore;
                  scoreIdx++)
              {
                (*currentProbs)[massIdx][scoreIdx] += aaProb
                    * (*currentProbs)[prevMass][scoreIdx - curPRMScore];
                /*
                 if((massIdx==97 and scoreIdx==1) or (massIdx==424 and scoreIdx==4)) {
                 cerr<<"prevMass = "<<prevMass<<", curPRMScore = "<<curPRMScore
                 <<", (*currentProbs)["<<"massIdx"<<"]["<<scoreIdx<<"] ("<<(*currentProbs)[massIdx][scoreIdx]
                 <<") = "<<aaProb<<" * "
                 <<(*currentProbs)[prevMass][scoreIdx - curPRMScore]<<" (in (*currentProbs)["<<prevMass<<"]["<<scoreIdx - curPRMScore<<"])\n";
                 }
                 */
              }
              /*if(massIdx==1035 and prevMass==872) {
               cerr<<"--- current score range was "<<curPRMScore<<" to "<<maxScore
               <<", (*currentProbs)[872][706] = "<<(*currentProbs)[872][706]
               <<", (*currentProbs)[1035][803] = "<<(*currentProbs)[1035][803]<<endl;
               }
               */
            }
          }

          // Iterate over predecessors at modification mass distances (predProbs contains
          //   spectral probabilities with one less mod than current stage)
          if (numMods > 0)
          {
            for (unsigned int deltaIdx = 0; deltaIdx < modMasses->size();
                deltaIdx++)
            {
              if ((*modMasses)[deltaIdx] <= massIdx)
              {
                unsigned int prevMass = massIdx - (*modMasses)[deltaIdx];
                for (unsigned int scoreIdx = curPRMScore; scoreIdx <= maxScore;
                    scoreIdx++)
                {
                  (*currentProbs)[massIdx][scoreIdx] += aaProb
                      * (*prevProbs)[prevMass][scoreIdx - curPRMScore];
                  /*
                   if((massIdx==97 and scoreIdx==1) or (massIdx==424 and scoreIdx==4)) {
                   cerr<<"prevMass = "<<prevMass<<", curPRMScore = "<<curPRMScore
                   <<", (*currentProbs)["<<"massIdx"<<"]["<<scoreIdx<<"] ("<<(*currentProbs)[massIdx][scoreIdx]
                   <<") = "<<aaProb<<" * "
                   <<(*prevProbs)[prevMass][scoreIdx - curPRMScore]<<" (in (*prevProbs)["<<prevMass<<"]["<<scoreIdx - curPRMScore<<"])\n";
                   }
                   */
                }
              }
            }
          }
        }

        // Accumulate probabilities with those of peptides with less mods
        for (unsigned int scoreIdx = 0; scoreIdx <= maxScore; scoreIdx++)
        {
          //DEBUG_VAR(scoreIdx);
          //DEBUG_VAR((*currentProbs)[maxMass][scoreIdx]);
          m_gfProbs[scoreIdx] += (*currentProbs)[maxMass][scoreIdx];
          //DEBUG_VAR(m_gfProbs[scoreIdx]);
        }

        if (outputGFTables)
        {
          stringstream s;
          s << "currentProbs_scan_" << scan << "_numMods_" << numMods << ".bin";
          Save_binArray(s.str().c_str(), (*currentProbs));
        }

      } // for(unsigned int numMods = 0; numMods <= maxNumAAMods; numMods++) {

      // Sum up cumulative probabilities to get p-values
      for (unsigned int scoreIdx = 1; scoreIdx <= maxScore; scoreIdx++)
      {
        //DEBUG_VAR(scoreIdx);
        //DEBUG_VAR(m_gfProbs[maxScore - scoreIdx + 1]);
        m_gfProbs[maxScore - scoreIdx] += m_gfProbs[maxScore - scoreIdx + 1];
        //DEBUG_VAR(m_gfProbs[maxScore - scoreIdx]);
      }

      delete currentProbs;
      delete prevProbs;

    } // if(aaMasses) {

    if (DEBUG_SPEC_PROB) DEBUG_VAR(score);

    int iScore = (int)score;

    if (DEBUG_SPEC_PROB) DEBUG_VAR(iScore);
    if (DEBUG_SPEC_PROB) {
      DEBUG_VAR(m_gfProbs.size());
      for (int g = 0; g < m_gfProbs.size(); g++) {
        DEBUG_MSG(g << "  " << m_gfProbs[g])
      }
    }
    if (iScore >= 0 and iScore < m_gfProbs.size()) {
      if (DEBUG_SPEC_PROB) DEBUG_VAR(m_gfProbs[iScore]);
      return m_gfProbs[iScore];
    } else {
      WARN_MSG("Spectrum::getGFProbability(): score "<<score<<" out of range [0,"<<m_gfProbs.size()-1<<"], returning prob 1.0\n");
    }

    return 1.0;
  }

  // -------------------------------------------------------------------------
  bool Spectrum::compare(Spectrum &toSpec)
  {
    if (peakList.size() != toSpec.peakList.size()
        || parentMass != toSpec.parentMass)
      return false;
    for (int i = 0; i < peakList.size(); i++)
      if (peakList[i][0] != toSpec.peakList[i][0]
          || peakList[i][1] != toSpec.peakList[i][1]
          || peakTols[i] != toSpec.getTolerance(i))
        return false;
    return true;
  }
  /*
   // -------------------------------------------------------------------------
   void Spectrum::mergeCommon(Spectrum &withSpectrum,
   Spectrum *toSpec,
   float shift,
   float peakTol,
   short mergeType)
   {
   WARN_MSG("mergeCommon has not been updated for per-peak tolerances");
   cerr << "mergeCommon has not been updated for per-peak tolerances\n";


   vector<int> idx1, idx2;
   unsigned int pivot;
   FindMatchPeaksAll2(*this, withSpectrum, shift, peakTol, idx1, idx2);
   Spectrum tmpToSpec;
   tmpToSpec.resize(idx1.size());
   tmpToSpec.copyNP(*this); // Temporary spectrum, in case toSpec==0

   float s, r1, r2;
   for (pivot = 0; pivot < idx1.size(); pivot++)
   {
   s = peakList[idx1[pivot]][1] + withSpectrum[idx2[pivot]][1];
   r1 = peakList[idx1[pivot]][1] / s;
   r2 = withSpectrum[idx2[pivot]][1] / s;
   tmpToSpec[pivot][0] = r1 * peakList[idx1[pivot]][0] + r2
   * (withSpectrum[idx2[pivot]][0] + shift);

   switch (mergeType)
   {
   case 0:
   tmpToSpec[pivot][1] = peakList[idx1[pivot]][1]
   + withSpectrum[idx2[pivot]][1];
   break;
   case 1:
   tmpToSpec[pivot][1] = max(peakList[idx1[pivot]][1],
   withSpectrum[idx2[pivot]][1]);
   break;
   case 2:
   tmpToSpec[pivot][1] = min(peakList[idx1[pivot]][1],
   withSpectrum[idx2[pivot]][1]);
   break;
   case 3:
   tmpToSpec[pivot][1] = peakList[idx1[pivot]][1];
   break;
   case 4:
   tmpToSpec[pivot][1] = withSpectrum[idx2[pivot]][1];
   break;
   }
   }
   sort(tmpToSpec.peakList.begin(), tmpToSpec.peakList.end());
   if (toSpec == (Spectrum *)0)
   (*this) = tmpToSpec;
   else
   {
   toSpec->copyNP(*this);
   (*toSpec) = tmpToSpec;
   }
   }
   */
  // -------------------------------------------------------------------------
  bool Spectrum::saveDTA(const char* outfile)
  {
    FILE* output = fopen(outfile, "wb");
    fprintf(output, "%.3f", parentMass);
    fprintf(output, " ");
    fprintf(output, "%d", parentCharge);
    fprintf(output, "\n");
    for (int idx = 0; idx < peakList.size(); idx++)
    {
      fprintf(output, "%.3f", peakList[idx][0]);
      fprintf(output, " ");
      fprintf(output, "%.3f", peakList[idx][1]);
      fprintf(output, "\n");
    }
    fclose(output);
    return true;
  }

  bool Spectrum::saveToBinaryStream(FILE* fp)
  {
    if (fp == 0)
    {
      return false;
    }

    unsigned int count;

    count = fwrite(&scan, sizeof(unsigned int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save scan #");
      return false;
    }

    count = fwrite(&msLevel, sizeof(short), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save msLevel");
      return false;
    }

    short fragT = (short)msFragType;
    count = fwrite(&fragT, sizeof(short), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save frag type");
      return false;
    }

    short massA = (short)msMassAnalyzerType;
    count = fwrite(&massA, sizeof(short), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save mass analyzer type");
      return false;
    }

    count = fwrite(&m_reversed, sizeof(bool), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save reversed flag");
      return false;
    }

    count = fwrite(&fileIndex, sizeof(int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save file index");
      return false;
    }

    vector<string> fName(1);
    fName[0] = fileName;
    if (!writeStringsToBinaryStream(fp, fName))
    {
      ERROR_MSG("Failed to save filename");
      return false;
    }

    count = fwrite(&parentMass, sizeof(float), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save parent mass");
      return false;
    }

    count = fwrite(&parentMassTol, sizeof(float), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save parent mass tolerance");
      return false;
    }

    count = fwrite(&parentCharge, sizeof(short), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save parent charge");
      return false;
    }

    count = fwrite(&parentMZ, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save parent m/z");
      return false;
    }

    count = fwrite(&collision_energy, sizeof(float), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save collision energy");
      return false;
    }

    unsigned int numPeaks = size();

    count = fwrite(&numPeaks, sizeof(unsigned int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save # of peaks");
      return false;
    }

    if (numPeaks == 0)
    {
      return true;
    }

    float* peakBuf = (float*)malloc(sizeof(float) * numPeaks);
    for (unsigned int i = 0; i < numPeaks; i++)
    {
      peakBuf[i] = peakList[i][0];
    }
    count = fwrite(peakBuf, sizeof(float), numPeaks, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save peak masses");
      free(peakBuf);
      return false;
    }

    for (unsigned int i = 0; i < numPeaks; i++)
    {
      peakBuf[i] = getTolerance(i);
    }
    count = fwrite(peakBuf, sizeof(float), numPeaks, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save peak tolerances");
      free(peakBuf);
      return false;
    }

    for (unsigned int i = 0; i < numPeaks; i++)
    {
      peakBuf[i] = peakList[i][1];
    }
    count = fwrite(peakBuf, sizeof(float), numPeaks, fp);
    free(peakBuf);
    if (count == 0)
    {
      ERROR_MSG("Failed to save peak intensities");
      return false;
    }

    //Writing Retention Time
    count = fwrite(&(retention_time), sizeof(float), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save retention time");
      return false;
    }

    //Writing Precursor intensity
    count = fwrite(&(precursor_intensity), sizeof(float), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save precursor intensity");
      return false;
    }

    return true;
  }

  bool Spectrum::loadFromBinaryStream(FILE* fp,
                                      map<string, unsigned short>& versions)
  {
    if (fp == 0)
    {
      return false;
    }
    unsigned short version = versions[BIN_VERSION_ID];
    unsigned short subVersion = versions[BIN_SUBVERSION_ID];

    if (version > BIN_VERSION
        || (version == BIN_VERSION && subVersion > BIN_SUBVERSION))
    {
      ERROR_MSG("Unsupported spectrum version " << version << "." << subVersion << " (this release supports up to " << BIN_VERSION << "." << BIN_SUBVERSION << "), you must obtain a later release of the code base to load this file.");
      return false;
    }

    unsigned int count;

    count = fread(&scan, sizeof(unsigned int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read scan #");
      return false;
    }

    count = fread(&msLevel, sizeof(short), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read msLevel");
      return false;
    }

    short fragT = 0;
    count = fread(&fragT, sizeof(short), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read fragmentation type");
      return false;
    }

    if (fragT == (short)FragType_CID)
    {
      msFragType = FragType_CID;
    }
    else if (fragT == (short)FragType_HCD)
    {
      msFragType = FragType_HCD;
    }
    else if (fragT == (short)FragType_ETD)
    {
      msFragType = FragType_ETD;
    }
    else if (fragT == (short)FragType_PRM)
    {
      msFragType = FragType_PRM;
    }
    else if (fragT == (short)FragType_PRM_ETD)
    {
      msFragType = FragType_PRM_ETD;
    }
    else
    {
      ERROR_MSG("Found unsupported fragmentation ID " << fragT);
      return false;
    }

    if (subVersion >= 4)
    {
      short massA = 0;
      count = fread(&massA, sizeof(short), 1, fp);
      if (count == 0)
      {
        ERROR_MSG("Failed to read mass analyzer type");
        return false;
      }

      if (massA == (short)MassAnalyzer_ION_TRAP)
      {
        msMassAnalyzerType = MassAnalyzer_ION_TRAP;
      }
      else if (massA == (short)MassAnalyzer_ORBI_TRAP)
      {
        msMassAnalyzerType = MassAnalyzer_ORBI_TRAP;
      }
      else if (massA == (short)MassAnalyzer_QTOF)
      {
        msMassAnalyzerType = MassAnalyzer_QTOF;
      }
      else
      {
        ERROR_MSG("Found unsupported mass analyzer ID " << massA);
        return false;
      }
    }

    count = fread(&m_reversed, sizeof(bool), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read reversed flag");
      return false;
    }

    if (version >= 1 && subVersion >= 3)
    {
      count = fread(&fileIndex, sizeof(int), 1, fp);
      if (count == 0)
      {
        ERROR_MSG("Failed to read file index");
        return false;
      }
    }
    else
    {
      fileIndex = -1;
    }

    vector<string> fName;
    if ((!readStringsFromBinaryStream(fp, fName)) || fName.size() == 0)
    {
      ERROR_MSG("Failed to read filename");
      return false;
    }
    fileName = fName[0];

    count = fread(&parentMass, sizeof(float), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read parent mass");
      return false;
    }

    count = fread(&parentMassTol, sizeof(float), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read parent mass tolerance");
      return false;
    }

    count = fread(&parentCharge, sizeof(short), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read parent charge");
      return false;
    }

    count = fread(&parentMZ, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read parent m/z");
      return false;
    }

    if (subVersion >= 4)
    {
      count = fread(&collision_energy, sizeof(float), 1, fp);
      if (count == 0)
      {
        ERROR_MSG("Failed to read collision energy");
        return false;
      }
    }

    unsigned int numPeaks = 0;
    count = fread(&numPeaks, sizeof(unsigned int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read # of peaks");
      return false;
    }
    resize(numPeaks);

    if (numPeaks == 0)
    {
      return true;
    }

    float* peakBuf = (float*)malloc(sizeof(float) * numPeaks);
    count = fread(peakBuf, sizeof(float), numPeaks, fp);
    if (numPeaks > 0 && count == 0)
    {
      ERROR_MSG("Failed to read peak masses");
      free(peakBuf);
      return false;
    }
    for (unsigned int i = 0; i < numPeaks; i++)
    {
      peakList[i][0] = peakBuf[i];
    }

    count = fread(peakBuf, sizeof(float), numPeaks, fp);
    if (numPeaks > 0 && count == 0)
    {
      ERROR_MSG("Failed to read peak tolerances");
      free(peakBuf);
      return false;
    }
    for (unsigned int i = 0; i < numPeaks; i++)
    {
      setTolerance(i, peakBuf[i]);
    }

    count = fread(peakBuf, sizeof(float), numPeaks, fp);
    if (numPeaks > 0 && count == 0)
    {
      ERROR_MSG("Failed to read peak intensities");
      free(peakBuf);
      return false;
    }
    for (unsigned int i = 0; i < numPeaks; i++)
    {
      peakList[i][1] = peakBuf[i];
    }
    free(peakBuf);

    //Version 1.2 or higher
    if (version >= 1 && subVersion >= 2)
    {
      count = fread(&retention_time, sizeof(float), 1, fp);
      if (count == 0)
      {
        ERROR_MSG("Failed to read retention time");
        return false;
      }

      count = fread(&precursor_intensity, sizeof(float), 1, fp);
      if (count == 0)
      {
        ERROR_MSG("Failed to read precursor intensity");
        return false;
      }

    }

    return true;
  }
  // -------------------------------------------------------------------------
  float Spectrum::getTotalIonCurrent(void) const
  {
    float total = 0;
    for (int i = 0; i < peakList.size(); i++)
    {
      total += peakList[i][1];
    }
    return total;
  }

  // -------------------------------------------------------------------------

  void Spectrum::computeSpectralProbabilities(vector<pair<unsigned int, bool> > & ntermMods,
                                              vector<unsigned int> & mods,
                                              float peakTol,
                                              bool useExistingRawScore)
  {
    if (psmList.empty())
    {
      return;
    }

    float matchScore;
    bool isReversed;
    Spectrum tmpSpec = *this;

    if (!useExistingRawScore)
    {
      // Need to scale the scores so that spectral prob doesn't run out of memory
      float maxScore = 0.0;
      for (unsigned int i = 0; i < tmpSpec.size(); i++)
      {
        maxScore = maxScore < tmpSpec[i][1] ? tmpSpec[i][1] : maxScore;
      }
      for (unsigned int i = 0; i < tmpSpec.size(); i++)
      {
        // GF uses rounded masses/scores
        tmpSpec[i][0] = round(0.9995 * tmpSpec[i][0]);
        float scaledPeak = ceil(tmpSpec[i][1] / maxScore * MAX_SPEC_PROB_PEAK); // Scale the peaks at the same time
        tmpSpec[i][1] = scaledPeak;
      }
      tmpSpec.parentMass = round(0.9995 * tmpSpec.parentMass);
      tmpSpec.setResolution(1, true); // Merge peaks with the same rounded mass
    }
    tmpSpec.filterLowMassPeaks(50);
    tmpSpec.filterHighMassPeaks(tmpSpec.parentMass - AAJumps::massMH - 50);

    if (DEBUG_SPEC_PROB) outputDebug();

    AAJumps aminoacids(1); // Used to calculate spectral probabilities
    for (unsigned int i = 0; i < aminoacids.size(); i++)
    {
      aminoacids.masses[i] = round(aminoacids.masses[i]);
    }
    aminoacids.forceTolerance(0, 1, true); // Enforce 1 Da resolution

    list<psmPtr>::iterator itr = psmList.begin();
    list<psmPtr>::iterator itr_end = psmList.end();
    for (; itr != itr_end; itr++)
    {
      psmPtr psm = *itr;

      if (DEBUG_SPEC_PROB)
        DEBUG_VAR(useExistingRawScore);
      if (useExistingRawScore)
      {
        if (DEBUG_SPEC_PROB)
          DEBUG_VAR(psm->m_score);
        if (psm->m_score < 0)
        {
          psm->m_score = 10.0;
          if (DEBUG_SPEC_PROB)
            DEBUG_MSG(scan << "\t" << psm->m_annotation << "\t" << psm->m_score);
          continue;
        }
        matchScore = psm->m_score;
        if (DEBUG_SPEC_PROB)
          DEBUG_VAR(matchScore);
      }
      else
      {
        matchScore = MatchSpecToPeptide(tmpSpec,
                                        psm->m_annotation.c_str(),
                                        peakTol,
                                        0,
                                        false,
                                        &isReversed,
                                        &aminoacids);
        if (DEBUG_SPEC_PROB)
          DEBUG_VAR(matchScore);
        if (isReversed)
        {
          if (DEBUG_SPEC_PROB)
            DEBUG_MSG("Reversed");
          tmpSpec.reverse(0.0);
          tmpSpec.filterLowMassPeaks(50);
          tmpSpec.filterHighMassPeaks(tmpSpec.parentMass - AAJumps::massMH
              - 50);
          tmpSpec.setResolution(1, true); // Merge peaks with the same rounded mass

          // matchScore needs to be recalculated because some PRMs may fall
          //  out of range after reversing the spectrum (e.g., if at mass 57 before reversing)
          matchScore = MatchSpecToPeptide(tmpSpec,
                                          psm->m_annotation.c_str(),
                                          peakTol,
                                          0,
                                          false,
                                          &isReversed,
                                          &aminoacids);
          if (DEBUG_SPEC_PROB)
            DEBUG_VAR(matchScore);
        } // if(isReversed) {

      }

      vector<float> pepMasses;
      aminoacids.getPRMMasses(psm->m_annotation.c_str(), pepMasses);

      vector<float> modifications;
      vector<unsigned int> positions;
      if (DEBUG_SPEC_PROB)
        DEBUG_VAR(psm->m_annotation);
      psm->getModificationsAndPositions(modifications, positions, false);
      if (DEBUG_SPEC_PROB)
        DEBUG_VAR(modifications.size());

      psm->m_pValue =
          tmpSpec.getGFProbability(matchScore,
                                   (unsigned int)round(pepMasses[pepMasses.size()
                                       - 1]),
                                   false,
                                   &aminoacids,
                                   modifications.size(), // max(2, modifications.size())
                                   &mods,
                                   &ntermMods,
                                   false,
                                   useExistingRawScore);

      if (DEBUG_SPEC_PROB)
        DEBUG_VAR(psm->m_pValue);

      // A little hack because 1e-38 is better than 0
      if (psm->m_pValue == 0.0)  {
        psm->m_pValue = 1e-38;
      }

      if (DEBUG_SPEC_PROB)
        DEBUG_MSG(scan << "\t" << psm->m_annotation << "\t" << psm->m_pValue);
    }

    return;
  }

  int Spectrum::compare(Spectrum &toSpec, vector<int> &diffIndexes)
  {
    int size1 = size();
    int size2 = toSpec.size();

    if (size1 != size2)
    {
      diffIndexes.push_back(-1);
      diffIndexes.push_back(size1);
      diffIndexes.push_back(size2);
      return -1;
    }

    for (unsigned i = 0; i < size1; i++)
    {
      bool mzCompResult = floatCompare(peakList[i][0],
                                       toSpec.peakList[i][0],
                                       4);
      bool intCompResult = floatCompare(peakList[i][1],
                                        toSpec.peakList[i][1],
                                        4);

      if (mzCompResult || intCompResult)
        diffIndexes.push_back(i);
    }

    return diffIndexes.size();
  }

  float Spectrum::scoreMatch(const Spectrum &spectrum_comparison,
                             float peakTol,
                             unsigned int & matchPeaks,
                             float & score1,
                             float & score2,
                             bool forceExact,
                             bool projectedCosine,
                             bool preprocess) const
  {

    vector<int> idxMatched1_zero, idxMatched1_other;
    idxMatched1_zero.reserve(1500);
    idxMatched1_other.reserve(1500);
    vector<int> idxMatched2_zero, idxMatched2_other;
    idxMatched2_zero.reserve(1500);
    idxMatched2_other.reserve(1500);

    vector<GPCAux> peakMatches;

    vector<char> peakUsed1, peakUsed2;

    Spectrum spec1 = *this;
    Spectrum spec2 = spectrum_comparison;

    //sqrt intensities
    if (preprocess == true)
    {
      for (unsigned int j = 0; j < spec1.size(); j++)
        spec1[j][1] = sqrt(spec1[j][1]);

      for (unsigned int j = 0; j < spec2.size(); j++)
        spec2[j][1] = sqrt(spec2[j][1]);

      //Noramlizing spectra
      spec1.normalize2();
      spec2.normalize2();
    }

    peakUsed1.resize(spec1.size());
    peakUsed2.resize(spec2.size());
    float pmDiff = spec1.parentMass - spec2.parentMass;

    if (forceExact)
    {
      pmDiff = 0.f;
    }

    // Greedy calculation of highest-possible cosine (maximum bipartite heuristic)
    FindMatchPeaksAll2_efficient(spec1,
                       spec2,
                       0,
                       peakTol,
                       idxMatched1_zero,
                       idxMatched2_zero);




    FindMatchPeaksAll2_efficient(spec1,
                       spec2,
                       pmDiff,
                       peakTol,
                       idxMatched1_other,
                       idxMatched2_other);


    //cerr<<"Got "<<idxMatched1_zero.size()<<"/"<<idxMatched2_zero.size()<<" (at 0) + "<<idxMatched1_other.size()<<"/"<<idxMatched2_other.size()<<" (at "<<pmDiff<<") matched peaks\n";

    peakMatches.resize(idxMatched1_zero.size() + idxMatched1_other.size());
    for (unsigned int i = 0; i < spec1.size(); i++)
      peakUsed1[i] = 0;
    for (unsigned int i = 0; i < spec2.size(); i++)
      peakUsed2[i] = 0;
    for (unsigned int i = 0; i < idxMatched1_zero.size(); i++)
    {
      peakMatches[i].i = idxMatched1_zero[i];
      peakMatches[i].j = idxMatched2_zero[i];
      peakMatches[i].score = spec1[idxMatched1_zero[i]][1]
          * spec2[idxMatched2_zero[i]][1];
    }
    for (unsigned int i = 0, j = idxMatched1_zero.size();
        i < idxMatched1_other.size(); i++, j++)
    {
      peakMatches[j].i = idxMatched1_other[i];
      peakMatches[j].j = idxMatched2_other[i];
      peakMatches[j].score = spec1[idxMatched1_other[i]][1]
          * spec2[idxMatched2_other[i]][1];
    }
    std::sort(peakMatches.begin(), peakMatches.end(), GPCAux_cmp);

    float cosine = 0;
    score1 = 0;
    score2 = 0;
    float matchedPeaksNorm1 = 0;
    float matchedPeaksNorm2 = 0;
    unsigned int numMatchedPeaks = 0;
    for (int i = (int)peakMatches.size() - 1; i >= 0; i--)
    {
      if (peakUsed1[peakMatches[i].i] == 0 and peakUsed2[peakMatches[i].j] == 0)
      {
        cosine += peakMatches[i].score;
        peakUsed1[peakMatches[i].i] = 1;
        peakUsed2[peakMatches[i].j] = 1;
        score1 += spec1[peakMatches[i].i][1];
        score2 += spec2[peakMatches[i].j][1];
        matchedPeaksNorm1 += spec1[peakMatches[i].i][1]
            * spec1[peakMatches[i].i][1];
        matchedPeaksNorm2 += spec2[peakMatches[i].j][1]
            * spec2[peakMatches[i].j][1];
        numMatchedPeaks++;
        //cerr<<" --- used ("<<peakMatches[i].i<<","<<specSet[spec1][peakMatches[i].i][0]<<") / ("<<peakMatches[i].j<<","<<specSet[spec2][peakMatches[i].j][0]<<")\n";
      }
      else
      {
        //cerr<<" --- skipped ("<<peakMatches[i].i<<","<<specSet[spec1][peakMatches[i].i][0]<<") / ("<<peakMatches[i].j<<","<<specSet[spec2][peakMatches[i].j][0]<<")\n";
      }
    }

    matchPeaks = numMatchedPeaks;

    if (not projectedCosine)
      cosine = 1 * cosine;
    else
      cosine = cosine / (matchedPeaksNorm1 * matchedPeaksNorm2);

    return cosine;
  }

  unsigned int Spectrum::filterLowSNR(float SNR_threshold)
  {
    float spectrum_noise_level = this->getNoiseLevel();
    this->filterLowIntensity(SNR_threshold * spectrum_noise_level);
    return 0;
  }

  float Spectrum::getNoiseLevel()
  {
    //Calculating Mean
    float total_int = 0.0;
    int total_count = 0;
    vector<float> peak_intensity_list;

    for (int peak_idx = 0; peak_idx < peakList.size(); peak_idx++)
    {
      peak_intensity_list.push_back(peakList[peak_idx][1]);
    }

    //sorting and finding bottom 25%
    sort(peak_intensity_list.begin(), peak_intensity_list.end());
    if (peak_intensity_list.size() < 10)
    {
      return 0.f;
    }

    //float bottom_threshold = peak_intensity_list[peak_intensity_list.size()/4];
    for (int peak_idx = 0; peak_idx < peak_intensity_list.size() / 4;
        peak_idx++)
    {
      total_int += peak_intensity_list[peak_idx];
      total_count++;
    }

    float mean = total_int / total_count;

    return mean;
  }

// MergeSpectra3 auxiliar variables:
//static vector<vector<int> > gaux_peaksIdx(0); // Indices of selected peaks per spectrum
//static vector<float> gaux_peaksScores(0); // Summed scores of the selected peaks per spectrum

// -------------------------------------------------------------------------

}
