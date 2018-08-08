/*
 * SpecSet.cpp
 *
 *  Created on: Aug 30, 2011
 *      Author: jsnedecor
 */

//
// **************************************************************
//    SpecSet methods
// **************************************************************
//
#include "SpecSet.h"
#include "tags.h"
//#include "mzxml.h"
//#include "PWizInterface.h"

namespace specnets
{

  // -------------------------------------------------------------------------
  SpecSet::SpecSet(unsigned int sz, bool enableIndexing) :
      m_index(0x0)
  {
    specs.resize(sz);

    if (enableIndexing)
    {
      m_index = new SpectrumMap(sz);
    }
  }
  // -------------------------------------------------------------------------
  SpecSet::SpecSet(char *filename) :
      m_index(0x0)
  {
    LoadSpecSet_pkl(filename);
  }

  SpecSet::~SpecSet()
  {
    if (m_index != 0x0)
    {
      delete m_index;
    }
  }

  // -------------------------------------------------------------------------
  void SpecSet::insert(vector<Spectrum>::iterator position,
                       vector<Spectrum>::iterator start,
                       vector<Spectrum>::iterator end)
  {
    specs.insert(position, start, end);

    if (m_index != 0x0)
    {
      unsigned curIdx = (*m_index)[position->getUniqueID()];
      for (vector<Spectrum>::iterator specIt = start; specIt != end; specIt++)
      {
        (*m_index)[specIt->getUniqueID()] = curIdx;
        curIdx++;
      }
    }
  }
  // -------------------------------------------------------------------------
  SpecSet &SpecSet::operator=(const SpecSet &other)
  {
    if (this == &other)
    {
      return *this;
    }
    specs.resize(other.specs.size());
    for (unsigned int i = 0; i < other.specs.size(); i++)
      specs[i] = other[i];

    if (m_index != 0x0)
    {
      computeIndex();
    }
    return *this;
  }
  // -------------------------------------------------------------------------
  void SpecSet::push_back(const Spectrum & x)
  {
    specs.push_back(x);

    if (m_index != 0x0)
    {
      (*m_index)[x.getUniqueID()] = size() - 1;
    }
  }

  void SpecSet::swap(SpecSet& other)
  {
    specs.swap(other.specs);
    if (m_index != 0x0)
    {
      computeIndex();
    }
  }

  void SpecSet::clear()
  {
    resize(0);
    vector<Spectrum>().swap(specs);
  }
  // -------------------------------------------------------------------------
  unsigned int SpecSet::resize(unsigned int newSize)
  {
    specs.resize(newSize);

    if (m_index != 0x0)
    {
      computeIndex();
    }

    return specs.size();
  }

  // -------------------------------------------------------------------------
  vector<Spectrum>::iterator SpecSet::begin(void)
  {
    return specs.begin();
  }
  // -------------------------------------------------------------------------
  vector<Spectrum>::iterator SpecSet::end(void)
  {
    return specs.end();
  }

  void SpecSet::setPeakTolerance(float tolerance, bool applyPPM)
  {
    for (int i = 0; i < size(); i++)
    {
      if (applyPPM)
      {
        specs[i].setPeakTolerancePPM(tolerance);
      }
      else
      {
        specs[i].setPeakTolerance(tolerance);
      }
    }
  }

  void SpecSet::setParentMassTolerance(float tolerance, bool applyPPM)
  {
    for (int i = 0; i < size(); i++)
    {
      if (applyPPM)
      {
        specs[i].setParentMassTolPPM(tolerance);
      }
      else
      {
        specs[i].setParentMassTol(tolerance);
      }
    }
  }

  void SpecSet::swapAppendSpecSet(SpecSet& other, bool updateScans)
  {
    unsigned int prevSize = specs.size();
    specs.resize(prevSize + other.specs.size());
    if (m_index != 0x0)
    {
      m_index->rehash(specs.size());
    }
    for (unsigned int i = 0; i < other.specs.size(); i++)
    {
      specs[prevSize + i].swap(other.specs[i]);
      if (updateScans)
      {
        specs[prevSize + i].scan = prevSize + other[i].scan;
      }
      if (m_index != 0x0)
      {
        (*m_index)[specs[prevSize + i].getUniqueID()] = prevSize + i;
      }
    }
    other.clear();
  }

  // -------------------------------------------------------------------------
  void SpecSet::appendSpecSet(const SpecSet& other, bool updateScans)
  {
    unsigned int prevSize = specs.size();
    specs.resize(prevSize + other.specs.size());
    if (m_index != 0x0)
    {
      m_index->rehash(specs.size());
    }
    for (unsigned int i = 0; i < other.specs.size(); i++)
    {
      specs[prevSize + i] = other.specs[i];
      if (updateScans)
      {
        specs[prevSize + i].scan = prevSize + other[i].scan;
      }
      if (m_index != 0x0)
      {
        (*m_index)[specs[prevSize + i].getUniqueID()] = prevSize + i;
      }
    }
  }
  // -------------------------------------------------------------------------
  vector<list<int> > &SpecSet::getMassesHistogram(vector<list<int> > &results,
                                                  float resolution) const
  {
    int i;
    double t, massIdx = 0; // Value is always integer but max is (double,double)

    for (i = 0; i < specs.size(); i++)
    {
      t = round(specs[i].parentMass / resolution);
      if (massIdx < t)
        massIdx = t;
    }
    results.resize((int)massIdx + 1);

    for (i = 0; i < specs.size(); i++)
      results[(int)(round(specs[i].parentMass / resolution))].push_front(i);
    return results;
  }
  // -------------------------------------------------------------------------
  void SpecSet::setResolution(float newResolution, bool enforceRounding)
  {
    for (unsigned int i = 0; i < specs.size(); i++)
      specs[i].setResolution(newResolution, enforceRounding);
  }
  // -------------------------------------------------------------------------
  void SpecSet::setFilename(string &filename)
  {
    for (unsigned int i = 0; i < specs.size(); i++)
    {
      specs[i].fileName = filename;
    }
  }
  // -------------------------------------------------------------------------
  void SpecSet::addZPMpeaks(float tolerance, float ionOffset, bool includeY0k)
  {
    for (unsigned int i = 0; i < specs.size(); i++)
      specs[i].addZPMpeaks(tolerance, ionOffset, includeY0k);
  }

  // -------------------------------------------------------------------------
  float SpecSet::averageIntensity()
  {
    float avgIntensity = 0.0;
    float totalPeaks = 0.0;
    for (int s = 0; s < specs.size(); s++)
    {
      totalPeaks += specs[s].size();
      for (int p = 0; p < specs[s].size(); p++)
      {
        avgIntensity += specs[s][p][1];
      }
    }

    if (totalPeaks == 0.0)
    {
      return 0.0;
    }

    avgIntensity /= totalPeaks;
    return avgIntensity;
  }

  // -------------------------------------------------------------------------
  void SpecSet::normalize(float newTotalInt)
  {
    for (unsigned int i = 0; i < specs.size(); i++)
      specs[i].normalize(newTotalInt);
  }
  // -------------------------------------------------------------------------
  void SpecSet::clearPsms(void)
  {
    for (int i = 0; i < specs.size(); i++)
    {
      specs[i].psmList.clear();
    }
  }
  // -------------------------------------------------------------------------
  template<class T> unsigned int SpecSet::extract(vector<T> &features,
                                                  T featureValue,
                                                  SpecSet &output)
  {
    output.resize(max(features.size(), specs.size()));
    unsigned int keptSpectra = 0, pivot;

    for (pivot = 0; pivot < output.size(); pivot++)
      if (features[pivot] == featureValue)
        output[keptSpectra++] = specs[pivot];

    output.resize(keptSpectra);
    return keptSpectra;
  }
  // -------------------------------------------------------------------------
  void instantiate_SpecSet_extract()
  {
    SpecSet specs, specsOut;

    vector<int> features;
    specs.extract(features, (int)0, specsOut);
  }
  // -------------------------------------------------------------------------
  bool SpecSet::maximumParsimony(void)
  {
    // Initialize assignment vector
    vector<psmPtr> assignments(specs.size());
    for (int pepIdx = 0; pepIdx < specs.size(); pepIdx++)
    {
      psmPtr p(new PeptideSpectrumMatch);
      assignments[pepIdx] = p;
      assignments[pepIdx]->m_dbIndex = -1;
    }

    //DEBUG_TRACE;
    int protIdx = 0;
    // Find maximum index of a matched protein
    list<psmPtr>::iterator iter;
    for (int pepIdx = 0; pepIdx < specs.size(); pepIdx++)
    {
      for (iter = specs[pepIdx].psmList.begin();
          iter != specs[pepIdx].psmList.end(); iter++)
      {
        protIdx = max((int)protIdx, (*iter)->m_dbIndex);
      }
    }
    vector<list<psmPtr> > pepsPerProt(protIdx + 1);
    //DEBUG_VAR(pepsPerProt.size());

    unsigned int maxHits, maxHitsIdx;
    float maxHitsMass;
    do
    {
      maxHits = 0;
      maxHitsIdx = 0;

      // Clear the peptides per protein vector
      for (int protIdx = 0; protIdx < pepsPerProt.size(); protIdx++)
        pepsPerProt[protIdx].clear();

      // Create the lists of peptides per protein and compute the max
      for (int pepIdx = 0; pepIdx < specs.size(); pepIdx++)
      {
        for (iter = specs[pepIdx].psmList.begin();
            iter != specs[pepIdx].psmList.end(); iter++)
        {
          int dbIndex = (*iter)->m_dbIndex;
          pepsPerProt[dbIndex].push_back(*iter);
          if (pepsPerProt[dbIndex].size() > maxHits)
          {
            maxHits = pepsPerProt[dbIndex].size();
            maxHitsIdx = dbIndex;
          }
          //cout << (*iter)->m_scanNum << ", " << dbIndex << endl;
        }
      }
      //DEBUG_VAR(maxHits);
      //DEBUG_VAR(maxHitsIdx);

      if (maxHits > 0)
      {
        for (iter = pepsPerProt[maxHitsIdx].begin();
            iter != pepsPerProt[maxHitsIdx].end(); iter++)
        {
          Spectrum * spectrum = (*iter)->m_spectrum;
          if (spectrum == 0)
          {
            ERROR_MSG("Unable to get scan [" << (*iter)->m_scanNum << "]");
            return false;
          }
          spectrum->psmList.clear();
          assignments[(*iter)->m_scanNum - 1] = (*iter);
        }
      }

    } while (maxHits > 0);

    //DEBUG_TRACE;
    // Copy the assignments back to the spectrum
    for (int pepIdx = 0; pepIdx < specs.size(); pepIdx++)
    {
      if (assignments[pepIdx]->m_dbIndex != -1)
      {
        specs[pepIdx].psmList.clear();
        specs[pepIdx].psmList.push_back(assignments[pepIdx]);
      }
    }

    //DEBUG_TRACE;
    return true;
  }

  unsigned int SpecSet::extractSpectra(SpecSet& outputSpecs,
                                       Spectrum::FragType fragType)
  {

    outputSpecs.resize(size());
    unsigned int idxUse = 0;

    for (unsigned int i = 0; i < size(); i++)
    {
      if (specs[i].msFragType == fragType)
      {
        outputSpecs[idxUse] = specs[i];
        ++idxUse;
      }
    }
    outputSpecs.resize(idxUse);
    return outputSpecs.size();
  }

  unsigned int SpecSet::swapExtractSpectra(SpecSet& outputSpecs,
                                           Spectrum::FragType fragType)
  {
    outputSpecs.resize(0);
    outputSpecs.resize(size());
    unsigned int idxUse = 0;

    for (unsigned int i = 0; i < size(); i++)
    {
      if (specs[i].msFragType == fragType)
      {
        outputSpecs.specs[idxUse++].swap(specs[i]);
      }
    }
    outputSpecs.resize(idxUse);
    return outputSpecs.size();
  }

  unsigned int SpecSet::extractPairedScans(vector<vector<unsigned int> >& outputPairs,
                                           int numConsec,
                                           int maxConsec)
  {
    outputPairs.resize(size());
    if (size() == 0)
    {
      return 0;
    }

    for (unsigned int i = 0; i < outputPairs.size(); i++)
    {
      outputPairs[i].resize(0);
    }

    unsigned int idxUse = 0;

    Spectrum* prevSpec = &specs[0];
    Spectrum* nextSpec;
    outputPairs[idxUse].push_back(prevSpec->scan);

    unsigned int numMerged = 0;

    int curConsec = 1;

    for (unsigned int i = 1; i < size(); i++)
    {
      nextSpec = &specs[i];

      if ((maxConsec <= 0 || curConsec < maxConsec)
          && prevSpec->parentMass == nextSpec->parentMass
          && (numConsec <= 0 || (numConsec > 0 && curConsec < numConsec)))
      {
        curConsec++;
        outputPairs[idxUse].push_back(nextSpec->scan);
        ++numMerged;
      }
      else
      {
        ++idxUse;
        outputPairs[idxUse].push_back(nextSpec->scan);
        curConsec = 1;
      }
      prevSpec = nextSpec;
    }
    outputPairs.resize(idxUse + 1);
    return size() - numMerged;
  }

  // -------------------------------------------------------------------------
  int SpecSet::saveMatchedProts(const char *filename)
  {
    vector<vector<int> > tempMatchedProts(specs.size());
    for (unsigned int i = 0; i < specs.size(); i++)
    {
      tempMatchedProts[i].resize(3);
      if (specs[i].psmList.size() == 0)
      {
        tempMatchedProts[i][0] = -1;
        tempMatchedProts[i][1] = 0;
        tempMatchedProts[i][2] = 0;
      }
      else
      {
        tempMatchedProts[i][0] = specs[i].psmList.front()->m_dbIndex;
        tempMatchedProts[i][1] = specs[i].psmList.front()->m_numMods;
        tempMatchedProts[i][2] = specs[i].psmList.front()->m_matchOrientation;
      }
    }
    return Save_binArray(filename, tempMatchedProts);
  }
  // -------------------------------------------------------------------------
  unsigned int SpecSet::Load(const char *filename, const char *ext)
  {
    // split the filename into components
    FilenameManager fm(filename);

    string extension;
    if (ext)
      extension = ext;
    else
      extension = fm.extension;

    fm.lowerCaseExtension();

    // select load method according to extension
    if (extension.compare("mgf") == 0)
      return LoadSpecSet_mgf(filename);

    if (extension.compare("pkl") == 0)
      return LoadSpecSet_pkl(filename);

    if (extension.compare("pklbin") == 0)
      return loadPklBin(filename);

    if (extension.compare("ms2") == 0)
      return LoadSpecSet_ms2(filename);

    if (extension.compare("prmsv3") == 0)
      return LoadSpecSet_prmsv3(filename);

    if (extension.compare("prms") == 0)
      return LoadSpecSet_prms(filename);

    //if(extension.compare("mzxml") == 0)  {
    //  vector<short> msLevel;
    //  try {
    //    return LoadMzxml(filename, *this, &msLevel, 2);
    //  } catch (...) {
    //  }
    //}

    //return loadDataUsingPWiz(fm.filenameFull, *this);

    return 0;
  }
  // -------------------------------------------------------------------------
  bool SpecSet::LoadSpecSet_pkl_mic(const char* filename)
  {
    BufferedLineReader blr;
    TwoValues<float> peak;
    Spectrum loadSpec;
    specs.resize(0);
    if (blr.Load(filename) <= 0 || blr.size() < 2)
    {
      cerr << "ERROR: Not enough lines loaded\n";
      return false;
    }
    bool inspec = false;
    for (unsigned int i = 0; i < blr.size(); i++)
    {

      vector<string> peak_vals;
      if (blr.getline(i) == NULL)
        break;
      const char* next_line = blr.getline(i);
      if (!splitText(next_line, peak_vals, "\t "))
      {
        cerr << "ATTENTION: stopped reading peaks before eof in " << filename
            << "\n";
        break;
      }
      if (peak_vals.size() == 0)
      {
        if (inspec)
        {
          loadSpec.sortPeaks();
          specs.push_back(loadSpec);
        }
        inspec = false;
        continue;
      }
      if (peak_vals.size() < 2)
      {
        cerr << "ATTENTION: found wierd line after header: <" << next_line
            << "> ---- skippping ...\n";
        continue;
      }

      if (!inspec && peak_vals.size() == 3)
      {
        loadSpec.resize(0);
        loadSpec.psmList.resize(0);
        loadSpec.parentMass = (float)atof(peak_vals[0].c_str());
        loadSpec.parentMZ = (double)strtod(peak_vals[0].c_str(), NULL);
        loadSpec.parentCharge = (short)atoi(peak_vals[2].c_str());
        if (loadSpec.parentCharge > 0)
        {
          loadSpec.parentMass = loadSpec.parentMass
              * ((float)loadSpec.parentCharge) - ((float)loadSpec.parentCharge)
              + AAJumps::massHion;
        }
        inspec = true;
        continue;
      }

      if (inspec && peak_vals.size() < 2)
      {
        cerr << "ATTENTION: not enough peak info after header: <" << next_line
            << "> ---- skippping ...\n";
        continue;
      }

      if (inspec)
      {
        peak[0] = (float)atof(peak_vals[0].c_str());
        peak[1] = (float)atof(peak_vals[1].c_str());
        loadSpec.insertPeak(peak[0], peak[1], 0);
      }
    }
    if (inspec)
    {
      loadSpec.sortPeaks();
      specs.push_back(loadSpec);
    }

    if (m_index != 0x0)
    {
      computeIndex();
    }
    return true;
  }
  // -------------------------------------------------------------------------
  Spectrum* SpecSet::getScan(int scan_num)
  {
    //cheating a bit, check to see whether scans are in order and we can just jump to the correct
    //scan
    if (scan_num < specs.size() && scan_num > 0)
    { // only try this if scan number is within range.
      Spectrum *curr_scan = &(specs[scan_num - 1]);

      if (curr_scan->scan == scan_num)
      {
        return curr_scan;
      }
    }

    vector<Spectrum>::iterator spec_iter;
    for (spec_iter = specs.begin(); spec_iter != specs.end(); spec_iter++)
    {
      if (scan_num == spec_iter->scan)
      {
        return &(*spec_iter);
      }
    }
    return (Spectrum*)NULL;
  }

  int SpecSet::getIndexFromID(const string& specID)
  {
    if (m_index == 0 || m_index->count(specID) == 0)
    {
      return -1;
    }
    else
    {
      return (*m_index)[specID];
    }
  }

  Spectrum* SpecSet::getIndex(const string& specID)
  {
    if (m_index == 0 || m_index->count(specID) == 0)
    {
      return (Spectrum*)0;
    }
    else
    {
      return &specs[(*m_index)[specID]];
    }
  }

  const Spectrum* SpecSet::getIndex(const string& specID) const
  {
    if (m_index == 0 || m_index->count(specID) == 0)
    {
      return (Spectrum*)0;
    }
    else
    {
      return &specs[(*m_index)[specID]];
    }
  }

  void SpecSet::index()
  {
    if (m_index == 0x0)
    {
      m_index = new SpectrumMap(size());
    }
    computeIndex();
  }

  unsigned int SpecSet::LoadSpecSet_mgf(const char *filename,
                                        unsigned int scan,
                                        int index)
  {
    BufferedLineReader blr;
    resize(0);
    //cout << "HERE" << endl;
    if (blr.Load(filename) <= 0 or blr.size() < 3)
      return 0; // A valid file must have at least 3 lines: BEGIN_IONS, END_IONS and one mass/intensity peak

    unsigned int lineIdx, specIdx, peakIdx;

    // Counts number of spectra and number of peaks per spectrum
    list<int> peaksPerSpec;
    int extract_index = index;
    int first;
    int numPeaks = 0; // Counts number of peaks in the spectrum
    for (lineIdx = 0; lineIdx < blr.size(); lineIdx++)
    {
      if (strncmp("END IONS", blr.getline(lineIdx), 8) == 0)
      {
        peaksPerSpec.push_back(numPeaks);
        numPeaks = -1;
        continue;
      }
      if (strncmp("BEGIN IONS", blr.getline(lineIdx), 10) == 0)
      {
        numPeaks = 0;
        continue;
      }
      //That means we need to determine an index
      if (scan > 0 && index == -1)
      {
        if (strncmp("SCANS=", blr.getline(lineIdx), sizeof("SCANS=") - 1) == 0)
        {
          //cout << blr.getline(lineIdx) << endl;
          unsigned int current_scan = (unsigned int)atoi(blr.getline(lineIdx)
              + sizeof("SCANS=") - 1);
          if (current_scan == scan)
          {
            extract_index = peaksPerSpec.size();
            //cout << (extract_index) << endl;
          }
        }
      }
      else
      {
        //Failure to specify the correct scan index combo
        return -1;
      }
      first = (int)blr.getline(lineIdx)[0];
      if (numPeaks >= 0 and first >= 48 and first <= 57)
        numPeaks++;
    }

    // Parse spectra
    resize(peaksPerSpec.size());
    lineIdx = 0;
    char *token, *line;
    for (specIdx = 0; specIdx < specs.size(); specIdx++)
    {
      // Skip empty lines
      while (lineIdx < blr.size()
          and (blr.getline(lineIdx)[0] == 0
              or strncmp("BEGIN IONS", blr.getline(lineIdx), 10) != 0))
        lineIdx++;

      if (lineIdx == blr.size())
      {
        cerr << "Error loading " << filename << " - " << specIdx
            << " spectra instead of " << specs.size() << "?\n";
        resize(0);
        return 0;
      }

      // Start of spectrum
      if (strncmp("BEGIN IONS", blr.getline(lineIdx), 10) != 0)
      {
        cerr << "ERROR: Expected BEGIN IONS, found '" << blr.getline(lineIdx)
            << "' (line " << lineIdx + 1 << ")\n";
        return 0;
      }
      else
        lineIdx++;

      // Read peaks/charge/parent mass
      if (specIdx == extract_index)
      {
        specs[specIdx].resize(peaksPerSpec.front());
      }
      peaksPerSpec.pop_front();
      peakIdx = 0;
      specs[specIdx].parentCharge = 0;
      specs[specIdx].parentMassTol = 0;
      specs[specIdx].scan = 0;
      while (lineIdx < blr.size()
          and strncmp("END IONS", blr.getline(lineIdx), 8) != 0)
      {

        line = blr.getline(lineIdx++);

        if (specIdx != extract_index)
        {
          if (strncmp("SCANS=", line, sizeof("SCANS=") - 1) == 0)
          {
            specs[specIdx].scan = (unsigned int)atoi(line + sizeof("SCANS=")
                - 1);
          }
          continue;
        }

        if (line[0] >= 48 and line[0] <= 57)
        {
          token = strtok(line, " \t");
          if (!token)
          {
            cerr << "Error loading " << filename
                << " - could not parse peak mass on line " << lineIdx << "!\n";
            resize(0);
            return 0;
          }
          specs[specIdx][peakIdx][0] = atof(token);
          token = strtok(NULL, " \t");
          if (!token)
          {
            cerr << "Error loading " << filename
                << " - could not parse peak intensity on line " << lineIdx
                << "!\n";
            resize(0);
            return 0;
          }
          specs[specIdx][peakIdx][1] = atof(token);
          peakIdx++;
          continue;
        }

        if (strncmp("CHARGE=+", line, 8) == 0)
        {
          specs[specIdx].parentCharge = (short)atof(&line[8]);
        }
        else if (strncmp("CHARGE=", line, 7) == 0)
        {
          specs[specIdx].parentCharge = (short)atof(&line[7]);
        }

        if (strncmp("TITLE=Scan Number: ", line, 19) == 0)
          specs[specIdx].scan = (unsigned int)atof(&line[19]);

        if (strncmp("SCANS=", line, sizeof("SCANS=") - 1) == 0)
          specs[specIdx].scan = (unsigned int)atoi(line + sizeof("SCANS=") - 1);
        if (strncmp("PEPMASS=", line, 8) == 0)
        {
          specs[specIdx].parentMass = (double)strtod(&line[8], NULL);
          specs[specIdx].parentMZ = (double)strtod(&line[8], NULL);
        }
        if (strncmp("PRECURSOR=", line, 10) == 0)
        {
          specs[specIdx].parentMZ = (double)strtod(&line[10], NULL);
        }
        if (strncmp("RTINSECONDS=", line, 10) == 0)
        {
          specs[specIdx].retention_time = atof(line + sizeof("RTINSECONDS=")
              - 1);
        }
        if (strncmp("PEPTIDE=", line, sizeof("PEPTIDE=") - 1) == 0)
        {
          psmPtr psm(new PeptideSpectrumMatch());
          psm->m_annotation = line + sizeof("PEPTIDE=") - 1;
          specs[specIdx].psmList.push_back(psm);
        }
        if (strncmp("SEQ=", line, sizeof("SEQ=") - 1) == 0)
        {
          psmPtr psm(new PeptideSpectrumMatch());
          psm->m_annotation = line + sizeof("SEQ=") - 1;
          specs[specIdx].psmList.push_back(psm);
          specs[specIdx].psmList.front()->m_spectrumFile = filename;
          specs[specIdx].psmList.front()->m_dbIndex = specIdx + 1;
        }
        if (strncmp("MSLEVEL=", line, sizeof("MSLEVEL=") - 1) == 0)
        {
          specs[specIdx].msLevel = atoi(line + sizeof("MSLEVEL=") - 1);
        }
        if (strncmp("SOURCE_INSTRUMENT=",
                    line,
                    sizeof("SOURCE_INSTRUMENT=") - 1) == 0)
        {
          specs[specIdx].instrument_name = line + sizeof("SOURCE_INSTRUMENT=")
              - 1;
        }
        if (strncmp("ITOL=", line, sizeof("ITOL=") - 1) == 0)
        {
          specs[specIdx].ITOL = atof(line + sizeof("ITOL=") - 1);
        }
        if (strncmp("ITOLU=", line, sizeof("ITOLU=") - 1) == 0)
        {
          specs[specIdx].ITOLU = line + sizeof("ITOLU=") - 1;
        }
        if (strncmp("TOL=", line, sizeof("TOL=") - 1) == 0)
        {
          specs[specIdx].TOL = atof(line + sizeof("TOL=") - 1);
        }
        if (strncmp("TOLU=", line, sizeof("TOLU=") - 1) == 0)
        {
          specs[specIdx].TOLU = line + sizeof("TOLU=") - 1;
        }
        if (strncmp("FILENAME=", line, sizeof("FILENAME=") - 1) == 0)
        {
          specs[specIdx].fileName = line + sizeof("FILENAME=") - 1;
        }
        if (strncmp("SPECTRUMQUALITY=", line, sizeof("SPECTRUMQUALITY=") - 1)
            == 0)
        {
          specs[specIdx].spectrum_quality = atoi(line
              + sizeof("SPECTRUMQUALITY=") - 1);
        }
        if (strncmp("SUBMISSION_METADATA=",
                    line,
                    sizeof("SUBMISSION_METADATA=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string submission_metadata = line
                + sizeof("SUBMISSION_METADATA=") - 1;
            specs[specIdx].psmList.front()->m_submission_metadata =
                (submission_metadata);
          }
        }
        if (strncmp("ORGANISM=", line, sizeof("ORGANISM=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string organism = line + sizeof("ORGANISM=") - 1;
            specs[specIdx].psmList.front()->m_organism = (organism);
          }
        }
        //ACCEPTING BOTH COMPOUNDNAME and NAME
        if (strncmp("COMPOUNDNAME=", line, sizeof("COMPOUNDNAME=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string compound = line + sizeof("COMPOUNDNAME=") - 1;
            specs[specIdx].psmList.front()->m_compound_name = (compound);
          }
        }
        if (strncmp("NAME=", line, sizeof("NAME=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string compound = line + sizeof("NAME=") - 1;
            specs[specIdx].psmList.front()->m_compound_name = (compound);
          }
        }
        if (strncmp("SMILES=", line, sizeof("SMILES=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string smiles = line + sizeof("SMILES=") - 1;
            specs[specIdx].psmList.front()->m_smiles = (smiles);
          }
        }
        if (strncmp("INCHI=", line, sizeof("INCHI=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string InChI = line + sizeof("INCHI=") - 1;
            specs[specIdx].psmList.front()->m_InChI = (InChI);
          }
        }
        if (strncmp("INCHIAUX=", line, sizeof("INCHIAUX=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string InChI_aux = line + sizeof("INCHIAUX=") - 1;
            specs[specIdx].psmList.front()->m_InChI_Aux = (InChI_aux);
          }
        }
        if (strncmp("NOTES=", line, sizeof("NOTES=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string notes = line + sizeof("NOTES=") - 1;
            specs[specIdx].psmList.front()->m_notes = notes;
          }
        }
        if (strncmp("IONMODE=", line, sizeof("IONMODE=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string ionmode = line + sizeof("IONMODE=") - 1;
            specs[specIdx].psmList.front()->m_ionmode = ionmode;
          }
        }

        if (strncmp("SLGF=", line, sizeof("SLGF=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string slgf_prob_pair = line + sizeof("SLGF=") - 1;

            float cosine_bin =
                atof(slgf_prob_pair.substr(0, slgf_prob_pair.find(":")).c_str());
            float cosine_prob =
                atof(slgf_prob_pair.substr(slgf_prob_pair.find(":") + 1).c_str());
            pair<float, float> cosine_prob_pair;
            cosine_prob_pair.first = cosine_bin;
            cosine_prob_pair.second = cosine_prob;
            specs[specIdx].psmList.front()->SLGF_distribution.push_back(cosine_prob_pair);
          }
        }

        if (strncmp("ACTIVATION=", line, sizeof("ACTIVATION=") - 1) == 0)
        {
          string activ = line + sizeof("ACTIVATION=") - 1;
          specs[specIdx].msFragType = Spectrum::parseActivation(activ);
        }

        if (strncmp("EXACTMASS=", line, sizeof("EXACTMASS=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string exactmass = line + sizeof("EXACTMASS=") - 1;
            specs[specIdx].psmList.front()->m_exactmass =
                atof(exactmass.c_str());
          }
        }

        if (strncmp("LIBRARYQUALITY=", line, sizeof("LIBRARYQUALITY=") - 1)
            == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string lib_quality = line + sizeof("LIBRARYQUALITY=") - 1;
            specs[specIdx].psmList.front()->m_library_quality =
                atoi(lib_quality.c_str());
          }
        }
        if (strncmp("SPECTRUMID=", line, sizeof("SPECTRUMID=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string spec_ID = line + sizeof("SPECTRUMID=") - 1;
            specs[specIdx].psmList.front()->m_spectrumID = spec_ID;
          }
        }

      }

      if (specs[specIdx].parentCharge > 0)
        specs[specIdx].parentMass = (specs[specIdx].parentMass
            * specs[specIdx].parentCharge)
            - (AAJumps::massHion * (specs[specIdx].parentCharge - 1));
      if (specs[specIdx].scan <= 0)
      {
        specs[specIdx].scan = specIdx + 1;
      }
      specs[specIdx].sortPeaks();
      lineIdx++;
    }

    if (m_index != 0x0)
    {
      computeIndex();
    }

    return specs.size();
  }

  unsigned SpecSet::FindSpectrumInMgfByScanOrIndex(const char * filename,
                                                   unsigned scan,
                                                   int index,
                                                   unsigned int &filestartposition)
  {
    //Testing Code
    std::fstream fs;
    fs.open(filename, std::fstream::in);

    string temp_line;
    unsigned int cur_index = 0;

    unsigned int start_spectrum_file_position = 0;
    unsigned int running_total_file_position = 0;
    while (fs.good())
    {
      getline(fs, temp_line);

      //Checking if begin ions
      if (temp_line.compare(0, sizeof("BEGIN IONS") - 1, "BEGIN IONS") == 0)
      {
        start_spectrum_file_position = running_total_file_position;
        cur_index++;

        if (cur_index == index && index > 0)
        {
          filestartposition = start_spectrum_file_position;
          return cur_index;
        }
      }

      //Checking if we are doing scan lookup
      if (temp_line.compare(0, sizeof("SCANS=") - 1, "SCANS=") == 0)
      {
        string scan_string = temp_line.substr(sizeof("SCANS=") - 1);
        unsigned int current_scan = (unsigned int)atoi(scan_string.c_str());

        if (current_scan == scan && scan > 0)
        {
          filestartposition = start_spectrum_file_position;
          return cur_index;
        }
      }

      //Updating the file position
      running_total_file_position += temp_line.length() + 1;

    }

    return -1;
  }

  unsigned SpecSet::FindSpectrumInMgfByScan(BufferedLineReader &blr,
                                            unsigned &lineIdx,
                                            unsigned scan)
  {
    unsigned index = 0;
    unsigned lastBeginLine = -1;
    // Cycle thu all read lines
    for (lineIdx = 0; lineIdx < blr.size(); lineIdx++)
    {
      // Continue in case an END or BEGIN is found
      //if (strncmp("END IONS", blr.getline(lineIdx), 8) == 0)
      //  continue;
      if (strncmp("BEGIN IONS", blr.getline(lineIdx), 10) == 0)
      {
        lastBeginLine = lineIdx;
        index++;
        continue;
      }

      // Check if it is a SCANS string
      if (strncmp("SCANS=", blr.getline(lineIdx), sizeof("SCANS=") - 1) == 0)
      {

        unsigned int current_scan = (unsigned int)atoi(blr.getline(lineIdx)
            + sizeof("SCANS=") - 1);
        if (current_scan == scan)
        {
          lineIdx = lastBeginLine;
          return index - 1;
        }
      }
    }

    // scan not found
    return -1;
  }

  bool SpecSet::FindSpectrumInMgfByIndex(BufferedLineReader &blr,
                                         unsigned &lineIdx,
                                         unsigned index)
  {
    unsigned currentIndex = 0;

    // Cycle thu all read lines
    for (lineIdx = 0; lineIdx < blr.size(); lineIdx++)
    {
      // Continue in case an END or BEGIN is found
      if (strncmp("END IONS", blr.getline(lineIdx), 8) == 0)
        continue;
      if (strncmp("BEGIN IONS", blr.getline(lineIdx), 10) == 0)
      {
        if (currentIndex == index)
          return true;
        currentIndex++;
      }
    }

    // scan not found
    return false;
  }

  unsigned int SpecSet::Calculate_MGF_index(const char *filename,
                                            vector<vector<unsigned int> > &file_index)
  {
    //Testing Code
    std::fstream fs;
    fs.open(filename, std::fstream::in);

    string temp_line;
    unsigned int cur_index = 0;

    unsigned int start_spectrum_file_position = 0;
    unsigned int running_total_file_position = 0;

    while (fs.good())
    {
      getline(fs, temp_line);

      //Checking if begin ions
      if (temp_line.compare(0, sizeof("BEGIN IONS") - 1, "BEGIN IONS") == 0)
      {
        start_spectrum_file_position = running_total_file_position;
        cur_index++;
      }

      //Checking if we are doing scan lookup
      if (temp_line.compare(0, sizeof("SCANS=") - 1, "SCANS=") == 0)
      {
        string scan_string = temp_line.substr(sizeof("SCANS=") - 1);
        unsigned int current_scan = (unsigned int)atoi(scan_string.c_str());

        vector<unsigned int> spectrum_index_object;
        spectrum_index_object.push_back(current_scan);
        spectrum_index_object.push_back(start_spectrum_file_position);
        file_index.push_back(spectrum_index_object);
      }

      //Updating the file position
      running_total_file_position += temp_line.length() + 1;

    }

    return 0;

  }

  unsigned int SpecSet::LoadSpectrum_mgf(const char *filename,
                                         unsigned int scan,
                                         int index,
                                         int use_index)
  {
    unsigned int filestartposition;
    if (use_index == 1)
    {
      //First check if the index exists
      string index_name(filename);
      index_name += ".specplotindex";

      //each row represents a spectrm in the file.
      //The vector represents in index 0 the spectrum scan, then the offset into the file
      vector<vector<unsigned int> > index_object;

      Load_binArray(index_name.c_str(), index_object);

      if (index_object.size() == 0)
      {
        cout << "Problems loading index, lets calculate it" << endl;
        Calculate_MGF_index(filename, index_object);
        string temp_name = index_name;
        srand(time(NULL));
        stringstream ss(stringstream::in | stringstream::out);
        ss << rand();
        temp_name += ss.str();
        cout << "Saving index: " << temp_name << endl;
        Save_binArray(temp_name.c_str(), index_object);
        cout << "Moving index to: " << index_name << endl;
        rename(temp_name.c_str(), index_name.c_str());
      }

      if (scan > 0)
      {
        for (int i = 0; i < index_object.size(); i++)
        {
          if (index_object[i][0] == scan)
          {
            filestartposition = index_object[i][1];
          }
        }
      }
      else
      {
        if (index > index_object.size() || index <= 0)
        {
          return 0;
        }
        else
        {
          //Index is valid
          filestartposition = index_object[index - 1][1];
        }
      }
    }
    else
    {
      int ret_val = FindSpectrumInMgfByScanOrIndex(filename,
                                                   scan,
                                                   index,
                                                   filestartposition);
      if (ret_val == -1)
      {
        return 0;
      }
    }

    std::fstream fs;
    fs.open(filename, std::fstream::in);
    fs.seekg(filestartposition);

    string temp_line;

    unsigned int specIdx = 0, peakIdx;

    // Parse spectra
    resize(1);
    specs[specIdx].resize(0);

    char *token;
    char * line = new char[8096];

    // Read peaks/charge/parent mass
    peakIdx = 0;
    specs[specIdx].parentCharge = 0;
    specs[specIdx].parentMassTol = 0;
    specs[specIdx].scan = 0;

    while (fs.good())
    {
      getline(fs, temp_line);
      //Checking if begin ions
      if (temp_line.compare(0, sizeof("END IONS") - 1, "END IONS") == 0)
      {
        break;
      }

      //line = temp_line.c_str();
      temp_line.copy(line, temp_line.length());
      line[temp_line.length()] = '\0';

      // SCANS
      if (strncmp("SCANS=", line, sizeof("SCANS=") - 1) == 0)
      {
        specs[specIdx].scan = (unsigned int)atoi(line + sizeof("SCANS=") - 1);
        continue;
      }

      if (line[0] >= 48 and line[0] <= 57)
      {
        specs[specIdx].resize(peakIdx + 1);
        token = strtok(line, " \t");
        if (!token)
        {
          ERROR_MSG("Error loading " << filename << " - could not parse peak mass on line " << temp_line << "!");
          resize(0);
          return 0;
        }
        specs[specIdx][peakIdx][0] = atof(token);
        token = strtok(NULL, " \t");
        if (!token)
        {
          ERROR_MSG("Error loading " << filename << " - could not parse peak mass on line " << temp_line << "!");
          resize(0);
          return 0;
        }
        specs[specIdx][peakIdx][1] = atof(token);
        peakIdx++;
        continue;
      }

      if (strncmp("CHARGE=+", line, 8) == 0)
      {
        specs[specIdx].parentCharge = (short)atof(&line[8]);
        continue;
      }
      else if (strncmp("CHARGE=", line, 7) == 0)
      {
        specs[specIdx].parentCharge = (short)atof(&line[7]);
        continue;
      }

      if (strncmp("TITLE=Scan Number: ", line, 19) == 0)
      {
        specs[specIdx].scan = (unsigned int)atof(&line[19]);
        continue;
      }

      if (strncmp("PEPMASS=", line, 8) == 0)
      {
        specs[specIdx].parentMass = (double)strtod(&line[8], NULL);
        specs[specIdx].parentMZ = (double)strtod(&line[8], NULL);
        continue;
      }

      if (strncmp("PRECURSOR=", line, 10) == 0)
      {
        specs[specIdx].parentMZ = (double)strtod(&line[10], NULL);
        continue;
      }

      if (strncmp("PEPTIDE=", line, sizeof("PEPTIDE=") - 1) == 0)
      {
        psmPtr psm(new PeptideSpectrumMatch());
        psm->m_annotation = line + sizeof("PEPTIDE=") - 1;
        specs[specIdx].psmList.push_back(psm);
        continue;
      }

      if (strncmp("SEQ=", line, sizeof("SEQ=") - 1) == 0)
      {
        psmPtr psm(new PeptideSpectrumMatch());
        psm->m_annotation = line + sizeof("SEQ=") - 1;
        specs[specIdx].psmList.push_back(psm);
        specs[specIdx].psmList.front()->m_spectrumFile = filename;
        specs[specIdx].psmList.front()->m_dbIndex = specIdx + 1;
        continue;
      }

      if (strncmp("MSLEVEL=", line, sizeof("MSLEVEL=") - 1) == 0)
      {
        specs[specIdx].msLevel = atoi(line + sizeof("MSLEVEL=") - 1);
        continue;
      }

      if (strncmp("SOURCE_INSTRUMENT=", line, sizeof("SOURCE_INSTRUMENT=") - 1)
          == 0)
      {
        specs[specIdx].instrument_name = line + sizeof("SOURCE_INSTRUMENT=")
            - 1;
        continue;
      }

      if (strncmp("ITOL=", line, sizeof("ITOL=") - 1) == 0)
      {
        specs[specIdx].ITOL = atof(line + sizeof("ITOL=") - 1);
        continue;
      }

      if (strncmp("ITOLU=", line, sizeof("ITOLU=") - 1) == 0)
      {
        specs[specIdx].ITOLU = line + sizeof("ITOLU=") - 1;
        continue;
      }

      if (strncmp("TOL=", line, sizeof("TOL=") - 1) == 0)
      {
        specs[specIdx].TOL = atof(line + sizeof("TOL=") - 1);
        continue;
      }

      if (strncmp("TOLU=", line, sizeof("TOLU=") - 1) == 0)
      {
        specs[specIdx].TOLU = line + sizeof("TOLU=") - 1;
        continue;
      }

      if (strncmp("FILENAME=", line, sizeof("FILENAME=") - 1) == 0)
      {
        specs[specIdx].fileName = line + sizeof("FILENAME=") - 1;
        continue;
      }

      if (strncmp("SPECTRUMQUALITY=", line, sizeof("SPECTRUMQUALITY=") - 1)
          == 0)
      {
        specs[specIdx].spectrum_quality = atoi(line + sizeof("SPECTRUMQUALITY=")
            - 1);
        continue;
      }

      if (strncmp("SUBMISSION_METADATA=",
                  line,
                  sizeof("SUBMISSION_METADATA=") - 1) == 0)
      {
        if (specs[specIdx].psmList.size() > 0)
        {
          std::string submission_metadata = line
              + sizeof("SUBMISSION_METADATA=") - 1;
          specs[specIdx].psmList.front()->m_submission_metadata =
              (submission_metadata);
        }
        continue;
      }

      if (strncmp("ORGANISM=", line, sizeof("ORGANISM=") - 1) == 0)
      {
        if (specs[specIdx].psmList.size() > 0)
        {
          std::string organism = line + sizeof("ORGANISM=") - 1;
          specs[specIdx].psmList.front()->m_organism = (organism);
        }
        continue;
      }

      //ACCEPTING BOTH COMPOUNDNAME and NAME
      if (strncmp("COMPOUNDNAME=", line, sizeof("COMPOUNDNAME=") - 1) == 0)
      {
        if (specs[specIdx].psmList.size() > 0)
        {
          std::string compound = line + sizeof("COMPOUNDNAME=") - 1;
          specs[specIdx].psmList.front()->m_compound_name = (compound);
        }
        continue;
      }

      if (strncmp("NAME=", line, sizeof("NAME=") - 1) == 0)
      {
        if (specs[specIdx].psmList.size() > 0)
        {
          std::string compound = line + sizeof("NAME=") - 1;
          specs[specIdx].psmList.front()->m_compound_name = (compound);
        }
        continue;
      }

      if (strncmp("SMILES=", line, sizeof("SMILES=") - 1) == 0)
      {
        if (specs[specIdx].psmList.size() > 0)
        {
          std::string smiles = line + sizeof("SMILES=") - 1;
          specs[specIdx].psmList.front()->m_smiles = (smiles);
        }
        continue;
      }

      if (strncmp("INCHI=", line, sizeof("INCHI=") - 1) == 0)
      {
        if (specs[specIdx].psmList.size() > 0)
        {
          std::string InChI = line + sizeof("INCHI=") - 1;
          specs[specIdx].psmList.front()->m_InChI = (InChI);
        }
        continue;
      }

      if (strncmp("INCHIAUX=", line, sizeof("INCHIAUX=") - 1) == 0)
      {
        if (specs[specIdx].psmList.size() > 0)
        {
          std::string InChI_aux = line + sizeof("INCHIAUX=") - 1;
          specs[specIdx].psmList.front()->m_InChI_Aux = (InChI_aux);
        }
        continue;
      }

      if (strncmp("NOTES=", line, sizeof("NOTES=") - 1) == 0)
      {
        if (specs[specIdx].psmList.size() > 0)
        {
          std::string notes = line + sizeof("NOTES=") - 1;
          specs[specIdx].psmList.front()->m_notes = notes;
        }
        continue;
      }

      if (strncmp("IONMODE=", line, sizeof("IONMODE=") - 1) == 0)
      {
        if (specs[specIdx].psmList.size() > 0)
        {
          std::string ionmode = line + sizeof("IONMODE=") - 1;
          specs[specIdx].psmList.front()->m_ionmode = ionmode;
        }
        continue;
      }

      if (strncmp("SLGF=", line, sizeof("SLGF=") - 1) == 0)
      {
        if (specs[specIdx].psmList.size() > 0)
        {
          std::string slgf_prob_pair = line + sizeof("SLGF=") - 1;

          float cosine_bin =
              atof(slgf_prob_pair.substr(0, slgf_prob_pair.find(":")).c_str());
          float cosine_prob =
              atof(slgf_prob_pair.substr(slgf_prob_pair.find(":") + 1).c_str());
          pair<float, float> cosine_prob_pair;
          cosine_prob_pair.first = cosine_bin;
          cosine_prob_pair.second = cosine_prob;
          specs[specIdx].psmList.front()->SLGF_distribution.push_back(cosine_prob_pair);
        }
        continue;
      }

      if (strncmp("ACTIVATION=", line, sizeof("ACTIVATION=") - 1) == 0)
      {
        string activ = line + sizeof("ACTIVATION=") - 1;
        specs[specIdx].msFragType = Spectrum::parseActivation(activ);
        continue;
      }

      if (strncmp("EXACTMASS=", line, sizeof("EXACTMASS=") - 1) == 0)
      {
        if (specs[specIdx].psmList.size() > 0)
        {
          std::string exactmass = line + sizeof("EXACTMASS=") - 1;
          specs[specIdx].psmList.front()->m_exactmass = atof(exactmass.c_str());
        }
        continue;
      }

      if (strncmp("LIBRARYQUALITY=", line, sizeof("LIBRARYQUALITY=") - 1) == 0)
      {
        if (specs[specIdx].psmList.size() > 0)
        {
          std::string lib_quality = line + sizeof("LIBRARYQUALITY=") - 1;
          specs[specIdx].psmList.front()->m_library_quality =
              atoi(lib_quality.c_str());
        }
        continue;
      }

      if (strncmp("SPECTRUMID=", line, sizeof("SPECTRUMID=") - 1) == 0)
      {
        if (specs[specIdx].psmList.size() > 0)
        {
          std::string spec_ID = line + sizeof("SPECTRUMID=") - 1;
          specs[specIdx].psmList.front()->m_spectrumID = spec_ID;
        }
        continue;
      }
    }

    /*
     BufferedLineReader blr;
     resize(0);

     // A valid file must have at least 3 lines: BEGIN_IONS, END_IONS and one mass/intensity peak
     if (blr.Load(filename) <= 0 or blr.size() < 3)
     return 0;

     unsigned int lineIdx, specIdx = 0, peakIdx;

     if(scan > 0) {
     index = FindSpectrumInMgfByScan(blr, lineIdx, scan);
     if(index == -1) return 0;
     } else {
     bool ret = FindSpectrumInMgfByIndex(blr, lineIdx, index);
     if(!ret) return 0;
     }


     while (lineIdx < blr.size() and strncmp("END IONS", blr.getline(lineIdx), 8) != 0) {

     line = blr.getline(lineIdx++);

     // SCANS
     if (strncmp("SCANS=", line, sizeof("SCANS=") - 1) == 0) {
     specs[specIdx].scan = (unsigned int)atoi(line + sizeof("SCANS=") - 1);
     continue;
     }

     if (line[0] >= 48 and line[0] <= 57) {
     specs[specIdx].resize(peakIdx+1);
     token = strtok(line, " \t");
     if (!token) {
     cerr << "Error loading " << filename
     << " - could not parse peak mass on line " << lineIdx << "!\n";
     resize(0);
     return 0;
     }
     specs[specIdx][peakIdx][0] = atof(token);
     token = strtok(NULL, " \t");
     if (!token)  {
     cerr << "Error loading " << filename
     << " - could not parse peak intensity on line " << lineIdx
     << "!\n";
     resize(0);
     return 0;
     }
     specs[specIdx][peakIdx][1] = atof(token);
     peakIdx++;
     continue;
     }

     if (strncmp("CHARGE=+", line, 8) == 0) {
     specs[specIdx].parentCharge = (short)atof(&line[8]);
     continue;
     } else if (strncmp("CHARGE=", line, 7) == 0) {
     specs[specIdx].parentCharge = (short)atof(&line[7]);
     continue;
     }

     if (strncmp("TITLE=Scan Number: ", line, 19) == 0) {
     specs[specIdx].scan = (unsigned int)atof(&line[19]);
     continue;
     }

     if (strncmp("PEPMASS=", line, 8) == 0) {
     specs[specIdx].parentMass = (double)strtod(&line[8], NULL);
     specs[specIdx].parentMZ = (double)strtod(&line[8], NULL);
     continue;
     }

     if (strncmp("PRECURSOR=", line, 10) == 0) {
     specs[specIdx].parentMZ = (double)strtod(&line[10], NULL);
     continue;
     }

     if (strncmp("PEPTIDE=", line, sizeof("PEPTIDE=") - 1) == 0) {
     psmPtr psm(new PeptideSpectrumMatch());
     psm->m_annotation = line + sizeof("PEPTIDE=") - 1;
     specs[specIdx].psmList.push_back(psm);
     continue;
     }

     if (strncmp("SEQ=", line, sizeof("SEQ=") - 1) == 0) {
     psmPtr psm(new PeptideSpectrumMatch());
     psm->m_annotation = line + sizeof("SEQ=") - 1;
     specs[specIdx].psmList.push_back(psm);
     specs[specIdx].psmList.front()->m_spectrumFile = filename;
     specs[specIdx].psmList.front()->m_dbIndex = specIdx + 1;
     continue;
     }

     if (strncmp("MSLEVEL=", line, sizeof("MSLEVEL=") - 1) == 0) {
     specs[specIdx].msLevel = atoi(line + sizeof("MSLEVEL=") - 1);
     continue;
     }

     if (strncmp("SOURCE_INSTRUMENT=", line, sizeof("SOURCE_INSTRUMENT=") - 1) == 0) {
     specs[specIdx].instrument_name = line + sizeof("SOURCE_INSTRUMENT=") - 1;
     continue;
     }

     if (strncmp("ITOL=", line, sizeof("ITOL=") - 1) == 0) {
     specs[specIdx].ITOL = atof(line + sizeof("ITOL=") - 1);
     continue;
     }

     if (strncmp("ITOLU=", line, sizeof("ITOLU=") - 1) == 0) {
     specs[specIdx].ITOLU = line + sizeof("ITOLU=") - 1;
     continue;
     }

     if (strncmp("TOL=", line, sizeof("TOL=") - 1) == 0) {
     specs[specIdx].TOL = atof(line + sizeof("TOL=") - 1);
     continue;
     }

     if (strncmp("TOLU=", line, sizeof("TOLU=") - 1) == 0) {
     specs[specIdx].TOLU = line + sizeof("TOLU=") - 1;
     continue;
     }

     if (strncmp("FILENAME=", line, sizeof("FILENAME=") - 1) == 0) {
     specs[specIdx].fileName = line + sizeof("FILENAME=") - 1;
     continue;
     }

     if (strncmp("SPECTRUMQUALITY=", line, sizeof("SPECTRUMQUALITY=") - 1) == 0) {
     specs[specIdx].spectrum_quality = atoi(line + sizeof("SPECTRUMQUALITY=") - 1);
     continue;
     }

     if (strncmp("SUBMISSION_METADATA=", line, sizeof("SUBMISSION_METADATA=") - 1) == 0) {
     if (specs[specIdx].psmList.size() > 0) {
     std::string submission_metadata = line + sizeof("SUBMISSION_METADATA=") - 1;
     specs[specIdx].psmList.front()->m_submission_metadata = (submission_metadata);
     }
     continue;
     }

     if (strncmp("ORGANISM=", line, sizeof("ORGANISM=") - 1) == 0) {
     if (specs[specIdx].psmList.size() > 0) {
     std::string organism = line + sizeof("ORGANISM=") - 1;
     specs[specIdx].psmList.front()->m_organism = (organism);
     }
     continue;
     }

     //ACCEPTING BOTH COMPOUNDNAME and NAME
     if (strncmp("COMPOUNDNAME=", line, sizeof("COMPOUNDNAME=") - 1) == 0) {
     if (specs[specIdx].psmList.size() > 0) {
     std::string compound = line + sizeof("COMPOUNDNAME=") - 1;
     specs[specIdx].psmList.front()->m_compound_name = (compound);
     }
     continue;
     }

     if (strncmp("NAME=", line, sizeof("NAME=") - 1) == 0) {
     if (specs[specIdx].psmList.size() > 0) {
     std::string compound = line + sizeof("NAME=") - 1;
     specs[specIdx].psmList.front()->m_compound_name = (compound);
     }
     continue;
     }

     if (strncmp("SMILES=", line, sizeof("SMILES=") - 1) == 0) {
     if (specs[specIdx].psmList.size() > 0) {
     std::string smiles = line + sizeof("SMILES=") - 1;
     specs[specIdx].psmList.front()->m_smiles = (smiles);
     }
     continue;
     }

     if (strncmp("INCHI=", line, sizeof("INCHI=") - 1) == 0) {
     if (specs[specIdx].psmList.size() > 0) {
     std::string InChI = line + sizeof("INCHI=") - 1;
     specs[specIdx].psmList.front()->m_InChI = (InChI);
     }
     continue;
     }

     if (strncmp("INCHIAUX=", line, sizeof("INCHIAUX=") - 1) == 0) {
     if (specs[specIdx].psmList.size() > 0) {
     std::string InChI_aux = line + sizeof("INCHIAUX=") - 1;
     specs[specIdx].psmList.front()->m_InChI_Aux = (InChI_aux);
     }
     continue;
     }

     if (strncmp("NOTES=", line, sizeof("NOTES=") - 1) == 0) {
     if (specs[specIdx].psmList.size() > 0) {
     std::string notes = line + sizeof("NOTES=") - 1;
     specs[specIdx].psmList.front()->m_notes = notes;
     }
     continue;
     }

     if (strncmp("IONMODE=", line, sizeof("IONMODE=") - 1) == 0) {
     if (specs[specIdx].psmList.size() > 0)  {
     std::string ionmode = line + sizeof("IONMODE=") - 1;
     specs[specIdx].psmList.front()->m_ionmode = ionmode;
     }
     continue;
     }

     if (strncmp("SLGF=", line, sizeof("SLGF=") - 1) == 0) {
     if (specs[specIdx].psmList.size() > 0) {
     std::string slgf_prob_pair = line + sizeof("SLGF=") - 1;

     float cosine_bin  = atof(slgf_prob_pair.substr(0, slgf_prob_pair.find(":")).c_str());
     float cosine_prob = atof(slgf_prob_pair.substr(slgf_prob_pair.find(":") + 1).c_str());
     pair<float, float> cosine_prob_pair;
     cosine_prob_pair.first = cosine_bin;
     cosine_prob_pair.second = cosine_prob;
     specs[specIdx].psmList.front()->SLGF_distribution.push_back(cosine_prob_pair);
     }
     continue;
     }

     if (strncmp("ACTIVATION=", line, sizeof("ACTIVATION=") - 1) == 0) {
     string activ = line + sizeof("ACTIVATION=") - 1;
     specs[specIdx].msFragType = Spectrum::parseActivation(activ);
     continue;
     }

     if (strncmp("EXACTMASS=", line, sizeof("EXACTMASS=") - 1) == 0) {
     if (specs[specIdx].psmList.size() > 0) {
     std::string exactmass = line + sizeof("EXACTMASS=") - 1;
     specs[specIdx].psmList.front()->m_exactmass = atof(exactmass.c_str());
     }
     continue;
     }

     if (strncmp("LIBRARYQUALITY=", line, sizeof("LIBRARYQUALITY=") - 1)   == 0) {
     if (specs[specIdx].psmList.size() > 0) {
     std::string lib_quality = line + sizeof("LIBRARYQUALITY=") - 1;
     specs[specIdx].psmList.front()->m_library_quality =
     atoi(lib_quality.c_str());
     }
     continue;
     }

     if (strncmp("SPECTRUMID=", line, sizeof("SPECTRUMID=") - 1) == 0) {
     if (specs[specIdx].psmList.size() > 0) {
     std::string spec_ID = line + sizeof("SPECTRUMID=") - 1;
     specs[specIdx].psmList.front()->m_spectrumID = spec_ID;
     }
     continue;
     }

     } // while
     */

    // TODO: UPDATED by dbeyter
//    if (specs[specIdx].parentCharge > 0)
//      specs[specIdx].parentMass = (specs[specIdx].parentMass
//          * specs[specIdx].parentCharge)
//          - specs[specIdx].parentCharge;

     if (specs[specIdx].parentCharge > 0)
          specs[specIdx].parentMass = (specs[specIdx].parentMass
              * specs[specIdx].parentCharge)
              - (AAJumps::massHion * (specs[specIdx].parentCharge - 1));

    if (specs[specIdx].scan <= 0)
    {
      specs[specIdx].scan = specIdx + 1;
    }

    //specs[specIdx].sortPeaks();
    //lineIdx++;
    //}

    //if (m_index != 0x0) {
    //  computeIndex();
    //}

    return specs.size();
  }

  // -------------------------------------------------------------------------
  unsigned int SpecSet::LoadSpecSet_mgf(const char *filename)
  {
    BufferedLineReader blr;
    resize(0);
    if (blr.Load(filename) <= 0 or blr.size() < 3)
      return 0; // A valid file must have at least 3 lines: BEGIN_IONS, END_IONS and one mass/intensity peak

    unsigned int lineIdx, specIdx, peakIdx;

    // Counts number of spectra and number of peaks per spectrum
    list<int> peaksPerSpec;
    int first;
    int numPeaks = 0; // Counts number of peaks in the spectrum
    for (lineIdx = 0; lineIdx < blr.size(); lineIdx++)
    {
      if (strncmp("END IONS", blr.getline(lineIdx), 8) == 0)
      {
        peaksPerSpec.push_back(numPeaks);
        numPeaks = -1;
        continue;
      }
      if (strncmp("BEGIN IONS", blr.getline(lineIdx), 10) == 0)
      {
        numPeaks = 0;
        continue;
      }
      first = (int)blr.getline(lineIdx)[0];
      if (numPeaks >= 0 and first >= 48 and first <= 57)
        numPeaks++;
    }

    // Parse spectra
    resize(peaksPerSpec.size());
    lineIdx = 0;
    char *token, *line;
    for (specIdx = 0; specIdx < specs.size(); specIdx++)
    {
      // Skip empty lines
      while (lineIdx < blr.size()
          and (blr.getline(lineIdx)[0] == 0
              or strncmp("BEGIN IONS", blr.getline(lineIdx), 10) != 0))
        lineIdx++;

      if (lineIdx == blr.size())
      {
        cerr << "Error loading " << filename << " - " << specIdx
            << " spectra instead of " << specs.size() << "?\n";
        resize(0);
        return 0;
      }

      // Start of spectrum
      if (strncmp("BEGIN IONS", blr.getline(lineIdx), 10) != 0)
      {
        cerr << "ERROR: Expected BEGIN IONS, found '" << blr.getline(lineIdx)
            << "' (line " << lineIdx + 1 << ")\n";
        return 0;
      }
      else
        lineIdx++;

      // Read peaks/charge/parent mass
      specs[specIdx].resize(peaksPerSpec.front());
      peaksPerSpec.pop_front();
      peakIdx = 0;
      specs[specIdx].parentCharge = 0;
      specs[specIdx].parentMassTol = 0;
      specs[specIdx].scan = 0;
      while (lineIdx < blr.size()
          and strncmp("END IONS", blr.getline(lineIdx), 8) != 0)
      {
        line = blr.getline(lineIdx++);
        if (line[0] >= 48 and line[0] <= 57)
        {
          token = strtok(line, " \t");
          if (!token)
          {
            cerr << "Error loading " << filename
                << " - could not parse peak mass on line " << lineIdx << "!\n";
            resize(0);
            return 0;
          }
          specs[specIdx][peakIdx][0] = atof(token);
          token = strtok(NULL, " \t");
          if (!token)
          {
            cerr << "Error loading " << filename
                << " - could not parse peak intensity on line " << lineIdx
                << "!\n";
            resize(0);
            return 0;
          }
          specs[specIdx][peakIdx][1] = atof(token);
          peakIdx++;
          continue;
        }

        if (strncmp("CHARGE=+", line, 8) == 0)
        {
          specs[specIdx].parentCharge = (short)atof(&line[8]);
        }
        else if (strncmp("CHARGE=", line, 7) == 0)
        {
          specs[specIdx].parentCharge = (short)atof(&line[7]);
        }

        if (strncmp("TITLE=Scan Number: ", line, 19) == 0)
          specs[specIdx].scan = (unsigned int)atof(&line[19]);

        if (strncmp("SCANS=", line, sizeof("SCANS=") - 1) == 0)
          specs[specIdx].scan = (unsigned int)atof(line + sizeof("SCANS=") - 1);
        if (strncmp("PEPMASS=", line, 8) == 0)
        {
          specs[specIdx].parentMass = (double)strtod(&line[8], NULL);
          specs[specIdx].parentMZ = (double)strtod(&line[8], NULL);
        }
        if (strncmp("PRECURSOR=", line, 10) == 0)
        {
          specs[specIdx].parentMZ = (double)strtod(&line[10], NULL);
        }
        if (strncmp("RTINSECONDS=", line, 10) == 0)
        {
          specs[specIdx].retention_time = atof(line + sizeof("RTINSECONDS=")
              - 1);
        }
        if (strncmp("PEPTIDE=", line, sizeof("PEPTIDE=") - 1) == 0)
        {
          psmPtr psm(new PeptideSpectrumMatch());
          psm->m_annotation = line + sizeof("PEPTIDE=") - 1;
          specs[specIdx].psmList.push_back(psm);
        }
        if (strncmp("SEQ=", line, sizeof("SEQ=") - 1) == 0)
        {
          psmPtr psm(new PeptideSpectrumMatch());
          psm->m_annotation = line + sizeof("SEQ=") - 1;
          specs[specIdx].psmList.push_back(psm);
          specs[specIdx].psmList.front()->m_spectrumFile = filename;
          specs[specIdx].psmList.front()->m_dbIndex = specIdx + 1;
        }
        if (strncmp("MSLEVEL=", line, sizeof("MSLEVEL=") - 1) == 0)
        {
          specs[specIdx].msLevel = atoi(line + sizeof("MSLEVEL=") - 1);
        }
        if (strncmp("SOURCE_INSTRUMENT=",
                    line,
                    sizeof("SOURCE_INSTRUMENT=") - 1) == 0)
        {
          specs[specIdx].instrument_name = line + sizeof("SOURCE_INSTRUMENT=")
              - 1;
        }
        if (strncmp("ITOL=", line, sizeof("ITOL=") - 1) == 0)
        {
          specs[specIdx].ITOL = atof(line + sizeof("ITOL=") - 1);
        }
        if (strncmp("ITOLU=", line, sizeof("ITOLU=") - 1) == 0)
        {
          specs[specIdx].ITOLU = line + sizeof("ITOLU=") - 1;
        }
        if (strncmp("TOL=", line, sizeof("TOL=") - 1) == 0)
        {
          specs[specIdx].TOL = atof(line + sizeof("TOL=") - 1);
        }
        if (strncmp("TOLU=", line, sizeof("TOLU=") - 1) == 0)
        {
          specs[specIdx].TOLU = line + sizeof("TOLU=") - 1;
        }
        if (strncmp("FILENAME=", line, sizeof("FILENAME=") - 1) == 0)
        {
          specs[specIdx].fileName = line + sizeof("FILENAME=") - 1;
        }
        if (strncmp("SPECTRUMQUALITY=", line, sizeof("SPECTRUMQUALITY=") - 1)
            == 0)
        {
          specs[specIdx].spectrum_quality = atoi(line
              + sizeof("SPECTRUMQUALITY=") - 1);
        }
        if (strncmp("SUBMISSION_METADATA=",
                    line,
                    sizeof("SUBMISSION_METADATA=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string submission_metadata = line
                + sizeof("SUBMISSION_METADATA=") - 1;
            specs[specIdx].psmList.front()->m_submission_metadata =
                (submission_metadata);
          }
        }
        if (strncmp("ORGANISM=", line, sizeof("ORGANISM=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string organism = line + sizeof("ORGANISM=") - 1;
            specs[specIdx].psmList.front()->m_organism = (organism);
          }
        }
        //ACCEPTING BOTH COMPOUNDNAME and NAME
        if (strncmp("COMPOUNDNAME=", line, sizeof("COMPOUNDNAME=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string compound = line + sizeof("COMPOUNDNAME=") - 1;
            specs[specIdx].psmList.front()->m_compound_name = (compound);
          }
        }
        if (strncmp("NAME=", line, sizeof("NAME=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string compound = line + sizeof("NAME=") - 1;
            specs[specIdx].psmList.front()->m_compound_name = (compound);
          }
        }
        if (strncmp("SMILES=", line, sizeof("SMILES=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string smiles = line + sizeof("SMILES=") - 1;
            specs[specIdx].psmList.front()->m_smiles = (smiles);
          }
        }
        if (strncmp("INCHI=", line, sizeof("INCHI=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string InChI = line + sizeof("INCHI=") - 1;
            specs[specIdx].psmList.front()->m_InChI = (InChI);
          }
        }
        if (strncmp("INCHIAUX=", line, sizeof("INCHIAUX=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string InChI_aux = line + sizeof("INCHIAUX=") - 1;
            specs[specIdx].psmList.front()->m_InChI_Aux = (InChI_aux);
          }
        }
        if (strncmp("NOTES=", line, sizeof("NOTES=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string notes = line + sizeof("NOTES=") - 1;
            specs[specIdx].psmList.front()->m_notes = notes;
          }
        }
        if (strncmp("IONMODE=", line, sizeof("IONMODE=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string ionmode = line + sizeof("IONMODE=") - 1;
            specs[specIdx].psmList.front()->m_ionmode = ionmode;
          }
        }

        if (strncmp("SLGF=", line, sizeof("SLGF=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string slgf_prob_pair = line + sizeof("SLGF=") - 1;

            float cosine_bin =
                atof(slgf_prob_pair.substr(0, slgf_prob_pair.find(":")).c_str());
            float cosine_prob =
                atof(slgf_prob_pair.substr(slgf_prob_pair.find(":") + 1).c_str());
            pair<float, float> cosine_prob_pair;
            cosine_prob_pair.first = cosine_bin;
            cosine_prob_pair.second = cosine_prob;
            specs[specIdx].psmList.front()->SLGF_distribution.push_back(cosine_prob_pair);
          }
        }

        if (strncmp("ACTIVATION=", line, sizeof("ACTIVATION=") - 1) == 0)
        {
          string activ = line + sizeof("ACTIVATION=") - 1;
          specs[specIdx].msFragType = Spectrum::parseActivation(activ);
        }

        if (strncmp("EXACTMASS=", line, sizeof("EXACTMASS=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string exactmass = line + sizeof("EXACTMASS=") - 1;
            specs[specIdx].psmList.front()->m_exactmass =
                atof(exactmass.c_str());
          }
        }

        if (strncmp("LIBRARYQUALITY=", line, sizeof("LIBRARYQUALITY=") - 1)
            == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string lib_quality = line + sizeof("LIBRARYQUALITY=") - 1;
            specs[specIdx].psmList.front()->m_library_quality =
                atoi(lib_quality.c_str());
          }
        }

        if (strncmp("SPECTRUMID=", line, sizeof("SPECTRUMID=") - 1) == 0)
        {
          if (specs[specIdx].psmList.size() > 0)
          {
            std::string spec_ID = line + sizeof("SPECTRUMID=") - 1;
            specs[specIdx].psmList.front()->m_spectrumID = spec_ID;
          }
        }

      }

      // TODO: UPDATED by dbeyter

	if (specs[specIdx].parentCharge > 0)
	  specs[specIdx].parentMass = (specs[specIdx].parentMass
		  * specs[specIdx].parentCharge)
		  - (AAJumps::massHion * specs[specIdx].parentCharge);



//      if (specs[specIdx].parentCharge > 0)
//        specs[specIdx].parentMass = (specs[specIdx].parentMass
//            * specs[specIdx].parentCharge)
//            - (AAJumps::massHion * (specs[specIdx].parentCharge - 1));



      if (specs[specIdx].scan <= 0)
      {
        specs[specIdx].scan = specIdx + 1;
      }
      specs[specIdx].sortPeaks();
      lineIdx++;
    }

    string fullPath = filename;
    FilenameManager mngr(fullPath);
    for (unsigned int i = 0; i < specs.size(); i++)
    {
      if (specs[i].fileName.length() == 0)
      {
        specs[i].fileName = mngr.getFilenameWithExtension();
      }
    }

    if (m_index != 0x0)
    {
      computeIndex();
    }

    return specs.size();
  }
  // -------------------------------------------------------------------------
  unsigned int SpecSet::LoadSpecSet_ms2(const char *filename)
  {
    BufferedLineReader blr;
    resize(0);
    if (blr.Load(filename) <= 0 or blr.size() < 3)
      return 0; // A valid file must have at least 3 lines: 1 header (2 lines) + 1 (m/z,intensity) pair

    unsigned int lineIdx, specIdx, peakIdx;

    // Counts number of spectra and number of peaks per spectrum
    list<int> peaksPerSpec;
    int numPeaks = 1; // Counts number of peaks in the spectrum
    for (lineIdx = 0; lineIdx < blr.size(); lineIdx++)
      if (blr.getline(lineIdx)[0] == ':')
      {
        if (numPeaks > 0)
          peaksPerSpec.push_back(numPeaks);
        numPeaks = -1;
      }
      else if (blr.getline(lineIdx)[0] != 0)
        numPeaks++;
    peaksPerSpec.push_back(numPeaks); // Number of peaks in the last spectrum
    peaksPerSpec.pop_front(); // First element is just the '1' used to initialize numPeaks

    // Parse spectra
    resize(peaksPerSpec.size());
    lineIdx = 0;
    char *token;
    for (specIdx = 0; specIdx < specs.size(); specIdx++)
    {
      // Skip empty lines
      while (blr.getline(lineIdx)[0] == 0 and lineIdx < blr.size())
        lineIdx++;
      if (lineIdx == blr.size())
      {
        cerr << "Error loading " << filename << " - " << specIdx
            << " spectra instead of " << specs.size() << "?\n";
        resize(0);
        return 0;
      }

      // Parse header(s)
      while (blr.getline(lineIdx)[0] == ':')
      {
        lineIdx++;
        token = strtok(blr.getline(lineIdx++), " \t");
        if (!token)
        {
          cerr << "Error loading " << filename
              << " - could not parse parent mass for spectrum " << specIdx + 1
              << "?\n";
          resize(0);
          return 0;
        }
        specs[specIdx].parentMass = (float)atof(token);
        specs[specIdx].parentMZ = (double)strtod(token, NULL);
        token = strtok(NULL, " \t");
        if (!token)
        {
          cerr << "Error loading " << filename
              << " - could not parse parent charge for spectrum " << specIdx + 1
              << "?\n";
          resize(0);
          return 0;
        }
        specs[specIdx].parentCharge = (short)atof(token);
        if (lineIdx == blr.size())
        {
          cerr << "Error loading " << filename << " - " << specIdx
              << " spectra instead of " << specs.size() << "?\n";
          resize(0);
          return 0;
        }
      }

      // Read spectrum peaks
      specs[specIdx].resize(peaksPerSpec.front());
      peaksPerSpec.pop_front();
      for (peakIdx = 0; peakIdx < specs[specIdx].size(); peakIdx++)
      {
        token = strtok(blr.getline(lineIdx++), " \t");
        if (!token)
        {
          cerr << "Error loading " << filename
              << " - could not parse peak mass for spectrum " << specIdx + 1
              << ", peak " << peakIdx + 1 << "?\n";
          resize(0);
          return 0;
        }
        specs[specIdx][peakIdx][0] = (float)atof(token);
        token = strtok(NULL, " \t");
        if (!token)
        {
          cerr << "Error loading " << filename
              << " - could not parse peak intensity for spectrum "
              << specIdx + 1 << ", peak " << peakIdx + 1 << "?\n";
          resize(0);
          return 0;
        }
        specs[specIdx][peakIdx][1] = (float)atof(token);
        if (lineIdx == blr.size() and peakIdx < specs[specIdx].size() - 1)
        {
          cerr << "Error loading " << filename
              << " - end of file before end of spectrum " << specIdx + 1
              << " ended, got only " << peakIdx + 1 << " peaks instead of "
              << specs.size() << "?\n";
          resize(0);
          return 0;
        }
      }
    }

    string fullPath = filename;
    FilenameManager mngr(fullPath);
    for (unsigned int i = 0; i < specs.size(); i++)
    {
      if (specs[i].fileName.length() == 0)
      {
        specs[i].fileName = mngr.getFilenameWithExtension();
      }
    }

    if (m_index != 0x0)
    {
      computeIndex();
    }

    return specs.size();
  }
  // -------------------------------------------------------------------------
  unsigned int SpecSet::LoadSpecSet_prms(const char *filename)
  {
    unsigned int specIdx, lineIdx, peakIdx;
    list<unsigned int> specSizes;
    char *line;
    BufferedLineReader blr;

    // Load whole file into memory
    if (blr.Load(filename) <= 0)
      return 0;

    // Count # spectra in the file
    bool inSpectrum = false;
    peakIdx = 0;
    for (lineIdx = 0; lineIdx < blr.size(); lineIdx++)
    {
      line = blr.getline(lineIdx);
      if (line[0] != 0 and line[0] == '>' and line[1] != 0 and line[1] == '>')
      {
        inSpectrum = true;
        peakIdx = 0;
        continue;
      }
      if (line[0] == 0)
      {
        if (inSpectrum)
          specSizes.push_back(peakIdx);
        inSpectrum = false;
      }
      else if (inSpectrum)
        peakIdx++;
    }
    if (inSpectrum)
      specSizes.push_back(peakIdx); // In case there is no empty line after the last spectrum
    specs.resize(specSizes.size());

    // Parse the text
    peakIdx = 0;
    specIdx = 0;
    list<unsigned int>::iterator sizesIter = specSizes.begin();
    char *token;
    for (lineIdx = 0; lineIdx < blr.size(); lineIdx++)
    {
      line = blr.getline(lineIdx);

      if (line[0] == 0)
      {
        if (inSpectrum)
        { // End of current spectrum
          if (peakIdx > 0 and peakIdx <= specs[specIdx].size())
            specs[specIdx].parentMass = specs[specIdx][peakIdx - 1][0]
                + AAJumps::massMH;
          specIdx++;
          inSpectrum = false;
        }
        continue;
      }

      // Check for start of a new spectrum
      if (line[0] != 0 and line[0] == '>' and line[1] != 0 and line[1] == '>')
      {
        //specs[specIdx].info = (char *)malloc(strlen(line) - 1);
        //strcpy(specs[specIdx].info, &line[2]);
        char *tok = strtok(line, " \t");
        tok = strtok(NULL, " \t"); // Skip ">> " and "<file_index> "
        tok = strtok(NULL, " \t");
        specs[specIdx].scan = (unsigned int)strtoul(tok, NULL, 10);
        specs[specIdx].resize(*sizesIter);
        specs[specIdx].psmList.resize(0);
        sizesIter++;
        peakIdx = 0;
        inSpectrum = true;
        continue;
      }

      // Peak <mass> <intensity> pair
      token = strtok(line, " \t");
      if (token[0] == 0)
      {
        //cerr << "ERROR reading peak mass for peak " << peakIdx
        //    << " for the spectrum entitled (" << specs[specIdx].info
        //    << ") in file " << filename << "!\n";
        specs.resize(0);
        return 0;
      }
      specs[specIdx][peakIdx][0] = atof(token);
      token = strtok(NULL, " \t");
      if (token[0] == 0)
      {
        //cerr << "ERROR reading peak intensity for peak " << peakIdx
        //    << " for the spectrum entitled (" << specs[specIdx].info
        //    << ") in file " << filename << "!\n";
        specs.resize(0);
        return 0;
      }
      specs[specIdx][peakIdx][1] = atof(token);
      peakIdx++;
    }

    string fullPath = filename;
    FilenameManager mngr(fullPath);
    for (unsigned int i = 0; i < specs.size(); i++)
    {
      if (specs[i].fileName.length() == 0)
      {
        specs[i].fileName = mngr.getFilenameWithExtension();
      }
    }

    if (m_index != 0x0)
    {
      computeIndex();
    }

    return specs.size();
  }
  // -------------------------------------------------------------------------
  unsigned int SpecSet::LoadSpecSet_prmsv3(const char *filename,
                                           vector<vector<string> >* prmOrigins)
  {
    unsigned int specIdx, lineIdx, peakIdx;
    list<unsigned int> specSizes;
    char *line;
    BufferedLineReader blr;

    // Load whole file into memory
    if (blr.Load(filename) <= 0)
      return 0;

    // Count # spectra in the file
    bool inSpectrum = false;
    peakIdx = 0;
    for (lineIdx = 0; lineIdx < blr.size(); lineIdx++)
    {
      line = blr.getline(lineIdx);
      if (line[0] != 0 and line[0] == '>' and line[1] != 0 and line[1] == '>')
      {
        inSpectrum = true;
        peakIdx = 0;
        continue;
      }

      if (line[0] == 0)
      {
        if (inSpectrum)
        {
          specSizes.push_back(peakIdx);
        }
        inSpectrum = false;
      }
      else if (line[0] == 'E' || line[0] == '#')
      {
        if (inSpectrum)
        {
          specSizes.push_back(peakIdx);
          inSpectrum = false;
        }
      }
      else if (inSpectrum and line[0] != 'C')
        peakIdx++;
    }
    if (inSpectrum)
      specSizes.push_back(peakIdx); // In case there is no empty line after the last spectrum
    specs.resize(specSizes.size());
    if (prmOrigins)
    {
      prmOrigins->resize(specs.size());
    }

    // Parse the text
    peakIdx = 0;
    specIdx = 0;
    list<unsigned int>::iterator sizesIter = specSizes.begin();
    char *token;
    for (lineIdx = 0; lineIdx < blr.size(); lineIdx++)
    { // Skip header lines
      line = blr.getline(lineIdx);
      if (line[0] != 0 and line[0] == '>' and line[1] != 0 and line[1] == '>')
        break;
    }

    for (; lineIdx < blr.size(); lineIdx++)
    {
      line = blr.getline(lineIdx);
      if (line[0] == 0 || line[0] == 'E' || line[0] == '#')
      {
        if (inSpectrum)
        { // End of current spectrum
          //  -- Not valid for Pepnovo v3
          //        if(peakIdx>0 and peakIdx<=specs[specIdx].size()) specs[specIdx].parentMass = specs[specIdx][peakIdx-1][0]+AAJumps::massMH;
          specIdx++;
          inSpectrum = false;
        }
        // DEBUG_VAR(lineIdx);
        continue;
      }

      // Check for start of a new spectrum
      if (line[0] != 0 and line[0] == '>' and line[1] != 0 and line[1] == '>')
      {
        if (inSpectrum)
        {
          specIdx++;
          inSpectrum = false;
        }
        //specs[specIdx].info = (char *)malloc(strlen(line) - 1);
        //strcpy(specs[specIdx].info, &line[2]);
        char *tok = strtok(line, " \t");
        tok = strtok(NULL, " \t"); // Skip ">> " and "<file_index> "
        tok = strtok(NULL, " \t");
        tok = strtok(NULL, " \t");
        if (strncmp("Scan", tok, 4) == 0)
        {
          tok = strtok(NULL, " \t");
          tok = strtok(NULL, " \t");
          specs[specIdx].scan = (unsigned int)atof(tok);
        }

        //specs[specIdx].scan = (unsigned int)strtoul(tok, NULL, 10);

        // Read charge/mass header
        line = blr.getline(++lineIdx);
        if (line[0] != 0 and (line[0] == '#' or line[0] == 'E'))
        {
          WARN_MSG("PepNovo failed to process spectrum " << specIdx << " (scan=" << specs[specIdx].scan << ") from file " << filename << ": \'" << line << "\'");
          inSpectrum = false;
          specIdx++;
          sizesIter++;
          continue;
        }
        specs[specIdx].parentCharge = (short)atoi(&line[8]);
        specs[specIdx].parentMass = (float)atof(&line[19]);
        peakIdx = 0;
        inSpectrum = true;

        specs[specIdx].resize(*sizesIter);
        if (prmOrigins)
        {
          (*prmOrigins)[specIdx].resize(specs[specIdx].size());
        }
        specs[specIdx].psmList.resize(0);
        sizesIter++;
        continue;
      }

      if (!inSpectrum)
      {
        continue;
      }

      // Peak <mass> <intensity> pair
      token = strtok(line, " \t");
      if (token[0] == 0)
      {
        //cerr << "ERROR reading peak mass for peak " << peakIdx
        //    << " for the spectrum entitled (" << specs[specIdx].info
        //    << ") in file " << filename << "!\n";
        specs.resize(0);
        return 0;
      }

      specs[specIdx][peakIdx][0] = atof(token);
      token = strtok(NULL, " \t");
      if (token[0] == 0)
      {
        //cerr << "ERROR reading peak intensity for peak " << peakIdx
        //    << " for the spectrum entitled (" << specs[specIdx].info
        //    << ") in file " << filename << "!\n";
        specs.resize(0);
        return 0;
      }

      specs[specIdx][peakIdx][1] = atof(token);

      if (prmOrigins)
      {
        token = strtok(NULL, " \t");
        if (token[0] == 0)
        {
          WARN_MSG("could not find MS/MS origin peak for PRM at index " << peakIdx << " in spectrum " << specIdx);
        }
        else
        {
          (*prmOrigins)[specIdx][peakIdx] = token;
        }
      }
      peakIdx++;
    }

    string fullPath = filename;
    FilenameManager mngr(fullPath);
    for (unsigned int i = 0; i < specs.size(); i++)
    {
      if (specs[i].fileName.length() == 0)
      {
        specs[i].fileName = mngr.getFilenameWithExtension();
      }
    }

    if (m_index != 0x0)
    {
      computeIndex();
    }

    return specs.size();
  }
  // -------------------------------------------------------------------------
  unsigned int SpecSet::LoadSpecSet_pklbin_with_annotation(const char * spectra_filename,
                                                           const char * annotation_filename)
  {
    SpecSet temp_specs;
    temp_specs.loadPklBin(spectra_filename);

    PeptideSpectrumMatchSetSpectralLibraryLoader psm_set;
    psm_set.loadSpecnetsResultsFile(annotation_filename);

    int original_specs_size = specs.size();
    //cout << "Original Size: " << original_specs_size << endl;
    //cout << temp_specs.size() << endl;
    specs.insert(specs.end(), temp_specs.specs.begin(), temp_specs.specs.end());

    for (int psm_idx = 0; psm_idx < psm_set.size(); psm_idx++)
    {
      psm_set[psm_idx]->m_spectrum = &(specs[psm_set[psm_idx]->m_scanNum - 1
          + original_specs_size]);
      psm_set[psm_idx]->m_spectrum->psmList.push_back(psm_set[psm_idx]);
    }

    return 0;
  }
  // -------------------------------------------------------------------------
  unsigned int SpecSet::LoadSpecSet_mgf_with_annotation(const char * spectra_filename,
                                                        const char * annotation_filename)
  {
    SpecSet temp_specs;
    temp_specs.LoadSpecSet_mgf(spectra_filename);

    PeptideSpectrumMatchSetSpectralLibraryLoader psm_set;
    psm_set.loadSpecnetsResultsFile(annotation_filename);

    int original_specs_size = specs.size();
    //cout << "Original Size: " << original_specs_size << endl;
    //cout << temp_specs.size() << endl;
    specs.insert(specs.end(), temp_specs.specs.begin(), temp_specs.specs.end());

    for (int psm_idx = 0; psm_idx < psm_set.size(); psm_idx++)
    {
      psm_set[psm_idx]->m_spectrum = &(specs[psm_set[psm_idx]->m_scanNum - 1
          + original_specs_size]);
      psm_set[psm_idx]->m_spectrum->psmList.push_back(psm_set[psm_idx]);
    }

    return 0;
  }
  // -------------------------------------------------------------------------
  unsigned int SpecSet::LoadSpecSet_pkl(const char *filename)
  {
    int numSpecs = 1; // Assume file contains at least one spectrum and count one additional spectrum per blank line
    ifstream input(filename, ios::binary);
    if (!input.is_open() || !input.good())
    {
      cerr << "ERROR: cannot open " << filename << "\n";
      specs.resize(0);
      return 0;
    }

    list<int> specSizes; // Register number of peaks per spectrum

    char *lineBuffer = (char *)malloc(1025);
    //    char *lineBuffer = new char [1024];
    int numPeaks = -1; // zero only after reading first tuple of (parent mass, intensity, charge)
    while (!input.eof() && !input.fail())
    {
      input.getline(lineBuffer, 1024, '\n');
      if (lineBuffer[0] == 0 or lineBuffer[0] == '\n' or lineBuffer[0] == '\r')
      {
        if (numPeaks >= 0)
        {
          numSpecs++;
          specSizes.push_back(numPeaks);
          numPeaks = -1;
        }
      }
      else
      {
        if (((int)lineBuffer[0] >= (int)'0' and (int)lineBuffer[0] <= (int)'9')
            or lineBuffer[0] == '-')
          numPeaks++;
      }
    }
    if (numPeaks >= 0)
    {
      specSizes.push_back(numPeaks);
    }
    else
    {
      numSpecs--; // In case there was no empty line at the end of the pkl file
    }

    specs.resize(numSpecs);
    //    input.seekg(0,ios_base::beg);
    //    input.close();  input.open(filename);
    //    Neither of the above worked, so ...
    input.close();
    ifstream input2(filename, ios::binary);

    float foo;
    for (int i = 0; i < numSpecs; i++)
    {
      input2 >> specs[i].parentMass >> foo >> specs[i].parentCharge; // Read precursor mass, total intensity, charge state
      specs[i].parentMZ = specs[i].parentMass;
      if (specs[i].parentCharge > 0)
        specs[i].parentMass = specs[i].parentMass * specs[i].parentCharge
            - specs[i].parentCharge + AAJumps::massHion; // Convert to PM+19

      //        specs[i].numPeaks = specSizes.front();   specSizes.pop_front();
      //        specs[i].peakList.resize(specs[i].numPeaks);
      specs[i].resize(specSizes.front());
      specs[i].psmList.resize(0);
      specSizes.pop_front();

      for (int j = 0; j < specs[i].size(); j++)
        input2 >> specs[i][j][0] >> specs[i][j][1];
      //            input2 >> specs[i].peakList[j][0] >> specs[i].peakList[j][1];
      specs[i].scan = i + 1;
      input2.getline(lineBuffer, 1024, '\n'); // Read intermediate newline
    }

    //    delete lineBuffer;
    input2.close();
    free(lineBuffer);

    string fullPath = filename;
    FilenameManager mngr(fullPath);
    for (unsigned int i = 0; i < specs.size(); i++)
    {
      if (specs[i].fileName.length() == 0)
      {
        specs[i].fileName = mngr.getFilenameWithExtension();
      }
    }

    if (m_index != 0x0)
    {
      computeIndex();
    }

    return 1;
  }

  // -------------------------------------------------------------------------
  short SpecSet::SaveSpecSet_pkl(const char *filename)
  {
    ofstream output(filename, ios::binary);
    if (!output)
    {
      cerr << "ERROR: cannot open " << filename << "\n";
      return -1;
    }

    for (int i = 0; i < specs.size(); i++)
    {
      output << specs[i].parentMass << " -1 " << specs[i].parentCharge << "\n";
      for (int j = 0; j < specs[i].size(); j++)
      {
        output << specs[i][j][0] << " " << specs[i][j][1] << endl;
        //            output << specs[i].peakList[j][0] << " " << specs[i].peakList[j][1] << endl;
      }
      if (i < specs.size() - 1)
      {
        output << endl;
      }
    }

    output.close();
    return 1;
  }
  // -------------------------------------------------------------------------
  bool SpecSet::SaveAnnotations(const char* filename)
  {
    FILE* out_buf = fopen(filename, "wb");

    if (out_buf == NULL)
      return false;

    for (int i = 0; i < size(); i++)
    {
      if (specs[i].psmList.size() == 0)
      {
        fprintf(out_buf, "\n");
      }
      else
      {
        list<psmPtr>::iterator it;
        for (it = specs[i].psmList.begin(); it != specs[i].psmList.end(); it++)
        {
          fprintf(out_buf,
                  "%s,",
                  (*it)->m_annotation.substr(2,
                                             (*it)->m_annotation.length() - 4).c_str());
        }
        fprintf(out_buf, "\n");
      }
    }
    fclose(out_buf);
    return true;
  }
  // -------------------------------------------------------------------------
  short SpecSet::SaveSpecSet(const char* filename, bool outputScans)
  {
    // split the filename into components
    FilenameManager fm(filename);

    fm.lowerCaseExtension();

    // select load method according to extension
    if (fm.extension.compare("mgf") == 0)
      return SaveSpecSet_mgf(filename, outputScans);

    if (fm.extension.compare("pkl") == 0)
      return SaveSpecSet_pkl(filename);

    if (fm.extension.compare("pklbin") == 0)
      return savePklBin(filename);

    if (fm.extension.compare("ms2") == 0)
      return SaveSpecSet_ms2(filename);

    if (fm.extension.compare("dta") == 0)
      return SaveSpecSet_dta(filename);

    return -2;
  }
  // -------------------------------------------------------------------------
  short SpecSet::SaveSpecSet_ms2(const char* filename)
  {
    for (int i = 0; i < specs.size(); i++)
    {
      ostringstream outs;
      outs << filename << i << ".ms2";
      ofstream output(outs.str().c_str(), ios::binary);
      if (!output)
      {
        cerr << "ERROR: cannot open " << filename << "\n";
        return -1;
      }
      specs[i].output_ms2(output);
      output.close();
    }
    return 1;
  }
  // -------------------------------------------------------------------------
  short SpecSet::SaveSpecSet_mgf(const char* filename, bool outputScans)
  {
    FILE* output = fopen(filename, "wb");
    if (!output)
    {
      cerr << "ERROR: cannot open " << filename << "\n";
      return -1;
    }

    float mz, z, pepmass;
    for (int i = 0; i < specs.size(); i++)
    {
      z = (float)specs[i].parentCharge;
      if (specs[i].psmList.size() > 0)
      {
        if ((*specs[i].psmList.begin())->m_charge != 0)
        {
          z = (*specs[i].psmList.begin())->m_charge;
        }
      }
      mz = specs[i].parentMZ;
      pepmass = specs[i].parentMass;
      if (z > 0 && mz == 0)
        mz = (pepmass + (z - 1) * AAJumps::massHion) / z;
      fprintf(output, "BEGIN IONS\nPEPMASS=%.5f\nCHARGE=%.0f+\n", mz, z);
      //output << "BEGIN IONS\nPEPMASS="<<mz<<"\nCHARGE=+"<<z<<"\n";

      if (specs[i].msLevel > 0)
        fprintf(output, "MSLEVEL=%d\n", specs[i].msLevel);
      if (specs[i].instrument_name.length() > 0)
        fprintf(output,
                "SOURCE_INSTRUMENT=%s\n",
                specs[i].instrument_name.c_str());
      if (specs[i].ITOL > 0)
        fprintf(output, "ITOL=%f\n", specs[i].ITOL);
      if (specs[i].ITOLU.length() > 0)
        fprintf(output, "ITOLU=%s\n", specs[i].ITOLU.c_str());
      if (specs[i].TOL > 0)
        fprintf(output, "TOL=%f\n", specs[i].TOL);
      if (specs[i].precursor_kl > 0)
        fprintf(output, "PRECURSORKL=%f\n", specs[i].precursor_kl);
      if (specs[i].precursor_intensity > 0)
        fprintf(output, "PRECURSORINTENSITY=%f\n", specs[i].precursor_intensity);
      if (specs[i].TOLU.length() > 0)
        fprintf(output, "TOLU=%s\n", specs[i].TOLU.c_str());
      if (specs[i].spectrum_quality > 0)
        fprintf(output, "SPECTRUMQUALITY=%d\n", specs[i].spectrum_quality);
      if (specs[i].fileName.length() > 0)
        fprintf(output, "FILENAME=%s\n", specs[i].fileName.c_str());

      if (specs[i].psmList.size() > 0)
      {
        list<psmPtr>::iterator it;
        for (it = specs[i].psmList.begin(); it != specs[i].psmList.end(); it++)
        {
          fprintf(output, "SEQ=%s\n", (*it)->m_annotation.c_str());
          fprintf(output, "NOTES=%s\n", (*it)->m_notes.c_str());
          fprintf(output, "IONMODE=%s\n", (*it)->m_ionmode.c_str());

          if ((*it)->m_exactmass > 0)
          {
            fprintf(output, "EXACTMASS=%f\n", (*it)->m_exactmass);
          }

          for (int SLGF_Distirbution_idx = 0;
              SLGF_Distirbution_idx < (*it)->SLGF_distribution.size();
              SLGF_Distirbution_idx++)
          {
            fprintf(output,
                    "SLGF=%f:%f\n",
                    (*it)->SLGF_distribution[SLGF_Distirbution_idx].first,
                    (*it)->SLGF_distribution[SLGF_Distirbution_idx].second);
          }

          if ((*it)->m_submission_metadata.length() > 0)
            fprintf(output,
                    "SUBMISSION_METADATA=%s\n",
                    (*it)->m_submission_metadata.c_str());
          if ((*it)->m_organism.length() > 0)
            fprintf(output, "ORGANISM=%s\n", (*it)->m_organism.c_str());
          if ((*it)->m_compound_name.length() > 0)
            fprintf(output, "NAME=%s\n", (*it)->m_compound_name.c_str());
          if ((*it)->m_smiles.length() > 0)
            fprintf(output, "SMILES=%s\n", (*it)->m_smiles.c_str());
          if ((*it)->m_InChI.length() > 0)
            fprintf(output, "INCHI=%s\n", (*it)->m_InChI.c_str());
          if ((*it)->m_InChI_Aux.length() > 0)
            fprintf(output, "INCHIAUX=%s\n", (*it)->m_InChI_Aux.c_str());

          if ((*it)->m_library_quality > 0)
            fprintf(output, "LIBRARYQUALITY=%d\n", (*it)->m_library_quality);

          if ((*it)->m_spectrumID.length() > 0)
            fprintf(output, "SPECTRUMID=%s\n", (*it)->m_spectrumID.c_str());

        }
      }

      string activID = "ACTIVATION=";
      activID += Spectrum::activationToString(specs[i].msFragType);
      activID += "\n";
      fprintf(output, "%s", activID.c_str());

      string analType = "INSTRUMENT=";
      analType += Spectrum::massAnalyzerToString(specs[i].msMassAnalyzerType);
      analType += "\n";
      fprintf(output, "%s", analType.c_str());

      if (specs[i].scan >= 0)
      {
        fprintf(output, "TITLE=Scan Number: %d\n", specs[i].scan);

        if (outputScans)
          fprintf(output, "SCANS=%d\n", specs[i].scan);

        //output<<"TITLE=Scan Number: "<<specs[i].scan<<"\n";
      }

      for (int j = 0; j < specs[i].size(); j++)
      {
        fprintf(output, "%.6f %.6f\n", specs[i][j][0], specs[i][j][1]);
        //output << specs[i][j][0] << " " << specs[i][j][1] << endl;
        //            output << specs[i].peakList[j][0] << " " << specs[i].peakList[j][1] << endl;
      }
      fprintf(output, "END IONS\n");
      //output << "END IONS\n";
      if (i < specs.size() - 1)
      {
        fprintf(output, "\n");
        //output << endl;
      }
    }
    fclose(output);
    //output.close();
    return 1;
  }

  bool SpecSet::SaveSpecSet_dta(const char *filename)
  {

    int maxScanLength = 0, maxIdxLength = 0;
    ;
    for (int i = 0; i < size(); i++)
    {
      string scanStr = parseInt(specs[i].scan);
      maxScanLength = max(maxScanLength, (int)scanStr.length());
      string idxStr = parseInt(i);
      maxIdxLength = max(maxIdxLength, (int)idxStr.length());
    }

    float mz, z, pepmass;
    for (int i = 0; i < size(); i++)
    {
      FilenameManager fm(filename);
      fm.filename += ".";
      fm.filename += parseInt(i, maxIdxLength);
      fm.filename += ".";
      fm.filename += parseInt(specs[i].scan, maxScanLength);
      fm.joinFilename();

      FILE* output = fopen(fm.filenameFull.c_str(), "wb");
      if (!output)
      {
        ERROR_MSG("ERROR: cannot open " << fm.filenameFull);
        return false;
      }

      mz = specs[i].parentMZ;
      pepmass = specs[i].parentMass;
      z = (float)specs[i].parentCharge;
      if (z > 0 && mz == 0)
        mz = (pepmass + (z - 1) * AAJumps::massHion) / z;
      fprintf(output, "%.5f %.0f\n", mz, z);

      for (int j = 0; j < specs[i].size(); j++)
      {
        fprintf(output, "%.6f %.6f\n", specs[i][j][0], specs[i][j][1]);
      }
      fclose(output);
    }
    return true;
  }
  // -------------------------------------------------------------------------
  short SpecSet::SaveSpecSet_bin(const char *filename)
  {
    WARN_MSG("SaveSpecSet_bin has been DEPRECATED!");
    WARN_MSG("Data previously saved in bin file is now save in pklbin");
    WARN_MSG("Use savePklBin() and loadPklBin().");
    // Matrix to store data
    vector<vector<int> > data;
    // get # of spectra
    unsigned int numSpecs = specs.size();
    // cycle thru all spectra
    for (int i = 0; i < numSpecs; i++)
    {
      // auxiliary vector to store a pair (scan #, msLevel)
      vector<int> aux;
      // store scan #
      aux.push_back(specs[i].scan);
      // store msLevel
      aux.push_back(specs[i].msLevel);
      // add spectrum data to matrix
      data.push_back(aux);
    }
    // call procedure to save .bin file and return it's return value
    return Save_binArray(filename, data);
  }
  //-------------------------------------------------------------------------
  //  save - saves a SpecSet to a file in binary format.
  int SpecSet::savePklBin(const char * spectrumFilename,
                          const char * psmFilename /* = 0x0 */,
                          const char * peaksFilename /* = 0x0 */)
  {
    FILE *fp;
    unsigned int numSpecs = specs.size();
    unsigned int i, p;

    fp = fopen(spectrumFilename, "wb");
    if (fp == 0)
    {
      cerr << "ERROR: cannot open " << spectrumFilename << "\n";
      return -1;
    }

    unsigned int dummy = 0;
    unsigned int count = fwrite(&dummy, sizeof(int), 1, fp); // Empty 4 bytes

    // Latest version
    char version = 3;
    count = fwrite(&version, sizeof(char), 1, fp); // Version of file
    char subversion = 1;
    count = fwrite(&subversion, sizeof(char), 1, fp); // Sub-version of file

    count = fwrite(&numSpecs, sizeof(int), 1, fp); // Number of spectra in the file

    map<string, unsigned short> versions;
    versions[Spectrum::BIN_VERSION_ID] = Spectrum::BIN_VERSION;
    versions[Spectrum::BIN_SUBVERSION_ID] = Spectrum::BIN_SUBVERSION;

    if (!writeStringMapToBinaryStream<unsigned short>(fp, versions))
    {
      ERROR_MSG("Error saving version info");
      fclose(fp);
      return false;
    }

    for (unsigned int i = 0; i < numSpecs; i++)
    {
      if (!specs[i].saveToBinaryStream(fp))
      {
        ERROR_MSG("Failed to save spectrum " << i);
        fclose(fp);
        return -1;
      }
    }

    fclose(fp);

    if (psmFilename != 0x0)
    {
      PeptideSpectrumMatchSet psmSetTemp;
      psmSetTemp.getPSMSet(this);
      psmSetTemp.saveToFile(psmFilename);
    }

    if (peaksFilename)
    {
      SpecSet tempMatchedPeaks(specs.size());
      for (int i = 0; i < specs.size(); i++)
      {
        list<psmPtr>::iterator litr = specs[i].psmList.begin();
        if (litr == specs[i].psmList.end())
        {
          continue;
        }
        int peakListSize = (*litr)->m_matchedPeaks.size();
        tempMatchedPeaks[i].resize(peakListSize);
        for (int j = 0; j < peakListSize; j++)
        {
          tempMatchedPeaks[i][j].set((*litr)->m_matchedPeaks[j][0],
                                     (*litr)->m_matchedPeaks[j][1]);
        }
      }
      tempMatchedPeaks.savePklBin(peaksFilename);
    }

    return 1;
  }
  bool SpecSet::loadPklBin_1(FILE* fp,
                             int numSpecs,
                             bool oldVersion,
                             char subversion)
  {
    float *data; // Pointer to array containing spectrum data
    unsigned int i, p, dataIdx, count, numValues;
    unsigned short *numPeaks; // Pointer to array containing 1+number of peaks per spectrum
    unsigned int dataSize = 0; // Size of the data buffer

    specs.resize(numSpecs);

    // Read in the scan numbers and MS Levels if this is a second generation file
    // First generation files didn't have the dummy field nor this data
    if (!oldVersion)
    {

      // Read the scan numbers
      unsigned int *scanNums = (unsigned int *)malloc(sizeof(unsigned int)
          * numSpecs);
      if (scanNums == (unsigned int *)0)
      {
        ERROR_MSG("Not enough memory for " << numSpecs << " scan numbers");
        fclose(fp);
        return false;
      }

      count = fread(scanNums, sizeof(unsigned int), numSpecs, fp);
      if (count != numSpecs)
      {
        free(scanNums);
        return false;
      }
      for (unsigned int i = 0; i < numSpecs; i++)
      {
        specs[i].scan = scanNums[i];
        if (specs[i].scan <= 0)
        {
          specs[i].scan = i + 1;
        }
      }
      free(scanNums);
      // Read the MS levels
      short *msLevels = (short *)malloc(sizeof(short) * numSpecs);
      if (msLevels == (short *)0)
      {
        ERROR_MSG("Not enough memory for " << numSpecs << " MS Levels.");
        fclose(fp);
        return false;
      }
      count = fread(msLevels, sizeof(short), numSpecs, fp);
      if (count != numSpecs)
      {
        free(msLevels);
        fclose(fp);
        return false;
      }
      for (unsigned int i = 0; i < numSpecs; i++)
      {
        specs[i].msLevel = msLevels[i];
      }
      free(msLevels);
    }
    else
    {
      for (unsigned int i = 0; i < numSpecs; i++)
      {
        specs[i].scan = i + 1;
      }
    }

    // Number of peaks per spectrum
    numPeaks = (unsigned short *)malloc(sizeof(unsigned short) * numSpecs);
    if (numPeaks == (unsigned short *)0)
    {
      ERROR_MSG("Not enough memory for " << numSpecs << " spectra");
      fclose(fp);
      return false;
    }
    count = fread(numPeaks, sizeof(unsigned short), numSpecs, fp);
    if (count != numSpecs)
    {
      free(numPeaks);
      fclose(fp);
      return false;
    }
    for (i = 0; i < numSpecs; i++)
    {
      if (numPeaks[i] > dataSize)
        dataSize = numPeaks[i];

      //cout << "numPeaks[" << i << "] = " << numPeaks[i] << endl;
    }
    dataSize = 2 * dataSize + 2; // Each spectrum 'line' has 2 values and there's one additional 'line' with parent mass/charge
    data = (float *)malloc(sizeof(float) * dataSize);
    if (data == (float *)0)
    {
      ERROR_MSG("Not enough memory for " << dataSize << " floats");
      fclose(fp);
      return false;
    }

    for (i = 0; i < numSpecs; i++)
    {
      unsigned int numValues = 2 * numPeaks[i] + 2;
      count = fread(data, sizeof(float), numValues, fp);
      if (count != numValues)
      {
        ERROR_MSG("Not enough memory for " << numSpecs << " numValues");
        free(numPeaks);
        free(data);
        fclose(fp);
        return false;
      }

      specs[i].parentMass = data[0];
      specs[i].parentCharge = (int)data[1];

      if (specs[i].parentCharge > 0)
      {
        specs[i].parentMZ = (specs[i].parentMass
            + ((specs[i].parentCharge - 1.0) * AAJumps::massHion))
            / specs[i].parentCharge;
      }
      else
      {
        specs[i].parentMZ = specs[i].parentMass;
      }

      //specs[i].scan = i + 1;
      specs[i].resize(numPeaks[i]);
      specs[i].psmList.resize(0);
      for (unsigned int p = 0, dataIdx = 2; dataIdx < numValues;
          dataIdx += 2, p++)
        specs[i][p].set(data[dataIdx], data[dataIdx + 1]);
    }

    free(numPeaks);
    free(data);
    fclose(fp);

    return true;
  }

  bool SpecSet::loadPklBin_2(FILE* fp, int numSpecs, char subversion)
  {
    float *data; // Pointer to array containing spectrum data
    unsigned int i, p, dataIdx, count, numValues;
    unsigned int dataSize = 0; // Size of the data buffer
    float *PMTols, *PMs;
    short *msLevels, *fragTypes, *precCharges;

    specs.resize(numSpecs);

    // Read in the scan numbers and MS Levels if this is a second generation file
    // First generation files didn't have the dummy field nor this data

    // Read the scan numbers
    unsigned int *scanNums = (unsigned int *)malloc(sizeof(unsigned int)
        * numSpecs);
    if (scanNums == (unsigned int *)0)
    {
      ERROR_MSG("Not enough memory for " << numSpecs << " scan numbers.");
      goto load_fail;
    }

    count = fread(scanNums, sizeof(unsigned int), numSpecs, fp);
    if (count != numSpecs)
    {
      goto load_fail;
    }
    for (unsigned int i = 0; i < numSpecs; i++)
    {
      specs[i].scan = scanNums[i];
      if (specs[i].scan <= 0)
      {
        specs[i].scan = i + 1;
      }
    }
    free(scanNums);

    // Read the MS levels
    msLevels = (short *)malloc(sizeof(short) * numSpecs);
    if (msLevels == (short *)0)
    {
      ERROR_MSG("Not enough memory for " << numSpecs << " MS Levels.");
      goto load_fail;
    }
    count = fread(msLevels, sizeof(short), numSpecs, fp);
    if (count != numSpecs)
    {
      free(msLevels);
      goto load_fail;
    }
    for (unsigned int i = 0; i < numSpecs; i++)
    {
      specs[i].msLevel = msLevels[i];
    }
    free(msLevels);

    // Read the fragmentation types
    fragTypes = (short *)malloc(sizeof(short) * numSpecs);
    if (fragTypes == (short *)0)
    {
      ERROR_MSG("Not enough memory for " << numSpecs << " MS frag types.");
      goto load_fail;
    }
    count = fread(fragTypes, sizeof(short), numSpecs, fp);
    if (count != numSpecs)
    {
      free(fragTypes);
      goto load_fail;
    }
    for (unsigned int i = 0; i < numSpecs; i++)
    {
      if (fragTypes[i] == (short)Spectrum::FragType_CID)
      {
        specs[i].msFragType = Spectrum::FragType_CID;
      }
      else if (fragTypes[i] == (short)Spectrum::FragType_HCD)
      {
        specs[i].msFragType = Spectrum::FragType_HCD;
      }
      else if (fragTypes[i] == (short)Spectrum::FragType_ETD)
      {
        specs[i].msFragType = Spectrum::FragType_ETD;
      }
      else
      {
        ERROR_MSG("Found unsupported fragmentation ID " << fragTypes[i]);
        free(fragTypes);
        goto load_fail;
      }
    }
    free(fragTypes);

    // Read the parent masses
    PMs = (float *)malloc(sizeof(float) * numSpecs);
    if (PMs == (float *)0)
    {
      ERROR_MSG("Not enough memory for " << numSpecs << " parent masses.");
      goto load_fail;
    }
    count = fread(PMs, sizeof(float), numSpecs, fp);
    if (count != numSpecs)
    {
      free(PMs);
      goto load_fail;
    }
    for (unsigned int i = 0; i < numSpecs; i++)
    {
      specs[i].parentMass = PMs[i];
    }
    free(PMs);

    // Read the MS precursor charges
    precCharges = (short *)malloc(sizeof(short) * numSpecs);
    if (precCharges == (short *)0)
    {
      ERROR_MSG("Not enough memory for " << numSpecs << " MS Charges.");
      goto load_fail;
    }
    count = fread(precCharges, sizeof(short), numSpecs, fp);
    if (count != numSpecs)
    {
      free(precCharges);
      goto load_fail;
    }
    for (unsigned int i = 0; i < numSpecs; i++)
    {
      specs[i].parentCharge = precCharges[i];
      if (specs[i].parentCharge > 0)
      {
        specs[i].parentMZ = (specs[i].parentMass
            + ((specs[i].parentCharge - 1.0) * AAJumps::massHion))
            / specs[i].parentCharge;
      }
      else
      {
        specs[i].parentMZ = specs[i].parentMass;
      }
    }
    free(precCharges);

    // Read the parent masses tolerances
    PMTols = (float *)malloc(sizeof(float) * numSpecs);
    if (PMTols == (float *)0)
    {
      ERROR_MSG("Not enough memory for " << numSpecs << " parent mass tolerances.");
      goto load_fail;
    }
    count = fread(PMTols, sizeof(float), numSpecs, fp);
    if (count != numSpecs)
    {
      free(PMTols);
      goto load_fail;
    }
    for (unsigned int i = 0; i < numSpecs; i++)
    {
      specs[i].parentMassTol = PMTols[i];
    }
    free(PMTols);

    unsigned int* numPeaks; // Pointer to array containing 1+number of peaks per spectrum
    numPeaks = (unsigned int *)malloc(sizeof(unsigned int) * numSpecs);
    if (numPeaks == (unsigned int *)0)
    {
      ERROR_MSG("Not enough memory for " << numSpecs << " spectra");
      goto load_fail;
    }

    if (subversion > 0)
    {
      // New in 2.1: numPeaks stored as integer to accommodate spectra w/ > 64k peaks
      count = fread(numPeaks, sizeof(unsigned int), numSpecs, fp);
      if (count != numSpecs)
      {
        free(numPeaks);
        goto load_fail;
      }
    }
    else
    {
      unsigned short* numPeaks_short =
          (unsigned short *)malloc(sizeof(unsigned short) * numSpecs);

      if (numPeaks_short == (unsigned short *)0)
      {
        ERROR_MSG("Not enough memory for " << numSpecs << " spectra");
        free(numPeaks);
        goto load_fail;
      }

      count = fread(numPeaks_short, sizeof(unsigned short), numSpecs, fp);
      if (count != numSpecs)
      {
        free(numPeaks);
        free(numPeaks_short);
        goto load_fail;
      }
      for (unsigned int i = 0; i < numSpecs; i++)
      {
        numPeaks[i] = (int)numPeaks_short[i];
      }
      free(numPeaks_short);
    }

    for (i = 0; i < numSpecs; i++)
    {
      if (numPeaks[i] > dataSize)
      {
        dataSize = numPeaks[i];
      }

      //cout << "numPeaks[" << i << "] = " << numPeaks[i] << endl;
    }
    dataSize = 3 * dataSize; // Each spectrum 'line' has 3 values
    data = (float *)malloc(sizeof(float) * dataSize);
    if (data == (float *)0)
    {
      ERROR_MSG("Not enough memory for " << dataSize << " floats.");
      goto load_fail;
    }

    for (i = 0; i < numSpecs; i++)
    {
      unsigned int numValues = 3 * numPeaks[i];
      count = fread(data, sizeof(float), numValues, fp);
      if (count != numValues)
      {
        free(numPeaks);
        free(data);
        goto load_fail;
      }

      specs[i].resize(numPeaks[i]);
      specs[i].psmList.resize(0);
      for (unsigned int p = 0, dataIdx = 0; dataIdx < numValues;
          dataIdx += 3, p++)
      {
        specs[i][p].set(data[dataIdx], data[dataIdx + 1]);
        specs[i].setTolerance(p, data[dataIdx + 2]);
      }
    }

    free(numPeaks);
    free(data);
    fclose(fp);

    return true;

    load_fail:

    fclose(fp);

    return false;

  }

  bool SpecSet::loadPklBin_3(FILE* fp, int numSpecs, char subversion)
  {
    specs.resize(numSpecs);

    map<string, unsigned short> versions;
    if (!readStringMapFromBinaryStream<unsigned short>(fp, versions))
    {
      ERROR_MSG("Error reading version info");
      fclose(fp);
      return false;
    }

    for (unsigned int i = 0; i < numSpecs; i++)
    {
      if (!specs[i].loadFromBinaryStream(fp, versions))
      {
        ERROR_MSG("Failed to load spectrum " << i);
        fclose(fp);
        return false;
      }
    }

    fclose(fp);
    return true;
  }

  //-------------------------------------------------------------------------
  //  load - loads a SpecSet from a binary format file.
  int SpecSet::loadPklBin(const char * filename,
                          const char * psmFilename /* = 0x0 */,
                          const char * peaksFilename /* = 0x0 */)
  {
    FILE *fp;
    unsigned int numSpecs = 0;
    unsigned int i, p, dataIdx, count, numValues;
    char version = 0;
    char subversion = 0;

    fp = fopen(filename, "rb");
    if (fp == 0)
    {
      ERROR_MSG("Opening " << filename);
      return 0;
    }

    count = fread(&numSpecs, sizeof(unsigned int), 1, fp); // Number of spectra in the file
    if (count != 1)
    {
      ERROR_MSG("Reading number of spectra from " << filename);
      fclose(fp);
      return 0;
    }

    // Read header information, including which version the file was saved in
    bool oldVersion = true;
    if (numSpecs != 0)
    {
      WARN_MSG("PKLBIN file non-zero first value.");
      WARN_MSG("Assuming first value is number of spectra [" << numSpecs << "]");
    }
    else
    {
      oldVersion = false;
      version = 0;
      count = fread(&version, sizeof(char), 1, fp); // Version of file
      if (count != 1)
      {
        ERROR_MSG("Reading version number from " << filename);
        fclose(fp);
        return 0;
      }
      subversion = 0;
      count = fread(&subversion, sizeof(char), 1, fp); // Sub-version of file
      if (count != 1)
      {
        ERROR_MSG("Reading sub-version number from " << filename);
        fclose(fp);
        return 0;
      }
      count = fread(&numSpecs, sizeof(unsigned int), 1, fp); // Number of spectra in the file
      if (count != 1)
      {
        ERROR_MSG("Reading number of spectra from " << filename);
        fclose(fp);
        return 0;
      }
    }

    // Load the SpecSet and Spectrum data fields here
    // Each version has its own load method that should be called here
    if (oldVersion or version == 1)
    {
      // Load version 1
      if (!loadPklBin_1(fp, numSpecs, oldVersion, subversion))
      {
        ERROR_MSG("Reading " << filename);
        return 0;
      }
    }
    else if (version == 2)
    {
      // Load version 2
      if (!loadPklBin_2(fp, numSpecs, subversion))
      {
        ERROR_MSG("Reading " << filename);
        return 0;
      }
    }
    else if (version == 3)
    {
      // Load version 2
      if (!loadPklBin_3(fp, numSpecs, subversion))
      {
        ERROR_MSG("Reading " << filename);
        return 0;
      }
    }
    else
    {
      ERROR_MSG("Found unsupported version " << (int)version);
      ERROR_MSG("Reading " << filename);
      return 0;
    }
    /*
     if (countSpectraOnly)
     return numSpecs;
     */

    if (psmFilename)
    {
      PeptideSpectrumMatchSet psmSetTemp;
      psmSetTemp.loadFromFile(psmFilename);
      //DEBUG_VAR(psmSetTemp.size());
      psmSetTemp.addSpectra(this, true);
    }

    if (peaksFilename)
    {
      //DEBUG_TRACE;
      SpecSet tempMatchedPeaks;
      if (tempMatchedPeaks.loadPklBin(peaksFilename) <= 0)
      {
        ERROR_MSG("Problem loading matched peaks from [" << peaksFilename << "]");
        return 0;
      }

      //DEBUG_VAR(tempMatchedPeaks.size());
      for (int i = 0; i < tempMatchedPeaks.size(); i++)
      {
        if (specs[i].psmList.size() == 0)
        {
          continue;
        }
        specs[i].psmList.front()->m_matchedPeaks.resize(tempMatchedPeaks[i].size());
        for (int j = 0; j < tempMatchedPeaks[i].size(); j++)
        {
          specs[i].psmList.front()->m_matchedPeaks[j].set((int)tempMatchedPeaks[i][j][0],
                                                          (int)tempMatchedPeaks[i][j][1]);
        }
      }
      //DEBUG_TRACE;
    }

    if (m_index != 0x0)
    {
      computeIndex();
    }

    string fullPath = filename;
    FilenameManager mngr(fullPath);
    for (unsigned int i = 0; i < specs.size(); i++)
    {
      if (specs[i].fileName.length() == 0)
      {
        specs[i].fileName = mngr.getFilenameWithExtension();
      }
    }

    return specs.size();
  }

  void SpecSet::saveTags(const std::string & filename, int tagLen, int gap, int tagNum, float fragTol, bool isPRMSpec)
  {
	FILE *fp = fopen(filename.c_str(), "w");
	if (fp == 0) {
	  ERROR_MSG("Can not open: " << filename);
	  return;
	}
	DEBUG_MSG("Generated tag file: " << filename);

	float ionOffset = (isPRMSpec)? 0 : AAJumps::massHion;

	AAJumps jumps(1);
	for (int i = 0; i < specs.size(); i++) {

	  specs[i].addZPMpeaks22(fragTol, ionOffset, true, 1);
	  for(int k=specs[i].size()-1; k>-1; k--){
			if( specs[i][k][0] > specs[i].parentMass ) specs[i].removePeak(k);
			else break;
	  }

	  fprintf(fp, ">>%d$\t%.2f\t%d\n", i+1, specs[i].parentMass, specs[i].parentCharge);

	  list<Tag> tags;
	  ExtractTags(specs[i], tags, fragTol, tagLen, gap, tagNum);

	  for (list<Tag>::iterator itr = tags.begin(); itr != tags.end(); itr++) {
		string seq;
		for (int c = 0; c < itr->sequence.size(); c++) {
		  seq += jumps.aaLetters[itr->sequence[c]];
		}
		fprintf(fp, "%.3f\t%s\t%.3f\t%.3f\n", itr->flankingPrefix-ionOffset, seq.c_str(), itr->flankingSuffix+ionOffset, itr->score);
	  }
	  fprintf(fp, "\n");
	}
	fclose(fp);
  }

  // -------------------------------------------------------------------------
  /* DEPRECATED
   *
   short SpecSet::SaveSpecSet_pklbin(const char *filename,
   const char *binFilename)
   {
   if (binFilename)
   {
   WARN_MSG("Bin filename [" << binFilename << "] passed to SaveSpecSet_pklbin");
   WARN_MSG("Bin file [" << binFilename << "] will NOT be written.");
   }
   return (short)savePklBin(filename);
   }
   // -------------------------------------------------------------------------
   unsigned int SpecSet::LoadSpecSet_pklbin(const char *filename,
   bool countSpectraOnly)
   {
   if (countSpectraOnly == true)
   {
   WARN_MSG("countSpectraOnly is no longer supported");
   }
   return loadPklBin(filename);
   }
   */
  // -------------------------------------------------------------------------
  void SpecSet::SaveScanNums(const char *filename)
  {
    vector<unsigned int> scanNums(specs.size());
    for (unsigned int i = 0; i < specs.size(); i++)
      scanNums[i] = specs[i].scan;
    Save_binArray(filename, scanNums);
  }
  // -------------------------------------------------------------------------
  /*
   short SaveSpecSet_pklbin(const char *filename, vector<Spectrum *> specs)
   {
   FILE *fp;
   unsigned int numSpecs = specs.size();
   unsigned int i, p;

   fp = fopen(filename, "wb");
   if (fp == 0)
   {
   cerr << "ERROR: cannot open " << filename << "\n";
   return -1;
   }

   fwrite(&numSpecs, sizeof(int), 1, fp); // Number of spectra in the file

   unsigned short *numPeaks = (unsigned short *)malloc(sizeof(short)
   * numSpecs);
   unsigned short maxNumPeaks = 0;
   for (i = 0; i < numSpecs; i++)
   {
   numPeaks[i] = specs[i]->size();
   maxNumPeaks = max(maxNumPeaks, numPeaks[i]);
   }
   fwrite(numPeaks, sizeof(short), numSpecs, fp); // Number of peaks per spectrum in the file

   float *peaksBuffer = (float *)malloc(2 * (maxNumPeaks + 1) * sizeof(float));
   unsigned int pbIdx;
   for (i = 0; i < numSpecs; i++)
   {
   peaksBuffer[0] = specs[i]->parentMass;
   peaksBuffer[1] = (float)specs[i]->parentCharge;
   for (pbIdx = 2, p = 0; p < numPeaks[i]; p++)
   {
   peaksBuffer[pbIdx++] = (*specs[i])[p][0];
   peaksBuffer[pbIdx++] = (*specs[i])[p][1];
   }
   fwrite(peaksBuffer, sizeof(float), 2 * (numPeaks[i] + 1), fp); // [parentMass charge] followed by [masses intensities]
   }

   free(peaksBuffer);
   free(numPeaks);
   fclose(fp);
   return 1;
   }
   */

  void SpecSet::computeIndex()
  {
    unsigned int numSpecs = size();
    m_index->rehash(numSpecs);

    for (unsigned int i = 0; i < numSpecs; i++)
    {
      (*m_index)[specs[i].getUniqueID()] = i;
    }
  }

  // typedef vector<pair<int,vector<int> > > SpecCompareData;

  int SpecSet::compare(SpecSet &toSet, SpecCompareData &cd)
  {
    int differentSpectra = 0;
    int differentPeaks = 0;
    int size1 = size();
    int size2 = toSet.size();

    if (size1 != size2)
    {
      vector<int> aux;
      aux.push_back(size1);
      aux.push_back(size2);
      cd.push_back(pair<int, vector<int> >(-1, aux));
      return -1;
    }

    for (unsigned i = 0; i < size1; i++)
    {
      vector<int> aux;
      int ret = specs[i].compare(toSet[i], aux);
      if (ret == -1)
      {
        differentSpectra++;
      }
      else if (ret)
      {
        differentSpectra++;
        differentPeaks += ret;
      }
      cd.push_back(pair<int, vector<int> >(i, aux));
    }
    return differentSpectra;
  }

  void SpecSet::removeSpectraBelowMinPeaks(int min_peaks, bool remove)
  {
    for (int i = 0; i < size(); i++)
    {
      if (specs[i].size() < min_peaks)
      {
        if (remove)
        {
          //TODO: Find efficient way to delete
        }
        else
        {
          specs[i].resize(0);
        }
      }
    }
  }

}
