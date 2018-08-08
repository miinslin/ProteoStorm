#include "PeptideSpectrumMatchSet.h"

#define DEBUG_PTM_TABLE 0
#define DEBUG_PTM_TABLE2 0

namespace specnets
{
  const unsigned int SCAN_INDEX = 0;
  const unsigned int SPECID_INDEX = 1;
  const unsigned int SPECFILE_INDEX = 2;
  const unsigned int ANNOTATION_INDEX = 3;
  const unsigned int ORIGANNOTATION_INDEX = 4;
  const unsigned int PROTEIN_INDEX = 5;
  const unsigned int DBINDEX_INDEX = 6;
  const unsigned int NUMMODS_INDEX = 7;
  const unsigned int MATCHORIENT_INDEX = 8;
  const unsigned int STARTMASS_INDEX = 9;
  const unsigned int CHARGE_INDEX = 10;
  const unsigned int SCORE_INDEX = 11;
  const unsigned int PVALUE_INDEX = 12;
  const unsigned int ISDECOY_INDEX = 13;
  const unsigned int FDR_INDEX = 14;
  const unsigned int STRICTENVELOPE_INDEX = 15;
  const unsigned int UNSTRICTENVELOPE_INDEX = 16;
  const unsigned int COMPOUND_INDEX = 17;
  const unsigned int ORGANISM_INDEX = 18;
  const unsigned int FILESCAN_INDEX = 19;
  const unsigned int LIBRARYNAME_INDEX = 20;
  const unsigned int LIBMETADATA_INDEX = 21;
  const unsigned int SMILES_INDEX = 22;
  const unsigned int INCHI_INDEX = 23;
  const unsigned int INCHIAUX_INDEX = 24;
  const unsigned int LIBSEARCHSHAREDPEAKS_INDEX = 25;
  const unsigned int MZERRORPPM_INDEX = 26;
  const unsigned int ABUNDANCE_INDEX = 27;
  const unsigned int PARENTMASSDIFF_INDEX = 28;
  const unsigned int EXACTMASS_INDEX = 29;
  const unsigned int LIBRARYSPECID_INDEX = 30;
  const unsigned int SPECMZ_INDEX = 31;
  const unsigned int DBSTART_INDEX = 32;
  const unsigned int DBEND_INDEX = 33;
  const unsigned int VARIANTGROUP_INDEX = 34;
  const unsigned int PEPTIDEREGION_INDEX = 35;
  const unsigned int PEPTIDEREGIONSTRING_INDEX = 36;
  const unsigned int COMPLEMENTARY_SCORE_INDEX = 37;
  const unsigned int LAST_INDEX = 38; // MAKE SURE TO INCREASE THIS IF YOU ADD MORE!!

  const float MINIMUM_AA_MASS = 57.0214637230;

  // -------------------------------------------------------------------------
  void debugExclusionList(vector<pair<string,float> > * pexclusionList)
  {
    for (int i = 0; i < pexclusionList->size(); i++) {
      DEBUG_MSG((*pexclusionList)[i].first << "  " << (*pexclusionList)[i].second)
    }
    return;
  }

  // -------------------------------------------------------------------------
  bool isExcluded(string cleanAnnotation,
                  int position,
                  int length,
                  float modification,
                  vector<pair<string,float> > * pexclusionList)
  {
    const float EXCLUDE_DELTA = 0.1;

    // If list is null then no exclusions
    if (pexclusionList == 0x0) {
      return false;
    }
    // We will not exclude any mod in a gap
    if (length > 1) {
      return false;
    }

    int pos = position - length;

    // If the position is 0 it is an N-term mod
    if (pos == 0) {
      return false;
    }
    string modChar("X");
    //DEBUG_MSG(cleanAnnotation << "  " << pos);
    modChar[0] = cleanAnnotation[pos];
    //DEBUG_MSG(modChar << "  " << modification);

    bool exclude = false;
    for (int i = 0; i < pexclusionList->size(); i++) {
      if (modChar != (*pexclusionList)[i].first) {
        continue;
      }
      if (fabs(modification - (*pexclusionList)[i].second) > EXCLUDE_DELTA) {
        continue;
      }
      exclude = true;
      break;
    }

    return exclude;
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchSet::createMapIndex(vector<string> & header, vector<int> & mapIndex)
  {
    mapIndex.resize(LAST_INDEX);
    // Initialize all the map entries to -1
    for (unsigned int i = 0; i < LAST_INDEX; i++)
    {
      mapIndex[i] = -1;
    }

    for (int i = 0; i < header.size(); i++)
    {

      if (header[i] == "ScanNum" || header[i] == "Scan#"
          || header[i] == "#Scan#" || header[i] == "ScanNo"
          || header[i] == "Cluster_index" || header[i] == "ClusterIndex")
      {
        mapIndex[SCAN_INDEX] = i;
      }

      if (header[i] == "SpecID" || header[i] == "SpecIndex" || header[i] == "Index")
      {
        mapIndex[SPECID_INDEX] = i;
      }

      if (header[i] == "SpectrumFile" || header[i] == "#SpectrumFile"
          || header[i] == "Original_filepath" || header[i] == "#SpecFile")
      {
        mapIndex[SPECFILE_INDEX] = i;
      }

      if (header[i] == "Annotation" || header[i] == "Peptide"||
          header[i] == "HomologSequence")
      {
        mapIndex[ANNOTATION_INDEX] = i;
      }

      if ((m_doubleLoad == 1 &&  header[i] == "Annotation1") ||
          (m_doubleLoad == 2 &&  header[i] == "Annotation2")) {
        mapIndex[ANNOTATION_INDEX] = i;
      }

      if (header[i] == "OrigAnnotation" || header[i] == "DeNovoSequence")
      {
        mapIndex[ORIGANNOTATION_INDEX] = i;
      }

      if (header[i] == "Protein" || header[i] == "ProteinName")
      {
        mapIndex[PROTEIN_INDEX] = i;
      }
      if ((m_doubleLoad == 1 &&  header[i] == "Protein1") ||
          (m_doubleLoad == 2 &&  header[i] == "Protein2")) {
        mapIndex[PROTEIN_INDEX] = i;
      }

      if (header[i] == "dbIndex" || header[i] == "ProteinIndex")
      {
        mapIndex[DBINDEX_INDEX] = i;
      }

      if (header[i] == "numMods")
      {
        mapIndex[NUMMODS_INDEX] = i;
      }

      if (header[i] == "matchOrientation")
      {
        mapIndex[MATCHORIENT_INDEX] = i;
      }

      if (header[i] == "startMass")
      {
        mapIndex[STARTMASS_INDEX] = i;
      }

      if (header[i] == "Charge" ||
          header[i] == "Charge1" ||
          header[i] == "Charge2")
      {
        mapIndex[CHARGE_INDEX] = i;
      }

      if (header[i] == "MQScore" || header[i] == "MSGFScore" ||
          header[i] == "SpecProb")
      {
        mapIndex[SCORE_INDEX] = i;
      }
      if ((m_doubleLoad == 1 &&  header[i] == "svm1-score") ||
          (m_doubleLoad == 2 &&  header[i] == "svm2-score")) {
        mapIndex[SCORE_INDEX] = i;
      }

      if (header[i] == "MSGFScore")
      {
        mapIndex[COMPLEMENTARY_SCORE_INDEX] = i;
       }
      
      if (header[i] == "p-value" || header[i] == "P-value" ||
          header[i] == "SpecEValue" || header[i] == "QValue" ||
          header[i] == "Probability")
      {
        mapIndex[PVALUE_INDEX] = i;
      }

      if (header[i] == "isDecoy")
      {
        mapIndex[ISDECOY_INDEX] = i;
      }

      if (header[i] == "FDR" || header[i] == "InspectFDR")
      {
        mapIndex[FDR_INDEX] = i;
      }

      if (header[i] == "StrictEnvelope" || header[i] == "StrictEnvelopeScore")
      {
        mapIndex[STRICTENVELOPE_INDEX] = i;
      }

      if (header[i] == "UnstrictEnvelope" || header[i] == "UnstrictEvelopeScore")
      {
        mapIndex[UNSTRICTENVELOPE_INDEX] = i;
      }

      if (header[i] == "CompoundName")
      {
        mapIndex[COMPOUND_INDEX] = i;
      }

      if (header[i] == "Organism")
      {
        mapIndex[ORGANISM_INDEX] = i;
      }

      if (header[i] == "FileScanUniqueID")
      {
        mapIndex[FILESCAN_INDEX] = i;
      }

      if (header[i] == "LibraryName")
      {
        mapIndex[LIBRARYNAME_INDEX] = i;
      }

      if (header[i] == "LibMetaData")
      {
        mapIndex[LIBMETADATA_INDEX] = i;
      }

      if (header[i] == "Smiles")
      {
        mapIndex[SMILES_INDEX] = i;
      }

      if (header[i] == "Inchi")
      {
        mapIndex[INCHI_INDEX] = i;
      }

      if (header[i] == "InchiAux")
      {
        mapIndex[INCHIAUX_INDEX] = i;
      }

      if (header[i] == "LibSearchSharedPeaks")
      {
        mapIndex[LIBSEARCHSHAREDPEAKS_INDEX] = i;
      }

      if (header[i] == "mzErrorPPM" || header[i] == "PrecursorMZError" 
          || header[i] == "PrecursorError(Da)")
      {
        mapIndex[MZERRORPPM_INDEX] = i;
      }

      if (header[i] == "Abundance")
      {
        mapIndex[ABUNDANCE_INDEX] = i;
      }

      if (header[i] == "ParentMassDiff")
      {
        mapIndex[PARENTMASSDIFF_INDEX] = i;
      }
      if (header[i] == "ExactMass")
      {
        mapIndex[EXACTMASS_INDEX] = i;
      }

      if (header[i] == "LibrarySpectrumID")
      {
        mapIndex[LIBRARYSPECID_INDEX] = i;
      }

      if (header[i] == "SpecMZ" || header[i] == "Precursor")
      {
        mapIndex[SPECMZ_INDEX] = i;
      }
      
      if (header[i] == "StartAA" || header[i] == "start")
      {
        mapIndex[DBSTART_INDEX] = i;
      }
      if (header[i] == "EndAA" || header[i] == "end")
      {
        mapIndex[DBEND_INDEX] = i;
      }
      if (header[i] == "variantGroup")
      {
        mapIndex[VARIANTGROUP_INDEX] = i;
      }
      if (header[i] == "peptideRegionNum")
      {
        mapIndex[PEPTIDEREGION_INDEX] = i;
      }
      if (header[i] == "peptideRegion")
      {
        mapIndex[PEPTIDEREGIONSTRING_INDEX] = i;
      }
    }

    return;
  }

  float computeTotalMass(psmPtr p, AAJumps & jumps)
  {
    string cleanAnno;
    PeptideSpectrumMatch::getUnmodifiedPeptide(p->m_annotation, cleanAnno);
    float totalMass = jumps.getPeptideMass(cleanAnno);
	    
    vector<float> modifications;
    vector<unsigned int> positions;
    p->getModificationsAndPositions(modifications, positions);
    
    for (int i = 0; i < modifications.size(); i++) {
      totalMass += modifications[i];
    }
    
    return totalMass;
    
  }

  inline int convertStringToInt(const char * stringValue)
  {
    int returnValue = -1;
    returnValue = atoi(stringValue);
    return returnValue;
  }

  inline float convertStringToFloat(const char * stringValue)
  {
    float returnValue = -1.0;
    returnValue = atof(stringValue);
    return returnValue;
  }

  inline int convertStringToInt(string & stringValue)
  {
    int returnValue = -1;
    returnValue = atoi(stringValue.c_str());
    return returnValue;
  }

  inline float convertStringToFloat(string & stringValue)
  {
    float returnValue = -1.0;
    returnValue = atof(stringValue.c_str());
    return returnValue;
  }

  // -------------------------------------------------------------------------
  template<typename T> bool parseLineT(psmPtr currMatch,
                                       vector<T> & line,
                                       vector<int> & mapIndex,
                                       bool zeroIndexed,
                                       bool isInspect)
  {
    bool haveScanNum = false;
    currMatch->m_scanNum = -1;
    int scanIndex = mapIndex[SCAN_INDEX];
    if (scanIndex != -1)
    {
      currMatch->m_scanNum = convertStringToInt(line[scanIndex]);
      if (currMatch->m_scanNum >= 0)
      {
        haveScanNum = true;
      }
    }

    int specIdxIndex = mapIndex[SPECID_INDEX];
    if (!haveScanNum && specIdxIndex != -1)
    {
      currMatch->m_scanNum = convertStringToInt(line[specIdxIndex]);
    }

    if (currMatch->m_scanNum >= 0 && zeroIndexed)
    {
      currMatch->m_scanNum += 1;
    }

    int spectrumFileIndex = mapIndex[SPECFILE_INDEX];
    if (spectrumFileIndex != -1)
    {
      currMatch->m_spectrumFile = line[spectrumFileIndex];
    }

    int annotationIndex = mapIndex[ANNOTATION_INDEX];
    if (annotationIndex != -1)
    {
      if (isInspect)
      {
        currMatch->inspectToSpecNets(line[annotationIndex],
                                     currMatch->m_annotation);
        currMatch->m_origAnnotation = line[annotationIndex];
      }
      else
      {
        currMatch->stripPrecedingAndFollowing(line[annotationIndex],
                                              currMatch->m_annotation);
      }
      //DEBUG_VAR(currMatch->m_annotation);
    }

    int origAnnotationIndex = mapIndex[ORIGANNOTATION_INDEX];
    if (origAnnotationIndex != -1)
    {
      if (isInspect)
      {
        currMatch->inspectToSpecNets(line[origAnnotationIndex],
                                     currMatch->m_origAnnotation);
      } else {
        currMatch->stripPrecedingAndFollowing(line[origAnnotationIndex],
                                     currMatch->m_origAnnotation);
      }
    }

    int proteinNameIndex = mapIndex[PROTEIN_INDEX];
    if (proteinNameIndex != -1)
    {
      currMatch->m_protein = line[proteinNameIndex];
      if (currMatch->m_protein.substr(0, 3).compare("XXX") == 0 ||
          currMatch->m_protein.substr(0, 3).compare("REV") == 0)
      {
        currMatch->m_isDecoy = true;
      }
      else
      {
        currMatch->m_isDecoy = false;
      }
    }

    int dbIndexIndex = mapIndex[DBINDEX_INDEX];
    if (dbIndexIndex != -1)
    {
      currMatch->m_dbIndex = convertStringToInt(line[dbIndexIndex]);
    }

    int numModsIndex = mapIndex[NUMMODS_INDEX];
    if (numModsIndex != -1)
    {
      currMatch->m_numMods = convertStringToInt(line[numModsIndex]);
    }

    int matchOrientationIndex = mapIndex[MATCHORIENT_INDEX];
    if (matchOrientationIndex != -1)
    {
      currMatch->m_matchOrientation =
          convertStringToInt(line[matchOrientationIndex]);
    }

    int startMassIndex = mapIndex[STARTMASS_INDEX];
    if (startMassIndex != -1)
    {
      currMatch->m_startMass = convertStringToFloat(line[startMassIndex]);
    }

    int chargeIndex = mapIndex[CHARGE_INDEX];
    if (chargeIndex != -1)
    {
      currMatch->m_charge = convertStringToInt(line[chargeIndex]);
    }

    int scoreIndex = mapIndex[SCORE_INDEX];
    if (scoreIndex != -1)
    {
      currMatch->m_score = convertStringToFloat(line[scoreIndex]);
    }

    int compScoreIndex = mapIndex[COMPLEMENTARY_SCORE_INDEX];
    if (compScoreIndex != -1)
    {
      currMatch->m_compScore = convertStringToFloat(line[compScoreIndex]);
    }
    
    int pvalueIndex = mapIndex[PVALUE_INDEX];
    if (pvalueIndex != -1)
    {
      currMatch->m_pValue = convertStringToFloat(line[pvalueIndex]);
    }

    int decoyIndex = mapIndex[ISDECOY_INDEX];
    if (decoyIndex != -1)
    {
      int isDecoy = false;
      isDecoy = convertStringToInt(line[decoyIndex]);
      currMatch->m_isDecoy = (bool)isDecoy;
    }

    int fdrIndex = mapIndex[FDR_INDEX];
    if (fdrIndex != -1)
    {
      currMatch->m_fdr = convertStringToFloat(line[fdrIndex]);
    }

    int strictenvIndex = mapIndex[STRICTENVELOPE_INDEX];
    if (strictenvIndex != -1)
    {
      currMatch->m_strict_envelope_score =
          convertStringToFloat(line[strictenvIndex]);
    }

    int unstrictenvIndex = mapIndex[UNSTRICTENVELOPE_INDEX];
    if (unstrictenvIndex != -1)
    {
      currMatch->m_unstrict_envelope_score =
          convertStringToFloat(line[unstrictenvIndex]);
    }

    int compoundNameIndex = mapIndex[COMPOUND_INDEX];
    if (compoundNameIndex != -1)
    {
      currMatch->m_compound_name = (line[compoundNameIndex]);
    }

    int organismIndex = mapIndex[ORGANISM_INDEX];
    if (organismIndex != -1)
    {
      currMatch->m_organism = line[organismIndex];
    }

    int fileScanUniqueIndex = mapIndex[FILESCAN_INDEX];
    if (fileScanUniqueIndex != -1)
    {
      //currMatch->m_spectrumFile = line[fileScanUniqueIndex];
    }

    int libraryNameIndex = mapIndex[LIBRARYNAME_INDEX];
    if (libraryNameIndex != -1)
    {
      currMatch->m_library_name = line[libraryNameIndex];
    }

    int libmetadataIndex = mapIndex[LIBMETADATA_INDEX];
    if (libmetadataIndex != -1)
    {
      currMatch->m_notes = (line[libmetadataIndex]);
    }

    int smilesIndex = mapIndex[SMILES_INDEX];
    if (smilesIndex != -1)
    {
      currMatch->m_smiles = (line[smilesIndex]);
    }

    int inchiIndex = mapIndex[INCHI_INDEX];
    if (inchiIndex != -1)
    {
      currMatch->m_InChI = (line[inchiIndex]);
    }

    int inchiAuxIndex = mapIndex[INCHIAUX_INDEX];
    if (inchiAuxIndex != -1)
    {
      currMatch->m_InChI_Aux = (line[inchiAuxIndex]);
    }

    int sharedpeakIndex = mapIndex[LIBSEARCHSHAREDPEAKS_INDEX];
    if (sharedpeakIndex != -1)
    {
      currMatch->m_shared_peaks = convertStringToInt(line[sharedpeakIndex]);
    }
    int ppmerrIndex = mapIndex[MZERRORPPM_INDEX];
    if (ppmerrIndex != -1)
    {
      currMatch->m_mz_error_ppm = convertStringToFloat(line[ppmerrIndex]);
    }

    int abundanceIndex = mapIndex[ABUNDANCE_INDEX];
    if (abundanceIndex != -1)
    {
      currMatch->m_abundance = convertStringToFloat(line[abundanceIndex]);
    }

    int massdiffIndex = mapIndex[PARENTMASSDIFF_INDEX];
    if (massdiffIndex != -1)
    {
      currMatch->m_parentmass_difference =
          convertStringToFloat(line[massdiffIndex]);
    }

    int specmzIndex = mapIndex[SPECMZ_INDEX];
    if (specmzIndex != -1)
    {
      currMatch->m_mz = convertStringToFloat(line[specmzIndex]);
    }

    int exactMassIndex = mapIndex[EXACTMASS_INDEX];
    if (exactMassIndex != -1)
    {
      currMatch->m_exactmass = convertStringToFloat(line[exactMassIndex]);
    }

    int spectrumidIndex = mapIndex[LIBRARYSPECID_INDEX];
    if (spectrumidIndex != -1)
    {
      currMatch->m_spectrumID = (line[spectrumidIndex]);
    }

    int startIndex = mapIndex[DBSTART_INDEX];
    if (startIndex != -1)
    {
      currMatch->m_startIndex = convertStringToInt(line[startIndex]);
    }
    int endIndex = mapIndex[DBEND_INDEX];
    if (endIndex != -1)
    {
      currMatch->m_endIndex = convertStringToInt(line[endIndex]);
    }
    int varGroupIndex = mapIndex[VARIANTGROUP_INDEX];
    if (varGroupIndex != -1)
    {
      currMatch->m_variantGroup = convertStringToInt(line[varGroupIndex]);
    }
    int pepRegIndex = mapIndex[PEPTIDEREGION_INDEX];
    if (pepRegIndex != -1)
    {
      currMatch->m_peptideRegionGroup = line[pepRegIndex];
    }
    int pepRegStringIndex = mapIndex[PEPTIDEREGIONSTRING_INDEX];
    if (pepRegIndex != -1)
    {
      currMatch->m_peptideRegion = line[pepRegIndex];
    }

    return true;
  }

  // -------------------------------------------------------------------------
  psmPtr & PeptideSpectrumMatchSet::operator[](unsigned int i)
  {
    return m_psmSet[i];
  }

  // -------------------------------------------------------------------------
  const psmPtr & PeptideSpectrumMatchSet::operator[](unsigned int i) const
  {
    return m_psmSet[i];
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchSet::removePsmSetItem(psmPtr item)
  {
    for (int i = 0; i < m_psmSet.size(); i++)
    {

      std::tr1::shared_ptr<PeptideSpectrumMatch> p =
          std::tr1::dynamic_pointer_cast<PeptideSpectrumMatch>(m_psmSet[i]); //I AM A PSMNETWORK POINTER

      // 	  		PeptideSpectrumMatch * p = (PeptideSpectrumMatch *)m_psmSet[i].get();

      if (p == item)
      {
        m_psmSet.erase(m_psmSet.begin() + i, m_psmSet.begin() + i + 1);
        return;
      }
    }
    DEBUG_MSG("PeptideSpectrumMatchSet::removePsmSetItem:: pointer not found!");
  }
  // -------------------------------------------------------------------------
  PeptideSpectrumMatchSet & PeptideSpectrumMatchSet::operator=(PeptideSpectrumMatchSet &other)
  {
    if (this == &other)
    {
      return *this;
    }
    m_psmSet.clear();
    for (int i = 0; i < other.m_psmSet.size(); i++)
    {
      psmPtr currMatch(new PeptideSpectrumMatch);
      *currMatch = *(other.m_psmSet[i]);
      m_psmSet.push_back(currMatch);
    }
    return *this;
  }
  // -------------------------------------------------------------------------
  unsigned int PeptideSpectrumMatchSet::push_back(const psmPtr &other)
  {
    psmPtr currMatch(new PeptideSpectrumMatch);
    *currMatch = *other;
    m_psmSet.push_back(currMatch);
    return m_psmSet.size();
  }

  // -------------------------------------------------------------------------
  unsigned int PeptideSpectrumMatchSet::resize(unsigned int newSize)
  {
    m_psmSet.resize(newSize);
    return m_psmSet.size();
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchSet::getPSMSet(SpecSet * spectra)
  {
    m_psmSet.clear(); //empty any current psm annotations

    for (int i = 0; i < spectra->size(); i++)
    {
      Spectrum spectrum = (*spectra)[i];
      if (spectrum.psmList.size() > 0)
      {
        list<psmPtr>::iterator it;

        for (it = spectrum.psmList.begin(); it != spectrum.psmList.end(); it++)
        {
          psmPtr currMatch(new PeptideSpectrumMatch);
          (*currMatch) = **it;
          currMatch->m_scanNum = spectrum.scan;
          m_psmSet.push_back(currMatch);
        }
      }
    }
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchSet::makePSMSetShared(SpecSet * spectra)
  {
    m_psmSet.clear(); //empty any current psm annotations
    for (int i = 0; i < spectra->size(); i++)
    {
      Spectrum spectrum = (*spectra)[i];
      if (spectrum.psmList.size() > 0)
      {
        list<psmPtr>::iterator it;

        for (it = spectrum.psmList.begin(); it != spectrum.psmList.end(); it++)
        {
          m_psmSet.push_back(*it);
        }
      }
    }
    return;
  }

  // -------------------------------------------------------------------------
  unsigned int PeptideSpectrumMatchSet::addSpectraByFilename(SpecSet * spectra,
                                                             bool addToSpectra,
                                                             bool suppressWarnings)
  {
    vector<psmPtr>::iterator it;
    spectra->index();
    unsigned int numMatched = 0;
    for (it = m_psmSet.begin(); it != m_psmSet.end(); it++)
    {
      Spectrum temp;
      temp.fileName = (**it).m_spectrumFile;
      temp.scan = (**it).m_scanNum;
      string key = temp.getUniqueID();
      Spectrum * spectrum = spectra->getIndex(key);
      if (spectrum == NULL)
      {
        if (!suppressWarnings)
          WARN_MSG("Unable to find spectrum for key " << key);
      }
      else
      {
        if (addToSpectra)
        {
          spectrum->psmList.push_back(*it);
        }
        (*it)->m_spectrum = spectrum;
        numMatched++;
      }
    }

    return numMatched;
  }
  // -------------------------------------------------------------------------
  unsigned int PeptideSpectrumMatchSet::addSpectra(SpecSet * spectra,
                                                   bool addToSpectra,
                                                   bool suppressWarnings)
  {
    vector<psmPtr>::iterator it;
    unsigned int numMatched = 0;
    for (it = m_psmSet.begin(); it != m_psmSet.end(); it++)
    {
      int scanNum = (*it)->m_scanNum;
      Spectrum * spectrum = spectra->getScan(scanNum); //getScan returns null if spectrum is not found

      if (spectrum == NULL)
      {
        if (!suppressWarnings)
          WARN_MSG("Unable to find spectrum for scan " << scanNum);
      }
      else
      {
        if (addToSpectra)
        {
          spectrum->psmList.push_back(*it);
        }
      }
      (*it)->m_spectrum = spectrum;
      numMatched++;
    }

    return numMatched;
  }
  // -------------------------------------------------------------------------
  unsigned int PeptideSpectrumMatchSet::addSpectra(SpecSet * spectra,
                                                   string filename,
                                                   bool addToSpectra,
                                                   bool suppressWarnings)
  {
    vector<psmPtr>::iterator it;
    unsigned int numMatched = 0;
    for (it = m_psmSet.begin(); it != m_psmSet.end(); it++)
    {
      string psmFile;
      extractFileName((*it)->m_spectrumFile, psmFile, false);

      if (filename.compare(psmFile) == 0)
      {
        int scanNum = (*it)->m_scanNum;

        Spectrum * spectrum = spectra->getScan(scanNum);
        if (spectrum == NULL)
        {
          if (!suppressWarnings)
            WARN_MSG("Unable to find spectrum for scan " << scanNum);
        }
        else
        {
          if (addToSpectra)
          {
            spectrum->psmList.push_back(*it);
          }
        }
        (*it)->m_spectrum = spectrum;
        numMatched++;
      }
    }
    return numMatched;
  }

  unsigned int PeptideSpectrumMatchSet::addMostSpectra(SpecSet * spectra,
                                                       bool addToSpectra,
                                                       bool suppressWarnings)
  {
    DEBUG_TRACE;
    unsigned int numMatchedFile = addSpectraByFilename(spectra, false, true);
    DEBUG_TRACE;
    unsigned int numMatchedScanOnly = 0; //addSpectra(spectra, false, true);

    if (numMatchedFile < size())
    {
      numMatchedScanOnly = addSpectra(spectra, false, true);
    }
    DEBUG_TRACE;

    if (numMatchedFile >= numMatchedScanOnly)
    {
      DEBUG_MSG("Adding spectra by filename + scan, matched " << numMatchedFile << " PSMs");
      addSpectraByFilename(spectra, addToSpectra, suppressWarnings);
      return numMatchedFile;
    }
    else
    {
      DEBUG_MSG("Adding spectra by scan only, matched " << numMatchedScanOnly << " PSMs");
      addSpectra(spectra, addToSpectra, suppressWarnings);
      return numMatchedScanOnly;
    }
  }

  int PeptideSpectrumMatchSet::cluster(map<int, list<int> >& clusterInfo,
                                       short mergeType)
  {
    WARN_MSG("Calling deprecated version of cluster(), assuming one spectrum file in PSM set...");
    tr1::unordered_set<string> filenames;
    for (unsigned int i = 0; i < m_psmSet.size(); i++)
    {
      FilenameManager mngr(m_psmSet[i]->m_spectrumFile);
      filenames.insert(mngr.filename);
    }
    if (filenames.size() > 1)
    {
      ERROR_MSG("Found multiple spectrum files in PSM set, must use non-deprecated cluster method");
      abort();
    }
    map<int, list<pair<int, string> > > newClustInfo;
    for (map<int, list<int> >::iterator clustIt = clusterInfo.begin();
        clustIt != clusterInfo.end(); clustIt++)
    {
      list<pair<int, string> > newL;
      for (list<int>::const_iterator childIt = clustIt->second.begin();
          childIt != clustIt->second.end(); childIt++)
      {
        newL.push_back(pair<int, string>(*childIt, *filenames.begin()));
      }
      newClustInfo[clustIt->first] = newL;
    }
    return cluster(newClustInfo, mergeType);
  }

  // -------------------------------------------------------------------------
  int PeptideSpectrumMatchSet::cluster(map<int, list<pair<int, string> > >& clusterInfo,
                                       short mergeType)
  {
    tr1::unordered_map<string, map<int, unsigned int> > scanToPSMIdx;
    for (unsigned int i = 0; i < m_psmSet.size(); i++)
    {
      FilenameManager mngr(m_psmSet[i]->m_spectrumFile);
      //DEBUG_MSG("Have PSM for (" << m_psmSet[i]->m_scanNum << ", \'" << mngr.filename << "\')");
      if (scanToPSMIdx.count(mngr.filename) > 0)
      {
        scanToPSMIdx[mngr.filename][m_psmSet[i]->m_scanNum] = i;
      }
      else
      {
        map<int, unsigned int> tempMap;
        tempMap[m_psmSet[i]->m_scanNum] = i;
        scanToPSMIdx[mngr.filename] = tempMap;
      }

    }

    //DEBUG_VAR(scanToPSMIdx.size());

    if (m_psmSet.size() < 2)
    {
      return 0;
    }

    vector<psmPtr> newPSMSet(0);
    set<int> clusteredScans;

    float numClusters = 0;
    float numAgreeClusters = 0;

    for (map<int, list<pair<int, string> > >::iterator clustIt =
        clusterInfo.begin(); clustIt != clusterInfo.end(); clustIt++)
    {
      int newScan = clustIt->first;

      //DEBUG_VAR(newScan);

      psmPtr clusteredPSM;
      bool foundPSM = false;
      map<string, pair<int, psmPtr> > peptideCounts;
      set<string> unModPeptides;
      if (mergeType == 0)
      {
        for (list<pair<int, string> >::const_iterator childIt =
            clustIt->second.begin(); childIt != clustIt->second.end();
            childIt++)
        {
          FilenameManager mngr(childIt->second);
          if (scanToPSMIdx.count(mngr.filename) == 0
              || scanToPSMIdx[mngr.filename].count(childIt->first) == 0)
          {
            //DEBUG_MSG("Failed to locate (" << childIt->first << ", \'" << mngr.filename << "\')");
            continue;
          }
          foundPSM = true;
          clusteredScans.insert(childIt->first);
          unsigned int psmIdx = scanToPSMIdx[mngr.filename][childIt->first];
          psmPtr nextPSM(m_psmSet[psmIdx]);
          unModPeptides.insert(AAJumps::stripMods(nextPSM->m_annotation));

          /*
           DEBUG_VAR(childIt->first);
           DEBUG_VAR(mngr.filename);
           DEBUG_VAR(nextPSM->m_annotation);
           */

          if (peptideCounts.count(nextPSM->m_annotation) == 0)
          {
            pair<int, psmPtr> newPair(1, nextPSM);
            peptideCounts[nextPSM->m_annotation] = newPair;
          }
          else
          {
            psmPtr prevPSM(peptideCounts[nextPSM->m_annotation].second);
            if (prevPSM->m_pValue > nextPSM->m_pValue)
            {
              peptideCounts[nextPSM->m_annotation].second = nextPSM;
            }
            peptideCounts[nextPSM->m_annotation].first++;
          }
        }

        bool initialized = false;
        int countSoFar = 0;

        for (map<string, pair<int, psmPtr> >::const_iterator countIt =
            peptideCounts.begin(); countIt != peptideCounts.end(); countIt++)
        {

          if (!initialized)
          {
            clusteredPSM = countIt->second.second;
            countSoFar = countIt->second.first;
            initialized = true;
            continue;
          }

          if (countSoFar < countIt->second.first)
          {
            clusteredPSM = countIt->second.second;
            countSoFar = countIt->second.first;
            continue;
          }

          if (clusteredPSM->m_pValue > countIt->second.second->m_pValue)
          {
            clusteredPSM = countIt->second.second;
            countSoFar = countIt->second.first;
            continue;
          }
        }
      }
      if (!foundPSM)
      {
        continue;
      }
      /*
       if (clustIt->second.size() > 3)
       {
       numClusters += 1.0;
       if (unModPeptides.size() == 1)
       {
       numAgreeClusters += 1.0;
       }
       else
       {
       DEBUG_VAR(peptideCounts.size());
       for (map<string, pair<int, psmPtr> >::const_iterator countIt =
       peptideCounts.begin(); countIt != peptideCounts.end(); countIt++)
       {
       DEBUG_VAR(countIt->first);
       }
       }
       }
       */
      //DEBUG_VAR(clusteredPSM->m_annotation);
      clusteredPSM->m_scanNum = newScan;
      newPSMSet.push_back(clusteredPSM);
    }
    /*
     for (unsigned int i = 0; i < m_psmSet.size(); i++)
     {
     if (clusteredScans.count(m_psmSet[i]->m_scanNum) > 0)
     {
     continue;
     }
     newPSMSet.push_back(m_psmSet[i]);
     }
     */

    m_psmSet = newPSMSet;

    //DEBUG_MSG((100.0 * numAgreeClusters)/numClusters << " percent of all clusters have agreeing PSMs (over " << numClusters << " clusters)");

    return clusteredScans.size();
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchSet::getModCounts(map<string, map<float, float> > & mapModCount,
                                             bool useOrig,
                                             vector<pair<string,float> > * pexclusionList)
  {
    //debugExclusionList(pexclusionList);

    AAJumps jumps(1);

    vector<psmPtr>::iterator it;
    for (it = m_psmSet.begin(); it != m_psmSet.end(); it++)
    {
      // Exclude decoys
      if ((*it)->m_isDecoy) {
        continue;
      }
      vector<float> modifications;
      vector<unsigned int> positions;
      vector<unsigned int> lengths;
      (*it)->getModificationsAndPositions(modifications, positions, lengths, useOrig);

      string annotation;
      if (useOrig) {
        annotation = (*it)->m_origAnnotation;
      } else {
        annotation = (*it)->m_annotation;
      }
      if (DEBUG_PTM_TABLE) DEBUG_VAR(annotation);

      string cleanAnnotation;
      PeptideSpectrumMatch::getUnmodifiedPeptide(annotation, cleanAnnotation);
      if (DEBUG_PTM_TABLE)
        DEBUG_VAR(cleanAnnotation);

      for (int i = 0; i < positions.size(); i++)  {
        if (DEBUG_PTM_TABLE) DEBUG_VAR(positions[i]);
        if (DEBUG_PTM_TABLE) DEBUG_VAR(lengths[i]);
        if (DEBUG_PTM_TABLE) DEBUG_VAR(modifications[i]);
        float length = (float)lengths[i];
        float modValue = (float)modifications[i];
        float intModValue = 0.0;
        if (modValue < 0) {
          intModValue = (float)int(modValue - 0.5);
        } else {
          intModValue = (float)int(modValue + 0.5);
        }
        if (DEBUG_PTM_TABLE)
          DEBUG_VAR(intModValue);

        if (isExcluded(cleanAnnotation,  positions[i], lengths[i], modifications[i], pexclusionList)) {
          continue;
        }

        float ntermAdd = 0.0;
        if (positions[i] == 0) {
          mapModCount["nTerm"][intModValue] = mapModCount["nTerm"][intModValue] + 1.0;
          if (DEBUG_PTM_TABLE) DEBUG_MSG("Add 1 to Nterm");
        } else if (positions[i] - lengths[i] == 0) {
          ntermAdd = 1.0;
          mapModCount["nTerm"][intModValue] += 1.0 / (length + ntermAdd);
          if (DEBUG_PTM_TABLE) DEBUG_MSG("Add "<< 1.0 / (length + ntermAdd) << " to Nterm");
        }
        if (DEBUG_PTM_TABLE) DEBUG_VAR(ntermAdd);

        float ctermAdd = 0.0;
        if (positions[i] == cleanAnnotation.length()) {
          if (lengths[i] == 0) {
            mapModCount["cTerm"][intModValue] += 1.0;
            if (DEBUG_PTM_TABLE) DEBUG_MSG("Add 1 to Cterm");
          } else {
            ctermAdd = 1.0;
            mapModCount["cTerm"][intModValue] += 1.0 / (length + ctermAdd);
            if (DEBUG_PTM_TABLE) DEBUG_MSG("Add "<< 1.0 / (length + ctermAdd) << " to Cterm");
          }
        }
        if (DEBUG_PTM_TABLE) DEBUG_VAR(ctermAdd);

#if 0
        int reducedLength = lengths[i];
        for (int k = 0; k < lengths[i]; k++) {
          int pos = positions[i] - k;
          if (DEBUG_PTM_TABLE) DEBUG_VAR(pos);
          float aaMass = getMass(cleanAnnotation[pos - 1]);
          if (DEBUG_PTM_TABLE) DEBUG_VAR(aaMass + modValue);
          if (aaMass + modValue < MINIMUM_AA_MASS) {
            reducedLength--;
            if (DEBUG_PTM_TABLE) DEBUG_VAR(reducedLength);
          }
        }

        for (int k = 0; k < lengths[i]; k++) {
          int pos = positions[i] - k;
          if (DEBUG_PTM_TABLE) DEBUG_VAR(pos);
          string modChar("X");
          modChar[0] = cleanAnnotation[pos - 1];
          if (DEBUG_PTM_TABLE) DEBUG_VAR(modChar);

          float aaMass = getMass(cleanAnnotation[pos - 1]);
          if (DEBUG_PTM_TABLE) DEBUG_VAR(aaMass + modValue);
          if (reducedLength == 0) {
            mapModCount[modChar][intModValue] += 1.0
                / (lengths[i] + ntermAdd + ctermAdd);
            if (DEBUG_PTM_TABLE) DEBUG_MSG("Add " <<
                       1.0 / (lengths[i] + ntermAdd + ctermAdd) <<
                       " " << modChar << " " << intModValue);
          } else {
            if (aaMass + modValue < MINIMUM_AA_MASS) {
              continue;
            }
            mapModCount[modChar][intModValue] += 1.0
                / (reducedLength + ntermAdd + ctermAdd);
            if (DEBUG_PTM_TABLE) DEBUG_MSG("Add " <<
                       1.0 / (reducedLength + ntermAdd + ctermAdd) <<
                       " " << modChar << " " << intModValue);
          }
        }
#endif
        for (int k = 0; k < lengths[i]; k++) {
          int pos = positions[i] - k;
          if (DEBUG_PTM_TABLE) DEBUG_VAR(pos);
          string modChar("X");
          modChar[0] = cleanAnnotation[pos - 1];
          if (DEBUG_PTM_TABLE) DEBUG_VAR(modChar);

          mapModCount[modChar][intModValue] += 1.0
              / (lengths[i] + ntermAdd + ctermAdd);
          if (DEBUG_PTM_TABLE) DEBUG_MSG("Add " <<
                     1.0 / (lengths[i] + ntermAdd + ctermAdd) <<
                     " " << modChar << " " << intModValue);
        }

      } // for (int iChar = 0; iChar < annotation.length(); iChar++) {

    } // for (it = m_psmSet.begin(); it != m_psmSet.end(); it++)

    return;
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchSet::getModCounts(map<string, map<float, float> > & mapModCount,
                                             map<float, vector<string> > mapPrefered,
                                             int nPref,
                                             bool useOrig,
                                             vector<pair<string,float> > * pexclusionList)
  {
    AAJumps jumps(1);

    //DEBUG_TRACE;
    vector<psmPtr>::iterator it;
    for (it = m_psmSet.begin(); it != m_psmSet.end(); it++)
    {
      // Exclude decoys
      if ((*it)->m_isDecoy) {
        continue;
      }
      vector<float> modifications;
      vector<unsigned int> positions;
      vector<unsigned int> lengths;
      (*it)->getModificationsAndPositions(modifications, positions, lengths, useOrig);

      string annotation;
      if (useOrig) {
        annotation = (*it)->m_origAnnotation;
      } else {
        annotation = (*it)->m_annotation;
      }
      if (DEBUG_PTM_TABLE2) DEBUG_VAR(annotation);

      string cleanAnnotation;
      PeptideSpectrumMatch::getUnmodifiedPeptide(annotation, cleanAnnotation);
      if (DEBUG_PTM_TABLE2) DEBUG_VAR(cleanAnnotation);

      for (int i = 0; i < positions.size(); i++) {
        if (DEBUG_PTM_TABLE2) DEBUG_VAR(positions[i]);
        int position = positions[i] - lengths[i];  // the START position
        if (DEBUG_PTM_TABLE2) DEBUG_VAR(position);
        if (DEBUG_PTM_TABLE2) DEBUG_VAR(lengths[i]);
        if (DEBUG_PTM_TABLE2) DEBUG_VAR(modifications[i]);
        int intLength = lengths[i];
        float floatLength = (float)intLength;

        float modValue = (float)modifications[i];
        float intModValue = 0.0;
        if (modValue < 0) {
          intModValue = (float)int(modValue - 0.5);
        } else {
          intModValue = (float)int(modValue + 0.5);
        }
        if (DEBUG_PTM_TABLE2) DEBUG_VAR(intModValue);

        if (isExcluded(cleanAnnotation,  positions[i], lengths[i], modifications[i], pexclusionList)) {
          continue;
        }

        if (intLength > 1 && mapPrefered.find(intModValue) != mapPrefered.end()) {
          vector<string> & vecPref = mapPrefered[intModValue];
          if (DEBUG_PTM_TABLE2) DEBUG_VAR(vecPref.size());
          for (int iPref = 0;
               iPref < vecPref.size() && ((nPref == -1) || iPref < nPref);
               iPref++) {
            // Lets check to see if we want to place the mod on a particular AA
            string preferedAA = vecPref[iPref];
            if (DEBUG_PTM_TABLE2) DEBUG_VAR(preferedAA);
            string gapString = cleanAnnotation.substr(position, intLength);
            if (DEBUG_PTM_TABLE2) DEBUG_VAR(gapString);
            size_t pos = gapString.find(preferedAA);
            if (DEBUG_PTM_TABLE2) DEBUG_VAR(pos);
            if (pos != string::npos) {
              position += pos;
              if (DEBUG_PTM_TABLE2) DEBUG_VAR(position);
              intLength = 1;
              floatLength = 1.0;
              break;  // We found the prefered AA
            }
          }
        }

        float ntermAdd = 0.0;
        if (position == 0 && intLength == 0) {
          mapModCount["nTerm"][intModValue] = mapModCount["nTerm"][intModValue] + 1.0;
          if (DEBUG_PTM_TABLE2) DEBUG_MSG("Add 1 to Nterm");
        } else if (position == 0) {
          ntermAdd = 1.0;
          mapModCount["nTerm"][intModValue] += 1.0 / (floatLength + ntermAdd);
          if (DEBUG_PTM_TABLE2) DEBUG_MSG("Add "<< 1.0 / (floatLength + ntermAdd) << " to Nterm");
        }
        if (DEBUG_PTM_TABLE2) DEBUG_VAR(ntermAdd);

        float ctermAdd = 0.0;
        if (position == cleanAnnotation.length()) {
          if (intLength == 0)  {
            mapModCount["cTerm"][intModValue] += 1.0;
            if (DEBUG_PTM_TABLE2) DEBUG_MSG("Add 1 to Cterm");
          } else {
            ctermAdd = 1.0;
            mapModCount["cTerm"][intModValue] += 1.0 / (floatLength + ctermAdd);
            if (DEBUG_PTM_TABLE2) DEBUG_MSG("Add "<< 1.0 / (floatLength + ctermAdd) << " to Cterm");
          }
        }
        if (DEBUG_PTM_TABLE2) DEBUG_VAR(ctermAdd);

#if 0
        float reducedLength = floatLength;
        vector<float> masses;
        for (int k = 0; k < intLength; k++) {
          int pos2 = position + k;
          if (DEBUG_PTM_TABLE2) DEBUG_VAR(pos2);
          float aaMass = getMass(cleanAnnotation[pos2]);
          if (DEBUG_PTM_TABLE2) DEBUG_VAR(aaMass + modValue);
          if (aaMass + modValue < MINIMUM_AA_MASS) {
            reducedLength--;
            if (DEBUG_PTM_TABLE2) DEBUG_VAR(reducedLength);
          }
        }

        for (int k = 0; k < intLength; k++) {
          int pos2 = position + k;
          if (DEBUG_PTM_TABLE2) DEBUG_VAR(pos2);
          string modChar("X");
          modChar[0] = cleanAnnotation[pos2];
          if (DEBUG_PTM_TABLE) DEBUG_VAR(modChar[0]);
          float aaMass = getMass(modChar[0]);
          if (DEBUG_PTM_TABLE2) DEBUG_VAR(aaMass + modValue);
          if (reducedLength == 0) {
            if (aaMass + modValue < MINIMUM_AA_MASS && reducedLength != 0.0) {
              continue;
            }
            mapModCount[modChar][intModValue] += 1.0 / (intLength + ntermAdd + ctermAdd);
            if (DEBUG_PTM_TABLE2) DEBUG_MSG("Add "<< 1.0 / (intLength + ntermAdd + ctermAdd) << " " << modChar << " " << intModValue);
          } else {
            mapModCount[modChar][intModValue] += 1.0 / (reducedLength + ntermAdd + ctermAdd);
            if (DEBUG_PTM_TABLE2) DEBUG_MSG("Add "<< 1.0 / (reducedLength + ntermAdd + ctermAdd) << " " << modChar << " " << intModValue);
          }
        }
#endif
        for (int k = 0; k < intLength; k++) {
          int pos2 = position + k;
          if (DEBUG_PTM_TABLE2) DEBUG_VAR(pos2);
          string modChar("X");
          modChar[0] = cleanAnnotation[pos2];
          if (DEBUG_PTM_TABLE) DEBUG_VAR(modChar[0]);
          mapModCount[modChar][intModValue] += 1.0 / (floatLength + ntermAdd + ctermAdd);
          if (DEBUG_PTM_TABLE2) DEBUG_MSG("Add "<< 1.0 / (floatLength + ntermAdd + ctermAdd) << " " << modChar << " " << intModValue);
        }

      } // for (int i = 0; i < positions.size(); i++) {

    } // for (it = m_psmSet.begin(); it != m_psmSet.end(); it++)

    //DEBUG_TRACE;
    return;
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchSet::recomputeGapsByPreferredPtm(
          map<float, vector<string> > & mapPrefered,
          int nPref,
          bool useOrig)
  {
    if (DEBUG_PTM_TABLE2) DEBUG_VAR(nPref);
    vector<psmPtr>::iterator it;
    for (it = m_psmSet.begin(); it != m_psmSet.end(); it++)
    {
      // Exclude decoys
      if ((*it)->m_isDecoy) {
        continue;
      }
      vector<float> modifications;
      vector<unsigned int> positions;
      vector<unsigned int> lengths;
      (*it)->getModificationsAndPositions(modifications, positions, lengths, useOrig);

      string annotation;
      if (useOrig) {
        annotation = (*it)->m_origAnnotation;
      } else {
        annotation = (*it)->m_annotation;
      }
      if (DEBUG_PTM_TABLE2) DEBUG_VAR(annotation);

      string cleanAnnotation;
      PeptideSpectrumMatch::getUnmodifiedPeptide(annotation, cleanAnnotation);
      if (DEBUG_PTM_TABLE2) DEBUG_VAR(cleanAnnotation);

      for (int i = 0; i < positions.size(); i++) {
        if (DEBUG_PTM_TABLE2) DEBUG_VAR(positions[i]);
        int position = positions[i] - lengths[i];  // the START position
        if (DEBUG_PTM_TABLE2) DEBUG_VAR(position);
        if (DEBUG_PTM_TABLE2) DEBUG_VAR(lengths[i]);
        if (DEBUG_PTM_TABLE2) DEBUG_VAR(modifications[i]);

        float modValue = (float)modifications[i];
        float intModValue = 0.0;
        if (modValue < 0) {
          intModValue = (float)int(modValue - 0.5);
        } else  {
          intModValue = (float)int(modValue + 0.5);
        }
        if (DEBUG_PTM_TABLE2) DEBUG_VAR(intModValue);

        if (lengths[i] > 1 && mapPrefered.find(intModValue) != mapPrefered.end()) {
          vector<string> & vecPref = mapPrefered[intModValue];
          if (DEBUG_PTM_TABLE2) DEBUG_VAR(vecPref.size());
          for (int iPref = 0;
              iPref < vecPref.size() && ((nPref == -1) || iPref < nPref);
              iPref++) {
            // Lets check to see if we want to place the mod on a particular AA
            string preferedAA = vecPref[iPref];
            if (DEBUG_PTM_TABLE2) DEBUG_VAR(preferedAA);
            string gapString = cleanAnnotation.substr(position, lengths[i]);
            if (DEBUG_PTM_TABLE2) DEBUG_VAR(gapString);
            size_t pos = gapString.find(preferedAA);
            if (DEBUG_PTM_TABLE2) DEBUG_VAR(pos);
            if (pos != string::npos) {
              position += pos;
              if (DEBUG_PTM_TABLE2) DEBUG_VAR(position);
              // Adjust the position and length
              positions[i] = position + 1;
              lengths[i] = 1;
              if (DEBUG_PTM_TABLE2) DEBUG_MSG(positions[i] << " " << lengths[i]);
              if (DEBUG_PTM_TABLE2) DEBUG_VAR(pos);
              break;// We found the prefered AA
            }
          }
        }

      } // for (int i = 0; i < positions.size(); i++) {

      // Reset the PTMs according to the new information
      string completeAnnotation;
      PeptideSpectrumMatch::makeAnnotationFromData(cleanAnnotation,
          modifications,
          positions,
          lengths,
          completeAnnotation);
      if (DEBUG_PTM_TABLE2) DEBUG_VAR(completeAnnotation);
      if (useOrig) {
        (*it)->m_origAnnotation = completeAnnotation;
      } else {
        (*it)->m_annotation = completeAnnotation;
      }

    } // for (it = m_psmSet.begin(); it != m_psmSet.end(); it++)

    //DEBUG_TRACE;
    return;
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchSet::getRecursiveModCounts(
                  map<string,map<float, float> > & mapModCount,
                  PeptideSpectrumMatchSet & copySet,
                  bool useOrig,
                  vector<pair<string,float> >* pexclusionList)
  {
    copySet = *this;

    for (int iPass = 0; iPass < 20; iPass++) {

      mapModCount.clear();  // Clear old results

      map<float, vector<string> > mapPrefered;  // mod = AA
      copySet.getModCounts(mapModCount, mapPrefered, iPass, useOrig, pexclusionList);

      //DEBUG_TRACE;

      map<float, vector<pair<float, string> > > ptmTable;  // mod = count, AA

      // First re-arrange the matrix for ease of use
      map<string, map<float, float> >::iterator itrc1 = mapModCount.begin();
      map<string, map<float, float> >::iterator itrc_end1 = mapModCount.end();
      for (; itrc1 != itrc_end1; itrc1++)  {
        // Skip "nterm" and "cterm" column names
        if (itrc1->first.length() > 1) {
          continue;
        }
        map<float, float>::iterator itrc2 = itrc1->second.begin();
        map<float, float>::iterator itrc_end2 = itrc1->second.end();
        for (; itrc2 != itrc_end2; itrc2++) {
          if (DEBUG_PTM_TABLE2) DEBUG_MSG(itrc2->first << "  " <<
                                          itrc2->second << ", " <<
                                          itrc1->first);
          ptmTable[itrc2->first].push_back(pair<float, string>(itrc2->second,
                                                                    itrc1->first));
        }
      }

      // Then figure out what are the preferred AA's for each mod mass
      map<float, vector<pair<float, string> > >::iterator itrb1 = ptmTable.begin();
      map<float, vector<pair<float, string> > >::iterator itrb_end1 = ptmTable.end();
      map<float, float> totals;          // mod = total
      for (; itrb1 != itrb_end1; itrb1++)  {
        sort(itrb1->second.begin(), itrb1->second.end());
        reverse(itrb1->second.begin(), itrb1->second.end());
        for (int i = 0; i < itrb1->second.size(); i++) {
          totals[itrb1->first] += itrb1->second[i].first;
          mapPrefered[itrb1->first].push_back(itrb1->second[i].second);
          if (DEBUG_PTM_TABLE2) DEBUG_MSG(itrb1->first << "  " <<
                                          itrb1->second[i].first << ", " <<
                                          itrb1->second[i].second);
        }
      }

      // Then recompute the gaps based on the the preferences
      copySet.recomputeGapsByPreferredPtm(mapPrefered, iPass, useOrig);
    }

    // Get the mod counts one last time for the final
    mapModCount.clear();
    copySet.getModCounts(mapModCount, useOrig, pexclusionList);

    return;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::saveModMatrix(const char * filename,
                                              bool useOrig,
                                              bool recurse,
                                              vector<pair<string,float> > * pexclusionList)
  {
    ofstream ofs(filename, ios::out);
    if (!ofs)
    {
      return false;
    }

    set<float> setAllMods;

    map<string, map<float, float> > mapModCount;   // AA = mod, count
    PeptideSpectrumMatchSet returnSet;
    if (DEBUG_PTM_TABLE) DEBUG_VAR(recurse);
    if (recurse) {
      getRecursiveModCounts(mapModCount, returnSet, useOrig, pexclusionList);
    } else {
      getModCounts(mapModCount, useOrig, pexclusionList);
    }

    map<string, map<float, float> >::iterator itrc1;
    map<string, map<float, float> >::iterator itrc_end1;
    if (DEBUG_PTM_TABLE)DEBUG_VAR(recurse);
    itrc1 = mapModCount.begin();
    itrc_end1 = mapModCount.end();

    if (DEBUG_PTM_TABLE) DEBUG_TRACE;
    ofs << "Mass\t";
    for (; itrc1 != itrc_end1; itrc1++) {
      ofs << itrc1->first << "\t";
      map<float, float>::iterator itrc2 = itrc1->second.begin();
      map<float, float>::iterator itrc_end2 = itrc1->second.end();
      for (; itrc2 != itrc_end2; itrc2++) {
        setAllMods.insert(itrc2->first);
      }
    }
    ofs << "Total" << endl;

    set<float>::iterator itrs = setAllMods.begin();
    set<float>::iterator itrs_end = setAllMods.end();
    for (; itrs != itrs_end; itrs++) {
      ofs << *itrs;
      float rowTotal = 0;

      map<string, map<float, float> >::iterator itrc2;
      map<string, map<float, float> >::iterator itrc_end2;
      itrc2 = mapModCount.begin();
      itrc_end2 = mapModCount.end();
      for (; itrc2 != itrc_end2; itrc2++) {
        ofs << "\t" << mapModCount[itrc2->first][*itrs];
        rowTotal += mapModCount[itrc2->first][*itrs];
      }
      ofs << "\t" << rowTotal << endl;
    }

    if (DEBUG_PTM_TABLE) DEBUG_TRACE;
    return true;
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchSet::maximumParsimony(void)
  {
    DEBUG_VAR(m_psmSet.size());

    // Initialize assignment vector
    int maxProtIdx = 0;
    int maxSpecIdx = 0;
    // Find maximum protein and maximum spectrum index of a matched protein
    for (size_t i = 0; i < m_psmSet.size(); i++)
    {
      maxProtIdx = max(maxProtIdx, m_psmSet[i]->m_dbIndex);
      maxSpecIdx = max(maxSpecIdx, m_psmSet[i]->m_scanNum);
    }

    //DEBUG_VAR(maxProtIdx);
    //DEBUG_VAR(maxSpecIdx);

    vector<psmPtr> assignments(maxSpecIdx);
    for (int pepIdx = 0; pepIdx < assignments.size(); pepIdx++)
    {
      psmPtr p(new PeptideSpectrumMatch);
      assignments[pepIdx] = p;
      assignments[pepIdx]->m_dbIndex = -1;
    }

    // Sort all the PSMs into scan based lists
    vector<list<psmPtr> > scanPsmLists(maxSpecIdx);
    for (size_t i = 0; i < m_psmSet.size(); i++)
    {
      scanPsmLists[m_psmSet[i]->m_scanNum - 1].push_back(m_psmSet[i]);
    }

    vector<list<psmPtr> > pepsPerProt(maxProtIdx + 1);
    //DEBUG_VAR(pepsPerProt.size());

    unsigned int maxHits, maxHitsIdx;
    float maxHitsMass;
    do
    {
      maxHits = 0;
      maxHitsIdx = 0;

      // Clear the peptides per protein vector
      for (int protIdx = 0; protIdx < pepsPerProt.size(); protIdx++)
      {
        pepsPerProt[protIdx].clear();
      }

      // Create the lists of peptides per protein and compute the max
      for (int scan = 0; scan < scanPsmLists.size(); scan++)
      {
        for (list<psmPtr>::iterator iter = scanPsmLists[scan].begin();
            iter != scanPsmLists[scan].end(); iter++)
        {
          int dbIndex = (*iter)->m_dbIndex;
          pepsPerProt[dbIndex].push_back(*iter);
          if (pepsPerProt[dbIndex].size() > maxHits)
          {
            maxHits = pepsPerProt[dbIndex].size();
            maxHitsIdx = dbIndex;
          }
        }
      }
      //DEBUG_VAR(maxHits);
      //DEBUG_VAR(maxHitsIdx);

      if (maxHits > 0)
      {
        for (list<psmPtr>::iterator iter = pepsPerProt[maxHitsIdx].begin();
            iter != pepsPerProt[maxHitsIdx].end(); iter++)
        {
          scanPsmLists[(*iter)->m_scanNum - 1].clear();
          assignments[(*iter)->m_scanNum - 1] = (*iter);
        }
      }

    } while (maxHits > 0);

    m_psmSet.clear();

    // Copy the assignments back to the vector
    for (int pepIdx = 0; pepIdx < assignments.size(); pepIdx++)
    {
      if (assignments[pepIdx]->m_dbIndex != -1)
      {
        m_psmSet.push_back(assignments[pepIdx]);
      }
    }

    DEBUG_VAR(m_psmSet.size());

    return;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::parseLine(psmPtr currMatch,
                                          vector<string> & line,
                                          vector<int> & mapIndex,
                                          bool zeroIndexed,
                                          bool isInspect)
  {
    return parseLineT(currMatch, line, mapIndex, zeroIndexed, isInspect);
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::parseLine(psmPtr currMatch,
                                          vector<char *> & line,
                                          vector<int> & mapIndex,
                                          bool zeroIndexed,
                                          bool isInspect)
  {
    return parseLineT(currMatch, line, mapIndex, zeroIndexed, isInspect);
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::readHeader(const char * filename,
                                           ifstream & ifs,
                                           string fieldDelim,
                                           string commentDelim,
                                           vector<int> & mapIndex)
  {
    vector<string> header;
    vector<string> requiredHeader;
    vector<int> requiredHeaderIndex;

    if (!DelimitedTextReader::loadHeader(filename,
                                         fieldDelim,
                                         commentDelim,
                                         header,
                                         requiredHeader,
                                         requiredHeaderIndex,
                                         ifs))
    {
      ERROR_MSG("Unable to read header from " << filename);
      return false;
    }

    createMapIndex(header, mapIndex);
    DEBUG_VAR(mapIndex.size());

    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::readNextPsm(ifstream & ifs,
                                            vector<int> & mapIndex,
                                            string fieldDelim,
                                            string commentDelim,
                                            bool zeroIndexed,
                                            bool isInspect,
                                            psmPtr nextPsm)
  {
    static int count = 0;
    count++;

    vector<char *> fields(mapIndex.size());
    for (int i = 0; i < fields.size(); i++)
    {
      fields[i] = 0x0;
    }

    if (!DelimitedTextReader::getNextLine(ifs,
                                          fieldDelim,
                                          commentDelim,
                                          fields))
    {
      return false;
    }

    if (!parseLine(nextPsm, fields, mapIndex, zeroIndexed, isInspect))
    {
      WARN_MSG("Unable to parse line!");
      return false;
    }

    // Free the memory used for the field strings
    for (int i = 0; i < fields.size(); i++)
    {
      delete[] fields[i];
      fields[i] = 0x0;
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::load(const char * fileName,
                                     vector<string> & requiredHeader,
                                     string fieldDelim,
                                     string commentDelim,
                                     bool zeroIndexed,
                                     bool isInspect,
                                     int firstScan,
                                     int lastScan)
  {
    vector<int> mapIndex;
    ifstream ifsTagFile;
    if (!readHeader(fileName, ifsTagFile, fieldDelim, commentDelim, mapIndex)) {
      return false;
    }
    psmPtr nextPsm;
    bool more = true;
    while (more)
    {
      psmPtr nextPsm(new PeptideSpectrumMatch);
      more = readNextPsm(ifsTagFile,
                         mapIndex,
                         fieldDelim,
                         commentDelim,
                         zeroIndexed,
                         isInspect,
                         nextPsm);
      if (more)
      {
        if (m_doubleLoad == 2 && nextPsm->m_annotation.empty()) {
          // If we are loading MSPLIT #2 then
          // Dump any PSMs that don't have an annotation #2
          continue;
        }
        if ((firstScan == -1 || nextPsm->m_scanNum >= firstScan)
            && (lastScan == -1 || nextPsm->m_scanNum <= lastScan))
        {
          m_psmSet.push_back(nextPsm);
        }
      }
    }
    ifsTagFile.close();

    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadFromFile(const char * filename,
                                             int firstScan, /* = -1 */
                                             int lastScan /* = -1 */)
  {
    // Clear any previously existing PSMs
    m_psmSet.clear();

    vector<string> requiredHeader;
    requiredHeader.push_back("#Scan#");
    requiredHeader.push_back("SpectrumFile");
    requiredHeader.push_back("Annotation");
    requiredHeader.push_back("OrigAnnotation");
    requiredHeader.push_back("Protein");
    requiredHeader.push_back("dbIndex");
    requiredHeader.push_back("numMods");
    requiredHeader.push_back("matchOrientation");
    requiredHeader.push_back("startMass");
    requiredHeader.push_back("Charge");
    requiredHeader.push_back("MQScore");
    requiredHeader.push_back("p-value");

    if (!load(filename,
              requiredHeader,
              "\t",
              "",
              false,
              false,
              firstScan,
              lastScan))
    {
      ERROR_MSG("Unable to read file! " << filename);
      return false;
    }

    DEBUG_VAR(m_psmSet.size());
    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadSizesFromFile(const char * filename,
                                                  map<int, int> & mapPsmCounts,
                                                  int & totalPsms,
                                                  int firstScan, /* = -1 */
                                                  int lastScan /* = -1 */)
  {
    vector<string> requiredHeader;
    requiredHeader.push_back("#Scan#");
    requiredHeader.push_back("SpectrumFile");
    requiredHeader.push_back("Annotation");
    requiredHeader.push_back("OrigAnnotation");
    requiredHeader.push_back("Protein");
    requiredHeader.push_back("dbIndex");
    requiredHeader.push_back("numMods");
    requiredHeader.push_back("matchOrientation");
    requiredHeader.push_back("startMass");
    requiredHeader.push_back("Charge");
    requiredHeader.push_back("MQScore");
    requiredHeader.push_back("p-value");

    vector<int> mapIndex;
    ifstream ifsTagFile;
    if (!readHeader(filename, ifsTagFile, "\t", "", mapIndex))
    {
      ERROR_MSG("Unable to read file! " << filename);
      return false;
    }
    psmPtr nextPsm;
    bool more = true;
    while (more)
    {
      psmPtr nextPsm(new PeptideSpectrumMatch);
      more = readNextPsm(ifsTagFile, mapIndex, "\t", "", false, false, nextPsm);
      if (more)
      {
        if ((firstScan == -1 || nextPsm->m_scanNum >= firstScan)
            && (lastScan == -1 || nextPsm->m_scanNum <= lastScan))
        {

          mapPsmCounts[nextPsm->m_scanNum]++;
        }
      }
    }
    ifsTagFile.close();

    totalPsms = 0;
    //DEBUG_VAR(mapPsmCounts.size());
    map<int, int>::iterator itrMap = mapPsmCounts.begin();
    map<int, int>::iterator itrMapEnd = mapPsmCounts.end();
    for (; itrMap != itrMapEnd; itrMap++)
    {
      //DEBUG_MSG(itrMap->first << "  " << itrMap->second);
      totalPsms += itrMap->second;
    }
    //DEBUG_VAR(totalPsms);

    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadFromFiles(const char * resultsFileList)
  {
    vector<vector<string> > fileList;
    map<string, unsigned int> fileListHeader;
    vector<string> requiredHeader;
    requiredHeader.push_back("Path");
    vector<int> requiredHeaderIndex;

    if (!DelimitedTextReader::loadDelimitedFile(resultsFileList,
                                                "\t",
                                                "#",
                                                fileListHeader,
                                                fileList,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to load files!" << resultsFileList);
      return false;
    }

    for (int i = 0; i < fileList.size(); i++)
    {
      //get path
      const char * path = fileList[i][requiredHeaderIndex[0]].c_str();
      DEBUG_VAR(path);

      if (!this->loadFromFile(path))
      {
        ERROR_MSG("Unable to load file! " << fileList[i][requiredHeaderIndex[0]]);
        return false;
      }
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadInspectResultsFile(const char * resultsFile,
                                                       bool zeroIndexed)
  {
    vector<string> requiredHeader;
    requiredHeader.push_back("Scan#");
    requiredHeader.push_back("Annotation");
    requiredHeader.push_back("Charge");
    requiredHeader.push_back("Protein");

    if (!load(resultsFile, requiredHeader, "\t", "", zeroIndexed, true))
    {
      ERROR_MSG("Unable to read file! " << resultsFile);
      return false;
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadModaResultsFile(const char * resultsFile,
                                                    bool zeroIndexed)
  {
    vector<string> requiredHeader;
    requiredHeader.push_back("ScanNo");
    requiredHeader.push_back("Peptide");
    requiredHeader.push_back("Charge");
    requiredHeader.push_back("Protein");
    requiredHeader.push_back("Probability");

    if (!load(resultsFile, requiredHeader, "\t", "", zeroIndexed, true))
    {
      ERROR_MSG("Unable to read file! " << resultsFile);
      return false;
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadMsplitResultsFile(const char * resultsFile,
                                                      bool zeroIndexed)
  {
    //LARS
    vector<string> requiredHeader;
    requiredHeader.push_back("Scan#");
    requiredHeader.push_back("Annotation1");
    requiredHeader.push_back("Charge1");
    requiredHeader.push_back("Protein1");
    requiredHeader.push_back("svm1-score");

    // This is a switch to control loading results 1 or results 2
    m_doubleLoad = 1;
    if (!load(resultsFile, requiredHeader, "\t", "", zeroIndexed, true))
    {
      ERROR_MSG("Unable to read file! " << resultsFile);
      return false;
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadMsplitResultsFile2(const char * resultsFile,
                                                      bool zeroIndexed)
  {
    vector<string> requiredHeader;
    requiredHeader.push_back("Scan#");
    requiredHeader.push_back("Annotation2");
    requiredHeader.push_back("Charge2");
    requiredHeader.push_back("Protein2");
    requiredHeader.push_back("svm2-score");

    // This is a switch to control loading results 1 or results 2
    m_doubleLoad = 2;
    if (!load(resultsFile, requiredHeader, "\t", "", zeroIndexed, true))
    {
      ERROR_MSG("Unable to read file! " << resultsFile);
      return false;
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadMultipassResultsFile(const char * resultsFile,
                                                         bool zeroIndexed)
  {
    vector<string> requiredHeader;
    requiredHeader.push_back("Cluster_index");
    requiredHeader.push_back("Peptide");
    requiredHeader.push_back("Charge");
    requiredHeader.push_back("Protein");

    if (!load(resultsFile, requiredHeader, "\t", "", zeroIndexed, false))
    {
      ERROR_MSG("Unable to read file! " << resultsFile);
      return false;
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadSpecnetsReportFile(const char * resultsFile)
  {
    vector<string> requiredHeader;
    requiredHeader.push_back("ClusterIndex");
    requiredHeader.push_back("HomologSequence");
    requiredHeader.push_back("Charge");

    if (!load(resultsFile, requiredHeader, ";", "", false, false))
    {
      ERROR_MSG("Unable to read file! " << resultsFile);
      return false;
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadSpecnetsResultsFile(const char * resultsFile,
                                                        bool zeroIndexed)
  {
    vector<string> requiredHeader;
    requiredHeader.push_back("Scan#");
    requiredHeader.push_back("Annotation");
    requiredHeader.push_back("Charge");
    requiredHeader.push_back("Protein");
    if (!load(resultsFile, requiredHeader, "\t", "", zeroIndexed, false))
    {
      ERROR_MSG("Unable to read file! " << resultsFile);
      return false;
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadInspectResultsFiles(const char * resultsFileList)
  {
    vector<vector<string> > inspectFileList;
    map<string, unsigned int> inspectFileListHeader;
    vector<string> requiredHeader;
    requiredHeader.push_back("Path");
    vector<int> requiredHeaderIndex;

    if (!DelimitedTextReader::loadDelimitedFile(resultsFileList,
                                                "\t",
                                                "#",
                                                inspectFileListHeader,
                                                inspectFileList,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to load inspect files!" << resultsFileList);
      return false;
    }

    //set up column numbers for optional parameters.
    int isZeroIndexedColumn = -1;
    map<string, unsigned int>::iterator it;
    it = inspectFileListHeader.find("isZeroIndexed");
    if (it != inspectFileListHeader.end())
    {
      isZeroIndexedColumn = it->second;
    }

    for (int i = 0; i < inspectFileList.size(); i++)
    {
      //get path
      const char * path = inspectFileList[i][requiredHeaderIndex[0]].c_str();
      DEBUG_VAR(path);

      bool isZeroIndexed = false;
      //make sure results are loaded with correct indices.
      if (isZeroIndexedColumn != -1)
      {
        int intZero = 0;
        intZero =
            convertStringToInt(inspectFileList[i][isZeroIndexedColumn].c_str());
        isZeroIndexed = (bool)intZero;
      }
      if (!this->loadInspectResultsFile(path, isZeroIndexed))
      {
        ERROR_MSG("Unable to load inspect file! " << inspectFileList[i][requiredHeaderIndex[0]]);
        return false;
      }
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadMSGFPlusResultsFile(const char * resultsFile,
                                                        bool zeroIndexed)
  {
    vector<string> header;
    vector<vector<string> > lines;
    vector<string> requiredHeader;
    vector<int> requiredHeaderIndex;

    //set up MSGFDB headers
    requiredHeader.push_back("#SpecFile");
    requiredHeader.push_back("SpecID");
    requiredHeader.push_back("ScanNum");
    requiredHeader.push_back("Charge");
    requiredHeader.push_back("Protein");
    requiredHeader.push_back("Peptide");
    //requiredHeader.push_back("QValue");

    if (!DelimitedTextReader::loadDelimitedFile(resultsFile,
                                                "\t",
                                                "",
                                                header,
                                                lines,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to open results! " << resultsFile);
      return false;
    }

    int scanIndex = -1;
    int specIdxIndex = -1;

    //map headers
    vector<int> mapIndex(LAST_INDEX);
    createMapIndex(header, mapIndex);

    // Need indices for this method
    for (int i = 0; i < header.size(); i++)
    {
      if (header[i].compare("ScanNum") == 0)
      {
        scanIndex = i;
      }
      if (header[i].compare("SpecID") == 0)
      {
        specIdxIndex = i;
      }
    }

    vector<string> specIDStrs;
    const char* specIDSep = "=";

    //parse results into m_psmSet
    for (int i = 0; i < lines.size(); i++)
    {
      psmPtr currMatch(new PeptideSpectrumMatch);

      if (!parseLine(currMatch, lines[i], mapIndex, zeroIndexed, true))
      {
        WARN_MSG("Unable to parse line" << i);
      }
      else
      {

        if (currMatch->m_scanNum <= 0)
        {
          splitText(lines[i][specIdxIndex].c_str(), specIDStrs, specIDSep);

          if (specIDStrs.size() > 1)
          {
            if (!sscanf(specIDStrs.back().c_str(),
                        "%d",
                        &(currMatch->m_scanNum)))
            {
              ERROR_MSG("Unable to get scan number from spectrum index!");
              return false;
            }
          }
        }
        /* cout << currMatch->m_spectrumFile << "\t" << currMatch->m_scanNum << "\t"
         << currMatch->m_annotation << endl; */
        m_psmSet.push_back(currMatch);
      }
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadMSGFDBResultsFile(const char * resultsFile,
                                                      bool zeroIndexed)
  {
    vector<string> header;
    vector<vector<string> > lines;
    vector<string> requiredHeader;
    vector<int> requiredHeaderIndex;

    //set up MSGFDB headers
    requiredHeader.push_back("SpecIndex");
    requiredHeader.push_back("Scan#");
    requiredHeader.push_back("Peptide");
    requiredHeader.push_back("Charge");
    requiredHeader.push_back("Protein");
    requiredHeader.push_back("P-value");
    requiredHeader.push_back("SpecProb");
    //requiredHeader.push_back("FDR");

    if (!DelimitedTextReader::loadDelimitedFile(resultsFile,
                                                "\t",
                                                "",
                                                header,
                                                lines,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to open results! " << resultsFile);
      return false;
    }

    int scanIndex = -1;
    int specIdxIndex = -1;

    vector<int> mapIndex(LAST_INDEX);
    createMapIndex(header, mapIndex);

    // Need indices for this method
    for (int i = 0; i < header.size(); i++)
    {
      if (header[i].compare("Scan#") == 0)
      {
        scanIndex = i;
      }
      if (header[i].compare("SpecIndex") == 0)
      {
        specIdxIndex = i;
      }
    }

    //parse results into m_psmSet
    string scanStr("");
    vector<string> clusteredScans(0);
    string specIdxStr("");
    vector<string> clusteredIdxs(0);
    // string fragStr("");
    //vector<string> clusteredFrags(0);

    const char* clustSep = "/";

    for (int i = 0; i < lines.size(); i++)
    {
      scanStr = lines[i][scanIndex];
      specIdxStr = lines[i][specIdxIndex];
      //fragStr = lines[i][fragIndex];

      splitText(scanStr.c_str(), clusteredScans, clustSep);
      splitText(specIdxStr.c_str(), clusteredIdxs, clustSep);

      if (clusteredScans.size() != clusteredIdxs.size())
      {
        ERROR_MSG("Number of clustered scans does not match the number of clustered indices for line " << i);
        DEBUG_MSG( lines[i][scanIndex]);
        DEBUG_MSG(lines[i][specIdxIndex]);
        return false;
      }

      vector<string> lineCopy(lines[i]);
      for (int j = 0; j < clusteredScans.size(); j++)
      {
        lineCopy[scanIndex] = clusteredScans[j];
        lineCopy[specIdxIndex] = clusteredIdxs[j];

        psmPtr currMatch(new PeptideSpectrumMatch);

        if (!parseLine(currMatch, lineCopy, mapIndex, zeroIndexed, true))
        {
          WARN_MSG("Unable to parse line" << i);
        }
        else
        {
          /* cout << currMatch->m_spectrumFile << "\t" << currMatch->m_scanNum << "\t"
           << currMatch->m_annotation << endl; */
          m_psmSet.push_back(currMatch);
        }
      }
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadMSGFDBResultsFiles(const char * resultsFileList)
  {
    vector<vector<string> > msgfFileList;
    map<string, unsigned int> msgfFileListHeader;
    vector<string> requiredHeader;
    requiredHeader.push_back("Path");
    vector<int> requiredHeaderIndex;

    if (!DelimitedTextReader::loadDelimitedFile(resultsFileList,
                                                "\t",
                                                "#",
                                                msgfFileListHeader,
                                                msgfFileList,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to load msgf files!" << resultsFileList);
      return false;
    }

    //set up column numbers for optional parameters.
    int isZeroIndexedColumn = -1;

    map<string, unsigned int>::iterator it;

    it = msgfFileListHeader.find("isZeroIndexed");

    if (it != msgfFileListHeader.end())
    {
      isZeroIndexedColumn = it->second;
    }

    for (int i = 0; i < msgfFileList.size(); i++)
    {
      //get path
      const char * path = msgfFileList[i][requiredHeaderIndex[0]].c_str();
      DEBUG_VAR(path);

      bool isZeroIndexed = false;
      if (isZeroIndexedColumn != -1)
      {
        int intZero = 0;
        intZero =
            convertStringToInt(msgfFileList[i][isZeroIndexedColumn].c_str());
        isZeroIndexed = (bool)intZero;
      }
      if (!this->loadMSGFDBResultsFile(path, isZeroIndexed))
      {
        ERROR_MSG("Unable to load MSGF file! " << msgfFileList[i][requiredHeaderIndex[0]]);
        return false;
      }
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSetSpectralLibraryLoader::loadSpecnetsResultsFile(const char * resultsFile,
                                                                             bool zeroIndexed)
  {
    vector<string> requiredHeader;
    requiredHeader.push_back("Scan#");
    requiredHeader.push_back("Charge");

    if (!load(resultsFile, requiredHeader, "\t", "", zeroIndexed, false))
    {
      ERROR_MSG("Unable to read file! " << resultsFile);
      return false;
    }
    return true;
  }

  bool PeptideSpectrumMatchSet::Load(const string& resultsFile,
                                     const string& fileType)
  {
    string fileTypeLower(fileType);
    std::transform(fileTypeLower.begin(),
                   fileTypeLower.end(),
                   fileTypeLower.begin(),
                   ::tolower);

    if (fileTypeLower == "")
    {
      DEBUG_TRACE;
      if (!loadFromFile(resultsFile.c_str()))
      {
        ERROR_MSG("Failed to load PSMs from file \'" << resultsFile << "\'!");
        return false;
      }
    }
    else if (fileTypeLower == "msgfdb")
    {
      if (!loadMSGFDBResultsFile(resultsFile.c_str()))
      {
        ERROR_MSG("Failed to load MSGFDB PSMs from file \'" << resultsFile << "\'!");
        return false;
      }
    }
    else if (fileTypeLower == "inspect")
    {
      if (!loadInspectResultsFile(resultsFile.c_str()))
      {
        ERROR_MSG("Failed to load InspecT PSMs from file \'" << resultsFile << "\'!");
        return false;
      }
    }
    else if (fileTypeLower == "moda")
    {
      if (!loadModaResultsFile(resultsFile.c_str()))
      {
        ERROR_MSG("Failed to load ModA PSMs from file \'" << resultsFile << "\'!");
        return false;
      }
    }
    else if (fileTypeLower == "msplit")
    {
      if (!loadMsplitResultsFile(resultsFile.c_str()))
      {
        ERROR_MSG("Failed to load MSplit PSMs from file \'" << resultsFile << "\'!");
        return false;
      }
    }
    else if (fileTypeLower == "msgfplus")
    {
      if (!loadMSGFPlusResultsFile(resultsFile.c_str()))
      {
        ERROR_MSG("Failed to load MSGF+ PSMs from file \'" << resultsFile << "\'!");
        return false;
      }
    }
    else
    {
      ERROR_MSG("Found unsupported PSM format \'" << fileTypeLower << "\'");
      return false;
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::saveToFile(const char * filename,
                                           bool includeHeader,
                                           bool variantCol)
  {
    ofstream ofs(filename);
    if (!ofs || !ofs.good())
    {
      ERROR_MSG("Unable to open file [" << filename << "]");
      return false;
    }

    if (includeHeader)
    {
      PeptideSpectrumMatch::saveHeaderToFile(ofs, variantCol);
    }

    for (int i = 0; i < size(); i++)
    {
      m_psmSet[i]->saveToFile(ofs, variantCol);
    }
    return true;
  }

  bool PSMComparatorBySpecIndex(psmPtr p1, psmPtr p2){
    return p1->m_scanNum < p2->m_scanNum;
  }
  void PeptideSpectrumMatchSet::sortBySpecIndex(){
  	sort(m_psmSet.begin(), m_psmSet.end(), PSMComparatorBySpecIndex);
  }

} // namespace specnets
