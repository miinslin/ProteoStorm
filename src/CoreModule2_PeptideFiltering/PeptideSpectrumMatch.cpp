#include "Logger.h"
#include "PeptideSpectrumMatch.h"
#include "utils.h"

#define DEBUG_GMPFA 0
#define DEBUG_GAFMP 0
#define DEBUG_MAKE_ANNO 0
#define DEBUG_CHANGE_ANNO_TO_SINGLE 0
#define DEBUG_CHANGE_CYSTEINES 0

namespace specnets {
// -------------------------------------------------------------------------
PeptideSpectrumMatch::PeptideSpectrumMatch() {
	m_spectrumFile = "";
	m_scanNum = -1;
	m_origAnnotation = "";
	m_annotation = "";
  m_precedingAA = 0;
  m_followingAA = 0;
	m_protein = "";
	m_charge = 0;
	m_score = -1;
	m_pValue = -1;
	m_isDecoy = false;
	m_spectrum = (Spectrum *) 0x0;
	m_peakAnnotations.resize(0);
	m_ionTypes.resize(0);
	m_strict_envelope_score = 0.f;
	m_unstrict_envelope_score = 0.f;
	m_protein = "";
	m_dbIndex = -1;
	m_numMods = -1;
	m_matchOrientation = 0;
	m_startMass = -1.0;
	m_fdr = -1.0;
	m_notes = "";
	m_ionmode = "";
	m_library_name = "";
	m_mz_error_ppm = -1;
	m_mz = 0.0;
	m_shared_peaks = 0;
	m_abundance = 0;
	m_parentmass_difference = 0.f;
	m_exactmass = 0.f;
	m_library_quality = 0;
	m_spectrumID = "";
	m_pepFdr = -1.0;
	m_useYendPts = 0;
	//m_submission_user = "";
	//m_submission_id = "";
	//m_submission_date = "";
	//m_organism = "";
	//m_molecule_mass = -1.0;
	//m_compound_name = "";
  m_variantGroup = -1;
  m_peptideRegionGroup = "-1";
  m_peptideRegion = "N/A";
  m_peptideRegionLength = -1;
  m_startIndex = -1;
  m_endIndex = -1;
	
}
// -------------------------------------------------------------------------
PeptideSpectrumMatch::~PeptideSpectrumMatch() {
	//EMPTY
}
// -------------------------------------------------------------------------
PeptideSpectrumMatch::PeptideSpectrumMatch(const PeptideSpectrumMatch &other) {
	internalCopy(other);
}

// -------------------------------------------------------------------------
PeptideSpectrumMatch &PeptideSpectrumMatch::operator=(
		const PeptideSpectrumMatch &other) {
	internalCopy(other);
	return (*this);
}

// -------------------------------------------------------------------------
void PeptideSpectrumMatch::internalCopy(const PeptideSpectrumMatch &other) {
	m_spectrumFile = other.m_spectrumFile;
	m_scanNum = other.m_scanNum;
	m_annotation = other.m_annotation;
	m_origAnnotation = other.m_origAnnotation;
  m_precedingAA = other.m_precedingAA;
  m_followingAA = other.m_followingAA;
	m_useYendPts = other.m_useYendPts;
	m_peakAnnotations.resize(other.m_peakAnnotations.size());
	m_ionTypes.resize(other.m_ionTypes.size());

	map<const ftIonFragment*, const ftIonFragment*> pointerMapping; //mapping between old pointer (pointer_mapping.first) and new pointer (pointer_mapping.second

	for (unsigned int i = 0; i < other.m_peakAnnotations.size(); i++) {
		//set old annotation fragment pointer to new annotation fragment pointer.
		const ftIonFragment* oldFragPtr = other.m_peakAnnotations[i].first;
		m_peakAnnotations[i].first = pointerMapping[oldFragPtr];
		m_peakAnnotations[i].second = other.m_peakAnnotations[i].second;
	}

	m_matchedPeaks = other.m_matchedPeaks;

	for (unsigned int i = 0; i < other.m_ionTypes.size(); i++) {
		m_ionTypes[i] = other.m_ionTypes[i];
		pointerMapping[&(other.m_ionTypes[i])] = &(m_ionTypes[i]);
	}

	m_protein = other.m_protein;
	m_charge = other.m_charge;
	m_score = other.m_score;
	m_pValue = other.m_pValue;
	m_isDecoy = other.m_isDecoy;
	m_spectrum = other.m_spectrum;
	m_strict_envelope_score = other.m_strict_envelope_score;
	m_unstrict_envelope_score = other.m_unstrict_envelope_score;

	m_protein = other.m_protein;
	m_dbIndex = other.m_dbIndex;
	m_numMods = other.m_numMods;
       m_variantGroup = other.m_variantGroup;
	m_matchOrientation = other.m_matchOrientation;
	m_startMass = other.m_startMass;
	m_organism = other.m_organism;
	m_compound_name = other.m_compound_name;
	m_smiles = other.m_smiles;
	m_InChI = other.m_InChI;
	m_InChI_Aux = other.m_InChI_Aux;
	m_notes = other.m_notes;
	m_ionmode = other.m_ionmode;
	m_fdr = other.m_fdr;
	m_library_name = other.m_library_name;
	m_mz_error_ppm = other.m_mz_error_ppm;
	m_shared_peaks = other.m_shared_peaks;
	m_abundance = other.m_abundance;
	m_exactmass = other.m_exactmass;
	m_library_quality = other.m_library_quality;
	m_spectrumID = other.m_spectrumID;
	m_pepFdr = other.m_pepFdr;

	m_ion_extraction = other.m_ion_extraction;
	SLGF_distribution = other.SLGF_distribution;
	m_mz = other.m_mz;
	m_parentmass_difference = other.m_parentmass_difference;

  m_variantGroup = other.m_variantGroup;
  m_peptideRegionGroup = other.m_peptideRegionGroup;
  m_peptideRegion = other.m_peptideRegion;
  m_peptideRegionLength = other.m_peptideRegionLength;
  m_startIndex = other.m_startIndex;
  m_endIndex = other.m_endIndex;

	return;
}

// -------------------------------------------------------------------------
/*
 * Helper function for annotate
 */
bool compareProbs(const ftIonFragment &i, const ftIonFragment &j) {
	return (i.prob > j.prob);
}

// -------------------------------------------------------------------------
void PeptideSpectrumMatch::setAnnotationToMatches(vector<int> &matches,
		vector<pair<const ftIonFragment*, short> > &annotation, int ionIdx,
		const ftIonFragment* currIonFrag) {
	for (int i = 0; i < matches.size(); i++) {
		short annot_index = ionIdx + 1;
		const ftIonFragment* last_ion = annotation[matches[i]].first;
		if (last_ion && last_ion->prob >= currIonFrag->prob) {
			// if the last annotation had greater probabiliy, ignore current annotation
		} else {
			annotation[matches[i]].first = currIonFrag;
			annotation[matches[i]].second = annot_index;
		}
	}
}

// -------------------------------------------------------------------------
void PeptideSpectrumMatch::setAnnotationToMatchesNoDups(vector<int> &matches,
		vector<pair<const ftIonFragment*, short> > &annotation, int ionIdx,
		const ftIonFragment* currIonFrag) {
	float maxIntensity = 0.0;
	int maxIntensityIndex = -1;

	for (int i = 0; i < matches.size(); i++) {
		float intensity = (*m_spectrum)[matches[i]][1];
		const ftIonFragment* last_ion = annotation[matches[i]].first;
		// Update max intensity
		if (last_ion && last_ion->prob >= currIonFrag->prob) {
			// if the last annotation had greater probabiliy, ignore current annotation
		} else {
			if (intensity > maxIntensity) {
				maxIntensity = intensity;
				maxIntensityIndex = matches[i];
			}
		}
	}

	short annot_index = ionIdx + 1;
	if (maxIntensityIndex > -1) {
		annotation[maxIntensityIndex].first = currIonFrag;
		annotation[maxIntensityIndex].second = annot_index;
	}

}

// -------------------------------------------------------------------------
void PeptideSpectrumMatch::setDbMatch(const string & protein, const int dbIndex,
		const float startMass) {
	m_protein = protein;
	m_dbIndex = dbIndex;
	m_startMass = startMass;
}

// -------------------------------------------------------------------------
void PeptideSpectrumMatch::getDbMatch(string & protein, int & dbIndex,
		float & startMass) const {
	protein = m_protein;
	dbIndex = m_dbIndex;
	startMass = m_startMass;
}

float PeptideSpectrumMatch::removeSILACMod(const bool modifySpectrum) {

	int pepLen = m_annotation.length();

	int i = pepLen - 1;
	if (pepLen <= 0 || m_annotation[i] != ')') {
		return 0;
	}

	while (i >= 0 && m_annotation[i] != ',') {
		i--;
	}

	string modSubstr = m_annotation.substr(i + 1, pepLen - i - 2);

	float modMass = getFloat(modSubstr.c_str());

	if (!MZRange::EqualWithinRange(modMass, (float) 8.0, (float) 0.5)
			&& !MZRange::EqualWithinRange(modMass, (float) 10.0, (float) 0.5)) {
		return 0;
	}

	int rightResIdx = i - 1;

	while (i >= 0 && m_annotation[i] != '(') {
		i--;
	}

	int leftResIdx = i + 1;

	string residues = m_annotation.substr(leftResIdx,
			rightResIdx - leftResIdx + 1);

	m_annotation = m_annotation.substr(0, i);
	m_annotation += residues;

	if (modifySpectrum) {
		if (m_spectrum == 0) {
			ERROR_MSG("Cannot modify spectrum, pointer is zero!!");
			abort();
		}
		m_spectrum->setParentMass(m_spectrum->parentMass - modMass);
	}

	return modMass;
}

// -------------------------------------------------------------------------

/**
 * helper function for annotate
 * sets m_ionTypes based on inclusion list
 * @param ionNamesInclude -  comma delimited string indicating which ions with which names should be annotated. If string is simply "all"
 * then just include all ions in MS2Model. ex. "y,b,y++,b++"
 * @param inputIonTypes - definition of ion types: mass offsets, ion probabilities, prefix/suffix.
 * is copied into m_ionTypes
 * @param ionTypes - output vector for fragments
 */
bool copyIonTypes(const string &_ionNamesInclude,
		const vector<ftIonFragment> &inputIonTypes,
		vector<ftIonFragment> &outputIonTypes) {
	string ionNamesInclude(_ionNamesInclude);
	std::transform(ionNamesInclude.begin(), ionNamesInclude.end(),
			ionNamesInclude.begin(), ::tolower);

	if (ionNamesInclude.compare("all") == 0) { //if we're including all ion types
		outputIonTypes.resize(inputIonTypes.size());
		copy(inputIonTypes.begin(), inputIonTypes.end(),
				outputIonTypes.begin());
		return inputIonTypes.size() == outputIonTypes.size(); // check to make sure our vector sizes match
	} else {
		vector < string > ionNames;

		splitText(ionNamesInclude.c_str(), ionNames, (const char*) ",");

		for (int i = 0; i < inputIonTypes.size(); i++) {
			string lowerCaseName(inputIonTypes[i].name);
			std::transform(lowerCaseName.begin(), lowerCaseName.end(),
					lowerCaseName.begin(), ::tolower);
			for (int j = 0; j < ionNames.size(); j++) {
				if (ionNames[j].compare(lowerCaseName) == 0) { //if current ion matches include list.
					ftIonFragment ion = inputIonTypes[i];
					outputIonTypes.push_back(ion); //copy ftIonFragment object.
					break;
				}
			}
		}
		return ionNames.size() == outputIonTypes.size();
	}
}
// -------------------------------------------------------------------------
bool PeptideSpectrumMatch::annotate(const string &peptide,
		const string &ionNamesInclude, const MS2ScoringModel &inputIonTypes,
		const float prmOffset, const float srmOffset, const float peakTol,
		const bool removeDuplicates, const bool ignoreParentCharge,
		const bool retainOldAnnotations) {
	AAJumps jumps(1);
	return annotate(peptide, ionNamesInclude, inputIonTypes, prmOffset,
			srmOffset, peakTol, jumps, removeDuplicates, ignoreParentCharge,
			retainOldAnnotations);
}
// -------------------------------------------------------------------------
bool PeptideSpectrumMatch::annotate(const string &peptide,
		const string &ionNamesInclude, const MS2ScoringModel &inputIonTypes,
		const float prmOffset, const float srmOffset, const float peakTol,
		const AAJumps &jumps, const bool removeDuplicates,
		const bool ignoreParentCharge, const bool retainOldAnnotations) {
	//check to make sure spectrum is defined:
	if (m_spectrum == 0x0) {
		WARN_MSG("Spectrum not defined for PSM!");
		return false;
	}
	m_spectrum->rememberTolerances();
	m_spectrum->setPeakTolerance(peakTol);
	bool res = annotate(peptide, ionNamesInclude, inputIonTypes, prmOffset,
			srmOffset, jumps, removeDuplicates, ignoreParentCharge,
			retainOldAnnotations);
	m_spectrum->revertTolerances();

	return res;
}

// ------------------------------------------------------------------------------
bool PeptideSpectrumMatch::annotate(const string &peptide,
		const string &ionNamesInclude, const MS2ScoringModelSet &inputIonTypes,
		const float prmOffset, const float srmOffset, const AAJumps &jumps,
		const bool removeDuplicates, const bool ignoreParentCharge,
		const bool retainOldAnnotations) {
	//check to make sure spectrum is defined:
	if (m_spectrum == 0x0) {
		WARN_MSG("Spectrum not defined for PSM!");
		return false;
	}

	string activID = Spectrum::activationToString(m_spectrum->msFragType);
	if (!inputIonTypes.containsModel(activID)) {
		WARN_MSG("Failed to locate fragmentation model \'" << activID << "\'");
		return false;
	}

	// Get spectrum-specific fragmentation model
	const MS2ScoringModel &fragModel = inputIonTypes.getModel(activID);
	bool res = annotate(peptide, ionNamesInclude, fragModel, prmOffset,
			srmOffset, jumps, removeDuplicates, ignoreParentCharge,
			retainOldAnnotations);
	return res;
}
// -------------------------------------------------------------------------
bool PeptideSpectrumMatch::annotate(const string &peptide,
		const string &ionNamesInclude, const MS2ScoringModel &inputIonTypes,
		const float prmOffset, const float srmOffset, const AAJumps &jumps,
		const bool removeDuplicates, const bool ignoreParentCharge,
		const bool retainOldAnnotations) {
	//check to make sure spectrum is defined:
	if (m_spectrum == 0x0) {
		WARN_MSG("Spectrum not defined for PSM!");
		return false;
	}
	m_annotation = peptide;

	if (!retainOldAnnotations || m_ionTypes.size() == 0) {
		//clear existing annotations
		m_peakAnnotations.clear();
		m_ionTypes.clear();

		//copy ftIonFragments into PSM ionTypes
		if (!copyIonTypes(ionNamesInclude, inputIonTypes.probs, m_ionTypes)) {
			cerr
					<< "Unable to copy ionNames from MS2ScoringModel! ionNamesInclude "
					<< ionNamesInclude << endl;
			return false;
		}

		//resize m_peakAnnotations to contain annotations for each peak
		m_peakAnnotations.resize(m_spectrum->size());

		//initialize annotate vector to be null values
		for (int i = 0; i < m_peakAnnotations.size(); i++) {
			m_peakAnnotations[i].first = (ftIonFragment*) NULL;
			m_peakAnnotations[i].second = 0;
		}

		sort(m_ionTypes.begin(), m_ionTypes.end(), compareProbs); //sort fragment types in descending order so that low probability
		//annotations will be overwritten.
	}
	//generate srm and prm masses, no offsets.
	vector<float> prm_masses;
	vector<float> srm_masses;

	//generate total peptide mass (summed value of aa masses)
	float peptide_mass;

	jumps.getPRMandSRMMasses(peptide, prm_masses, srm_masses, peptide_mass);

	short peptide_length = prm_masses.size(); //length of peptide

	short ionIdx;
	int last_srm_index = -1; //keep track of last srm index so other srm searches start from there
	int last_prm_index = -1; //keep track of last prm index so other prm searches start from there

	for (ionIdx = 0; ionIdx < peptide_length - 1; ionIdx++) {

		vector<ftIonFragment>::const_iterator currIonFrag;

		for (currIonFrag = m_ionTypes.begin(); currIonFrag != m_ionTypes.end();
				currIonFrag++) {
			vector<int> matches;
			if (currIonFrag->isIF) //we're looking at internal fragment
			{
				if (currIonFrag->isNTerm) {
					float curr_mass = (prm_masses[ionIdx]
							+ currIonFrag->massOffset + prmOffset)
							/ currIonFrag->charge;

#ifdef DEBUG
					std::cout << "mass " << curr_mass << endl;
					std::cout << "name " << ionIdx << " " << currIonFrag->name << endl;
#endif

					if (m_charge > 0 && !ignoreParentCharge) {
						if (m_charge >= currIonFrag->charge)
							m_spectrum->findMatches(curr_mass, -1.0, matches,
									-1);
					} else {
						m_spectrum->findMatches(curr_mass, -1.0, matches, -1);
					}

					if (matches.size() > 0)
						last_prm_index = matches[matches.size() - 1];

				} else {
					float curr_mass = (srm_masses[ionIdx]
							+ currIonFrag->massOffset + srmOffset)
							/ currIonFrag->charge;
#ifdef DEBUG
					std::cout << "mass " << curr_mass << endl;
					std::cout << "name " << ionIdx << " " << currIonFrag->name << endl;
#endif

					if (m_charge > 0 && !ignoreParentCharge) {
						if (m_charge >= currIonFrag->charge)
							m_spectrum->findMatches(curr_mass, -1.0, matches,
									-1);
					} else {
						m_spectrum->findMatches(curr_mass, -1.0, matches, -1);
					}

					if (matches.size() > 0)
						last_srm_index = matches[matches.size() - 1];
				}
			}
			if (removeDuplicates)
				setAnnotationToMatchesNoDups(matches, m_peakAnnotations, ionIdx,
						&(*currIonFrag));
			else
				setAnnotationToMatches(matches, m_peakAnnotations, ionIdx,
						&(*currIonFrag));
		}
	}

	int last_pm_index = -1; //keep track of last pm index so other searches start from there.

	vector<ftIonFragment>::const_iterator currIonFrag;

	// looking at peptide mass shift (i.e., total peptide, total peptide - water, etc.)
	for (currIonFrag = m_ionTypes.begin(); currIonFrag != m_ionTypes.end();
			currIonFrag++) {
		vector<int> matches;

		if (!currIonFrag->isIF) { // we're looking at peptide mass shift
			float curr_mass = (peptide_mass + currIonFrag->massOffset
					+ prmOffset + srmOffset) / currIonFrag->charge;
#ifdef DEBUG
			std::cout << "mass " << curr_mass << endl;
			std::cout << "name " << currIonFrag->name << endl;
#endif

			if (m_charge > 0 && !ignoreParentCharge) {
				if (m_charge >= currIonFrag->charge)
					m_spectrum->findMatches(curr_mass, -1.0, matches, -1);
			} else {
				m_spectrum->findMatches(curr_mass, -1.0, matches, -1);
			}
#ifdef DEBUG
			cout << "matches.size() " << matches.size() << endl;
#endif
			if (matches.size() > 0) {
				last_pm_index = matches[matches.size() - 1];
			}
		}
		if (removeDuplicates)
			setAnnotationToMatchesNoDups(matches, m_peakAnnotations, ionIdx,
					&(*currIonFrag));
		else
			setAnnotationToMatches(matches, m_peakAnnotations, ionIdx,
					&(*currIonFrag));
	}
	return true;
}
// -------------------------------------------------------------------------
bool matchIon(const string &ionName, const vector<ftIonFragment> &ionTypes,
		unsigned int &outputIndex) {
	bool found = false;
	for (int i = 0; i < ionTypes.size(); i++) {
		if (ionName.compare(ionTypes[i].name) == 0) {
			outputIndex = i;
			found = true;
		}
	}
	return found;
}
// -------------------------------------------------------------------------

bool PeptideSpectrumMatch::annotate(const vector<IonMass> &ionInputMasses,
		const MS2ScoringModel &model, const float peakTol,
		const bool removeDuplicates, const bool ignoreParentCharge,
		const bool retainOldAnnotations) {
	AAJumps jumps(1);
	return annotate(ionInputMasses, model, peakTol, jumps, removeDuplicates,
			ignoreParentCharge, retainOldAnnotations);
}
// -------------------------------------------------------------------------

bool PeptideSpectrumMatch::annotate(const vector<IonMass> &ionInputMasses,
		const MS2ScoringModel &model, const float peakTol, const AAJumps &jumps,
		const bool removeDuplicates, const bool ignoreParentCharge,
		const bool retainOldAnnotations) {
	//check to make sure spectrum is defined:
	if (m_spectrum == 0x0) {
		WARN_MSG("Spectrum not defined for PSM!");
		return false;
	}

	if (!retainOldAnnotations || m_ionTypes.size() == 0) {
		//clear existing annotations
		m_peakAnnotations.clear();
		m_ionTypes.clear();

		//copy ftIonFragments into PSM ionTypes
		if (!copyIonTypes("all", model.probs, m_ionTypes)) {
			cerr << "Unable to copy ionNames from ionInputFragments! " << endl;
			return false;
		}

		//resize m_peakAnnotations to contain annotations for each peak
		m_peakAnnotations.resize(m_spectrum->size());

		//initialize annotate vector to be null values
		for (int i = 0; i < m_peakAnnotations.size(); i++) {
			m_peakAnnotations[i].first = (ftIonFragment*) NULL;
			m_peakAnnotations[i].second = 0;
		}

		sort(m_ionTypes.begin(), m_ionTypes.end(), compareProbs); //sort fragment types in descending order so that low probability
		//annotations will be overwritten.
	}

	for (short ionIdx = 0; ionIdx < ionInputMasses.size(); ionIdx++) {
		vector<int> matches;

#ifdef DEBUG
		std::cout << "mass " << ionInputMasses[ionIdx].m_mass << endl;
		std::cout << "name " << ionInputMasses[ionIdx].m_name << endl;
#endif

		unsigned int ionFragmentIndex = -1;

		if (!matchIon(ionInputMasses[ionIdx].m_name, m_ionTypes,
				ionFragmentIndex)) {
			WARN_MSG("Unable to find ion ! " << ionInputMasses[ionIdx].m_name);
			continue;
		} else {

			m_spectrum->findMatches(ionInputMasses[ionIdx].m_mass, -1.0,
					matches);

			if (removeDuplicates) {
				setAnnotationToMatchesNoDups(matches, m_peakAnnotations,
						ionInputMasses[ionIdx].m_index,
						&(m_ionTypes[ionFragmentIndex]));
			} else {
				setAnnotationToMatches(matches, m_peakAnnotations,
						ionInputMasses[ionIdx].m_index,
						&(m_ionTypes[ionFragmentIndex]));
			}
		}
	}
	return true;
}
// -------------------------------------------------------------------------
/*pair<int, float> PeptideSpectrumMatch::countAnnotatedPeaks(string& _ionNamesInclude,
 map<string, pair<int, float> >* outputIonCounts) {
 if (outputIonCounts != 0) {
 outputIonCounts->clear();
 }

 int numMatched = 0;
 float matchedIntensity = 0;

 string ionNamesInclude(_ionNamesInclude);
 std::transform(ionNamesInclude.begin(), ionNamesInclude.end(),
 ionNamesInclude.begin(), ::tolower);

 vector < string > ionNames;
 splitText(ionNamesInclude.c_str(), ionNames, (const char*)",");
 set < string > ionNamesSet;
 for (int i = 0; i < ionNames.size(); i++) {
 ionNamesSet.insert(ionNames[i]);
 }

 for (int i = 0; i < m_spectrum->size(); i++) {
 if (m_peakAnnotations[i].first == (ftIonFragment*)NULL) {
 continue;
 }
 string ionName(m_peakAnnotations[i].first->name);
 string ionNameLower(ionName);
 std::transform(ionNameLower.begin(), ionNameLower.end(),
 ionNameLower.begin(), ::tolower);

 if (ionNamesSet.count(ionNameLower) > 0) {
 numMatched++;
 matchedIntensity += (*m_spectrum)[i][1];

 if (outputIonCounts) {
 if (outputIonCounts->count(ionName) == 0) {
 pair<int, float> pairCount(1, (*m_spectrum)[i][1]);
 (*outputIonCounts)[ionName] = pairCount;
 } else {
 (*outputIonCounts)[ionName].first += 1;
 (*outputIonCounts)[ionName].second += (*m_spectrum)[i][1];
 }
 }
 }
 }

 return pair<int, float>(numMatched, matchedIntensity);
 }*/
 
// -------------------------------------------------------------------------
bool PeptideSpectrumMatch::isTryptic(bool useOrig)
{
  string anno = m_annotation;
  if (useOrig) {
    anno = m_origAnnotation;
  }
  
  //DEBUG_VAR(anno);
  // If we have leading AA lets check it
  if (anno.length() > 1 && anno[1] == '.') {
    if (anno[0] != '_' && anno[0] != '_' && anno[0] != 'K' && anno[0] != 'R') {
      //DEBUG_MSG("  Non-tryptic nterm [" << anno[0] << "]");
      return false;
    }
  }

  string annoWithoutLeadTrail;
  stripPrecedingAndFollowing(anno, annoWithoutLeadTrail);
  //DEBUG_VAR(annoWithoutLeadTrail);
  string annoClean;
  getUnmodifiedPeptide(annoWithoutLeadTrail, annoClean);
  //DEBUG_VAR(annoClean);

  // Check final AA
  if (annoClean.length() > 0) {
    char last = annoClean[annoClean.length() - 1];
    //DEBUG_VAR(last);
    if (last != 'K' && last != 'R') {
      //DEBUG_MSG("  Non-tryptic cterm [" << last << "]");
      return false;
    }
  }

#if 0
  // Check for missed internal cleavage ??
  for (int i = 0; i < annoClean.length() - 1; i++) {
    if (annoClean[i] = 'K' || annoClean[i] == 'R') {
      DEBUG_MSG("  Missed cleavage [" << annoClean[i] << "]");
      return false;
    }
  }
#endif
  
  return true;
}

// -------------------------------------------------------------------------
bool PeptideSpectrumMatch::isModified(void) const {
	size_t found = m_annotation.find_first_of("([");

	return found != string::npos;
}
// -------------------------------------------------------------------------
bool PeptideSpectrumMatch::setChargeByAnnotation(void) {
	AAJumps jumps(1);
	return setChargeByAnnotation(jumps);
}
// -------------------------------------------------------------------------
bool PeptideSpectrumMatch::setChargeByAnnotation(const AAJumps &jumps) {
	if (m_spectrum == 0x0) {
		WARN_MSG("Spectrum not defined for PSM!");
		return false;
	}

	if (m_spectrum->parentMZ == 0) {
		WARN_MSG("Parent mass to charge not defined for spectrum!");
		return false;
	}

	if (m_annotation.empty()) {
		WARN_MSG("Annotation not defined for this peptide!");
		return false;
	}

	double charge = floor(
			(jumps.getPeptideMass(m_annotation) / m_spectrum->parentMZ) + .5);

	stringstream ss;
	ss << charge;
	ss >> m_charge;
	return true;
}

// -------------------------------------------------------------------------
void PeptideSpectrumMatch::getMatchedPeaksFromAnnotation(Spectrum & dbSpec,
		AAJumps & aaJumps, bool useOriginal) {
	if (!m_spectrum || m_spectrum->size() == 0) {
		return;
	}

	vector<float> massValues;
	if (DEBUG_GMPFA)
		DEBUG_VAR(useOriginal);
	string annotation;
	if (useOriginal) {
		annotation = m_origAnnotation;
	} else {
		annotation = m_annotation;
	}
	if (annotation.empty()) {
		return;
	}
	if (DEBUG_GMPFA)
		DEBUG_VAR(annotation);
	vector < string > tokens;
	aaJumps.getPRMMasses(annotation, massValues, 0.0, &tokens, true);

	if (DEBUG_GMPFA)
		DEBUG_TRACE;
	vector<int> tokenLengths;
	for (int x = 0; x < tokens.size(); x++) {
		if (DEBUG_GMPFA)
			DEBUG_VAR(tokens[x]);
		if (tokens[x][0] == '[') {
			tokenLengths.push_back(-1);
		} else {
			string stripped = aaJumps.stripMods(tokens[x]);
			if (DEBUG_GMPFA)
				DEBUG_VAR(stripped);
			tokenLengths.push_back(stripped.length());
		}
		if (DEBUG_GMPFA)
			DEBUG_VAR(tokenLengths[tokenLengths.size() - 1]);
	}
	bool startWithGap = tokenLengths[0] == -1;
	bool endWithGap = tokenLengths[tokenLengths.size() - 1] == -1;

	if (DEBUG_GMPFA)
		DEBUG_TRACE;
	// Find DB start index
	if (DEBUG_GMPFA)
		DEBUG_VAR(m_startMass);
	int dbStartIndex = 0;
	for (int idb = 1; idb < dbSpec.size(); idb++) {
		if (fabs(dbSpec[idb][0] - m_startMass) < 1.0) {
			dbStartIndex = idb;
			break;
		}
	}
	if (DEBUG_GMPFA)
		DEBUG_VAR(dbStartIndex);

	if (DEBUG_GMPFA)
		DEBUG_TRACE;
	if (DEBUG_GMPFA)
		DEBUG_VAR(massValues[0]);
	for (int i = 1; i < massValues.size(); i++) {
		if (DEBUG_GMPFA)
			DEBUG_VAR(massValues[i]);
		float massDiff = massValues[i] - massValues[i - 1];
		if (DEBUG_GMPFA)
			DEBUG_VAR(massDiff);
	}
	if (DEBUG_GMPFA)
		DEBUG_VAR(massValues.size());

	if (DEBUG_GMPFA)
		DEBUG_TRACE;
	if (DEBUG_GMPFA)
		DEBUG_VAR((*m_spectrum)[0][0]);
	for (int i = 1; i < m_spectrum->size(); i++) {
		if (DEBUG_GMPFA)
			DEBUG_VAR((*m_spectrum)[i][0]);
		float massDiff = (*m_spectrum)[i][0] - (*m_spectrum)[i - 1][0];
		if (DEBUG_GMPFA)
			DEBUG_VAR(massDiff);
	}
	if (DEBUG_GMPFA)
		DEBUG_VAR(m_spectrum->size());

	if (DEBUG_GMPFA)
		DEBUG_TRACE;
	m_matchedPeaks.clear();
	float peakTolerance = 1.0;
	vector<int> annoPeaks;
	vector<int> specPeaks;
	float firstSpecMass = startWithGap ? 0.0 : (*m_spectrum)[0][0];
	int is = startWithGap ? 1 : 0;
	int ia = startWithGap ? 1 : 0;
	while (ia < massValues.size() && is < m_spectrum->size()) {
		float diffMass;
		do {
			if (DEBUG_GMPFA)
				DEBUG_VAR(is);
			if (DEBUG_GMPFA)
				DEBUG_VAR(ia);
			if (DEBUG_GMPFA)
				DEBUG_VAR(massValues[ia]);
			float specMass = (*m_spectrum)[is][0] - firstSpecMass;
			if (DEBUG_GMPFA)
				DEBUG_VAR(specMass);
			diffMass = fabs(massValues[ia] - specMass);
			if (DEBUG_GMPFA)
				DEBUG_VAR(diffMass);
			if (diffMass < peakTolerance) {
				break;
			}
			if (specMass < massValues[ia]) {
				is++;
			} else {
				ia++;
			}

		} while (diffMass > peakTolerance && ia < massValues.size()
				&& is < m_spectrum->size());

		if (ia < massValues.size() && is < m_spectrum->size()) {
			annoPeaks.push_back(ia);
			specPeaks.push_back(is);
			if (DEBUG_GMPFA)
				DEBUG_VAR(is);
			if (DEBUG_GMPFA)
				DEBUG_VAR(ia);
			is++;
			ia++;
		}
	}

	if (DEBUG_GMPFA)
		DEBUG_TRACE;

	int id = dbStartIndex;
	int it = startWithGap ? 1 : 0;
	if (DEBUG_GMPFA)
		DEBUG_MSG("MATCHING PEAKS ARE:  " << specPeaks[0] << "  " << id);
	m_matchedPeaks.push_back(TwoValues<int>(specPeaks[0], id));

	int lengthMod = (endWithGap ? -1 : 0);
	for (int ip = 1; ip < annoPeaks.size() + lengthMod; ip++) {
		if (DEBUG_GMPFA)
			DEBUG_VAR(ip);
		if (DEBUG_GMPFA)
			DEBUG_VAR(annoPeaks[ip]);
		for (int x = 0; x < annoPeaks[ip] - annoPeaks[ip - 1]; x++) {
			id += tokenLengths[it];
			it++;
		}
		if (DEBUG_GMPFA)
			DEBUG_MSG("MATCHING PEAKS ARE:  " << specPeaks[ip] << "  " << id);
		m_matchedPeaks.push_back(TwoValues<int>(specPeaks[ip], id));
	}

	return;
}

// -------------------------------------------------------------------------
void PeptideSpectrumMatch::getAnnotationFromMatchedPeaks(Spectrum & dbSpec,
		const string & dbSeqStr, string & annotation) {
	annotation.clear(); // Make sure string is clear since we will be appending to it
	// Annotate gap at beginning
	if (m_matchedPeaks.size() == 0) {
		WARN_MSG("No matched peaks!");
		annotation = m_annotation;
		return;
	}

	if (m_matchedPeaks[0][0] != 0 && m_spectrum) {
		char modString[32];
		sprintf(modString, "%.3f", (*m_spectrum)[m_matchedPeaks[0][0]][0]);
		annotation = "[";
		annotation += modString;
		annotation += "]";
	}

	int firstDbIndex = m_matchedPeaks[0][1];
	if (DEBUG_GAFMP)
		DEBUG_VAR(m_matchedPeaks.size());
	for (int k = 0; k < m_matchedPeaks.size() - 1; k++) {
		if (DEBUG_GAFMP)
			DEBUG_VAR(k);
		int thisDbIndex = m_matchedPeaks[k][1];
		int nextDbIndex = m_matchedPeaks[k + 1][1];
		int peakWidth = nextDbIndex - thisDbIndex;
		int thisSpecIndex = m_matchedPeaks[k][0];
		int nextSpecIndex = m_matchedPeaks[k + 1][0];
		if (DEBUG_GAFMP)
			DEBUG_VAR(thisDbIndex);
		if (DEBUG_GAFMP)
			DEBUG_VAR(nextDbIndex);
		if (DEBUG_GAFMP)
			DEBUG_VAR(peakWidth);
		if (DEBUG_GAFMP)
			DEBUG_VAR(thisSpecIndex);
		if (DEBUG_GAFMP)
			DEBUG_VAR(nextSpecIndex);

		float dbMassDiff = dbSpec[nextDbIndex][0] - dbSpec[thisDbIndex][0];
		float specMassDiff = (*m_spectrum)[nextSpecIndex][0]
				- (*m_spectrum)[thisSpecIndex][0];
		if (DEBUG_GAFMP)
			DEBUG_VAR(dbMassDiff);
		if (DEBUG_GAFMP)
			DEBUG_VAR(specMassDiff);

		float gapMassDiff = specMassDiff - dbMassDiff;

		if (fabs(gapMassDiff) < 0.5) {
			for (int a = 0; a < peakWidth; a++) {
				annotation += dbSeqStr[a + thisDbIndex];
			}
		} else {
			annotation += "(";
			for (int a = 0; a < peakWidth; a++) {
				annotation += dbSeqStr[a + thisDbIndex];
			}
			annotation += ",";
			char modString[32];
			sprintf(modString, "%.3f", gapMassDiff);
			annotation += modString;
			annotation += ")";
		}
		if (DEBUG_GAFMP)
			DEBUG_VAR(annotation);
	}
	if (DEBUG_GAFMP)
		DEBUG_VAR(annotation);

	// Annotate gap at end
	int lastMatchedPeak = m_matchedPeaks[m_matchedPeaks.size() - 1][0];
	if (m_spectrum && (lastMatchedPeak != m_spectrum->size() - 1)) {
		char modString[32];
		sprintf(modString, "%.3f",
				(*m_spectrum)[m_spectrum->size() - 1][0]
						- (*m_spectrum)[lastMatchedPeak][0]);
		annotation += "[";
		annotation += modString;
		annotation += "]";
	}
	if (DEBUG_GAFMP)
		DEBUG_VAR(annotation);

	return;
}

// -------------------------------------------------------------------------
void PeptideSpectrumMatch::getUnmodifiedPeptide(const string &inputPeptide,
		string &outputPeptide) {

#if 1
  // Simply remove everything that isn't a letter
  outputPeptide.clear();
  string letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
  for (int i = 0; i < inputPeptide.size(); i++) {
    if (letters.find(inputPeptide[i]) != string::npos) {
      outputPeptide += inputPeptide[i];
    }
  }
#else
	outputPeptide = inputPeptide;
	//erase inspect syle mods
	size_t start = outputPeptide.find_first_of("(");

	//erase nterm mods
	if (outputPeptide.substr(start + 1, 1).compare("(") == 0) {
		outputPeptide.erase(start, 1);
	}
	while (start != string::npos) {
		outputPeptide.erase(start, 1);
		size_t modStart = outputPeptide.find_first_of(",", start);
		size_t modEnd = outputPeptide.find_first_of(")", modStart);
		outputPeptide.erase(modStart, modEnd - modStart + 1);

		start = outputPeptide.find_first_of("(", modStart);
	}
	//erase specnets style mods
	start = outputPeptide.find_first_of("[");

	while (start != string::npos) {
		outputPeptide.erase(start, 1);
		size_t modStart = start;
		size_t modEnd = outputPeptide.find_first_of("]", modStart);
		outputPeptide.erase(modStart, modEnd - modStart + 1);

		start = outputPeptide.find_first_of("[", modStart);
	}
#endif

}
// -------------------------------------------------------------------------
int PeptideSpectrumMatch::getPeptideLength() {
	int len = 0;
	for(int i=0; i<m_annotation.length(); i++){
		if( isalpha(m_annotation[i]) ) len++;
	}
	return len;
}
// -------------------------------------------------------------------------
bool PeptideSpectrumMatch::generateTheoreticalMasses(const string &peptide,
		const string &ionNamesInclude, const MS2ScoringModel &inputIonTypes,
		const float prmOffset, const float srmOffset, vector<string> &ionNames,
		vector<float> &theoreticalMasses) {
	AAJumps jumps(1);
	return generateTheoreticalMasses(peptide, ionNamesInclude, inputIonTypes,
			prmOffset, srmOffset, jumps, ionNames, theoreticalMasses);

}
// -------------------------------------------------------------------------
bool PeptideSpectrumMatch::generateTheoreticalMasses(const string &peptide,
		const string &ionNamesInclude, const MS2ScoringModel &inputIonTypes,
		const float prmOffset, const float srmOffset, const AAJumps &jumps,
		vector<string> &ionNames, vector<float> &theoreticalMasses) {
	//generate srm and prm masses, no offsets.
	vector<float> prmMasses;
	vector<float> srmMasses;

	//generate total peptide mass (summed value of aa masses)
	float peptideMass;

	jumps.getPRMandSRMMasses(peptide, prmMasses, srmMasses, peptideMass);

	short peptideLength = prmMasses.size(); //length of peptide

	vector<ftIonFragment> ionTypes;
	//copy fragments from m_model.

	//copy ftIonFragments into PSM ionTypes
	if (!copyIonTypes(ionNamesInclude, inputIonTypes.probs, ionTypes)) {
		cerr << "Unable to copy ionNames from MS2ScoringModel! ionNamesInclude "
				<< ionNamesInclude << endl;
		return false;
	}
	sort(ionTypes.begin(), ionTypes.end(), compareProbs); //sort fragment types in ascending order so that low probability
	//annotations will be overwritten.

	short ionIdx;

	for (ionIdx = 0; ionIdx < peptideLength - 1; ionIdx++) {
		vector<ftIonFragment>::const_iterator currIonFrag;

		for (currIonFrag = ionTypes.begin(); currIonFrag != ionTypes.end();
				currIonFrag++) {
			if (currIonFrag->isIF) //we're looking at internal fragment
			{
				float currMass;
				if (currIonFrag->isNTerm) {
					currMass = (prmMasses[ionIdx] + currIonFrag->massOffset
							+ prmOffset) / currIonFrag->charge;
				} else {
					currMass = (srmMasses[ionIdx] + currIonFrag->massOffset
							+ srmOffset) / currIonFrag->charge;
				}
				theoreticalMasses.push_back(currMass);
				stringstream ss;
				ss << currIonFrag->name << ionIdx + 1;
				string key = ss.str();
				ionNames.push_back(key);
			}
		}
	}

	vector<ftIonFragment>::const_iterator currIonFrag;

	// looking at peptide mass shift (i.e., total peptide, total peptide - water, etc.)
	for (currIonFrag = ionTypes.begin(); currIonFrag != ionTypes.end();
			currIonFrag++) {
		if (!currIonFrag->isIF) { // we're looking at peptide mass shift
			float currMass = (peptideMass + currIonFrag->massOffset + prmOffset
					+ srmOffset) / currIonFrag->charge;
			theoreticalMasses.push_back(currMass);
			stringstream ss;
			ss << currIonFrag->name << ionIdx + 1;
			string key = ss.str();
			ionNames.push_back(key);
		}
	}
	return true;
}
// -------------------------------------------------------------------------
float totalModMass(vector<float> &modMass, vector<int> &modIndices) {
	float totalMass = 0;
	for (int i = 0; i < modIndices.size(); i++) {
		totalMass += modMass[modIndices[i]];
	}
	return totalMass;
}
// -------------------------------------------------------------------------
bool PeptideSpectrumMatch::generateTheoreticalMassesWithMods(
		const string &peptide, const string &ionNamesInclude,
		const MS2ScoringModel &inputIonTypes, const float prmOffset,
		const float srmOffset, const AAJumps &jumps,
		const vector<vector<float> > &massShifts,
		const vector<float> &fixedMods,
		const vector<unsigned int> &fixedPositions,
		vector<IonMass> &ionOutputMasses) {
	vector<float> newMods;
	for (int i = 0; i < massShifts.size(); i++) {
		for (int j = 0; j < massShifts[i].size(); j++) {
			newMods.push_back(massShifts[i][j]);
		}
	}

	ionOutputMasses.resize(0);

	//generate mod combos.
	vector < vector<int> > modCombos;

	for (int k = 1; k <= newMods.size(); k++) {
		MathUtils::combinations(newMods.size(), k, modCombos);
	}

	//generate srm and prm masses, no offsets.
	vector<float> prmMasses;
	vector<float> srmMasses;

	//generate total peptide mass (summed value of aa masses)
	float peptideMass;

	string unmodPeptide;
	PeptideSpectrumMatch::getUnmodifiedPeptide(peptide, unmodPeptide);
	DEBUG_VAR(unmodPeptide);

	jumps.getPRMandSRMMasses(unmodPeptide, prmMasses, srmMasses, peptideMass);

	short peptideLength = prmMasses.size(); //length of peptide

	//copy fragments from m_model.
	vector<ftIonFragment> ionOutputFragments;
	copyIonTypes(ionNamesInclude, inputIonTypes.probs, ionOutputFragments);

	sort(ionOutputFragments.begin(), ionOutputFragments.end(), compareProbs); //sort fragment types in ascending order so that low probability
	//annotations will be overwritten.
#ifdef DEBUG
	cout << "ionOutputFragments.size()" << ionOutputFragments.size() << endl;
#endif

	short ionIdx;

	for (ionIdx = 0; ionIdx < peptideLength - 1; ionIdx++) {
		vector<ftIonFragment>::const_iterator currIonFrag;

		for (currIonFrag = ionOutputFragments.begin();
				currIonFrag != ionOutputFragments.end(); currIonFrag++) {
			if (currIonFrag->isIF) {
				IonMass currIonMass;
				currIonMass.m_name = currIonFrag->name;
				currIonMass.m_index = (unsigned int) ionIdx;
				for (int i = 0; i < modCombos.size(); i++) {
					//check to make sure ion length is enough
					if (ionIdx + 1 >= modCombos[i].size()) {
						float currMass;

						if (currIonFrag->isNTerm) {
							float totalMass = totalModMass(newMods,
									modCombos[i]);

							//add fixed mods

							for (int position = 0;
									position < fixedPositions.size();
									position++) {
								if (fixedPositions[position] == 0) {
									totalMass += fixedMods[position];
								} else if (ionIdx + 1
										>= fixedPositions[position]) //ionIdx indexed from 0, fixedPositions indexed from 1
										{
									totalMass += fixedMods[position];
								}
							}
							currMass = (prmMasses[ionIdx]
									+ currIonFrag->massOffset + totalMass)
									/ currIonFrag->charge;
#ifdef DEBUG
							std::cout << "mass " << currMass << endl;
							std::cout << "name " << ionIdx << " " << currIonFrag->name << endl;
#endif
						} else {
							float totalMass = totalModMass(newMods,
									modCombos[i]);

#ifdef DEBUG
							cout << "totalMass" << totalMass << endl;
#endif

							//add fixed mods
							for (int position = 0;
									position < fixedPositions.size();
									position++) {
								if (fixedPositions[position]
										== peptideLength + 1) {
									totalMass += fixedMods[position];
								} else if (ionIdx + 1
										< fixedPositions[position]) //ionIdx indexed from 0, fixedPositions indexed from 1
										{
#ifdef DEBUG
									cout << "add Cterm " << endl;
									cout << "index " << ionIdx << endl;
									cout << "fixedPositions[position]" << fixedPositions[position] << endl;
#endif
									totalMass += fixedMods[position];
								}
							}
							currMass = (srmMasses[ionIdx]
									+ currIonFrag->massOffset + totalMass)
									/ currIonFrag->charge;
#ifdef DEBUG
							cout << "srmMasses[ionIdx] " << srmMasses[ionIdx] << endl;
							cout << "currIonFrag->massOffset" << currIonFrag->massOffset << endl;
							cout << "totalMass" << totalMass << endl;

							std::cout << "mass " << currMass << endl;
							std::cout << "name " << ionIdx << " " << currIonFrag->name << endl;
#endif
						}
						currIonMass.m_mass = currMass;
						ionOutputMasses.push_back(currIonMass);
					}
				}
				//include unmodified ion
				float currMass;
				if (currIonFrag->isNTerm) {
					float totalMass = 0;

					//add fixed mods
					for (int position = 0; position < fixedPositions.size();
							position++) {
						if (fixedPositions[position] == 0) {
							totalMass += fixedMods[position];
						} else if (ionIdx + 1 >= fixedPositions[position]) //ionIdx indexed from 0, fixedPositions indexed from 1
								{
							totalMass += fixedMods[position];
						}
					}
					currMass = (prmMasses[ionIdx] + currIonFrag->massOffset
							+ totalMass) / currIonFrag->charge;
#ifdef DEBUG
					std::cout << "mass " << currMass << endl;
					std::cout << "name " << ionIdx << " " << currIonFrag->name << endl;
#endif
				} else {
					float totalMass = 0;

					//add fixed mods
					for (int position = 0; position < fixedPositions.size();
							position++) {
						if (fixedPositions[position] == peptideLength + 1) {
							totalMass += fixedMods[position];
						} else if (ionIdx + 1 < fixedPositions[position]) //ionIdx indexed from 0, fixedPositions indexed from 1
								{
							totalMass += fixedMods[position];
						}
					}
					currMass = (srmMasses[ionIdx] + currIonFrag->massOffset
							+ totalMass) / currIonFrag->charge;
#ifdef DEBUG
					std::cout << "mass " << currMass << endl;
					std::cout << "name " << ionIdx << " " << currIonFrag->name << endl;
#endif
				}
				currIonMass.m_mass = currMass;
				ionOutputMasses.push_back(currIonMass);
			}
		}
	}
	sort(ionOutputMasses.begin(), ionOutputMasses.end());
	return true;
}
// -------------------------------------------------------------------------
void PeptideSpectrumMatch::getUnmodifiedPeptide(string &outputPeptide) const {
	getUnmodifiedPeptide(m_annotation, outputPeptide);
}
// -------------------------------------------------------------------------
int PeptideSpectrumMatch::countGaps(void) {
	if (m_spectrum == (Spectrum *) NULL) {
		WARN_MSG("Warning: Spectrum is not defined!");
		return -1.0;
	}

	Spectrum * spectrum = m_spectrum;

	AAJumps jumps(1);

	size_t start = m_annotation.find_first_of("[");

	int gapCount = 0;

	while (start != string::npos) {
		if (start == 1) {
			//not necessarily a gap, we need to check to make sure we're not looking at an n-term mod
			size_t end = m_annotation.find_first_of("]", start);
			float modification;
			stringstream ss;
			ss << m_annotation.substr(start + 1, end - start - 1);
			ss >> modification;
			if (modification > -100 && modification < 100) {
				gapCount++;
			}
		} else {
			gapCount++;
			start = m_annotation.find_first_of("[", start + 1);
		}
	}
	return gapCount;
}

// -------------------------------------------------------------------------
unsigned int getAACount(const string &peptide) {
	unsigned int count = 0;
	for (unsigned int i = 0; i < peptide.length(); i++) {
		string currString = peptide.substr(i, 1);
		size_t found = currString.find_first_of("ARNDCEQGHILKMFPSTWYV");
		if (found != string::npos) {
			count++;
		}
	}
	return count;
}

// -------------------------------------------------------------------------
bool PeptideSpectrumMatch::getModifications(vector<float> &modifications,
		bool useOriginal) const {
	vector<unsigned int> positions;
	vector<unsigned int> lengths;
	return getModificationsAndPositions(modifications, positions, lengths,
			useOriginal);
}

// -------------------------------------------------------------------------
bool PeptideSpectrumMatch::getModificationsAndPositions(
		vector<float> & modifications, vector<unsigned int> & positions,
		bool useOriginal) const {
	vector<unsigned int> lengths;
	return getModificationsAndPositions(modifications, positions, lengths,
			useOriginal);
}

// -------------------------------------------------------------------------
bool PeptideSpectrumMatch::getModificationsAndPositions(
		vector<float> & modifications, vector<unsigned int> & positions,
		vector<unsigned int> & lengths, bool useOriginal) const {
	string annotation;
	if (useOriginal) {
		annotation = m_origAnnotation;
	} else {
		annotation = m_annotation;
	}

	size_t start = annotation.find_first_of("([");

	while (start != string::npos) {
		size_t mod;
		size_t end;

		if (annotation.substr(start, 1).compare("(") == 0) {
			mod = annotation.find_first_of(",", start);
			end = annotation.find_first_of(")", mod);
			float modification;
			stringstream ss;
			ss << annotation.substr(mod + 1, end - mod - 1);
			ss >> modification;
			if (ss.fail()) {
				ss.clear();
				WARN_MSG(
						"Unable to convert " << annotation.substr(mod+1, end-mod-1) << " to float!");
				return false;
			} else {
				unsigned int position = getAACount(annotation.substr(0, end));
				positions.push_back(position);
				modifications.push_back(modification);
				lengths.push_back(mod - start == 0 ? 0 : mod - start - 1);
			}
		} else {
			mod = start;
			end = annotation.find_first_of("]", mod);
			float modification;
			stringstream ss;
			ss << annotation.substr(mod + 1, end - mod - 1);
			ss >> modification;
			if (ss.fail()) {
				ss.clear();
				WARN_MSG(
						"Unable to convert " << annotation.substr(mod+1, end-mod-1) << " to float!");
				return false;
			} else {
				unsigned int position = getAACount(annotation.substr(0, end));
				positions.push_back(position);
				modifications.push_back(modification);
				lengths.push_back(mod - start == 0 ? 0 : mod - start - 1);
			}
		}
		start = annotation.find_first_of("([", end);
	}
	return true;
}

// -------------------------------------------------------------------------
void PeptideSpectrumMatch::makeAnnotationFromData(string & cleanAnnotation,
		vector<float> & modifications, vector<unsigned int> & positions,
		vector<unsigned int> & lengths, string & completeAnnotation) {

	if (DEBUG_MAKE_ANNO) DEBUG_VAR(cleanAnnotation);
	char modBuf[100];
	int iMod = 0;
	int iLength = 0;

	if (modifications.size() > 0 && positions[0] == 0) {
    int roundedMod = (int)(fabs(modifications[0]) + 0.5) * (modifications[0] < 0 ? -1 : 1);
		if (DEBUG_MAKE_ANNO) DEBUG_VAR(roundedMod);
		completeAnnotation += "[";
    sprintf(modBuf, "%d", roundedMod);
		completeAnnotation += modBuf;
		completeAnnotation += "]";
		iMod++;
	}

	for (int i = 0; i < cleanAnnotation.length(); i++) {
		if (DEBUG_MAKE_ANNO) DEBUG_VAR(i);
		if (DEBUG_MAKE_ANNO) DEBUG_VAR(iMod);
		if (DEBUG_MAKE_ANNO) DEBUG_VAR(positions[iMod]);
		if (DEBUG_MAKE_ANNO) DEBUG_VAR(lengths[iMod]);
		if (iMod < modifications.size()
				&& i == positions[iMod] - lengths[iMod]) {
			if (lengths[iMod] == 0) {
        int roundedMod = (int)(fabs(modifications[iMod]) + 0.5) * (modifications[iMod] < 0 ? -1 : 1);
				completeAnnotation += "[";
        sprintf(modBuf, "%d", roundedMod);
				completeAnnotation += modBuf;
				completeAnnotation += "]";
				iMod++;
				continue;
			} else {
				completeAnnotation += "(";
			}
			iLength = 0;
		}
		if (DEBUG_MAKE_ANNO) DEBUG_VAR(completeAnnotation);
		completeAnnotation += cleanAnnotation[i];
		if (DEBUG_MAKE_ANNO) DEBUG_VAR(completeAnnotation);
		iLength++;
		if (DEBUG_MAKE_ANNO) DEBUG_VAR(iLength);
		if (iMod < modifications.size() && i == positions[iMod] - 1) {
			if (lengths[iMod] == 0) {
        int roundedMod = (int)(fabs(modifications[iMod]) + 0.5) * (modifications[iMod] < 0 ? -1 : 1);
				completeAnnotation += "[";
        sprintf(modBuf, "%d", roundedMod);
				completeAnnotation += modBuf;
				completeAnnotation += "]";
			} else {
        int roundedMod = (int)(fabs(modifications[iMod]) + 0.5) * (modifications[iMod] < 0 ? -1 : 1);
        sprintf(modBuf, "%d", roundedMod);
				completeAnnotation += ",";
				completeAnnotation += modBuf;
				completeAnnotation += ")";
			}
			iMod++;
		}
		if (DEBUG_MAKE_ANNO) DEBUG_VAR(completeAnnotation);
		if (DEBUG_MAKE_ANNO) DEBUG_VAR(iMod);
	}

	if (iMod < modifications.size()) {
    int roundedMod = (int)(fabs(modifications[iMod]) + 0.5) * (modifications[iMod] < 0 ? -1 : 1);
		if (DEBUG_MAKE_ANNO) DEBUG_VAR(roundedMod);
    sprintf(modBuf, "%d", roundedMod);
		if (lengths[iMod] == 0) {
			completeAnnotation += "[";
			completeAnnotation += modBuf;
			completeAnnotation += "]";
		} else {
			completeAnnotation += ",";
			completeAnnotation += modBuf;
			completeAnnotation += ")";
		}
	}

	if (DEBUG_MAKE_ANNO) DEBUG_VAR(completeAnnotation);

	return;
}

// -------------------------------------------------------------------------
bool PeptideSpectrumMatch::insertModifications(const string &unmodifiedPeptide,
		const vector<float> &modifications,
		const vector<unsigned int> &positions, string &outputPeptide) {
	stringstream ss;
	unsigned int peptidePosition = 1; //current peptide position, one indexed to allow for nterm mods

	if (modifications.size() != positions.size()) {
		ERROR_MSG("Unequal modifications and positions array sizes!");
		return false;
	}

	for (unsigned int i = 0; i < positions.size(); i++) {
		if (positions[i] == 0) {
			ss << "((," << modifications[i] << ")";
			continue;
		}

		if (positions[i] > unmodifiedPeptide.length() + 1) {
			ERROR_MSG(
					"Modification position outside peptide! positions " << positions[i] << " length " << unmodifiedPeptide.length()+1);
			return false;
		}

		if (positions[i] < peptidePosition) {

			ERROR_MSG("Positions vector not sorted!");
			return false;
		}

		bool positionFound = false;
		while (!positionFound) {
			// assert(peptidePosition <= unmodifiedPeptide.length() && peptidePosition > 0);
			if (peptidePosition == positions[i]) {
				ss << '(' << unmodifiedPeptide[peptidePosition - 1] << ','
						<< modifications[i] << ')';
				positionFound = true;
			} else {
				ss << unmodifiedPeptide[peptidePosition - 1];
			}
			peptidePosition++;
		}
	}
	if (peptidePosition <= unmodifiedPeptide.length()) {
		ss << unmodifiedPeptide.substr(peptidePosition - 1);
	}
	outputPeptide = ss.str();
	return true;
}
// -------------------------------------------------------------------------
bool PeptideSpectrumMatch::insertModifications(
		const vector<float> &modifications,
		const vector<unsigned int> &positions, string &outputPeptide) const {
	return insertModifications(m_annotation, modifications, positions,
			outputPeptide);
}

// -------------------------------------------------------------------------
void PeptideSpectrumMatch::changeGapAnnosToSingle(void)
{
  AAJumps jumps(1);

  if (DEBUG_CHANGE_ANNO_TO_SINGLE) DEBUG_VAR(m_annotation);
  string cleanAnnotation;
  PeptideSpectrumMatch::getUnmodifiedPeptide(m_annotation, cleanAnnotation);
  if (DEBUG_CHANGE_ANNO_TO_SINGLE) DEBUG_VAR(cleanAnnotation);

  vector<float> modifications;
  vector<unsigned int> positions;
  vector<unsigned int> lengths;
  getModificationsAndPositions(modifications, positions, lengths);

  if (DEBUG_CHANGE_ANNO_TO_SINGLE){
    for (int i = 0; i < modifications.size(); i++) {
      DEBUG_MSG(i << "  " << modifications[i] << "  " << positions[i] << "  " << lengths[i]);
    }
  }

  for (int i = 0; i < modifications.size(); i++) {
    if (lengths[i] > 1) {
      // Move the mod to the first AA
      positions[i] = positions[i] - lengths[i] + 1;
      lengths[i] = 1;
    }
  }

  if (DEBUG_CHANGE_ANNO_TO_SINGLE){
    DEBUG_TRACE;
    for (int i = 0; i < positions.size(); i++) {
      DEBUG_MSG(i << "  " << modifications[i] << "  " << positions[i] << "  " << lengths[i]);
    }
  }
  string annotationOut;
  PeptideSpectrumMatch::makeAnnotationFromData(cleanAnnotation,
                         modifications,
                         positions,
                         lengths,
                         annotationOut);

  if (DEBUG_CHANGE_ANNO_TO_SINGLE) DEBUG_VAR(annotationOut);
  m_annotation = annotationOut;

  return;
}

// -------------------------------------------------------------------------
void PeptideSpectrumMatch::addFixedCysteineMods(void)
{
  if (DEBUG_CHANGE_CYSTEINES) DEBUG_VAR(m_annotation);
  string cleanAnnotation;
  PeptideSpectrumMatch::getUnmodifiedPeptide(m_annotation, cleanAnnotation);
  if (DEBUG_CHANGE_CYSTEINES) DEBUG_VAR(cleanAnnotation);

  vector<float> modifications;
  vector<unsigned int> positions;
  vector<unsigned int> lengths;
  getModificationsAndPositions(modifications, positions, lengths);

  if (DEBUG_CHANGE_CYSTEINES){
    for (int i = 0; i < modifications.size(); i++) {
      DEBUG_MSG(i << "  " << modifications[i] << "  " << positions[i] << "  " << lengths[i]);
    }
  }

  const double NOMINAL_MASS_CYSTEINE = 103.009184477;
  double massC = NOMINAL_MASS_CYSTEINE;
  AAJumps jumps(1);
  jumps.getAAref('C', massC); // Find the mass of C according to AA masses in use
  if (DEBUG_CHANGE_CYSTEINES) DEBUG_VAR(massC);

  double cysteineDiff = massC - NOMINAL_MASS_CYSTEINE;
  if (DEBUG_CHANGE_CYSTEINES) DEBUG_VAR(cysteineDiff);

  if (abs(cysteineDiff) > 1.0) {
    // We need put fixed mods on all the cysteines
    vector<float> newModifications;
    vector<unsigned int> newPositions;
    vector<unsigned int> newLengths;

    vector<int> cysteines;
    for (int i = 0; i < cleanAnnotation.size(); i++) {
      if (cleanAnnotation[i] == 'C') {
        if (DEBUG_CHANGE_CYSTEINES) DEBUG_MSG("C at " << i);
        cysteines.push_back(i+1);
      }
    }

    if (modifications.size() > 0 && positions[0] == 0) {
      newModifications.push_back(modifications[0]);
      newPositions.push_back(positions[0]);
      newLengths.push_back(lengths[0]);
    }

    // Combine the cysteine mods with the other mods
    int currCysteine = 0;
    int otherMod = 0;
    for (int i = 1; i <= cleanAnnotation.size(); i++) {
      float mod = 0;

      for (int j = 0; j < positions.size(); j++) {
        if (DEBUG_CHANGE_CYSTEINES) DEBUG_MSG(positions[j] << "  " << i);
        if (positions[j] == i && lengths[j] != 0) {
           mod += modifications[j];
           if (DEBUG_CHANGE_CYSTEINES) DEBUG_MSG(mod << " at " << i);
        }
      }
      for (int k = 0; k < cysteines.size(); k++) {
        if (DEBUG_CHANGE_CYSTEINES) DEBUG_MSG(cysteines[k] << "  " << i);
        if (cysteines[k] == i) {
           if (DEBUG_CHANGE_CYSTEINES) DEBUG_MSG("C at " << i);
           mod += cysteineDiff;
        }
      }
      if (DEBUG_CHANGE_CYSTEINES) DEBUG_MSG(mod << " at " << i);
      if (abs(mod) > 0) {
        newModifications.push_back(mod);
        newPositions.push_back(i);
        newLengths.push_back(1);  // All mods are length 1 now
      }
    }

    if (positions.size() > 0 &&
        positions[positions.size() - 1] == cleanAnnotation.size() &&
        lengths[lengths.size() - 1] == 0) {
       newModifications.push_back(modifications[modifications.size() - 1]);
       newPositions.push_back(positions[positions.size() - 1]);
       newLengths.push_back(lengths[lengths.size() - 1]);
    }

    if (DEBUG_CHANGE_CYSTEINES){
      for (int i = 0; i < newPositions.size(); i++) {
        DEBUG_MSG(i << "  " << newModifications[i] << "  " << newPositions[i] << "  " << newLengths[i]);
      }
    }

    string annotationOut;
    PeptideSpectrumMatch::makeAnnotationFromData(cleanAnnotation,
                                                 newModifications,
                                                 newPositions,
                                                 newLengths,
                                                 annotationOut);
    m_annotation = annotationOut;
  }

  if (DEBUG_CHANGE_CYSTEINES) DEBUG_VAR(m_annotation);

  return;
}


// -------------------------------------------------------------------------
void PeptideSpectrumMatch::mapIons(vector<vector<int> > &outputMatches,
		std::tr1::unordered_map<string, int> &ionMap) const {
	outputMatches.resize(0);

	for (int i = 0; i < m_peakAnnotations.size(); i++) {
		if (m_peakAnnotations[i].first != NULL) {
			stringstream ss;
			ss << m_peakAnnotations[i].first->name << " "
					<< m_peakAnnotations[i].second;
			string key = ss.str();

			if (ionMap.find(key) != ionMap.end()) {
				outputMatches[ionMap[key]].push_back(i);
			} else {
				vector<int> newVector;
				newVector.push_back(i);
				outputMatches.push_back(newVector);
				ionMap[key] = outputMatches.size() - 1;
			}
		}
	}
}

// -------------------------------------------------------------------------
void PeptideSpectrumMatch::stripPrecedingAndFollowing(const string & original,
		string & stripped) {

	if (original.length() == 0)
		return;

	stripped = original;
	if (stripped.length() > 2) {
		if (stripped[1] == '.'
				&& string("ABCDEFGHIJKLMNOPQRSTUVWXYZ-_").find(stripped[0])
						!= string::npos) {
      m_precedingAA = stripped[0];
			stripped = stripped.substr(2);
		}
	}
	if (stripped.length() > 2) {
		if (stripped[stripped.length() - 2] == '.'
				&& string("ABCDEFGHIJKLMNOPQRSTUVWXYZ-_").find(
						stripped[stripped.length() - 1]) != string::npos) {
      m_followingAA = stripped[stripped.length() - 1];
			stripped = stripped.substr(0, stripped.length() - 2);
		}
	}

	return;
}

// -------------------------------------------------------------------------

void PeptideSpectrumMatch::inspectToSpecNets(const string &inspectOriginal,
		string &specnets) {

	//specnets = "*.";

	//check for trailing characters

	if (inspectOriginal.length() == 0)
		return;

	string inspect = inspectOriginal;
	size_t test = inspect.find_first_of(".");
	char firstChar = inspect.at(0);

	if (test == string::npos
			|| (test != 1 && !(firstChar >= '0' && firstChar <= '9'))) {
		stringstream ss;
		ss << "*." << inspect << ".*";
		inspect = ss.str();
	}

	int inspectLength = inspect.length();
	bool modification = false;
	bool startMass = false;

	const char * ptrInspect = inspect.c_str();
	char prevChar = 0x0;

	int i;
	//note: inspect_length - 3 means we ignore the trailing characters around annotation
	for (i = 2; i < inspectLength - 2; i++) {
		char currChar = ptrInspect[i];
		//handle modifications
		if ('p' == currChar) {
			//phosphorylation
			specnets += "(";
			specnets += prevChar;
			specnets += ",80)";
		}

		// look for [] mass gaps
		else if ('[' == currChar) {

			// add previous char
			if (prevChar != 0x00)
				specnets += prevChar;

			// find next ] char
			int found = inspect.find_first_of(']', i + 1);

			// check for unexistent matching ]
			if (found == std::string::npos) {
				// ignore [
				continue;
			}

			// copy as is
			int move = found - i + 1;
			specnets += inspect.substr(i, move);

			// prevchar invalid
			prevChar = 0x00;
			currChar = (i + move >= inspectLength ? 0 : inspect[i + move]);

			// advance
			i += move;
		}

		// look for starting + ou -
		else if ((i == 2) && ('+' == currChar || '-' == currChar)) {
			startMass = true;
			if ('-' == currChar)
				specnets += "[-";
			else
				specnets += "[";

		}

		// look for a mass in the sequence (modification)
		else if ('+' == currChar || '-' == currChar) {
			modification = true;

			if ('-' == currChar) {
				specnets += "(";
				specnets += prevChar;
				specnets += ",-";
			} else {
				specnets += "(";
				specnets += prevChar;
				specnets += ",";
			}

		} else if (modification) {
			if (0 == currChar) // check for null
					{
				// do nothing
			} else if ('9' >= currChar || '.' == currChar) { //this is a number (or decimal)
				specnets += currChar;
			} else {
				specnets += ')';
				modification = false;
			}
		} else if (startMass) {
			if (0 == currChar) // check for null
					{
				// do nothing
			} else if ('9' >= currChar || '.' == currChar) { //this is a number (or decimal)
				specnets += currChar;
			} else {
				specnets += ']';
				startMass = false;
			}
		}
		//handle normal AAs
		else if ('a' > currChar) //this is a capital letter
				{
			if (prevChar && 'a' > prevChar) {
				specnets += prevChar;
			}
		}
		prevChar = currChar;
	}
	if ('9' <= prevChar && '.' != prevChar) {
		specnets += prevChar;
	}

	if (modification) {
		specnets += ')';
	}

	if (startMass) {
		specnets += ')';
	}
	//cout << specnets << endl;

	//specnets += ".*";
}

// -------------------------------------------------------------------------

bool PeptideSpectrumMatch::loadFromFile(std::ifstream & ifs) {
	ifs >> m_spectrumFile;
	ifs >> m_scanNum;
	ifs >> m_annotation;
	ifs >> m_origAnnotation;
	ifs >> m_protein;
	ifs >> m_dbIndex;
	ifs >> m_numMods;
  ifs >> m_variantGroup;
	ifs >> m_matchOrientation;
	ifs >> m_startMass;
	ifs >> m_charge;
	ifs >> m_score;
	ifs >> m_pValue;
	ifs >> m_isDecoy;
	ifs >> m_strict_envelope_score;
	ifs >> m_unstrict_envelope_score;

	std::string compound_name;
	ifs >> compound_name;
	m_compound_name = compound_name;

	return true;
}

// -------------------------------------------------------------------------

bool PeptideSpectrumMatch::saveHeaderToFile(std::ofstream & ofs, bool variantCol) {
  if (variantCol) {
	ofs << "#Scan#\tSpectrumFile\tAnnotation\tOrigAnnotation\tProtein\t"
			<< "dbIndex\tStartAA\tEndAA\tnumMods\tvariantGroup\tpeptideRegionNum\tpeptideRegion\tpepRegLength\t"
			<< "matchOrientation\tstartMass\tCharge\tMQScore\t"
			<< "p-value\tisDecoy\tStrictEnvelopeScore\tUnstrictEvelopeScore\t"
			<< "CompoundName\tOrganism\tFileScanUniqueID\tFDR\tLibraryName\t"
			<< "mzErrorPPM\tLibMetaData\tSmiles\tInchi\tLibSearchSharedPeaks\t"
			<< "Abundance\tParentMassDiff\tSpecMZ\tExactMass\tLibrarySpectrumID"
			<< endl;
  } else {
	ofs << "#Scan#\tSpectrumFile\tAnnotation\tOrigAnnotation\tProtein\t"
			<< "dbIndex\tnumMods\tmatchOrientation\tstartMass\tCharge\tMQScore\t"
			<< "p-value\tisDecoy\tStrictEnvelopeScore\tUnstrictEvelopeScore\t"
			<< "CompoundName\tOrganism\tFileScanUniqueID\tFDR\tLibraryName\t"
			<< "mzErrorPPM\tLibMetaData\tSmiles\tInchi\tLibSearchSharedPeaks\t"
			<< "Abundance\tParentMassDiff\tSpecMZ\tExactMass\tLibrarySpectrumID"
			<< endl;
  }
	return true;
}

// -------------------------------------------------------------------------

bool PeptideSpectrumMatch::saveToFile(std::ofstream & ofs, bool variantCol) {
	ofs << m_scanNum << "\t";
	if (m_spectrumFile.empty()) {
		ofs << "N/A" << "\t";
	} else {
		ofs << m_spectrumFile << "\t";
	}
	if (m_annotation.empty()) {
		ofs << "N/A" << "\t";
	} else {
		ofs << m_annotation << "\t";
	}
	if (m_origAnnotation.empty()) {
		ofs << "N/A" << "\t";
	} else {
		ofs << m_origAnnotation << "\t";
	}
	if (m_protein.empty()) {
		ofs << "N/A" << "\t";
	} else {
		ofs << m_protein << "\t";
	}
	ofs << m_dbIndex << "\t";
  if (variantCol) {
    ofs << m_startIndex << "\t";
    ofs << m_endIndex << "\t";
  }
	ofs << m_numMods << "\t";

  if (variantCol) {
    ofs << m_variantGroup << "\t";
    ofs << m_peptideRegionGroup << "\t";
	  if (m_peptideRegion.empty()) {
		  ofs << "N/A" << "\t";
	  } else {
      ofs << m_peptideRegion << "\t";
	  }
	  ofs << m_peptideRegionLength << "\t";
  }

	ofs << m_matchOrientation << "\t";
	ofs << m_startMass << "\t";
	ofs << m_charge << "\t";
	ofs << m_score << "\t";
	ofs << m_pValue << "\t";
	ofs << m_isDecoy << "\t";
	ofs << m_strict_envelope_score << "\t";
	ofs << m_unstrict_envelope_score << "\t";

	if (m_compound_name.size() == 0) {
		ofs << "N/A";
	} else {
		if (m_compound_name.length() == 0) {
			ofs << "N/A";
		} else {
			ofs << string_replace_all(m_compound_name, "\t", "");
		}
	}

	ofs << "\t";

	if (m_organism.size() == 0) {
		ofs << "N/A";
	} else {
		if (m_organism.length() == 0) {
			ofs << "N/A";
		} else {
			ofs << string_replace_all(m_organism, "\t", "");
		}
	}

	ofs << "\t";

	if (m_spectrumFile.empty()) {
		ofs << "N/A";
	} else {
		ofs << m_spectrumFile << m_scanNum;
	}

	ofs << "\t";

	ofs << m_fdr << "\t";

	if (m_library_name.empty()) {
		ofs << "N/A";
	} else {
		ofs << string_replace_all(m_library_name, "\t", "");
	}
	ofs << "\t";

	ofs << m_mz_error_ppm << "\t";

	if (m_notes.empty()) {
		ofs << "N/A";
	} else {
		ofs << m_notes;
	}

	ofs << "\t";

	if (m_smiles.size() == 0) {
		ofs << "N/A";
	} else {
		if (m_smiles.length() == 0) {
			ofs << "N/A";
		} else {
			ofs << string_replace_all(m_smiles, "\t", "");
		}
	}

	ofs << "\t";

	if (m_InChI.size() == 0) {
		ofs << "N/A";
	} else {
		if (m_InChI.length() == 0) {
			ofs << "N/A";
		} else {
			ofs << string_replace_all(m_InChI, "\t", "");
		}
	}

	ofs << "\t";

	ofs << m_shared_peaks;

	ofs << "\t";

	ofs << m_abundance;

	ofs << "\t";

	ofs << m_parentmass_difference;

	ofs << "\t";

	ofs << m_mz;

	ofs << "\t";

	ofs << m_exactmass;

	ofs << "\t";

	if (m_spectrumID.size() == 0) {
		ofs << "N/A";
	} else {
		ofs << string_replace_all(m_spectrumID, "\t", "");
	}

	ofs << endl;

	return true;
}

}
