// Header Include
#include "Logger.h"
#include "SpectrumPairSet.h"
#include "SpecnetsGraph.h"
#include "projectionutils.h"
#include "alignment_scoring.h"

#include <algorithm>
#include <math.h>
#include <iostream>
#include <set>

static bool DEBUG_FILTER = false;

using namespace std;

namespace specnets
{

  //---------------------------------------------------------------------------
  bool SpectrumPairComparator (SpectrumPair sp1, SpectrumPair sp2)
  {
    float val1 = sp1.score1 + sp1.score2;
    float val2 = sp2.score1 + sp2.score2;
    return (val1 >= val2);
  }

  //---------------------------------------------------------------------------
  bool SpectrumPairComparatorByScore (SpectrumPair sp1, SpectrumPair sp2)
  {//here, the bigger score is better
    float val1 = min(sp1.score1, sp1.score2);//in a pair, representative score is min
    float val2 = min(sp2.score1, sp2.score2);

    if( val1 < val2 ) return false;
    else if( val1 > val2 ) return true;

    val1 = max(sp1.score1, sp1.score2);
    val2 = max(sp2.score1, sp2.score2);

    if( val1 < val2 ) return false;
    else if( val1 > val2 ) return true;

    if( sp1.spec1 > sp2.spec1 ) return false;
    else if( sp1.spec1 < sp2.spec1 ) return true;

    if( sp1.spec2 > sp2.spec2 ) return false;
    else if( sp1.spec2 < sp2.spec2 ) return true;

    return true;
  }

  //---------------------------------------------------------------------------
  bool SpectrumPairComparatorByIndex (SpectrumPair sp1, SpectrumPair sp2)
  {
    /*	float val1 = sp1.spec1;
	float val2 = sp2.spec1;
	return (val1 < val2);//*/
    if( sp1.spec1 > sp2.spec1 ) return false;
    else if( sp1.spec1 < sp2.spec1 ) return true;

    if( sp1.spec2 > sp2.spec2 ) return false;
    else if( sp1.spec2 < sp2.spec2 ) return true;

    return true;
  }

  //---------------------------------------------------------------------------
  SpectrumPairSet::SpectrumPairSet(void)
  {
    resize(0);
  }

  //---------------------------------------------------------------------------
  SpectrumPairSet::SpectrumPairSet(unsigned int nsize)
  {
    thePairs.resize(nsize);
  }

  //---------------------------------------------------------------------------
  void SpectrumPairSet::copy(SpectrumPairSet & that)
  {
    thePairs.resize(that.size());
    for (size_t i = 0; i < that.size(); i++) {
      thePairs[i] = that.thePairs[i];
    }
    return;
  }

  //---------------------------------------------------------------------------
  unsigned int SpectrumPairSet::size(void) const
  {
    return thePairs.size();
  }

  //---------------------------------------------------------------------------
  void SpectrumPairSet::resize(unsigned int newSize)
  {
    thePairs.resize(newSize);
  }

  //---------------------------------------------------------------------------
  SpectrumPair const & SpectrumPairSet::operator[](unsigned int index) const
  {
    return thePairs[index];
  }

  //---------------------------------------------------------------------------
  SpectrumPair & SpectrumPairSet::operator[](unsigned int index)
  {
    return thePairs[index];
  }

  //---------------------------------------------------------------------------
  void SpectrumPairSet::sort_pairs()
  {
    // using function as comp
    sort (thePairs.begin(), thePairs.end(), SpectrumPairComparator);
  }
  void SpectrumPairSet::sort_descending_by_score()
  {
    // using function as comp
    sort (thePairs.begin(), thePairs.end(), SpectrumPairComparatorByScore);
  }
  //---------------------------------------------------------------------------
  void SpectrumPairSet::sort_pairs_by_index()
  {
    // using function as comp
    sort (thePairs.begin(), thePairs.end(), SpectrumPairComparatorByIndex);
  }
  //---------------------------------------------------------------------------
  int SpectrumPairSet::loadFromBinaryFile(const std::string & filename)
  {
    FILE *fp;
    unsigned int numEntries, vCount, resIdx, vIdx;
    fp = fopen(filename.c_str(), "rb");
    if (fp == 0)
    {
      ERROR_MSG("Can not open: " << filename);
      return -1;
    }

    fread(&numEntries, sizeof(unsigned int), 1, fp); // Number of entries in the file
    resize(numEntries);
    if (numEntries == 0)
    {
      fclose(fp);
      ERROR_MSG("No entries in file: " << filename);
      return 0;
    }
    vCount = operator[](0).loadSz();
    float *data = (float *)new float[vCount];
    int read_result;
    vector<float> dataV(vCount);
    for (unsigned int resIdx = 0; resIdx < numEntries; resIdx++)
    {
      read_result = fread(data, sizeof(float), vCount, fp);
      if (read_result != vCount)
      {
        resize(0);
        ERROR_MSG("Can not read " << filename);
        return -2;
      }
      for (vIdx = 0; vIdx < 2; vIdx++)
      {
        dataV[vIdx] = data[vIdx] - 1; // Spectrum indices are 1-based
      }
      for (vIdx = 2; vIdx < vCount; vIdx++)
      {
        dataV[vIdx] = data[vIdx];
      }
      operator[](resIdx).load(dataV);
      operator[](resIdx).spec2rev = operator[](resIdx).shift2 > -0.01;
    }
    delete[] data;
    fclose(fp);
    return numEntries;
  }

  //---------------------------------------------------------------------------
  bool SpectrumPairSet::saveToBinaryFile(const std::string & filename)
  {
    FILE *fp;
    vector<float> dataV;
    float data[6];

    fp = fopen(filename.c_str(), "wb");
    if (fp == 0)
    {
      ERROR_MSG("Can not open: " << filename);
      return false;
    }
    unsigned int numEntries = thePairs.size();
    fwrite(&numEntries, sizeof(unsigned int), 1, fp); // Number of entries in the file
    for (unsigned int i = 0; i < numEntries; i++)
    {
      thePairs[i].serialize(dataV);
      for (unsigned int j = 0; j < 2; j++)
      {
        data[j] = dataV[j] + 1; // Spectrum indices are 1-based
      }
      for (unsigned int j = 2; j < dataV.size(); j++)
      {
        data[j] = dataV[j];
      }
      fwrite(data, sizeof(float), dataV.size(), fp);
    }
    fclose(fp);

    return true;
  }

  //---------------------------------------------------------------------------
  void SpectrumPairSet::push_back(const SpectrumPair & newPair)
  {
    thePairs.push_back(newPair);
  }

  //-----------------------------------------------------------------------------
  bool SpectrumPairSet::getModificationFrequencies(float resolution,
                                                   float maxDiffMass,
                                                   map<float, float> & modFreqs)
  {
    size_t fpsize = thePairs.size();
    DEBUG_VAR(fpsize);

    float totalCounts = 0.0;
    for (size_t i = 0; i < fpsize; i++)
    {
      SpectrumPair sppair = thePairs[i];
      //DEBUG_VAR(sppair.shift1);
      //DEBUG_VAR(sppair.shift2);
      float shift1 = sppair.shift1;
      float shift2 = sppair.shift2;
      // Round to the desired precision
      shift1 = fabs((float)((int)(shift1 / resolution)) * resolution);
      shift2 = fabs((float)((int)(shift2 / resolution)) * resolution);
      //DEBUG_VAR(shift1);
      //DEBUG_VAR(shift2);
      if (shift1 != 0) {
        modFreqs[shift1] += 1.0;
        totalCounts += 1.0;
      }

      if (shift2 != 0) {
        modFreqs[shift2] += 1.0;
        totalCounts += 1.0;
      }
    } // for (size_t i = 0; i < fpsize; i++)

#if 0
    map<int, int> modHist;
    // Count the number of masses that have 'count' number of modifications
    map<float, float>::iterator itr = modFreqs.begin();
    map<float, float>::iterator itrEnd = modFreqs.end();
    for (; itr != itrEnd; itr++) {
      //DEBUG_MSG("modCount: " << itr->first << "  " << itr->second);
      if (itr->first <= 100.0) {
        modHist[itr->second] += 1;
      }
    }

    map<int, int>::iterator itrH = modHist.begin();
    map<int, int>::iterator itrHEnd = modHist.end();

    //for (; itrH != itrHEnd; itrH++) {
    //  DEBUG_MSG("modHist: " << itrH->first << "  " << itrH->second);
    //}
#endif

    map<float, float>::iterator itr = modFreqs.begin();
    map<float, float>::iterator itrEnd = modFreqs.end();
    for (; itr != itrEnd; itr++) {
      //DEBUG_MSG("modCount: " << itr->first << "  " << itr->second);
      if (itr->second <= 1 || itr->first > maxDiffMass) {
        itr->second = 0;
      } else {
        //DEBUG_MSG("modFreqs: " << itr->first << "  " << itr->second);
        itr->second /= totalCounts;
        //DEBUG_MSG("modFreqs: " << itr->first << "  " << itr->second);
      }
    }
    return true;
  }

  //---------------------------------------------------------------------------
  bool SpectrumPairSet::filter_by_component_size(unsigned int max_component_size)
  {
    if (max_component_size < 1) {
      return false;
    }

    SpectrumPairSet backup_pairs;
    for(int i = 0; i < thePairs.size(); i++){
      backup_pairs.push_back(thePairs[i]);
    }

    SpecnetsGraph spec_graph(backup_pairs);

    spec_graph.filter_graph_component_size(max_component_size);

    thePairs.resize(0);

    std::vector<unsigned int> deleted_edges = spec_graph.get_pairs_deleted();

    std::set<unsigned int> deleted_edgs_set;

    for(int i = 0; i < deleted_edges.size(); i++){
      deleted_edgs_set.insert(deleted_edges[i]);
    }

    for(int i = 0; i < backup_pairs.size(); i++){
      if(deleted_edgs_set.find(i) == deleted_edgs_set.end()){
        thePairs.push_back(backup_pairs[i]);
      }
    }
    if (thePairs.size() == 0) {
      ERROR_MSG("All pairs removed!");
      return false;
    }
    backup_pairs.resize(0);
    return true;
  }

  //---------------------------------------------------------------------------
  /**
   * Filter pairs by limiting the number of unique precursor masses per component
   * SpecSet & specSet         : spectrum set to get precursor mass info of MS2 spectra
   * int max_unique_mass_number: specify max number of precursor masses
   */
  //---------------------------------------------------------------------------
  bool SpectrumPairSet::filter_by_unique_mass_number_in_component (
      SpecSet & specSet,
      unsigned int max_unique_mass_number)
  {
    static float factor_to_nominal =  0.9995;

    if (max_unique_mass_number < 1) {
      return false;
    }

    SpectrumPairSet backup_pairs;
    for (int i = 0; i < thePairs.size(); i++) {
      backup_pairs.push_back(thePairs[i]);
    }
    backup_pairs.sort_descending_by_score();
    thePairs.resize(0);

    vector<int> spec_masses;
    for (int i=0; i<specSet.size(); i++) {
      spec_masses.push_back( round(specSet[i].parentMass*factor_to_nominal) );
    }

    vector<int> refer_to_comp(specSet.size(), -1); //which comps
    vector<vector<int> > specs_in_comp; //spectra in comps
    vector<set<int> > 	 pmass_in_comp; //unique masses in comps

    int removedEdges = 0;
    for(int i = 0; i < backup_pairs.size(); i++) {

      bool added = false;

      int spec1 = backup_pairs[i].spec1;
      int spec2 = backup_pairs[i].spec2;

      int spec1_comp = refer_to_comp[spec1];
      int spec2_comp = refer_to_comp[spec2];

      // two specs newly introduced into specnets
      if (spec1_comp == -1 && spec2_comp == -1) {
        //create new component
        vector<int> comp;
        comp.push_back(spec1);
        comp.push_back(spec2);
        specs_in_comp.push_back(comp);

        set<int> vart;
        vart.insert(spec_masses[spec1]);
        vart.insert(spec_masses[spec2]);
        pmass_in_comp.push_back(vart);

        refer_to_comp[spec1] = refer_to_comp[spec2] = specs_in_comp.size()-1;
        added = true;
      } else if (spec1_comp != -1 && spec2_comp == -1) { // spec1 in specnets, new spec2
        if( pmass_in_comp[spec1_comp].size() < max_unique_mass_number ||
            pmass_in_comp[spec1_comp].find(spec_masses[spec2]) != pmass_in_comp[spec1_comp].end() ) {
          specs_in_comp[spec1_comp].push_back(spec2);
          pmass_in_comp[spec1_comp].insert(spec_masses[spec2]);
          refer_to_comp[spec2] = spec1_comp;
          added = true;
        }
      }	else if (spec1_comp == -1 && spec2_comp != -1) { // new spec1 , spec2 in specnets
        if( pmass_in_comp[spec2_comp].size() < max_unique_mass_number ||
            pmass_in_comp[spec2_comp].find(spec_masses[spec1]) != pmass_in_comp[spec2_comp].end() ) {
          specs_in_comp[spec2_comp].push_back(spec1);
          pmass_in_comp[spec2_comp].insert(spec_masses[spec1]);
          refer_to_comp[spec1] = spec2_comp;
          added = true;
        }
      } else if (spec1_comp == spec2_comp) { // already in the same component
        added = true;
      }	else { // merging two diff components
        set<int> merged;
        merged.insert(pmass_in_comp[spec1_comp].begin(), pmass_in_comp[spec1_comp].end());
        merged.insert(pmass_in_comp[spec2_comp].begin(), pmass_in_comp[spec2_comp].end());

        if (merged.size() <= max_unique_mass_number) {
          for(int k=0; k<specs_in_comp[spec2_comp].size(); k++)
            specs_in_comp[spec1_comp].push_back(specs_in_comp[spec2_comp][k]);
          pmass_in_comp[spec1_comp]= merged;

          for (int k=0; k<specs_in_comp[spec1_comp].size(); k++) {
            refer_to_comp[specs_in_comp[spec1_comp][k]] = spec1_comp;
          }
          vector<int>().swap(specs_in_comp[spec2_comp]);
          set<int>().swap(pmass_in_comp[spec2_comp]);
          added = true;
        }
      }

      if (added) {
        thePairs.push_back(backup_pairs[i]);
      } else {
        removedEdges++;
      }
    }
    sort_pairs_by_index();

    DEBUG_MSG("[filter_by_unique_mass_number_in_component] In one component, maximum# of unique precursor masses was limited to " << max_unique_mass_number );
    DEBUG_MSG("[filter_by_unique_mass_number_in_component] Of " << backup_pairs.size() << " pairs, " << removedEdges << " were filtered out. -> " << thePairs.size() << " retained pairs.");
    backup_pairs.resize(0);
    return true;
  }

  bool SpectrumPairSet::filter_by_max_spec_per_variant(
						  SpecSet & specSet,
						  unsigned int max_per_variant)
  {
    static float factor_to_nominal =  0.9995;

    if (max_per_variant < 1) {
      return false;
    }

    SpectrumPairSet backup_pairs;
    for (int i = 0; i < thePairs.size(); i++) {
      backup_pairs.push_back(thePairs[i]);
    }
    backup_pairs.sort_descending_by_score();
    thePairs.resize(0);

    vector<int> spec_masses;
    for (int i=0; i<specSet.size(); i++) {
      spec_masses.push_back( round(specSet[i].parentMass*factor_to_nominal) );
    }

    vector<map<int, int> > varSpecCount(specSet.size());

    int removedEdges = 0;
    for(int i = 0; i < backup_pairs.size(); i++) {

    	bool added = false;

    	int spec1 = backup_pairs[i].spec1;
    	int spec2 = backup_pairs[i].spec2;

    	if( varSpecCount[spec1][spec_masses[spec2]] < max_per_variant &&
    			varSpecCount[spec2][spec_masses[spec1]] < max_per_variant ){
    //	if( varSpecCount[spec1].size() < max_per_variant &&
	//			varSpecCount[spec2].size() < max_per_variant ){
    		varSpecCount[spec1][spec_masses[spec2]]++;
    		varSpecCount[spec2][spec_masses[spec1]]++;
    		added = true;
    	}

    	if (added) {
			thePairs.push_back(backup_pairs[i]);
		} else {
			removedEdges++;
		}
    }
    sort_pairs_by_index();

    DEBUG_MSG("[filter_by_max_spec_per_variant] In spectral networks, #maximum spectra of a precursor mass connected with a spectrum was limited to " << max_per_variant );
    DEBUG_MSG("[filter_by_max_spec_per_variant] Of " << backup_pairs.size() << " pairs, " << removedEdges << " were filtered out. -> " << thePairs.size() << " retained pairs.");
    backup_pairs.resize(0);
    return true;
  }

  bool SpectrumPairSet::getSpecComponentID(SpecSet & specSet,
		  	  	  	  	  	  	  	  	   vector<int> & specCompID,
		  	  	  	  	  	  	  	  	   vector<int> & compSize)
  {
    specCompID.resize(specSet.size(), -1);

    vector<vector<int> > specs_in_comp; //spectra in comps

    for(int i = 0; i < thePairs.size(); i++) {

      int spec1 = thePairs[i].spec1;
      int spec2 = thePairs[i].spec2;

      int spec1_comp = specCompID[spec1];
      int spec2_comp = specCompID[spec2];

      // two specs newly introduced into specnets
      if (spec1_comp == -1 && spec2_comp == -1) {
        //create new component
        vector<int> comp;
        comp.push_back(spec1);
        comp.push_back(spec2);
        specs_in_comp.push_back(comp);
        specCompID[spec1] = specCompID[spec2] = specs_in_comp.size()-1;

      }
      else if (spec1_comp != -1 && spec2_comp == -1) { // spec1 in specnets, new spec2

          specs_in_comp[spec1_comp].push_back(spec2);
          specCompID[spec2] = spec1_comp;
      }
      else if (spec1_comp == -1 && spec2_comp != -1) { // new spec1 , spec2 in specnets

          specs_in_comp[spec2_comp].push_back(spec1);
          specCompID[spec1] = spec2_comp;
      }
      else if (spec1_comp == spec2_comp) { // already in the same component

      }
      else { // merging two diff components

          for(int k=0; k<specs_in_comp[spec2_comp].size(); k++)
            specs_in_comp[spec1_comp].push_back(specs_in_comp[spec2_comp][k]);

          for (int k=0; k<specs_in_comp[spec1_comp].size(); k++) {
        	  specCompID[specs_in_comp[spec1_comp][k]] = spec1_comp;
          }
          vector<int>().swap(specs_in_comp[spec2_comp]);
      }
    }

    compSize.resize(specs_in_comp.size());
    for(int i=0; i<specs_in_comp.size(); i++){
    	compSize[i] = specs_in_comp[i].size();
    }

    return true;
  }

  //---------------------------------------------------------------------------
  // filter_by_fdr helpers
  //---------------------------------------------------------------------------
  static void getPRM(string annotation, vector<float> & prm, const float aa_masses[]);
  static void matchTwoPRMs(vector<float> &AA, vector<float> &BB, pair<float, int> &matchQual, float fragTol);
  static void matchPartiallyTwoPRMs(vector<float> &AA, vector<float> &BB, float shift1, float shift2, pair<float, int> &matchQual, float pmTol, float fragTol);
  static void getMatchQuality(vector<float> &AA, vector<float> &BB, float shift1, float shift2, pair<float, int> &matchQual, float pmTol, float fragTol);
  static pair<float, float> fdr_thresholding(vector<float> &target, vector<float> &decoy, float FDRCut); //first:thr, second:fdr
  static int isCrossMatched(Spectrum &specAA, float toScore, vector<float> &AA, vector<float> &BB, float modDelta, float fragTol);
  //---------------------------------------------------------------------------
  /**
   * Calculate FDR using pairs annotated by MSGF-DB IDs and Filter pairs
   * SpecSet & specSet: spectrum set to get precursor mass info of MS2 spectra
   * PeptideSpectrumMatchSet & psmSet : Confident PSM set to annotate pairs. 'm_scanNum' should be matched to indice of specSet
   * float fdrCut     : specify edge FDR
   * float maxDelta   : MAX_MOD_MASS value in param file
   * float fragTol    : TOLERANCE_PEAK value in param file
   */
  //---------------------------------------------------------------------------
  bool SpectrumPairSet::cosine_precision_recall(SpecSet & specSet,
		  	  	  	  	  	  	  	  	  	  	PeptideSpectrumMatchSet & psmSet,
		  	  	  	  	  	  	  	  	  	  	const std::string & filename,
                                           	    float maxDelta,
                                           	    float fragTol,
                                           	    bool seperateCharge,
                                           	    bool canonicalPSMForm)
  {
    static float coPRMCut = 0.6;
    static int   coLenCut = 12;
    static float coElutionWindow = 3;

    float aa_masses[]={
						71.03711, 0, 103.00919, 115.02694, 129.04259,
						147.06841, 57.02146, 137.05891, 113.08406, 0,
						128.09496, 113.08406, 131.04049, 114.04293, 0,
						97.05276, 128.05858, 156.10111, 87.03203, 101.04768,
						0, 99.06841, 186.07931, 0, 163.06333, 0};

	if( !canonicalPSMForm ){
		AAJumps jumps(1);
		for(int i=0; i<jumps.size(); i++) aa_masses[jumps.aaLetters[i]-'A'] = jumps[i];
		DEBUG_VAR(aa_masses[2]);
	}

	vector<map<int, float> > pair_score_map(specSet.size()); //pair-score map
	for(int i = 0; i < thePairs.size(); i++){ // 2 lines fixed
	  pair_score_map[thePairs[i].spec1][thePairs[i].spec2] = thePairs[i].score1;//cosine score
	}

    vector<vector<float> > prmSpecs(psmSet.size());
    psmSet.sortBySpecIndex(); //must
	for (int i = 0; i < psmSet.size(); i++) {
	  if (psmSet[i]->m_annotation.find("!") != string::npos) {
		continue;
	  }

	  string massaged_annotation = specnets_annotation_to_msgf(msp_annotation_to_specnets(psmSet[i]->m_annotation));

	  vector<float> prm;
	  getPRM(massaged_annotation, prm, aa_masses);
	  prmSpecs[i] = prm;
	}

	FILE *fp = fopen(filename.c_str(), "w");
	if (fp == 0) {
	  ERROR_MSG("Can not open: " << filename);
	  return false;
	}
	DEBUG_MSG("Annotated pair file: " << filename);
//	fprintf(fp, "spec1\tspec2\tpept1\tpept2\tcorrect\tspspair\tcosine\tmzdiff\n" );
	fprintf(fp, "overlap\tdistance\t" );
	fprintf(fp, "spec1\tpept1\tcs1\tmsgfscore1\t-log(specprob1)\tprotein1\t" );
	fprintf(fp, "spec2\tpept2\tcs2\tmsgfscore2\t-log(specprob2)\tprotein2\t" );
	fprintf(fp, "correct\tspspair\tcosine\tmzdiff\n" );

	int numTruePairs= 0;
    vector<float> targetScore;
    vector<float> decoyScore;
    pair<float, int> matchQual; //first: %common PRM Masses, second: longest consecutive ions

    for (int i = 0; i < psmSet.size(); i++) {
      int spec1 = psmSet[i]->m_scanNum-1;
      if (spec1 >= specSet.size()) {
          ERROR_MSG("Invalid scan number [" << spec1 << "] Spec size is [" << specSet.size() << "]");
        return false;
      }

      for (int k = i+1; k < psmSet.size(); k++) {
        int spec2 = psmSet[k]->m_scanNum-1;
        if (spec2 >= specSet.size()) {
          ERROR_MSG("Invalid scan number [" << spec2 << "] Spec size is [" << specSet.size() << "]");
          return false;
        }

        if ( seperateCharge && psmSet[i]->m_charge != psmSet[k]->m_charge ) continue;
        if ( fabs(specSet[spec1].parentMass - specSet[spec2].parentMass) > maxDelta ) continue;

        int paired = ( pair_score_map[spec1].find(spec2) == pair_score_map[spec1].end() )? 0 : 1;

        matchTwoPRMs(prmSpecs[i], prmSpecs[k], matchQual, fragTol); //match two ID

        int correct = 2;//ambiguous
        if (matchQual.first > 0.999) correct= 1;
		else if (matchQual.first < coPRMCut && matchQual.second < coLenCut)  correct= 0;

        if( correct == 0 && paired == 0 ) continue;

  /*      fprintf(fp, "%d\t%d\t%s\t%s\t%d\t%d\t%.4f\t%.4f\n", spec1+1, spec2+1,
        					psmSet[i]->m_origAnnotation.c_str(), psmSet[k]->m_origAnnotation.c_str(), correct, paired,
                				pair_score_map[spec1][spec2], (specSet[spec1].parentMass-specSet[spec2].parentMass) );//*/

        int smallLen = prmSpecs[i].size(), largeLen = prmSpecs[k].size();
        if( smallLen > largeLen ){
        	smallLen = prmSpecs[k].size(),
        	largeLen = prmSpecs[i].size();
        }
        fprintf(fp, "%.4f\t%d\t", (float)smallLen/largeLen, (int)round(smallLen*(1-matchQual.first)) );
        fprintf(fp, "%d\t%s\t%d\t%d\t%f\t%s\t",
        		spec1+1, psmSet[i]->m_annotation.c_str(), psmSet[i]->m_charge, (int)psmSet[i]->m_score, -log10(psmSet[i]->m_score), psmSet[i]->m_protein.c_str());
        fprintf(fp, "%d\t%s\t%d\t%d\t%f\t%s\t",
        		spec2+1, psmSet[k]->m_annotation.c_str(), psmSet[k]->m_charge, (int)psmSet[k]->m_score, -log10(psmSet[i]->m_score), psmSet[k]->m_protein.c_str());
        fprintf(fp, "%d\t%d\t%.4f\t%.4f\n", correct, paired, pair_score_map[spec1][spec2], (specSet[spec1].parentMass-specSet[spec2].parentMass));


        if( fabs(specSet[spec1].parentMZ - specSet[spec2].parentMZ) < coElutionWindow ) continue;
        //now check accuracy, removed similar masses
		if( correct == 1 ) numTruePairs++; //#theoretical pairs

		if( paired ){
			if( correct == 1 ) targetScore.push_back(pair_score_map[spec1][spec2]);
			else if( correct == 0 ) decoyScore.push_back(pair_score_map[spec1][spec2]);
		}

      }//loop spec2
    }//loop spec1
    fclose(fp);

    //get precision, recall
    int tarSize = targetScore.size();
    int decSize = decoyScore.size();
    sort( targetScore.begin(), targetScore.end() );
    sort( decoyScore.begin(), decoyScore.end() );
    DEBUG_MSG("#Detected TRUE: " << tarSize << ", #FALSE: " << decSize);

    int tar=0, dec=0;
    std::cout<<"Cosine\tPrecision\tRecall"<<std::endl;
    for(float thr = 0.f; thr < 1.f; thr+= 0.10){

    	while( tar < tarSize && targetScore[tar] < thr ){
    		tar++;
    	}
    	while( dec < decSize && decoyScore[dec] < thr ){
    		dec++;
    	}

    	int tarAlive = tarSize-tar;
    	if( tarAlive == 0 ) break;

    	std::cout<< thr << "\t" <<
    					((float)tarAlive/(tarAlive+decSize-dec)*100) << "\t" <<
    							((float)tarAlive/numTruePairs*100) << std::endl;
    }

    return true;
  }

  bool SpectrumPairSet::filter_by_edge_fdr22(SpecSet & specSet,
                                             PeptideSpectrumMatchSet & psmSet,
                                             float fdrCut,
                                             float pmTol,
                                             float fragTol,
                                             bool canonicalPSMForm)
  {
	static float coPRMCut = 0.6;//0.5;
	static int   coLenCut = 12;//9;
	static float coElutionWindow = 3;

	float aa_masses[]={
	      				71.03711, 0, 103.00919, 115.02694, 129.04259,
	      				147.06841, 57.02146, 137.05891, 113.08406, 0,
	      				128.09496, 113.08406, 131.04049, 114.04293, 0,
	      				97.05276, 128.05858, 156.10111, 87.03203, 101.04768,
	      				0, 99.06841, 186.07931, 0, 163.06333, 0};

	if( !canonicalPSMForm ){
		AAJumps jumps(1);
		for(int i=0; i<jumps.size(); i++) aa_masses[jumps.aaLetters[i]-'A'] = jumps[i];
		DEBUG_VAR(aa_masses[2]);
	}

	psmSet.sortBySpecIndex(); //must
	map<int, int > psmMap;
	map<int, vector<float> > prmSpecs;
	map<int, float > prmDtScores;
	for (int i = 0; i < psmSet.size(); i++) {
	  if (psmSet[i]->m_annotation.find("!") != string::npos) {
		continue;
	  }

	  int sIndex = psmSet[i]->m_scanNum-1;

	  vector<float> prm;
	  getPRM(psmSet[i]->m_origAnnotation, prm, aa_masses);
	  prmSpecs[sIndex] = prm;

	  float toScore = 0;
	  for(int i=0; i<specSet[sIndex].size(); i++){
		  toScore += specSet[sIndex][i][1];
	  }
	  prmDtScores[sIndex] = toScore*0.1;

	  psmMap[sIndex] = i;
	}

	vector<float> targetScore;
	vector<float> decoyScore;
	pair<float, int> matchQual; //first: %common PRM Masses, second: longest consecutive ions

	for(int i = 0; i < thePairs.size(); i++){
	  if( prmSpecs.find(thePairs[i].spec1) == prmSpecs.end() || prmSpecs.find(thePairs[i].spec2) == prmSpecs.end() )
		  continue;//both spectra have to be identified

	  int spec1 = thePairs[i].spec1, spec2 = thePairs[i].spec2;
	  if( fabs(specSet[spec1].parentMZ - specSet[spec2].parentMZ) < coElutionWindow ) continue;

	  float alignGFScore = min(thePairs[i].score1, thePairs[i].score2);
	  getMatchQuality(prmSpecs[spec1], prmSpecs[spec2], thePairs[i].shift1, thePairs[i].shift2, matchQual, pmTol*2, fragTol);

	  if (DEBUG_FILTER) DEBUG_MSG("    " << matchQual.first << "  " << matchQual.second << "  " );
	  if (DEBUG_FILTER) DEBUG_MSG("    " << coPRMCut << "  " << coLenCut << "  " );

	  if(matchQual.first > 0.999) {
		  if (DEBUG_FILTER) DEBUG_MSG("    TARGET");
		  targetScore.push_back( alignGFScore );
	  } else if ( matchQual.first < coPRMCut && matchQual.second < coLenCut ) {

		  if( isCrossMatched(specSet[spec1], prmDtScores[spec1], prmSpecs[spec1], prmSpecs[spec2], thePairs[i].shift1, fragTol) || //match spec2 PSM to spec1 spectrum
				  	  isCrossMatched(specSet[spec2], prmDtScores[spec2], prmSpecs[spec2], prmSpecs[spec1], -thePairs[i].shift1, fragTol) ) { //match spec1 PSM to spec2 spectrum
		  }
		  else {
			  if (DEBUG_FILTER) DEBUG_MSG("    DECOY");
			  decoyScore.push_back( alignGFScore );
		  }
	  } else {
	  }
	}

	pair<float, float> thr_fdr = fdr_thresholding(targetScore, decoyScore, fdrCut);
	float threshold = thr_fdr.first;
	m_fdr = thr_fdr.second;

	long org_num = thePairs.size(), cutIndex = 0;
	sort_descending_by_score();

	for (int i = 0; i < thePairs.size(); i++) {
	  if (DEBUG_FILTER) DEBUG_VAR(min(thePairs[i].score1, thePairs[i].score2));
	  if (min(thePairs[i].score1, thePairs[i].score2) < threshold) {
		  thePairs.erase(thePairs.begin() + i, thePairs.end());
		  break;
	  }
	}

	sort_pairs_by_index();

	if (DEBUG_FILTER) DEBUG_VAR(targetScore.size());
	if (DEBUG_FILTER) DEBUG_VAR(decoyScore.size());
	if (DEBUG_FILTER) DEBUG_VAR(fdrCut);
	if (DEBUG_FILTER) DEBUG_VAR(threshold);

	DEBUG_MSG("Spectral pairs were evaluated based on " << psmSet.size() << " input PSMs.");
	DEBUG_MSG("The pairs were filtered out at edge FDR " << fdrCut);
	DEBUG_MSG("Align-GF p-value threshold was fixed to " << threshold << " (-log value).");
	DEBUG_MSG("Of " << org_num << " pairs, " << (org_num-thePairs.size())
			  << " were filtered out. Retained " << thePairs.size() << " pairs.");
	return true;
  }

  bool SpectrumPairSet::filter_by_edge_fdr33(SpecSet & specSet,
		  	  	  	  	  	  	  	  	  	 PeptideSpectrumMatchSet & psmSet,
		  	  	  	  	  	  	  	  	  	 float fdrCut,
		  	  	  	  	  	  	  	  	  	 float pmTol,
										   	 float fragTol,
										   	 const std::string & filename,
										   	 bool canonicalPSMForm)
  {
	static float coPRMCut = 0.6;//0.5;
	static int   coLenCut = 12;//9;
	static float coElutionWindow = 3;

	float aa_masses[]={
	      				71.03711, 0, 103.00919, 115.02694, 129.04259,
	      				147.06841, 57.02146, 137.05891, 113.08406, 0,
	      				128.09496, 113.08406, 131.04049, 114.04293, 0,
	      				97.05276, 128.05858, 156.10111, 87.03203, 101.04768,
	      				0, 99.06841, 186.07931, 0, 163.06333, 0};

	if( !canonicalPSMForm ){
		AAJumps jumps(1);
		for(int i=0; i<jumps.size(); i++) aa_masses[jumps.aaLetters[i]-'A'] = jumps[i];
		DEBUG_VAR(aa_masses[2]);
	}

	psmSet.sortBySpecIndex(); //must
	map<int, int > psmMap;
	map<int, vector<float> > prmSpecs;
	map<int, float > prmDtScores;

	for (int i = 0; i < psmSet.size(); i++) {
	  if (psmSet[i]->m_annotation.find("!") != string::npos) {
		continue;
	  }

	  int sIndex = psmSet[i]->m_scanNum-1;

	  vector<float> prm;
	  getPRM(psmSet[i]->m_origAnnotation, prm, aa_masses);
	  prmSpecs[sIndex] = prm;

	  float toScore = 0;
	  for(int i=0; i<specSet[sIndex].size(); i++){
		  toScore += specSet[sIndex][i][1];
	  }
	  prmDtScores[sIndex] = toScore*0.1;

	  psmMap[sIndex] = i;
	}

	FILE *fp = fopen(filename.c_str(), "w");
	if (fp == 0) {
	  ERROR_MSG("Can not open: " << filename);
	  return false;
	}
	DEBUG_MSG("Annotated pair file: " << filename);

	fprintf(fp, "spec1\tspec2\tscore1\tscore2\tshift1\tshift2\t");
	fprintf(fp, "pept1\tpept2\tcs1\tcs2\tsprob1\tsprob2\tcommonPRM\tconsecutiveAAs\talignGFScore\tannotation\n");

	vector<float> targetScore;
	vector<float> decoyScore;
	pair<float, int> matchQual; //first: %common PRM Masses, second: longest consecutive ions
	vector<int> xCount(specSet.size());

	for(int i = 0; i < thePairs.size(); i++){
	  if( prmSpecs.find(thePairs[i].spec1) == prmSpecs.end() || prmSpecs.find(thePairs[i].spec2) == prmSpecs.end() )
		  continue;//both spectra have to be identified

	  int spec1 = thePairs[i].spec1, spec2 = thePairs[i].spec2;
	  if( fabs(specSet[spec1].parentMZ - specSet[spec2].parentMZ) < coElutionWindow ) continue;

	  float alignGFScore = min(thePairs[i].score1, thePairs[i].score2);
	  fprintf(fp, "%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t",
	  				spec1+1, spec2+1, thePairs[i].score1, thePairs[i].score2, thePairs[i].shift1, thePairs[i].shift2);

	  getMatchQuality(prmSpecs[spec1], prmSpecs[spec2], thePairs[i].shift1, thePairs[i].shift2, matchQual, pmTol*2, fragTol);

	  fprintf(fp, "%s\t%s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%d\t%.4f\t",
			  	  	  	  psmSet[psmMap[spec1]]->m_origAnnotation.c_str(), psmSet[psmMap[spec2]]->m_origAnnotation.c_str(),
			  	  	  	  	  specSet[spec1].parentCharge, specSet[spec2].parentCharge,
			  	  	  	  	  -log10(psmSet[psmMap[spec1]]->m_score), -log10(psmSet[psmMap[spec2]]->m_score),
							  	  matchQual.first, matchQual.second, alignGFScore);

	  if (DEBUG_FILTER) DEBUG_MSG("    " << matchQual.first << "  " << matchQual.second << "  " );
	  if (DEBUG_FILTER) DEBUG_MSG("    " << coPRMCut << "  " << coLenCut << "  " );

	  if(matchQual.first > 0.999) {
		  if (DEBUG_FILTER) DEBUG_MSG("    TARGET");
		  targetScore.push_back( alignGFScore );
		  fprintf(fp, "T");
	  }
	  else if ( matchQual.first < coPRMCut && matchQual.second < coLenCut ) {
		  int id2_to_spec1= isCrossMatched(specSet[spec1], prmDtScores[spec1],
				  	  	  prmSpecs[spec1], prmSpecs[spec2], thePairs[i].shift1, fragTol);//match spec2 PSM to spec1 spectrum
		  int id1_to_spec2= isCrossMatched(specSet[spec2], prmDtScores[spec2],
				  	  	  prmSpecs[spec2], prmSpecs[spec1], -thePairs[i].shift1, fragTol);//match spec1 PSM to spec2 spectrum

		  if( id2_to_spec1 || id1_to_spec2 ) {
	//	  if( id2_to_spec1 == 2 || id1_to_spec2 == 2 ) {
			  if( id2_to_spec1 == 2 ) xCount[spec1]++;//probably wrong id1
			  if( id1_to_spec2 == 2 ) xCount[spec2]++;
			  fprintf(fp, "Z%d%d", id2_to_spec1, id1_to_spec2);//ambiguous pairs due to the error of ID
		  }
		  else {
			  if (DEBUG_FILTER) DEBUG_MSG("    DECOY");
			  decoyScore.push_back( alignGFScore );
			  fprintf(fp, "F");
		  }
	  }
	  else {
		  fprintf(fp, "X");//ambiguous pairs //ambiguous pairs due to similar ID
	  }
	  fprintf(fp, "\n");
	}
	fclose(fp);

	int numAmbiguousPSMs = 0;
	for(int i=0; i<xCount.size(); i++){
		if( xCount[i] ) numAmbiguousPSMs++;
	}
	DEBUG_VAR(numAmbiguousPSMs);

	pair<float, float> thr_fdr = fdr_thresholding(targetScore, decoyScore, fdrCut);
	float threshold = thr_fdr.first;
	m_fdr = thr_fdr.second;

	long org_num = thePairs.size();
	SpectrumPairSet filtered_pairs;

	sort_descending_by_score();
	for (int i = 0; i < thePairs.size(); i++) {
	  if (min(thePairs[i].score1, thePairs[i].score2) < threshold) break;
	  if( xCount[thePairs[i].spec1] || xCount[thePairs[i].spec2] ) continue;

	  filtered_pairs.push_back(thePairs[i]);
	}

	thePairs.resize(0);
	for (int i = 0; i < filtered_pairs.size(); i++) {
		thePairs.push_back(filtered_pairs[i]);
	}
	filtered_pairs.resize(0);
	sort_pairs_by_index();

	if (DEBUG_FILTER) DEBUG_VAR(targetScore.size());
	if (DEBUG_FILTER) DEBUG_VAR(decoyScore.size());
	if (DEBUG_FILTER) DEBUG_VAR(fdrCut);
	if (DEBUG_FILTER) DEBUG_VAR(threshold);

	DEBUG_MSG("Spectral pairs were evaluated based on " << psmSet.size() << " input PSMs.");
	DEBUG_MSG("The pairs were filtered out at edge FDR " << fdrCut);
	DEBUG_MSG("Align-GF p-value threshold was fixed to " << threshold << " (-log value).");
	DEBUG_MSG("Of " << org_num << " pairs, " << (org_num-thePairs.size())
			  << " were filtered out. Retained " << thePairs.size() << " pairs.");
	return true;
  }

  //---------------------------------------------------------------------------
  /*** helpers of SpectrumPairSet::filter_by_fdr *****************************/
  /*** getPRM, matchTwoPRM, getThreshold         *****************************/
  //should be merged to PSM class
  //---------------------------------------------------------------------------
  static void getPRM(string annotation, vector<float> & prm, const float aa_masses[])
  {
    string  		stripSeq = "";
    int 			modCount = 0;
    string 			delta;
    vector<int> 	msite;
    vector<string> 	msdelta;

    bool prevAA = true;
    for(int i=2; i<annotation.length()-2; i++){
      if( isalpha(annotation[i]) ){
        stripSeq += annotation[i];
        if(!prevAA){
          msdelta.push_back(delta);
          delta.clear();
        }
        prevAA = true;
      }
      else if(prevAA) {
        modCount++;
        if( stripSeq.length() == 0 ) msite.push_back(0);
        else msite.push_back(stripSeq.length()-1);
        prevAA = false;
        delta += annotation[i];
      }
      else if(!prevAA) {
    	if( annotation[i] == '+' || annotation[i] == '-' ){
    		msdelta.push_back(delta);
    		delta.clear();
    		msite.push_back(msite[msite.size()-1]);
    	}
    	delta += annotation[i];
      }
    }
    if(!prevAA){
      msdelta.push_back(delta);
    }

    prm.resize(stripSeq.length());
    vector<float> mods(stripSeq.length());
    if( modCount > 0 ){
      for(int k=0; k<msite.size(); k++){
        mods[msite[k]] += atof(msdelta[k].c_str());
      }
    }

    float prm_mass = 0;
    for(int i=0; i<stripSeq.length(); i++){
      prm_mass += aa_masses[stripSeq[i]-'A']+mods[i];
      prm[i] = prm_mass;
    }
  }

  static void getMatchQuality(vector<float> &AA,
						  	  vector<float> &BB,
						  	  float shift1,
						  	  float shift2,
						  	  pair<float, int> &matchQual,
						  	  float pmTol,
						  	  float fragTol)
  {
//	  if( fabs(shift2) < pmTol || fabs(shift1) < pmTol ){
	  if( fabs(shift2) < 1.001 || fabs(shift1) < 1 ){
		  matchTwoPRMs(AA, BB, matchQual, fragTol);
	//	  matchTwoPRMs(AA, BB, matchQual, 1.1);
		  return;
	  }

	  matchPartiallyTwoPRMs(AA, BB, shift1, shift2, matchQual, pmTol, fragTol);

	  if( matchQual.first < 0.999 ){
		  pair<float, int> srmMatch;
		  matchPartiallyTwoPRMs(AA, BB, shift2, shift1, srmMatch, pmTol, fragTol);
		  if( matchQual.first < srmMatch.first ) matchQual = srmMatch;
	  }
  }

  //---------------------------------------------------------------------------
  static void matchTwoPRMs(vector<float> &AA,
                           vector<float> &BB,
                           pair<float, int> &matchQual,
                           float fragTol)
  {
    int aaLen = AA.size()-1;
    int bbLen = BB.size()-1;
    double delta = BB[bbLen] - AA[aaLen];
    vector<int> aaMat(aaLen);
    vector<int> bbMat(bbLen);

    int start = 0;
    for(int i=0; i<aaLen; i++){
      for(int k=start; k<bbLen; k++){
        if( BB[k] < AA[i]-fragTol ) continue;
        else if( BB[k] > AA[i]+fragTol ) {
          start = k;
          break;
        }
        aaMat[i]=1;
        bbMat[k]=1;
        break;
      }
    }

    start = 0;
    for(int i=0; i<aaLen; i++){
      float shiftedAA = AA[i]+delta;
      for(int k=start; k<bbLen; k++){
        if( BB[k] < shiftedAA-fragTol ) continue;
        else if( BB[k] > shiftedAA+fragTol ) {
          start = k;
          break;
        }
        aaMat[i]=1;
        bbMat[k]=1;
        break;
      }
    }

    int aaHit = 0, bbHit = 0;
    int aaLongest = 0, bbLongest = 0;
    int consec = 0;

    for(int i=0; i<aaMat.size(); i++) {
      aaHit += aaMat[i];
      if( aaMat[i] != 0 ) consec++;
      else {
        if( aaLongest < consec ) aaLongest = consec;
        consec = 0;
      }
    }
    if( aaLongest < consec ) aaLongest = consec;


    consec = 0;
    for(int i=0; i<bbMat.size(); i++) {
      bbHit += bbMat[i];
      if( bbMat[i] != 0 ) consec++;
      else {
        if( bbLongest < consec ) bbLongest = consec;
        consec = 0;
      }
    }
    if( bbLongest < consec ) bbLongest = consec;

    //first: %common PRM Masses, second: longest consecutive ions
    matchQual.first = max((float)aaHit/aaMat.size(), (float)bbHit/bbMat.size());
    matchQual.second = min(aaLongest, bbLongest);
  }

  //---------------------------------------------------------------------------
    static void matchPartiallyTwoPRMs(vector<float> &AA,
                                      vector<float> &BB,
                                      float shift1,
                                      float shift2,
                                      pair<float, int> &matchQual,
                                      float pmTol,
                                      float fragTol)
    {
    	float buf= 18, commonPMError= pmTol+3, corrUnit = fragTol/4;

		int aaLen= AA.size()-1, bbLen= BB.size()-1;
		float aaMW= AA[aaLen], bbMW= BB[bbLen];

		int aaStart= 0, aaEnd= aaLen;
		int bbStart= 0, bbEnd= bbLen;

		if( shift2 > 0 ){
			//cutting spec1 Cterm
			for(int i=aaLen-1; i>-1; i--){
				if( AA[i] < aaMW-shift2-buf ){
					aaEnd = i+1;
					break;
				}
			}

		} else {
			//spec2
			for(int i=bbLen-1; i>-1; i--){
				if( BB[i] < bbMW+shift2-buf ){
					bbEnd = i+1;
					break;
				}
			}
		}

		if( shift1 > 0 ){
			//cutting spec1 Nterm
			for(int i=0; i<aaLen; i++){
				if( shift1+buf < AA[i] ){
					aaStart = i;
					break;
				}
			}
		} else {
			for(int i=0; i<bbLen; i++){
				if( -shift1+buf < BB[i] ){
					bbStart = i;
					break;
				}
			}
	  }

	  int maxMat = 0;
	  vector<int> aaMat(aaLen);
	  vector<int> bbMat(bbLen);

	  for(float adjshift=shift1-commonPMError; adjshift<=shift1+commonPMError; adjshift+=corrUnit) {

		  int localMat = 0;
		  vector<int> locAA(aaLen);
		  vector<int> locBB(bbLen);

		  int start = bbStart;
		  for(int i=aaStart; i<aaEnd; i++){
			for(int k=start; k<bbEnd; k++){
			  float shiftedMass= BB[k]+adjshift;
			  if( shiftedMass < AA[i]-fragTol ) continue;
			  else if( shiftedMass > AA[i]+fragTol ) {
				start = k;
				break;
			  }
			  locAA[i]=1;
			  locBB[k]=1;
			  localMat++;
			  break;
			}
		  }

		  if( maxMat < localMat ){
			  maxMat = localMat;
			  aaMat = locAA;
			  bbMat = locBB;
		  }
	  }

	  int aaHit = 0, bbHit = 0;
	  int aaLongest = 0, bbLongest = 0;
	  int consec = 0;

      for(int i=aaStart; i<aaEnd; i++) {
        aaHit += aaMat[i];
        if( aaMat[i] != 0 ) consec++;
        else {
          if( aaLongest < consec ) aaLongest = consec;
          consec = 0;
        }
      }
      if( aaLongest < consec ) aaLongest = consec;

      consec = 0;
      for(int i=bbStart; i<bbEnd; i++) {
        bbHit += bbMat[i];
        if( bbMat[i] != 0 ) consec++;
        else {
          if( bbLongest < consec ) bbLongest = consec;
          consec = 0;
        }
      }
      if( bbLongest < consec ) bbLongest = consec;

      //first: %common PRM Masses, second: longest consecutive ions
      matchQual.first = max((float)aaHit/(aaEnd-aaStart), (float)bbHit/(bbEnd-bbStart));
      matchQual.second = min(aaLongest, bbLongest);
    }


  //---------------------------------------------------------------------------
  static pair<float, float> fdr_thresholding(vector<float> &target, vector<float> &decoy, float FDRCut)
  {
    sort( target.begin(), target.end() );
    sort( decoy.begin(), decoy.end() );
    int tarSize = target.size();
    int decSize = decoy.size();

 //   if (DEBUG_FILTER)
    	DEBUG_MSG("[Thresholding] #Annotated Pairs: " << (tarSize + decSize)
                  << " | Target = " << tarSize << " | Decoy= " << decSize);

    if( decSize == 0 || 1 <= FDRCut ) return pair<float, float>(0, 0);

    int t_i=0, d_i=0;
    int targetHit=0, decoyHit=0;
    double threshold = 0;

    for( d_i=0; d_i < decSize; d_i++ ){
      while( t_i < tarSize && target[t_i] < decoy[d_i] ){
        t_i++;
      }
      int T = tarSize - t_i;
      if( T == 0 ) break;
      int D = decSize - d_i - 1;
      if( (double)D/(T+D) <= FDRCut ){
        targetHit= T;
        decoyHit = D;
        threshold= decoy[d_i];
        break;
      }
    }

//    if (DEBUG_FILTER)
    	DEBUG_MSG("[Thresholding] Threshold =  " << threshold
                  << " | FDR = " << (double)(decoyHit)/(targetHit+decoyHit)*100
                  << " | Target = " << targetHit << " | Decoy= " << decoyHit);

    return pair<float, float>(threshold, (double)(decoyHit)/(targetHit+decoyHit));
  }

  //---------------------------------------------------------------------------
  static int isCrossMatched(Spectrum &specAA,
		  	  	  	  	  	float deltaScore,
				   	   	    vector<float> &aaPRM,
				   	   	   	vector<float> &bbPRM,
				   	   	   	float modDelta,
				   	   	   	float fragTol)
  {
	int specLen = specAA.size();

	vector<float> AASpec;
	int aaLen = aaPRM.size()-1;
	float aaPM= aaPRM[aaLen];
	for( int a=0; a<aaLen; a++ ){
		AASpec.push_back(aaPRM[a]);//b-series
		AASpec.push_back(aaPM-aaPRM[a]+AAJumps::massH2O);//y-series
	}
	sort( AASpec.begin(), AASpec.end() );

	int start = 0;
	float aaScore = 0;
	for( int a=0; a<AASpec.size(); a++ ){
		for (int i=start; i<specLen; i++){
			if( specAA[i][0] < AASpec[a]-fragTol ) continue;
			else if( specAA[i][0] > AASpec[a]+fragTol ) {
				start = i;
				break;
			}
			aaScore += specAA[i][1];
		}
	}//get pept1 score

	int numModSite = bbPRM.size();
	vector<vector<float> > modBBSpec(numModSite);

	int bbLen = bbPRM.size()-1;
	float bbPM= bbPRM[bbLen]+modDelta;//to match to AA spcetrum
	for(int site=0; site<numModSite; site++){
		for(int b=0; b<site; b++){
			modBBSpec[site].push_back(bbPRM[b]);
			modBBSpec[site].push_back(bbPM-bbPRM[b]+AAJumps::massH2O);
		}
		for(int b=site; b<bbLen; b++){
			modBBSpec[site].push_back(bbPRM[b]+modDelta);
			modBBSpec[site].push_back(bbPM-bbPRM[b]-modDelta+AAJumps::massH2O);
		}
		sort(modBBSpec[site].begin(), modBBSpec[site].end());
	}

	int mPos = -1;
	float bbScore = 0;
	for(int mod=0; mod<numModSite; mod++){

		start = 0;
		float score = 0;
		for( int b=0; b<modBBSpec[mod].size(); b++ ){
			for (int i=start; i<specLen; i++){
				if( specAA[i][0] < modBBSpec[mod][b]-fragTol ) continue;
				else if( specAA[i][0] > modBBSpec[mod][b]+fragTol ) {
					start = i;
					break;
				}
				score += specAA[i][1];
			}
		}
		if( bbScore < score ) {
			bbScore = score;
			mPos = mod;
		}
	}

	if( aaScore < bbScore ) return 2;
	if( aaScore-bbScore < deltaScore || aaScore*0.8 < bbScore ) return 1;
	return 0;
  }

  bool SpectrumPairSet::cosine_precision_recall(PeptideSpectrumMatchSet & psmSet,
		  	  	  	  	  	  	  	  	  	  	const std::string & filename,
                                           	    float minOverlap,
                                           	    bool seperateCharge,
                                           	    bool canonicalPSMForm)
  {
    static float coPRMCut = 0.6;
    static int   coLenCut = 12;
    static float coElutionWindow = 3;

    float aa_masses[]={
						71.03711, 0, 103.00919, 115.02694, 129.04259,
						147.06841, 57.02146, 137.05891, 113.08406, 0,
						128.09496, 113.08406, 131.04049, 114.04293, 0,
						97.05276, 128.05858, 156.10111, 87.03203, 101.04768,
						0, 99.06841, 186.07931, 0, 163.06333, 0};

	if( !canonicalPSMForm ){
		AAJumps jumps(1);
		for(int i=0; i<jumps.size(); i++) aa_masses[jumps.aaLetters[i]-'A'] = jumps[i];
		DEBUG_VAR(aa_masses[2]);
	}

	int specSetSize= -1;

    vector<vector<float> > prmSpecs(psmSet.size());
    psmSet.sortBySpecIndex(); //must
	for (int i = 0; i < psmSet.size(); i++) {

	  if( specSetSize < psmSet[i]->m_scanNum-1 ) specSetSize = psmSet[i]->m_scanNum-1;
	  if (psmSet[i]->m_annotation.find("!") != string::npos) {
		continue;
	  }

	  string massaged_annotation = specnets_annotation_to_msgf(msp_annotation_to_specnets(psmSet[i]->m_annotation));

	  vector<float> prm;
	  getPRM(massaged_annotation, prm, aa_masses);
	  prmSpecs[i] = prm;
	}

	for(int i = 0; i < thePairs.size(); i++){
		if( specSetSize < thePairs[i].spec2 ) specSetSize = thePairs[i].spec2;
	}
	specSetSize += 1;

	vector<map<int, float> > pair_score_map(specSetSize); //pair-score map
	for(int i = 0; i < thePairs.size(); i++){ // 2 lines fixed
	  pair_score_map[thePairs[i].spec1][thePairs[i].spec2] = thePairs[i].score1;//cosine score
	}

	FILE *fp = fopen(filename.c_str(), "w");
	if (fp == 0) {
	  ERROR_MSG("Can not open: " << filename);
	  return false;
	}
	DEBUG_MSG("Annotated pair file: " << filename);

	fprintf(fp, "overlap\tdistance\t" );
	fprintf(fp, "spec1\tpept1\tcs1\tmsgfscore1\t-log(specprob1)\tprotein1\t" );
	fprintf(fp, "spec2\tpept2\tcs2\tmsgfscore2\t-log(specprob2)\tprotein2\t" );
	fprintf(fp, "correct\tspspair\tcosine\tmzdiff\n" );

	int numTruePairs= 0;
    vector<float> targetScore;
    vector<float> decoyScore;
    pair<float, int> matchQual; //first: %common PRM Masses, second: longest consecutive ions

    for (int i = 0; i < psmSet.size(); i++) {
      if ( prmSpecs[i].size() == 0 ) continue;

      int spec1 = psmSet[i]->m_scanNum-1;
      float MW1 = prmSpecs[i][prmSpecs[i].size()-1];

      for (int k = i+1; k < psmSet.size(); k++) {
	if ( prmSpecs[k].size() == 0 ) continue;

        int spec2 = psmSet[k]->m_scanNum-1;
        float MW2 = prmSpecs[k][prmSpecs[k].size()-1];

        if ( seperateCharge && psmSet[i]->m_charge != psmSet[k]->m_charge ) continue;
//        if ( fabs(MW1-MW2) > 200 ) continue;

        int smallLen = prmSpecs[i].size(), largeLen = prmSpecs[k].size();
		if( smallLen > largeLen ){
			smallLen = prmSpecs[k].size(),
			largeLen = prmSpecs[i].size();
		}
        float overplap = (float)smallLen/largeLen;

        if ( overplap < minOverlap ) continue;

        int paired = ( pair_score_map[spec1].find(spec2) == pair_score_map[spec1].end() )? 0 : 1;
        float cosine = ( paired == 1 )? pair_score_map[spec1][spec2] : 0.;

        matchTwoPRMs(prmSpecs[i], prmSpecs[k], matchQual, 0.5); //match two ID

        int correct= 2;//ambiguous
        if (matchQual.first > 0.999) correct= 1;
		else if (matchQual.first < coPRMCut && matchQual.second < coLenCut) correct= 0;//*/

//        int correct = (matchQual.first > 0.999)? 1 : 0;

        if( correct == 0 && paired == 0 ) continue;

        float deltaMass = MW1-MW2;

        fprintf(fp, "%.4f\t%d\t", overplap, (int)round(smallLen*(1-matchQual.first)) );
        fprintf(fp, "%d\t%s\t%d\t%d\t%.4f\t%s\t",
        		spec1+1, psmSet[i]->m_annotation.c_str(), psmSet[i]->m_charge, (int)psmSet[i]->m_compScore, -log10(psmSet[i]->m_score), psmSet[i]->m_protein.c_str());
        fprintf(fp, "%d\t%s\t%d\t%d\t%.4f\t%s\t",
        		spec2+1, psmSet[k]->m_annotation.c_str(), psmSet[k]->m_charge, (int)psmSet[k]->m_compScore, -log10(psmSet[k]->m_score), psmSet[k]->m_protein.c_str());
        fprintf(fp, "%d\t%d\t%.4f\t%.4f\n", correct, paired, cosine, deltaMass);

        if( overplap < 0.99 || fabs(deltaMass) < 0.5 ) continue;
        //now check accuracy, removed similar masses
		if( correct == 1 ) numTruePairs++; //#theoretical pairs

		if( paired == 1 ){
			if( correct == 1 ) {
				targetScore.push_back(cosine);
			}
			else if( correct == 0 ) decoyScore.push_back(cosine);
		}

      }//loop spec2
    }//loop spec1
    fclose(fp);

    //get precision, recall
    int tarSize = targetScore.size();
    int decSize = decoyScore.size();
    sort( targetScore.begin(), targetScore.end() );
    sort( decoyScore.begin(), decoyScore.end() );
    DEBUG_MSG("#TURE Pairs: " << numTruePairs);
    DEBUG_MSG("#Detected TRUE: " << tarSize << ", #FALSE: " << decSize);

    int tar=0, dec=0;
    std::cout<<"Cosine\tPrecision\tRecall"<<std::endl;
    for(float thr = 0.f; thr < 1.f; thr+= 0.10){

    	while( tar < tarSize && targetScore[tar] < thr ){
    		tar++;
    	}
    	while( dec < decSize && decoyScore[dec] < thr ){
    		dec++;
    	}

    	int tarAlive = tarSize-tar;
    	if( tarAlive == 0 ) break;

    	std::cout<< thr << "\t" <<
    					((float)tarAlive/(tarAlive+decSize-dec)*100) << "\t" <<
    							((float)tarAlive/numTruePairs*100) << std::endl;
    }

    return true;
  }

  bool SpectrumPairSet::cosine_precision_recall_curve(PeptideSpectrumMatchSet & psmSet,
                                           	          float minOverlap,
                                           	          bool seperateCharge,
                                           	          bool canonicalPSMForm)
  {
    static float coPRMCut = 0.6;
    static int   coLenCut = 12;
    static float coElutionWindow = 3;

    float aa_masses[]={
						71.03711, 0, 103.00919, 115.02694, 129.04259,
						147.06841, 57.02146, 137.05891, 113.08406, 0,
						128.09496, 113.08406, 131.04049, 114.04293, 0,
						97.05276, 128.05858, 156.10111, 87.03203, 101.04768,
						0, 99.06841, 186.07931, 0, 163.06333, 0};

	if( !canonicalPSMForm ){
		AAJumps jumps(1);
		for(int i=0; i<jumps.size(); i++) aa_masses[jumps.aaLetters[i]-'A'] = jumps[i];
		DEBUG_VAR(aa_masses[2]);
	}

	int specSetSize= -1;

    vector<vector<float> > prmSpecs(psmSet.size());
    psmSet.sortBySpecIndex(); //must
	for (int i = 0; i < psmSet.size(); i++) {

	  if( specSetSize < psmSet[i]->m_scanNum-1 ) specSetSize = psmSet[i]->m_scanNum-1;
	  if (psmSet[i]->m_annotation.find("!") != string::npos) {
		continue;
	  }

	  string massaged_annotation = specnets_annotation_to_msgf(msp_annotation_to_specnets(psmSet[i]->m_annotation));

	  vector<float> prm;
	  getPRM(massaged_annotation, prm, aa_masses);
	  prmSpecs[i] = prm;
	}

	for(int i = 0; i < thePairs.size(); i++){
		if( specSetSize < thePairs[i].spec2 ) specSetSize = thePairs[i].spec2;
	}
	specSetSize += 1;

	vector<map<int, float> > pair_score_map(specSetSize); //pair-score map
	for(int i = 0; i < thePairs.size(); i++){ // 2 lines fixed
	  pair_score_map[thePairs[i].spec1][thePairs[i].spec2] = thePairs[i].score1;//cosine score
	}

	int numTruePairs= 0;
    vector<float> targetScore;
    vector<float> decoyScore;
    pair<float, int> matchQual; //first: %common PRM Masses, second: longest consecutive ions

    for (int i = 0; i < psmSet.size(); i++) {
      if ( prmSpecs[i].size() == 0 ) continue;

      int spec1 = psmSet[i]->m_scanNum-1;
      float MW1 = prmSpecs[i][prmSpecs[i].size()-1];

      for (int k = i+1; k < psmSet.size(); k++) {
	if ( prmSpecs[k].size() == 0 ) continue;

        int spec2 = psmSet[k]->m_scanNum-1;
        float MW2 = prmSpecs[k][prmSpecs[k].size()-1];

        if ( seperateCharge && psmSet[i]->m_charge != psmSet[k]->m_charge ) continue;

        int smallLen = prmSpecs[i].size(), largeLen = prmSpecs[k].size();
		if( smallLen > largeLen ){
			smallLen = prmSpecs[k].size(),
			largeLen = prmSpecs[i].size();
		}
        float overplap = (float)smallLen/largeLen;

        if ( overplap < minOverlap ) continue;

        int paired = ( pair_score_map[spec1].find(spec2) == pair_score_map[spec1].end() )? 0 : 1;
        float cosine = ( paired == 1 )? pair_score_map[spec1][spec2] : 0.;

        matchTwoPRMs(prmSpecs[i], prmSpecs[k], matchQual, 0.5); //match two ID

        int correct= 2;//ambiguous
		if (matchQual.first > 0.999) correct= 1;
		else if (matchQual.first < coPRMCut && matchQual.second < coLenCut) correct= 0;

//        int correct = (matchQual.first > 0.999)? 1 : 0;

        if( correct == 0 && paired == 0 ) continue;

        float deltaMass = MW1-MW2;

        if( overplap < 0.99 || fabs(deltaMass) < 0.5 ) continue;
        //now check accuracy, removed similar masses
		if( correct == 1 ) numTruePairs++; //#theoretical pairs

		if( paired == 1 ){
			if( correct == 1 ) {
				targetScore.push_back(cosine);
			}
			else if( correct == 0 ) decoyScore.push_back(cosine);
		}

      }//loop spec2
    }//loop spec1

    //get precision, recall
    int tarSize = targetScore.size();
    int decSize = decoyScore.size();
    sort( targetScore.begin(), targetScore.end() );
    sort( decoyScore.begin(), decoyScore.end() );
    DEBUG_MSG("#TURE Pairs: " << numTruePairs);
    DEBUG_MSG("#Detected TRUE: " << tarSize << ", #FALSE: " << decSize);

    int tar=0, dec=0;
    std::cout<<"Cosine\tPrecision\tRecall"<<std::endl;
    for(float thr = 0.f; thr < 1.f; thr+= 0.10){

    	while( tar < tarSize && targetScore[tar] < thr ){
    		tar++;
    	}
    	while( dec < decSize && decoyScore[dec] < thr ){
    		dec++;
    	}

    	int tarAlive = tarSize-tar;
    	if( tarAlive == 0 ) break;

    	std::cout<< thr << "\t" <<
    					((float)tarAlive/(tarAlive+decSize-dec)*100) << "\t" <<
    							((float)tarAlive/numTruePairs*100) << std::endl;
    }

    return true;
  }



} // namespace specnets





