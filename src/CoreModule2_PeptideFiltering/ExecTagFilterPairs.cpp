/*
 * ExecTagFilterPairs.cpp
 *
 *  Created on: Mar 11, 2016
 *      Author: isna
 */

#include "ExecTagFilterPairs.h"
#include "ExecDeconvoluteMS2.h"

namespace specnets {

	static int 		maxTagSeqKey  = 8000;
	static int 		pm_err_reward = 3;
	static float	iso_err_reward = 1.5;
	static float	flanking_mass_limit = 10000;
	static int 		maxNtermFlankingMass = 0;
	static int 		maxCtermFlankingMass = 0;

	static int aaIndex(char aa){

		switch(aa)
		{
			case 'A' : return 0;
			case 'C' : return 1;
			case 'D' : return 2;
			case 'E' : return 3;
			case 'F' : return 4;
			case 'G' : return 5;
			case 'H' : return 6;
			case 'I' : return 9;
			case 'K' : return 8;
			case 'L' : return 9;
			case 'M' : return 10;
			case 'N' : return 11;
			case 'P' : return 12;
			case 'Q' : return 8;
			case 'R' : return 14;
			case 'S' : return 15;
			case 'T' : return 16;
			case 'V' : return 17;
			case 'W' : return 18;
			case 'Y' : return 19;
			default :  return -1;
		}
	}

	static bool tag_score_comparator(tag_T t1, tag_T t2){
		return t1.score > t2.score;
	} // Used to sort tags by descending tag score

	static bool tag_seq_comparator(tag_T t1, tag_T t2){
		if( t1.seqKey < t2.seqKey ) return true;
		else if( t1.seqKey > t2.seqKey ) return false;
		return t1.score > t2.score;
	}

	ExecTagFilterPairs::ExecTagFilterPairs(void) :
		ExecBase(), m_prmSpectra(0x0), isPrm(false), m_maxTagSize(50),
		m_spec_tags(0x0), nterm_prefix_fmMap(0x0), cterm_suffix_fmMap(0x0),
		ownInput(true), ownOutput(true)
	{
		m_name = "ExecTagFilterPairs";
		m_type = "ExecTagFilterPairs";
	}

	ExecTagFilterPairs::ExecTagFilterPairs(const ParameterList & inputParams) :
		ExecBase(inputParams), m_prmSpectra(0x0), isPrm(false), m_maxTagSize(50),
		m_spec_tags(0x0), nterm_prefix_fmMap(0x0), cterm_suffix_fmMap(0x0),
		ownInput(true), ownOutput(true)
	{
		m_name = "ExecTagFilterPairs";
		m_type = "ExecTagFilterPairs";
	}

	ExecTagFilterPairs::ExecTagFilterPairs(const ParameterList & inputParams,
										   SpecSet * prmSpectra) :
		ExecBase(inputParams), m_prmSpectra(prmSpectra), isPrm(false), m_maxTagSize(50),
		m_spec_tags(0x0), nterm_prefix_fmMap(0x0), cterm_suffix_fmMap(0x0),
		ownInput(false), ownOutput(false)
	{
		m_name = "ExecTagFilterPairs";
		m_type = "ExecTagFilterPairs";
	}

	ExecTagFilterPairs::~ExecTagFilterPairs() {

		if( ownInput ){
			delete m_prmSpectra;
		}
		if( m_spec_tags != 0x0 ) delete m_spec_tags;
		if( nterm_prefix_fmMap != 0x0 ) delete nterm_prefix_fmMap;
		if( cterm_suffix_fmMap != 0x0 ) delete cterm_suffix_fmMap;
	}

	// -------------------------------------------------------------------------
	ExecBase * ExecTagFilterPairs::clone(const ParameterList & inputParams) const {
		return new ExecTagFilterPairs(inputParams);
	}

	// -------------------------------------------------------------------------
	bool ExecTagFilterPairs::invoke(void) {

		m_maxTagSize 		= m_params.getValueInt("MAX_TAG_SIZE", 50);

		float fragTol 		= m_params.getValueDouble("TOLERANCE_PEAK", 0.5);
		bool alignPA 		= m_params.getValueBool("PARTIAL_OVERLAPS", false);

		int generateTag		= m_params.getValueInt("TAG_GEN", 0);
		string tags_file    = m_params.getValue("INPUT_TAGS", "");

		DEBUG_VAR(m_maxTagSize);
		DEBUG_VAR(alignPA);
		DEBUG_VAR(generateTag);
		DEBUG_VAR(tags_file);

		//tag generation in file
		if( tags_file != "" ) {
			if( m_params.getValueBool("TAG_FILTER", false) ) {
				m_prmSpectra->saveTags(tags_file, 3, 1, m_maxTagSize*2, fragTol, isPrm);
			}
			else {
				FILE *fp = fopen(tags_file.c_str(), "w");
				fclose(fp);
			}
		}

		return true;
	}

	// -------------------------------------------------------------------------
	bool ExecTagFilterPairs::prepare(string tags_file, int startSpec, bool matchingFlankingMass) {

		m_maxTagSize 		= m_params.getValueInt("MAX_TAG_SIZE", 50);

		DEBUG_VAR(m_maxTagSize);

		//get tags
		if( !readTagsFile(tags_file) ) return false;

		//construct index
		if( matchingFlankingMass ) indexTags( startSpec );

		return true;
	}


	// -------------------------------------------------------------------------
	bool ExecTagFilterPairs::loadInputData(void) {

		if( m_prmSpectra == 0x0 ) {
			ownInput = true;
			m_prmSpectra = new SpecSet;
		}

		if (m_params.exists("INPUT_PRM_SPECS")){
			m_prmSpectra->loadPklBin(m_params.getValue("INPUT_PRM_SPECS").c_str());
			DEBUG_MSG( "#PRM Spectra: " << m_prmSpectra->size() );
			isPrm = true;
		}
		else if (m_params.exists("INPUT_SPECS")){

			//non-prm spectra for cosine workflow
			isPrm = false;
			if( !m_params.getValueBool("DECONV_MS2", false) ){
				m_prmSpectra->loadPklBin(m_params.getValue("INPUT_SPECS").c_str());
				DEBUG_MSG( "#MS Spectra: " << m_prmSpectra->size() );
			}
			else {
				SpecSet inputSpecs;
				inputSpecs.loadPklBin(m_params.getValue("INPUT_SPECS").c_str());
				DEBUG_MSG( "#MS Spectra: " << inputSpecs.size() );

				m_params.setValue("MIN_PEAK_INT", "0");
				m_params.setValue("FILTER_STDDEV_PEAK_INT", "0");
				DEBUG_VAR(m_params.getValue("MIN_PEAK_INT", "0"));
				DEBUG_VAR(m_params.getValue("FILTER_STDDEV_PEAK_INT", "0"));

				IsoEnvelope isoModel;
				string modelPath = getPath(m_params.getValue("ISO_ENV"), "model_isoenv.bin", false);

				DEBUG_MSG("Loading input iso envelope from <" << modelPath << "> ...");
				if( !isoModel.LoadModel(modelPath.c_str()) ){
					ERROR_MSG("Could not load " << modelPath);
				}

				ExecDeconvoluteMS2 module(m_params,
										  &inputSpecs,
										  &isoModel,
										  m_prmSpectra);

				if (!module.invoke()) return false;
			}
			DEBUG_MSG( "#Deconvoluted Spectra: " << m_prmSpectra->size() );
		}
		return true;
	}

	// -------------------------------------------------------------------------
	bool ExecTagFilterPairs::saveOutputData(void) {
		return true;
	}

	// -------------------------------------------------------------------------
	bool ExecTagFilterPairs::saveInputData(std::vector<std::string> & filenames) {
		return true;
	}

	// -------------------------------------------------------------------------
	bool ExecTagFilterPairs::loadOutputData(void) {
		return true;
	}

	// -------------------------------------------------------------------------
	vector<ExecBase*> const & ExecTagFilterPairs::split(int numSplit) {
		return m_subModules;
	}

	// -------------------------------------------------------------------------
	bool ExecTagFilterPairs::merge(void) {
		return true;
	}

	// -------------------------------------------------------------------------
	bool ExecTagFilterPairs::validateParams(std::string & error) {

	  m_isValid = false;

//	  VALIDATE_PARAM_EXIST("INPUT_PRM_SPECS");

	  m_isValid = true;
	  return true;
	}

	bool ExecTagFilterPairs::readTagsFile(const string &tag_file) {

		char buf[512];
		FILE *fp = fopen(tag_file.c_str(), "r");
		if (fp == 0) {
		  ERROR_MSG("Can not open: " << tag_file);
		  return false;
		}

		m_spec_tags = new vector<vector<tag_T> >(m_prmSpectra->size());

		while( fgets(buf, 512, fp) != NULL ){

			if( strstr(buf, ">>") ){

				int spec = atoi(buf+2)-1; //specIndex is zero-based

				vector<tag_T> st;
				char tseq[4]="";
				float pm=0, sm=0, sc=0;
				while( fgets(buf, 512, fp) != NULL ){//read lines
					if( sscanf(buf, "%f\t%s\t%f\t%f", &pm, tseq, &sm, &sc) != 4 ) break;
					if( flanking_mass_limit < pm || flanking_mass_limit < sm ) continue;

					int seqID = aaIndex(tseq[0])*400 + aaIndex(tseq[1])*20 + aaIndex(tseq[2]);
					tag_T tag = {(int)round(pm), seqID, (int)round(sm+AAJumps::massH2O), sc};
					st.push_back( tag );
				}

				sort(st.begin(), st.end(), tag_seq_comparator);
				for(vector<tag_T>::iterator itr=st.begin(); itr!=st.end(); itr++) {
					for(vector<tag_T>::iterator ktr=itr+1; ktr!=st.end(); ktr++) {
						if( itr->seqKey == ktr->seqKey ){
							if( abs(itr->ntMass - ktr->ntMass ) < iso_err_reward ){
								st.erase(ktr);
								ktr--;
							}
						}
						else break;
					}
				}

				sort(st.begin(), st.end(), tag_score_comparator);
				int size = 0;
				for(vector<tag_T>::iterator itr=st.begin(); itr!=st.end(); itr++){
					if( size++ == m_maxTagSize ){
						st.erase(itr, st.end());
						break;
					}
					if( maxNtermFlankingMass < itr->ntMass ) maxNtermFlankingMass = itr->ntMass;
					if( maxCtermFlankingMass < itr->ctMass ) maxCtermFlankingMass = itr->ctMass;
				}
				(*m_spec_tags)[spec] = st;
			}
		}
		fclose(fp);

		maxNtermFlankingMass += pm_err_reward+1;
		maxCtermFlankingMass += pm_err_reward+1;

		DEBUG_MSG("Read tags: " << m_maxTagSize << " tags/spectrum selected.");
		DEBUG_MSG("Max flanking mass at Nterm: " << maxNtermFlankingMass << " Cterm: " << maxCtermFlankingMass);

		return true;
	}

	void ExecTagFilterPairs::indexTags(int startSpec) {

		int setSize = m_prmSpectra->size();

		nterm_prefix_fmMap= new vector<vector<list<int> > >(maxTagSeqKey, vector<list<int> >(maxNtermFlankingMass));
		cterm_suffix_fmMap= new vector<vector<list<int> > >(maxTagSeqKey, vector<list<int> >(maxCtermFlankingMass));
		for(int i=startSpec; i<setSize; i++){
			for(vector<tag_T>::iterator itr=(*m_spec_tags)[i].begin(); itr!=(*m_spec_tags)[i].end(); itr++){
				(*nterm_prefix_fmMap)[itr->seqKey][itr->ntMass].push_back(i);
				(*cterm_suffix_fmMap)[itr->seqKey][itr->ctMass].push_back(i);
			}
		}
		DEBUG_MSG("Tag-Indexing Done.");
	}

	//for ASP Alignment, i.e. finding modifications
	bool ExecTagFilterPairs::getPairs(int spec1, list<int> &candidates, list<int> &subpairs) {

		set<int> matS;
		for(vector<tag_T>::iterator itr=(*m_spec_tags)[spec1].begin(); itr!=(*m_spec_tags)[spec1].end(); itr++){
			int sid = itr->seqKey;

			int start = ( itr->ntMass > pm_err_reward )? itr->ntMass-pm_err_reward : 0;
			int end = itr->ntMass+pm_err_reward;
			for(int k=start; k <=end; k++){
				for(list<int>::iterator ktr = (*nterm_prefix_fmMap)[sid][k].begin(); ktr!=(*nterm_prefix_fmMap)[sid][k].end(); ktr++)
					if( spec1 < *ktr ) matS.insert(*ktr);
			}

			start = ( itr->ctMass > pm_err_reward )? itr->ctMass-pm_err_reward : 0;
			end = itr->ctMass+pm_err_reward;
			for(int k=start; k <=end; k++){
				for(list<int>::iterator ktr = (*cterm_suffix_fmMap)[sid][k].begin(); ktr!=(*cterm_suffix_fmMap)[sid][k].end(); ktr++)
					 if( spec1 < *ktr ) matS.insert(*ktr);
			}
		}

		vector<int> matList(matS.size());
		for(set<int>::iterator itr = matS.begin(); itr!=matS.end(); itr++){
			matList.push_back(*itr);
		}
		sort(matList.begin(), matList.end());

		//matList & candidates should be sorted;
		int ASize= matList.size(), BSize= candidates.size();

		int a = 0;
		list<int>::iterator b = candidates.begin();
		while( a < ASize && b != candidates.end() ){
			int cur = matList[a];
			if( cur > *b ) b++;
			else if( cur < *b ) a++;
			else { // matched
				subpairs.push_back(cur);
				a++;
				b++;
			}
		}
		return true;
	}

	//for Partial Overlap Alignment
	void ExecTagFilterPairs::getPossibleDeltas(set<int> &deltas, int spec1, int spec2, float spec1MW, float spec2MW, float minOvlpMass) {

		float leftM = minOvlpMass - pm_err_reward;
		float rightM = spec1MW + spec2MW - leftM;

		for(vector<tag_T>::iterator itr=(*m_spec_tags)[spec1].begin(); itr!=(*m_spec_tags)[spec1].end(); itr++){
			for(vector<tag_T>::iterator ktr=(*m_spec_tags)[spec2].begin(); ktr!=(*m_spec_tags)[spec2].end(); ktr++){

				if( itr->seqKey == ktr->seqKey ){

					float tpM = spec2MW + (itr->ntMass-ktr->ntMass);
					if( leftM <= tpM && tpM <= rightM ) deltas.insert( itr->ntMass-ktr->ntMass );

					tpM = spec2MW + (itr->ctMass-ktr->ctMass);
					if( leftM <= tpM && tpM <= rightM ) deltas.insert( itr->ctMass-ktr->ctMass );

				/*	deltas.insert( itr->ntMass-ktr->ntMass );
					deltas.insert( itr->ctMass-ktr->ctMass );//*/
				}
			}
		}
	}

	bool ExecTagFilterPairs::hasTags(){
		return m_spec_tags != 0x0;
	}


/*	//deprecated
	void ExecTagFilterPairs::getPossibleDeltas(int spec1, vector<set<int> > &pairDeltas) {

		pairDeltas.resize(m_prmSpectra->size());

		for(vector<tag_T>::iterator itr=(*m_spec_tags)[spec1].begin(); itr!=(*m_spec_tags)[spec1].end(); itr++){

			int sid = itr->seqKey;

			for(int k=0; k <maxNtermFlankingMass; k++){

				list<int>::iterator start = (*nterm_prefix_fmMap)[sid][k].begin(), end =(*nterm_prefix_fmMap)[sid][k].end();
				for(list<int>::iterator ktr = start; ktr != end; ktr++)
					if( spec1 < *ktr ) pairDeltas[*ktr].insert( itr->ntMass-k );
			}

			for(int k=0; k <maxCtermFlankingMass; k++){

				list<int>::iterator start = (*cterm_suffix_fmMap)[sid][k].begin(), end =(*cterm_suffix_fmMap)[sid][k].end();
				for(list<int>::iterator ktr = start; ktr != end; ktr++)
					if( spec1 < *ktr ) pairDeltas[*ktr].insert( itr->ctMass-k );
			}
		}
	}//*/

} /* namespace specnets */





