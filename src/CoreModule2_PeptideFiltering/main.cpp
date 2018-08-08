#include <stdio.h>
#include <iostream>	// std::cout
#include <iterator>
#include <string>
#include <regex>

#include <algorithm>	// std::sort
#include <vector>	// std::vector
#include <set>
#include <fstream>	// std::ifstream
#include <map>
#include <math.h>
#include <time.h>	//time_t, struct tm, difftime, time, mktime

#include <libgen.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <sstream>

#include <stdint.h>

#include "Logger.h"
#include "spectrum.h"
#include "SpecSet.h"
#include "SpectrumPairSet.h"
#include <stdlib.h>

#include "spectrum_window_filter.h"

//#include "SpectrumPeptideScoring.h"
typedef vector <uint32_t> A_spectrum_simple;
typedef vector<A_spectrum_simple> All_spectra_simple;

#include "DatabaseIndex.h"

#include <boost/algorithm/string.hpp>
#include <boost/range/algorithm/remove_if.hpp>
#include <boost/range/algorithm_ext/erase.hpp>
#include <boost/algorithm/string/replace.hpp>

#include <set>

#include <unordered_set>
#include <unordered_map>
#include <queue>

#include "tags.h"

using namespace specnets;

//=== CONSTANTS === 
const string DECOY_TAG = "XXX";
const int ISOTOPE_NUM = 1;
const float ISOTOPE_EXTRA_MASS = 1.003355;

//=== USAGE === 
static void show_usage(std::string name)
{
	std::cerr << "Usage: ProteoStorm.exe \n" 
		<< " -i <spectral_partition_input> \n"
		<< "-p <database_partition_input> \n"
		<< "-o <spectrum_peptide_pair_output> \n"
		<< "-s <SPC_cutoff> \n"
		<< "-ms1tol <precursor_mass_tolerance_ppm> \n"
		<< "-ms2tol <fragment_mass_tolerance_ppm> \n"
		<< "-tmt <TMT_labeling_boolean> \n"
		<< std::endl;
}

//=== DECLARE FUNCTIONS === 

// parse arguments function
bool parseArguments(int argc, const char* argv[], map<string, string>& args_map);

// read peptides from database function
double read_peptides_allpeaks(const char* peptides_fn, vector<string> &pep_sequences, vector<string> &protein_ids, vector<uint32_t> &pep_masses, int TMT_labeling);
void remove_fixedCmod_from_peptides(vector<string> &pep_sequences, int TMT_labeling);

// read spectra
void read_spectra_allpeaks(const char* mgf_fn, All_spectra_simple &spectra_vec, vector<uint32_t> &ionmass_vec, vector<uint32_t> &sorted_spectrum_indices, int TMT_labeling);

// counting SPC for spectra
void spectra_against_peptides_allpeaks(All_spectra_simple const& spectra_vec, vector<uint32_t> const& ionmass_vec, vector<uint32_t> const& sorted_spectrum_indices, vector<string> &protein_ids, vector<string> const& peptides, vector<uint32_t> const& peptide_masses, AllPeaks const& db_data, int ms1_tolerance, float ms2_tolerance, const char* output_fn, int score_cutoff);

void find_matched_peptide_indices(AllPeaks const& db_data, uint32_t spectrum_hit_ind, uint32_t ion_mass_min, uint32_t ion_mass_max, uint32_t theoretical_mass_ind_begin, uint32_t theoretical_mass_ind_end, vector <uint32_t> const& peptide_masses, vector <short> & spectrum_scores );


//=== MAIN FUNCTION === 
int main(int argc, const char* argv[]) // or const char** argv
{
	//check number of paramters
	if (argc < 8){
		//show usage text
		show_usage(argv[0]);
		return 1;
	}
	map<string, string> args_map;
	if (!parseArguments(argc, argv, args_map)) return 1;
	
	const char* mgf_fn = args_map["spectral_partition"].c_str();
	const char* database_peptides_fn = args_map["database_partition"].c_str();
	const char* output_fn = args_map["sp_pair_out"].c_str();
	const int score_cutoff = std::stoi(args_map["spc_score"].c_str());
	const int MS1_PPM_TOLERANCE = std::stoi(args_map["precursor_tol"].c_str()); // 10ppm for high res
	const float MS2_PPM_TOLERANCE = std::stof(args_map["fragment_tol"].c_str()); // 0.015 for high res, 0.6 for low res
	const int TMT_labeling = std::stoi(args_map["add_tmt_monomass"].c_str()); //0(no) or 1(yes)

	std::clock_t clock_t_t1;
	clock_t_t1 = std::clock();
	
	// read peptides
	vector<string> peptides;
    vector <string> protein_ids;
	vector <uint32_t> peptide_masses;
    // peptide masses include modifications - masses read from partition
    // because partition is sorted by mass, max_mass_float should be the highest mass in partition
	double max_mass_float = read_peptides_allpeaks(database_peptides_fn, peptides, protein_ids, peptide_masses, TMT_labeling);

	std::cout << "Read peptides from database partitions." << std::endl;
	
	// Compute theoretical b/y ions (create mass-ion index)
    // peptides are strings with mod annotations
    // peptide_masses include modifications in mass
	DatabaseIndex DBI;
	AllPeaks db_data_allpeaks = DBI.retrieve_database_index(peptide_masses, peptides, max_mass_float, MS2_PPM_TOLERANCE);
	
	// read spectra & create ds
	All_spectra_simple spectra_vec;
	vector<uint32_t> ionmass_vec;
	vector<uint32_t> sorted_spectrum_indices;
	std::cout << "Begin reading spectra: " << std::endl;
	read_spectra_allpeaks(mgf_fn, spectra_vec, ionmass_vec, sorted_spectrum_indices, TMT_labeling);

    // remove mods from peptides
    remove_fixedCmod_from_peptides(peptides, TMT_labeling);

	// counting SPC for spectra
	spectra_against_peptides_allpeaks(spectra_vec, ionmass_vec, sorted_spectrum_indices, protein_ids, peptides, peptide_masses, db_data_allpeaks, MS1_PPM_TOLERANCE, MS2_PPM_TOLERANCE, output_fn, score_cutoff);
	
	std::cout << "Core Module 1 (peptide filtering) completed in: " << (float)(std::clock() - clock_t_t1)/CLOCKS_PER_SEC << " seconds. " << std::endl;
	
	return 0;
}

//=== DEFINE FUNCTIONS ===

bool cmp(MassIndexTup lhs, MassIndexTup rhs)
{
	return lhs.first < rhs.first;
}

// parse arguments function
bool parseArguments(int argc, const char* argv[], map<string, string>& args_map)
{
	for (int i = 1; i < argc; i+=2) {
		std::string arg = argv[i];
		if ((arg == "-h") || (arg == "--help")) {
			show_usage(argv[0]);
			return 0;
		}
		else if (arg == "-i") {
			args_map["spectral_partition"] = std::string(argv[i+1]);
			std::cout << "-i "<<args_map["spectral_partition"]<< std::endl;
		}
		else if (arg == "-p") {
			args_map["database_partition"] = std::string(argv[i+1]);
			std::cout << "-p "<<args_map["database_partition"]<< std::endl;
		}
		else if (arg == "-o") {
			args_map["sp_pair_out"] = std::string(argv[i+1]);
			std::cout << "-o "<<args_map["sp_pair_out"]<< std::endl;
		}
		else if (arg == "-s") {
			args_map["spc_score"] = std::string(argv[i+1]);
			std::cout << "-s "<<args_map["spc_score"]<< std::endl;
		}
		else if (arg == "-ms1tol") {
			args_map["precursor_tol"] = std::string(argv[i+1]);
			std::cout << "-ms1tol "<<args_map["precursor_tol"]<< std::endl;
		}
		else if (arg == "-ms2tol") {
			args_map["fragment_tol"] = std::string(argv[i+1]);
			std::cout << "-ms2tol "<<args_map["fragment_tol"]<< std::endl;
		}
		// add tmt parameter
		else if (arg == "-tmt") {
			args_map["add_tmt_monomass"] = std::string(argv[i+1]);
			std::cout << "-tmt "<<args_map["add_tmt_monomass"]<< std::endl;
		}
		else{
			std::cerr << "Did not recognize parameter: "<< arg << std::endl;
			return false;
		}
	}
	return true;
}

// read peptides from database function
double read_peptides_allpeaks(const char* peptides_fn, vector<string> &pep_sequences, vector<string> &protein_ids, vector<uint32_t> &pep_masses, int TMT_labeling)
{
	std::ifstream infile(peptides_fn);
	std::string line;
	double mass;

	// core module 2 (peptide filter) input format:
	//774.3368778	R.SGGRNNNG.-	4790
	//774.35079657	R.APSGDSNK.A	XXX_4791

	while (std::getline(infile, line))
	{
		vector<string> v;
		boost::split(v, line, boost::is_any_of("\t"));
		mass = atof(v[0].c_str());
		pep_masses.push_back(round(mass*MULTIPLIER));
		// v[1] is peptide string with pre post
		string preaa = v[1].substr(0,2);
		string postaa = v[1].substr(v[1].length()-2, v[1].length());
		string peptide_remove_prepost = v[1].substr(2,v[1].length()-4);
		string modified_pep = boost::replace_all_copy(peptide_remove_prepost, "C", "C+57.021");
		// if TMT, add +229.162932 to all K and all n-term
		if(TMT_labeling==1) {
			// modify all K (lysine)
			modified_pep = boost::replace_all_copy(modified_pep, "K", "K+229.163");
			// modify all n-term
			modified_pep = "+229.163"+modified_pep;
		}
		modified_pep = preaa+modified_pep+postaa;

		pep_sequences.push_back(modified_pep);
		protein_ids.push_back(v[2]);
	}
	return mass;
}

//
void remove_fixedCmod_from_peptides(vector<string> &pep_sequences, int TMT_labeling)
{
    for(uint32_t i = 0; i < pep_sequences.size(); i++)
    {
        pep_sequences[i] = boost::replace_all_copy(pep_sequences[i], "C+57.021", "C");
		if(TMT_labeling==1) {
			pep_sequences[i] = boost::replace_all_copy(pep_sequences[i], "+229.163", "");
		}
    }

    return;
}

// read spectra
	/* filtering peaks
	unsigned int filterLowMassPeaks(float minMass);
	* Removes peaks with mass smaller than minMass
	* @param minMass Minimum mass to keep
	* @return Number of removed peaks
	
	int removeChargeReducedPrecursors(const float reducedPrecursorTol = 3.1,
		const bool removeNeutralLoss = true,
		const float minKMedianIntensity = 5.0);
	* Removes charge reduced precursors via Good et. al.′s Algorithm
	* @param reducedPrecursorTol Da tolerance around charge-reduced precursors to look for peaks
	* @param removeNeutralLoss if true also look for neutral-losses from charge-reduced peaks.
	* @param minKMedianIntensity only remove peaks whose intensity is greater than this value times
		the median intensity in the spectrum
	* @return number of peaks that were removed
	
	unsigned int filterLowSNR(float SNR_threshold);
	* Removes peaks with low SNR ratios
	* Noise level is calculated as mean of bottom 25% of peaks
	* @param SNR_threshold Maximum SNR to Keep
	* @return Number of removed peaks
	
	void rankFilterPeaks(int maxRank,
		float windowRadius = AAJumps::minAAmass - 1.0);
	* Removes peaks that have an intensity ranked worse than maxRank
	* compared to all neighboring peaks +/- windowRadius
	* @param maxRank maximum allowable rank of each peak
	* @param windowRadius radius of peak comparison in the spectrum
	* @return
	*/

void read_spectra_allpeaks(const char* mgf_fn, All_spectra_simple &spectra_vec, vector<uint32_t> &ionmass_vec, vector<uint32_t> &sorted_spectrum_indices, int TMT_labeling)
{
	SpecSet spectra;
	spectra.LoadSpecSet_mgf(mgf_fn);

    //filtering peaks in spectra
	for (unsigned int s =0; s < spectra.size(); s++){
        if(TMT_labeling==1) {
            // HCD is a collisional fragmentation method that generates ten unique reporter ions from 126 to 131Da (MS3)
            // some reporter ions may still be generated after CID in MS2
            spectra[s].filterLowMassPeaks(132);
        }

        spectra[s].filterLowSNR(1);
		filter_window(&spectra[s], 75, 7);
		//spectra[s].filterLowSNR(0.8);
	}
	std::cout << "spectra read from spectral partition."<< std::endl;
	
	vector< pair<uint32_t, uint32_t> > mass_spec_index_tuple;
	mass_spec_index_tuple.resize(spectra.size()*(ISOTOPE_NUM+1));
	
	// set up the mass and spectrum index elements.
	spectra_vec.resize(spectra.size()); // peaks in spectra
	ionmass_vec.resize(spectra.size()*(ISOTOPE_NUM+1)); //masses of parent
	sorted_spectrum_indices.resize(spectra.size()*(ISOTOPE_NUM+1)); // spectrum index
	
	for (unsigned int i = 0; i< spectra.size(); i++){
		
		for (int iso= 0; iso < (ISOTOPE_NUM+1); iso++){
			uint32_t index = i*(ISOTOPE_NUM+1)+iso;
			uint32_t mass = round((spectra[i].parentMass - ISOTOPE_EXTRA_MASS*(iso)) * MULTIPLIER);
			mass_spec_index_tuple[index] = std::make_pair(mass, i);
		}

		int peak_list_size = spectra[i].peakList.size();
		A_spectrum_simple spec;

		for (int j = 0; j < peak_list_size; j++){
			// .first is m/z value, .second is intensity
			float p_mass_raw = spectra[i].peakList[j].values.first;
			// DO NOT ACCEPT AN MS2 PEAK IF IT BREACHES THE MASS OF THE PARENT MASS.
			if (p_mass_raw <= spectra[i].parentMass){
				uint32_t p_mass = round(spectra[i].peakList[j].values.first* MULTIPLIER);
				//spec[j] = p_mass;
				spec.push_back(p_mass);
			}
		}
		spectra_vec[i] = spec;
	}
	
	sort(mass_spec_index_tuple.begin(), mass_spec_index_tuple.end());
	for (unsigned int i = 0; i< mass_spec_index_tuple.size(); i++){
		ionmass_vec[i] = mass_spec_index_tuple[i].first; // mass parent
		sorted_spectrum_indices[i] = mass_spec_index_tuple[i].second; // index
	}
	
	return;
	
}

// counting SPC for spectra
void spectra_against_peptides_allpeaks(All_spectra_simple const& spectra_vec, vector<uint32_t> const& ionmass_vec, vector<uint32_t> const& sorted_spectrum_indices, vector<string> &protein_ids, vector<string> const& peptides, vector<uint32_t> const& peptide_masses, AllPeaks const& db_data, int ms1_tolerance, float ms2_tolerance, const char* output_fn, int score_cutoff)
{

	// output file
	std::ofstream output_file;
	output_file.open(output_fn);
	output_file << "#spec_id" << '\t' << "score"<< '\t' << "peptide" << '\t' << "protein" << std::endl;
	// start time
	std::clock_t clock_t_t2;
	clock_t_t2 = std::clock();
	//
	
	vector< vector<short> > picked_scores;
	picked_scores.resize(spectra_vec.size());
	vector< vector<uint32_t> > picked_pep_indices;
	picked_pep_indices.resize(spectra_vec.size());

	uint32_t theoretical_mass_ind_begin = 0;
	uint32_t theoretical_mass_ind_end = 0;

	uint32_t peptide_index = 0;
    // ionmass_vec is spectra parent masses
    // for each spectrum, find SPC
	for (unsigned int i = 0; i < ionmass_vec.size(); i++){
		
		if (i % 2000 == 0){
			std::cout << i << std::endl;
			std::cout<< "took: " << (float)(std::clock() - clock_t_t2)/CLOCKS_PER_SEC << " seconds. " << std::endl;
		}

		uint32_t spec_ms1_tolerance = (ionmass_vec[i]* ms1_tolerance) / pow(10,6);

		uint32_t ion_mass_min = ionmass_vec[i] - spec_ms1_tolerance;
		uint32_t ion_mass_max = ionmass_vec[i] + spec_ms1_tolerance;
		uint32_t spectrum_index = sorted_spectrum_indices[i];

		A_spectrum_simple s1 = spectra_vec[spectrum_index]; // all peaks for spectrum

        // peptide masses include modifications
        // while peptide mass is within tolerance of spectrum's parent mass
        // determine theoretical_mass_ind_begin, theoretical_mass_ind_end
		while ((peptide_masses[peptide_index] <= ion_mass_max) && (peptide_index < peptide_masses.size())){
			if (peptide_masses[peptide_index] < ion_mass_min){
				theoretical_mass_ind_begin = peptide_index;
				theoretical_mass_ind_end = peptide_index;
			}
			else{
				theoretical_mass_ind_end = peptide_index;
			}
			peptide_index++;
		}

        // spectrum_scores is vector of SPCs, where spectrum_scores[i] is the SPC for the ith peptide in the partition
		vector <short> spectrum_scores(theoretical_mass_ind_end-theoretical_mass_ind_begin+1, 0);

		if (peptide_masses[theoretical_mass_ind_begin] < ion_mass_min){
			theoretical_mass_ind_begin += 1;
		}

		if (peptide_masses[theoretical_mass_ind_begin] <= ion_mass_max){
			//Determine the index of spectrum.
            // for each peak in spectrum, find indexes of peptides that have a matching b/y ion
			for (uint32_t S = 0; S < s1.size(); S++){
				uint32_t spectrum_hit_ind =	s1[S]/int(2*ms2_tolerance*MULTIPLIER);
				find_matched_peptide_indices(db_data, spectrum_hit_ind, ion_mass_min, ion_mass_max, theoretical_mass_ind_begin, theoretical_mass_ind_end, peptide_masses, spectrum_scores);
			}

			//Find max of the scores.
			int max_index = max_element(spectrum_scores.begin(), spectrum_scores.end()) - spectrum_scores.begin();
			short max_score = spectrum_scores[max_index];

            // output candidate as long as peptide has SPC >= max score-1 and >= score cutoff (e.g. 7)
            vector<pair<string, uint32_t>> unsorted_candidates;
			for (uint32_t ss = 0; ss < spectrum_scores.size(); ss++){
				if ((max_score > 0) && (spectrum_scores[ss] >= max_score-1) && (spectrum_scores[ss] >= score_cutoff)){
                    unsorted_candidates.push_back(make_pair(peptides[ss + theoretical_mass_ind_begin],ss));
                    //output_file << spectrum_index << '\t' << spectrum_scores[ss] << '\t'
                    //<< peptides[ss + theoretical_mass_ind_begin] << '\t'
                    //<< protein_ids[ss + theoretical_mass_ind_begin] << std::endl;
				}
				//spectrum_scores[ss] = 0;
			}
            //sort by peptide alphabetically, then write to file
            sort(unsorted_candidates.begin(), unsorted_candidates.end());

            for (uint32_t ss_p = 0; ss_p < unsorted_candidates.size(); ss_p++){
                uint32_t ss = unsorted_candidates[ss_p].second;
                output_file << spectrum_index << '\t' << spectrum_scores[ss] << '\t'
                            << peptides[ss + theoretical_mass_ind_begin] << '\t'
                            << protein_ids[ss + theoretical_mass_ind_begin] << std::endl;
                spectrum_scores[ss] = 0;
            }

		}

		// now take the max and initialize to zero.
		if ((theoretical_mass_ind_begin != 0) && (peptide_masses[theoretical_mass_ind_begin-1] < ion_mass_min)){
			theoretical_mass_ind_begin -= 1;
			peptide_index = theoretical_mass_ind_begin;
		}
	}

	std::cout<< "Search finished in: " << (float)(std::clock() - clock_t_t2)/CLOCKS_PER_SEC << " seconds. " << std::endl;
}

void find_matched_peptide_indices(AllPeaks const& db_data, uint32_t spectrum_hit_ind, uint32_t ion_mass_min, uint32_t ion_mass_max, uint32_t theoretical_mass_ind_begin, uint32_t theoretical_mass_ind_end, vector <uint32_t> const& peptide_masses, vector <short> & spectrum_scores )
{
	uint32_t id  = std::lower_bound(db_data[spectrum_hit_ind].begin(),
									db_data[spectrum_hit_ind].end(),
									std::make_pair(ion_mass_min, 0),cmp) - db_data[spectrum_hit_ind].begin();

    // db_data is mass-ion index.
    // all_peaks_data[int_m].push_back(std::make_pair(peptide_masses[i], i))
	while((id < db_data[spectrum_hit_ind].size()) && (db_data[spectrum_hit_ind][id].first <= ion_mass_max) ){
		uint32_t index_to_add = db_data[spectrum_hit_ind][id].second - theoretical_mass_ind_begin;
		uint32_t xx = theoretical_mass_ind_end-theoretical_mass_ind_begin+1;
		uint32_t yy = index_to_add;
		if (xx <= yy){
			std::cout << " SPACE ALLOWED " << xx << std::endl;
			std::cout << " INDEX: " << yy << std::endl;
			std::cout << " BEGIN: " << theoretical_mass_ind_begin << std::endl;
			std::cout << " END: " << theoretical_mass_ind_end << std::endl;
			std::cout << " min cand. mass: " << peptide_masses[theoretical_mass_ind_begin] << std::endl;
			std::cout << " max cand. mass: " << peptide_masses[theoretical_mass_ind_end] << std::endl;
			std::cout << " max cand. mass +1 : " << peptide_masses[theoretical_mass_ind_end+1] << std::endl;
			std::cout << " min actual mass: " << ion_mass_min << std::endl;
			std::cout << " max actual mass: " << ion_mass_max << std::endl;
			std::cout << " current mass: " << db_data[spectrum_hit_ind][id].first << std::endl;
			std::cout << " current mass index: " << db_data[spectrum_hit_ind][id].second << std::endl;
			std::cout << " WHY? : " << peptide_masses[4158950] << std::endl;
		}
		spectrum_scores[index_to_add]++;
		id++;
	}
}
