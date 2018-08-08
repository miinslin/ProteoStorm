/*
 * SpectrumPeptideScoring.cpp
 *
 *  Created on: Feb 5, 2017
 *      Author: dbeyter
 */

#include "SpectrumPeptideScoring.h"
using namespace std;


const int TAG_ARRAY_SIZE = 19*3;

namespace specnets
{

SpectrumPeptideScoring::SpectrumPeptideScoring() {
	// TODO Auto-generated constructor stub

}

SpectrumPeptideScoring::~SpectrumPeptideScoring() {
	// TODO Auto-generated destructor stub
}


short SpectrumPeptideScoring::merge_peak_count_tagged(A_spectrum const& s1, vector<unsigned int> const& s1_tags, uint32_t* theoretical_peaks, uint32_t begin_index, uint32_t end_index, bool* peptide_tag_presence, uint32_t pep_index, float ms2_tolerance)
{

	short score = 0;
//	bool tag_exists = false;
//	for (int t = 0; t < s1_tags.size(); t++)
//	{
//		if (peptide_tag_presence[TAG_ARRAY_SIZE*pep_index + s1_tags[t]])
//		{
//			tag_exists = true;
//			break;
//		}
//	}
//
//	if (tag_exists)
//	{
//
//		uint32_t* loc_lower = theoretical_peaks+begin_index;
//		uint32_t* last_element = theoretical_peaks+end_index;
//		for (int i = 0; i < s1.size(); i++)
//		{
//			loc_lower = std::lower_bound(loc_lower, theoretical_peaks+end_index, s1[i].first.first);
//			//loc_upper = std::upper_bound(loc_lower, loc_lower+1, s1[i].first.second);
//
//			if ((loc_lower != last_element) && (*loc_lower <= s1[i].first.second))
//			{
//				score += 1;
//			}
//
//		}
//	}
	return score;






}


short SpectrumPeptideScoring::merge_peak_count(A_spectrum const& s1, uint32_t* theoretical_peaks, uint32_t begin_index, uint32_t end_index, float ms2_tolerance)
{
	short score = 0;
//	uint32_t* loc_lower = theoretical_peaks+begin_index;
//	uint32_t* last_element = theoretical_peaks+end_index;
//	for (int i = 0; i < s1.size(); i++)
//	{
//		loc_lower = std::lower_bound(loc_lower, theoretical_peaks+end_index, s1[i].first.first);
//		//loc_upper = std::upper_bound(loc_lower, loc_lower+1, s1[i].first.second);
//
//		if ((loc_lower != last_element) && (*loc_lower <= s1[i].first.second))
//		{
//			score += 1;
//		}
//
//	}

	return score;
}



float SpectrumPeptideScoring::merge_peak_match(A_spectrum const& s1, uint32_t* theoretical_peaks, uint32_t begin_index, uint32_t end_index, float ms2_tolerance)
{
	float score = 0;
//	uint32_t* loc_lower = theoretical_peaks+begin_index;
//	uint32_t* last_element = theoretical_peaks+end_index;
//	for (int i = 0; i < s1.size(); i++)
//	{
//		loc_lower = std::lower_bound(loc_lower, theoretical_peaks+end_index, s1[i].first.first);
//		if ((loc_lower != last_element) && (*loc_lower <= s1[i].first.second))
//		{
//			score += s1[i].second;
//		}
//	}

	return score;
}





float SpectrumPeptideScoring::bisect_peak_match(A_spectrum s1, uint32_t* theoretical_peaks, uint32_t begin_index, uint32_t end_index, float ms2_tolerance){

//	cout << s1.peakList[0].values.first << endl;
//	cout << theoretical_peaks[0] << endl;

	float score = 0;
//	for (int i = 0; i < s1.size(); i++)
//	{
//
//		//cout << "begin: " << begin_index << " end: " << end_index << endl;
////		vector<uint32_t> theoretical_s(end_index-begin_index+1);
////		for(int p =0; p< theoretical_s.size(); p++)
////		{
////			theoretical_s[p] = theoretical_peaks[p];
////		}
////		std::vector<uint32_t>::iterator loc_lower,loc_upper;
////		loc_lower = std::lower_bound(theoretical_s.begin(), theoretical_s.end(), s1[i].first.first);
////		loc_upper = std::upper_bound(theoretical_s.begin(), theoretical_s.end(), s1[i].first.second);
//
//
//		// min and max for the peak itself.
//		uint32_t loc_lower = *std::lower_bound(theoretical_peaks+begin_index, theoretical_peaks+end_index, s1[i].first.first);
//		uint32_t loc_upper = *std::upper_bound(theoretical_peaks+begin_index, theoretical_peaks+end_index, s1[i].first.second);
//
//		if (loc_upper > loc_lower)
//		{
//			score += s1[i].second;
//		}
//
//	}

	return score;
}





}
