/*
 * SpectrumPeptideScoring.h
 *
 *  Created on: Feb 5, 2017
 *      Author: dbeyter
 */

#ifndef SPECTRUMPEPTIDESCORING_H_
#define SPECTRUMPEPTIDESCORING_H_

#include <algorithm>    // std::lower_bound, std::upper_bound, std::sort
#include "spectrum.h"
#include <unistd.h>
#include <string.h>
#include <stdint.h>


typedef pair<int, int> A_peak_tols;
typedef vector<A_peak_tols> A_spectrum;
typedef vector<A_spectrum> All_spectra;



typedef vector <uint32_t> A_spectrum_simple;
typedef vector<A_spectrum_simple> All_spectra_simple;

typedef pair<uint32_t, float> A_peak;
typedef vector <A_peak> A_spectrum_nonsimple;
typedef vector<A_spectrum_nonsimple> All_spectra_nonsimple;

namespace specnets
{


//class Spectrum;

class SpectrumPeptideScoring {
public:
	SpectrumPeptideScoring();
	virtual ~SpectrumPeptideScoring();


	// bisect_left and bisect_right idea.
	float bisect_peak_match(A_spectrum s1, uint32_t* theoretical_peaks, uint32_t begin_index, uint32_t end_index, float ms2_tolerance);
	float merge_peak_match(A_spectrum const& s1, uint32_t* theoretical_peaks, uint32_t begin_index, uint32_t end_index, float ms2_tolerance);
	short merge_peak_count(A_spectrum const& s1, uint32_t* theoretical_peaks, uint32_t begin_index, uint32_t end_index, float ms2_tolerance);
	short merge_peak_count_tagged(A_spectrum const& s1, vector<unsigned int> const& s1_tags, uint32_t* theoretical_peaks, uint32_t begin_index, uint32_t end_index, bool* peptide_tag_presence, uint32_t pep_index, float ms2_tolerance);


};

}

#endif /* SPECTRUMPEPTIDESCORING_H_ */
