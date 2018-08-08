/*
 * DatabaseIndex.h
 *
 *  Created on: May 2, 2017
 *      Author: dbeyter
 *  Revised by mslin
 *  Includes code from Seungjin Na - theoretical_spectrum.cpp, created April 7, 2016
 */

#ifndef DATABASEINDEX_H_
#define DATABASEINDEX_H_

#include <stdio.h>
#include <iostream>     // std::cout
#include <algorithm>    // std::sort
#include <vector>       // std::vector
#include <set>
#include <fstream>
#include <map>
#include <math.h>
#include <time.h>       /* time_t, struct tm, difftime, time, mktime */

#include <libgen.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <sstream>

#include <boost/algorithm/string.hpp>
#include <boost/range/algorithm/remove_if.hpp>
#include <boost/range/algorithm_ext/erase.hpp>


using namespace std;

static int DECIMAL = 4;
static int MULTIPLIER = pow(10, DECIMAL);

static int CHARGE_NUM = 2;

static double	Hydrogen = 1.007825035;
static double	Oxygen = 15.99491463;
static double	H2O = Hydrogen*2 + Oxygen;
static double	Proton = 1.00727649; 
static double	H2Oplus = H2O+Proton;

//a through z
// original aa masses.
static double org_aa_masses[]={
	71.037113805, 0.0, 103.009184505, 115.026943065, 129.042593135,
	147.068413945, 57.021463735, 137.058911875, 113.084064015, 0.0,
	128.09496305, 113.084064015, 131.040484645, 114.04292747, 0.0,
	97.052763875, 128.05857754, 156.10111105, 87.032028435, 101.047678505,
	0.0, 99.068413945, 186.07931298, 0.0, 163.063328575, 0.0};

/*// block cysteine + TMT
static double org_aa_masses[]={
		71.037113805, 0.0, 160.03064824, 115.026943065, 129.042593135,
		147.068413945, 57.021463735, 137.058911875, 113.084064015, 0.0,
		357.257895228, 113.084064015, 131.040484645, 114.04292747, 0.0,
		97.052763875, 128.05857754, 156.10111105, 87.032028435, 101.047678505,
		0.0, 99.068413945, 186.07931298, 0.0, 163.063328575, 0.0};*/

typedef pair <uint32_t, uint32_t> MassIndexTup;
typedef vector <MassIndexTup> PeakSharingMasses;
typedef vector<PeakSharingMasses> AllPeaks;


namespace specnets {

class DatabaseIndex {
public:
	DatabaseIndex();
	virtual ~DatabaseIndex();

	void getTheoreticalIons(string annotation, int charge, vector<float> & bions, vector<float> & yions);
	AllPeaks retrieve_database_index(vector <uint32_t> const&  peptide_masses, vector<string> const& peptides, double max_mass_float, float ms2_tolerance);

};

} /* namespace specnets */

#endif /* DATABASEINDEX_H_ */
