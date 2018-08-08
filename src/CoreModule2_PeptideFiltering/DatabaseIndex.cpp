/*
 * DatabaseIndex.cpp
 *
 *  Created on: May 2, 2017
 *      Author: dbeyter
  *  Revised by mslin
  *  Includes code from Seungjin Na - theoretical_spectrum.cpp, created April 7, 2016
 */

#include "DatabaseIndex.h"

namespace specnets {

DatabaseIndex::DatabaseIndex() {
	// TODO Auto-generated constructor stub

}

DatabaseIndex::~DatabaseIndex() {
	// TODO Auto-generated destructor stub
}

AllPeaks DatabaseIndex::retrieve_database_index(vector <uint32_t> const&  peptide_masses, vector<string> const& peptides, double max_mass_float,float ms2_tolerance)
{
	//vector <uint32_t> peptide_masses;
//	double max_mass = -1;

//	vector<string> peptide_seqs;
//	vector <uint32_t> peptide_masses;
//	read_peptides_boost(peptides_fn, peptide_seqs, peptide_masses);
//	max_mass = peptide_masses[peptide_masses.size()-1];

	//vector <string> peptide_seqs = read_peptides_boost(peptides_fn, peptide_masses, max_mass);

	AllPeaks all_peaks_data;
	cout << max_mass_float << endl;
	uint32_t max_mass_ind = round((max_mass_float+Proton)*MULTIPLIER)/(2*ms2_tolerance*MULTIPLIER);
	cout << max_mass_ind << endl;
	// initialize
	for (uint32_t p = 0; p <= max_mass_ind; p++)
	{
		all_peaks_data.push_back(PeakSharingMasses({}));
	}
	cout << "Initialization done." << endl;

	vector<float> bions;
	vector<float> yions;
	bions.resize(0);
	yions.resize(0);

    // for each peptide in partition
	for (int i = 0; i < peptides.size(); i ++)
	{
		vector <float> all_peptide_ions;
        string stripped_peptide = peptides[i].substr(2, (peptides[i].length() - 4)); // peptide with modifications and without pre post aa

        //cout << "stripped peptide: " << stripped_peptide << endl;

		getTheoreticalIons(stripped_peptide, 3, bions, yions);
        // getTheoreticalIons(peptides[i], 3, bions, yions);
		all_peptide_ions.reserve( bions.size() + yions.size() );
		all_peptide_ions.insert( all_peptide_ions.end(), bions.begin(), bions.end() );
		all_peptide_ions.insert( all_peptide_ions.end(), yions.begin(), yions.end() );


		if ((i % 100000) == 0)
		{
			cout << i << endl;
		}

		for (int f=0; f < all_peptide_ions.size(); f++)
		{
			uint32_t int_m = round(all_peptide_ions[f]*MULTIPLIER)/int(2*ms2_tolerance*MULTIPLIER);
			all_peaks_data[int_m].push_back(std::make_pair(peptide_masses[i], i));

		}

	}

	return(all_peaks_data);


}
//
//void DatabaseIndex::read_peptides_boost(const char* peptides_fn, vector<string> &pep_sequences, vector<uint32_t> &pep_masses)
//{
//
//	//vector<string> return_peps;
//
//	ifstream infile(peptides_fn);
//	string line;
//	while (std::getline(infile, line))
//	{
//
//		vector<string> v;
//		boost::split(v, line, boost::is_any_of("\t"));
//		//boost::remove_erase_if(v[1], !boost::is_alpha());
//		pep_masses.push_back(round(atof(v[0].c_str())*MULTIPLIER));
//		pep_sequences.push_back(v[1]);
//		//cout << v[1] << endl;
//	}
//
//	return;
//
//}

//
//vector <string> DatabaseIndex::read_peptides(const char* peptides_fn, vector <uint32_t> & masses, double & max_mass)
//{
//
//	vector <string> peptides;
//	char * pch;
//	ifstream infile(peptides_fn);
//	string line;
//	while (std::getline(infile, line))
//	{
//
//		stringstream ss(line);
//		string item;
//		vector<string> tokens;
//		int index = 0;
//		while (getline(ss, item, '\t'))
//		{
//			//tokens.push_back(item);
//			if (index == 0)
//			{
//				double orig_mass = strtod(item.c_str(), 0);
//				uint32_t int_mass =  round(orig_mass*MULTIPLIER);
//				masses.push_back(int_mass);
//
//				if (orig_mass > max_mass)
//				{
//					max_mass = orig_mass;
//				}
//
//
//			}
//			else if (index == 1)
//			{
//				peptides.push_back(item);
//			}
//			index++;
//		}
//
//
//
//	}
//
//	return peptides;
//
//}


//---------------------------------------------------------------------------
/**
 * get theoretical b-and y-ion m/z values
 * string annotation     : msgf+ peptide id (e.g. "DC+57.021M+15.995RTC+57.021GGAK")
 * int charge            : charge state of ms2 spectrum
 * vector<float> & bions : saved b-ion m/z
 * vector<float> & yions : saved y-ion m/z
 */
//---------------------------------------------------------------------------
void DatabaseIndex::getTheoreticalIons(string annotation, int charge, vector<float> & bions, vector<float> & yions)
{
    string  		stripSeq = "";
    int 			modCount = 0;
    string 			delta;
    vector<int> 	msite;
    vector<string> 	msdelta;

    bool prevAA = true;
    for(int i=0; i<annotation.length(); i++){
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
        delta += annotation[i];
      }
    }
    if(!prevAA){
      msdelta.push_back(delta);
    }

    vector<float> gti_mods;
    gti_mods.resize(stripSeq.length());
    if( modCount > 0 ){
      for(int k=0; k<msite.size(); k++){
    	  gti_mods[msite[k]] += atof(msdelta[k].c_str());
      }
    }

    bions.resize(0);
    yions.resize(0);

    //get b-ion series
    float theo_mz = Proton;
    for(int i=0; i<stripSeq.length(); i++){
      theo_mz += org_aa_masses[stripSeq[i]-'A']+gti_mods[i];
      bions.push_back(theo_mz);
    }

    //get y-ion series
    theo_mz = H2Oplus;
    for(int i=stripSeq.length()-1; i>-1; i--){
	  theo_mz += org_aa_masses[stripSeq[i]-'A']+gti_mods[i];
	  yions.push_back(theo_mz);
	}

    if( charge > 2 ){//add multiply-charged ions
      for(int cs=2; cs<charge; cs++){
    	for(int i=0; i<stripSeq.length(); i++){
    	  bions.push_back( (bions[i]+(Proton*cs-1))/cs );
    	  yions.push_back( (yions[i]+(Proton*cs-1))/cs );
    	}
      }
    }


    // sort(bions.begin(), bions.end());
    // sort(yions.begin(), yions.end());

}






} /* namespace specnets */
