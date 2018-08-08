#include "db_fasta.h"
#include "aminoacid.h"
#include "Logger.h"
#include "alignment_scoring.h"

#include <string.h>

#include <utility>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <limits.h>

#define DEBUG_TAG 0

namespace specnets
{

  char DB_fasta::emptyStr[] = "(empty)";

  DB_fasta::~DB_fasta()
  {
    reset();
  }

  char *DB_fasta::operator[](char *id) const
  {
    for (unsigned int i = 0; i < IDs.size(); i++)
      if (strcmp(id, IDs[i]) == 0)
        return sequences[i];
    return (char *) 0;
  }

  /**
      * Added by Christina Boucher, cboucher@ucsd.edu, Jul 2011
      */
  void DB_fasta::addDecoyReversed()
  {
    unsigned int sizeDb = size();

    // Double size of the DB
    desc.resize(sizeDb * 2);
    sequences.resize(sizeDb * 2);
    IDs.resize(sizeDb * 2);
    masses.resize(sizeDb * 2);

    // Copy all the existing entries
    for (unsigned int i = 0; i < sizeDb; i++) {
      if (desc[i]) {
        desc[i + sizeDb] = new char[strlen(desc[i]) + 1];
        strcpy(desc[i + sizeDb], desc[i]);
      }
      else
        desc[i + sizeDb] = (char *) 0x0;
      if (sequences[i]) {
        int lenStr = strlen(sequences[i]);
        sequences[i + sizeDb] = new char[lenStr + 1];
        // Decoy is reverse of protein
        for (int j = 0; j < lenStr; j++) {
          sequences[i + sizeDb][lenStr - j - 1] = sequences[i][j];
        }
        sequences[i + sizeDb][lenStr] = '\0'; // Null terminate it
      }
      else
        sequences[i] = (char *) 0x0;
      if (IDs[i]) {
        // New ID is old ID + "XXX" at the beginning
        IDs[i + sizeDb] = new char[strlen(IDs[i]) + 4]; // 1 for null, 3 for XXX
        strcpy(IDs[i + sizeDb], "XXX");
        strcpy(&IDs[i + sizeDb][3], IDs[i]);
      }
      else
        IDs[i] = (char *) 0x0;
      masses[i + sizeDb] = masses[i];
    }
    return;
  }

  void DB_fasta::shuffleSequence(char * original, char * shuffled)
  {
    //DEBUG_MSG("original [" << original << "]");
    int lenStr = strlen(original);
    strcpy(shuffled, original);
    srand ( 100 );
    for (int j = 0; j < lenStr; j++) {
      int rndChar = rand() % (lenStr - j) + 1;
      char temp = shuffled[lenStr - j - 1];
      shuffled[lenStr - j - 1] = shuffled[rndChar];
      shuffled[rndChar] = temp;
    }
    shuffled[lenStr] = '\0'; // Null terminate it
    //DEBUG_MSG("shuffled [" << shuffled << "]");
  }

  void DB_fasta::replaceDecoyShuffled()
  {
    unsigned int sizeDb = size();

    // Replace all the existing entries
    for (unsigned int i = 0; i < sizeDb; i++) {
      if (sequences[i]) {
        int lenStr = strlen(sequences[i]);
        char * temp = new char[lenStr + 1];
        shuffleSequence(sequences[i], temp);
        strcpy(sequences[i], temp);
        sequences[i][lenStr] = '\0'; // Null terminate it
        delete[] temp;
      }
      if (IDs[i]) {
        // New ID is old ID + "XXX" at the beginning
        char * decoy = new char[strlen(IDs[i]) + 4]; // 1 for null, 3 for XXX
        strcpy(decoy, "XXX");
        strcpy(&decoy[3], IDs[i]);
        char * old = IDs[i]; // save the old string
        IDs[i] = decoy;      // replace with new
        delete[] old;        // free memory from old
      }
    }
    return;
  }
  
  void DB_fasta::addDecoyShuffled()
  {
    unsigned int sizeDb = size();

    // Double size of the DB
    desc.resize(sizeDb * 2);
    sequences.resize(sizeDb * 2);
    IDs.resize(sizeDb * 2);
    masses.resize(sizeDb * 2);

    // Copy all the existing entries
    for (unsigned int i = 0; i < sizeDb; i++) {
      if (desc[i]) {
        desc[i + sizeDb] = new char[strlen(desc[i]) + 1];
        strcpy(desc[i + sizeDb], desc[i]);
      } else {
        desc[i + sizeDb] = (char *) 0x0;
      }
      if (sequences[i]) {
        int lenStr = strlen(sequences[i]);
        sequences[i + sizeDb] = new char[lenStr + 1];
        shuffleSequence(sequences[i], sequences[i + sizeDb]);
      } else {
        sequences[i] = (char *) 0x0;
      }
      if (IDs[i]) {
        // New ID is old ID + "XXX" at the beginning
        IDs[i + sizeDb] = new char[strlen(IDs[i]) + 4]; // 1 for null, 3 for XXX
        strcpy(IDs[i + sizeDb], "XXX");
        strcpy(&IDs[i + sizeDb][3], IDs[i]);
      } else {
        IDs[i] = (char *) 0x0;
        masses[i + sizeDb] = masses[i];
      }
    }
    return;
  }

  /**
    * Added by Christina Boucher, cboucher@ucsd.edu, Jul 2011
    */
  void DB_fasta::populateIndex(bool printout)
  {
	  massesIndex.resize(3000);

	   for(int i = 0; i < desc.size(); i++)
	   {
		//   DEBUG_VAR(i);
		//   DEBUG_VAR(desc[i]);
		//   DEBUG_VAR(sequences[i]);

		   Spectrum sp = getMassesSpec(i);
		   string peptide = string(sequences[i]);

		   for(int startPosition = 0; startPosition <= masses[i].size(); startPosition++)
		   {
	 			  for(int stopPosition = startPosition + 1; stopPosition < masses[i].size() && stopPosition < startPosition + 30; stopPosition++)
	 			  {
	 				  //int mass = (int) round ((masses[i].peakList[stopPosition][0] - masses[i].peakList[startPosition ][0] )*0.9995);
	 				  string partial_peptide = peptide.substr(startPosition, stopPosition - startPosition);
	 				  string str2("X");
	 				  size_t found=partial_peptide.find(str2);

	 				  if(found==string::npos)
	 				  {
	 					  AAJumps myjumps(1);
	 					  vector<float> mymasses;
	 					  float mass_float = myjumps.getPeptideMass(partial_peptide);
	 					  int mass_int = (int)round(mass_float*0.9995);

	 					  if(mass_int < 3000)
	 					  {
	 						  string str2("X");
	 						  size_t found=partial_peptide.find(str2);
	 						  pair<int, int> p = make_pair(i, startPosition );
	 						  massesIndex[mass_int].push_back( p );
	 					  }
	 				  }
	 			  }
	 		  }
	 	  }
  }

  /**
    * Added by Christina Boucher, cboucher@ucsd.edu, Jul 2011
    */
  vector<string> DB_fasta::findPeptideForSingleMass(int m)
  {
	  if(!peptideIndices.empty())
		  peptideIndices.clear();

	  vector<string> retvec;
	  if(m >= 3000 || m == 0)
		  return retvec;

	  if(massesIndex[m].empty())
		  return retvec;

	  vector< pair< int, int> > final = massesIndex[m];

	  for(int i = 0; i < final.size(); i++)
	  {
		  string peptide = sequences[final[i].first];

		  // find stop position
		  for(int p = final[i].second + 1; p < final[i].second + 30; p++ )
		  {
			  AAJumps myjumps(1);
			  vector<float> mymasses;
			  string partial_peptide = peptide.substr(final[i].second, p - final[i].second );
			  float mass = myjumps.getPeptideMass(partial_peptide);
			  int mass_int = (int) round(mass * 0.9995);
			  if(mass_int == (m))
			  {
				  peptideIndices.push_back(final[i].first);
				  retvec.push_back(partial_peptide);
			  }
		  }
	  }
	  std::sort(retvec.begin(), retvec.end());
	  retvec.erase(unique(retvec.begin(), retvec.end()), retvec.end());
	  return retvec;
  }

  /**
    * Added by Christina Boucher, cboucher@ucsd.edu, Jul 2011
    */
  vector<string> DB_fasta::findPeptides(int m1, int m2, int m3)
  {

	  if(!peptideIndices.empty())
		  peptideIndices.clear();
	  vector<string> retvec;

	  // check if all the masses are equal to zero.
	  if((m1 == 0 && m2 == 0 && m3 == 0) || m1 >= 3000 || m2 >= 3000 || m3 >= 3000)
		  return retvec;


	  // check if two of the masses are equal to zero.
	  if(checkSingleton(m1, m2, m3))
	  {

		  int m = returnSingleton(m1, m2, m3);


		  if(massesIndex[m].empty())
			  return retvec;

		  vector< pair< int, int> > final = massesIndex[m];

		  for(int i = 0; i < final.size(); i++)
		  {
			  string peptide = sequences[final[i].first];

			  // find stop position
			  for(int p = final[i].second + 1; p < final[i].second + 30; p++ )
			  {
				  AAJumps myjumps(1);
				  vector<float> mymasses;
				  string partial_peptide = peptide.substr(final[i].second, p - final[i].second );
				  float mass = myjumps.getPeptideMass(partial_peptide);
			 	  int mass_int = (int) round(mass * 0.9995);
			 	  if(mass_int == (m))
			 	  {
			 		  peptideIndices.push_back(final[i].first);
			 		  retvec.push_back(partial_peptide);
			 	  }
			  }
		  }

		  std::sort(retvec.begin(), retvec.end());
		  retvec.erase(unique(retvec.begin(), retvec.end()), retvec.end());

		  return retvec;
	  }

	  // check if one of the masses is equal to zero.
	  if(m1 == 0 || m2 == 0 || m3 == 0)
	  {
		  if(m1 == 0)
		  {
			  m1 = m2;
			  m2 = m3;
		  }
		  else if(m2 == 0)
		  {
			  m2 = m3;
		  }

		  vector< pair< int, int> > final = intersectSets(massesIndex[m1], massesIndex[m1 + m2]);

		  for(int i = 0; i < final.size(); i++)
		  {
			  string peptide = sequences[final[i].first];

			  // find stop position
			  for(int p = final[i].second + 1; p < final[i].second + 30; p++ )
			  {
				  AAJumps myjumps(1);
				  vector<float> mymasses;
				  string partial_peptide = peptide.substr(final[i].second, p - final[i].second );
				  float mass = myjumps.getPeptideMass(partial_peptide);
			 	  int mass_int = (int) round(mass * 0.9995);
			 	  if(mass_int == (m1 + m2))
			 	  {
			 		  retvec.push_back(partial_peptide);
			 		 peptideIndices.push_back(final[i].first);
			 	  }
			  }
		  }

		  std::sort(retvec.begin(), retvec.end());
		  retvec.erase(unique(retvec.begin(), retvec.end()), retvec.end());

		  return retvec;
	  }

	  if(massesIndex[m1].empty() || massesIndex[m1 + m2].empty()  || massesIndex[m1 + m2 + m3].empty())
	  {
	//	  DEBUG_MSG("One set is empty");
		  return retvec;
	  }
	  vector< pair< int, int> > intSet_m1_m1m2 = intersectSets(massesIndex[m1], massesIndex[m1 + m2]);
	  vector< pair< int, int> > intSet_m1_m1m2m3 = intersectSets(massesIndex[m1], massesIndex[m1 + m2 + m3]);
	  vector< pair< int, int> > intSet_m1m2_m1m2m3 = intersectSets(massesIndex[m1 + m2], massesIndex[m1 + m2 + m3]);

	  if(intSet_m1_m1m2.empty() || intSet_m1_m1m2m3.empty() || intSet_m1m2_m1m2m3.empty())
	  {
		//  DEBUG_MSG("One of the intersection sets is empty");
		  return retvec;
	  }
	  vector< pair< int, int> > semiFinalSet = intersectSets(intSet_m1_m1m2, intSet_m1_m1m2m3);
	  vector< pair< int, int> > final = intersectSets(semiFinalSet, intSet_m1m2_m1m2m3);

	  if(semiFinalSet.empty() || final.empty() )
	  {
	// 		  DEBUG_MSG("One of the second last intersection sets is empty");
	 		  return retvec;
	 	  }
	  for(int i = 0; i < final.size(); i++)
	  {
		  string peptide = sequences[final[i].first];

		  // find stop position
		  for(int p = final[i].second + 1; p < final[i].second + 30; p++ )
		  {
			  AAJumps myjumps(1);
			  vector<float> mymasses;
			  string partial_peptide = peptide.substr(final[i].second, p - final[i].second );
			  float mass = myjumps.getPeptideMass(partial_peptide);
		 	  int mass_int = (int) round(mass * 0.9995);
		 	  if(mass_int == (m1 + m2 + m3))
		 	  {
		 		 peptideIndices.push_back(final[i].first);
		 		 retvec.push_back(partial_peptide);
		 	  }
		  }
	  }

	  std::sort(retvec.begin(), retvec.end());
	  retvec.erase(unique(retvec.begin(), retvec.end()), retvec.end());

	  return retvec;
  }

  /**
    * Added by Christina Boucher, cboucher@ucsd.edu, Jul 2011
    */
  vector<string> DB_fasta::findPeptidesForTriple(int m11, int m12, int m13, int m21, int m22, int m23)
  {
	  vector<string> retvec;

	  if(m11 >= 3000 || m12 >= 3000 || m13 >= 3000 || m21 >= 3000 || m22 >= 3000 || m23 >= 3000 )
		  return retvec;

	  vector<string> set1 = findPeptides(m11, m12, m13);
	  vector<int> s1Indices = peptideIndices;
	  vector<string> set2 = findPeptides(m21, m22, m23);
	  vector<int> s2Indices = peptideIndices;

	  if(!peptideIndices.empty())
	  		  peptideIndices.clear();

	  for(int s1 = 0; s1 < set1.size(); s1++ )
	  {
		  string peptide = set1[s1];
		  if(std::find(set2.begin(), set2.end(), peptide) != set2.end() )
		  {
			  peptideIndices.push_back(s1Indices[s1]);
			  retvec.push_back(peptide);
		  }
//		  DEBUG_MSG(peptide);
	  }

	  return retvec;
  }


  /**
   * Added by Christina Boucher, cboucher@ucsd.edu, Jul 2011
   */
  vector<string> DB_fasta::findPeptideMultiDimIndex(int firstShift, int secondShift, int thirdShift)
   {

	  vector<string> retvec;

	  if(firstShift >= 3000 || secondShift >= 3000 || thirdShift >= 3000)
	  {
		DEBUG_MSG("Array out of bounds.");
	  }

	  if(massOffsets[firstShift].empty() || massOffsets[secondShift].empty())
		  return retvec;

	  it = massOffsets[firstShift].find(secondShift);
	  vector<pair<int, int> > intvec1 = it->second;
	  map<int, vector<pair<int, int> > >::iterator it2 = massOffsets[secondShift].find(thirdShift);
	  vector<pair<int, int> > intvec2 = it2->second;


	  for(int vecindex1 = 0; vecindex1 < intvec1.size(); vecindex1++)
	  {

		  for(int vecindex2 = 0; vecindex2 < intvec2.size(); vecindex2++)
		  {


			  if(intvec1[vecindex1].first == intvec2[vecindex2].first )
			  {
					  string s = sequences[intvec1[vecindex1].first];
					  string str = s.substr(intvec1[vecindex1].second, 30);

					  // increment through mass 1
					  vector<char> sequence = vector<char> (s.begin(), s.end());
					  vector<float> masses;
					  getMasses(sequence, masses);
					  int mass = 0;
					  int k = 0;
					  for (k = 0; k < masses.size(); k++)
					  {
						  if(mass == firstShift)
							  break;
						  mass += (int) round(masses[k] * 0.9995);
					  }

					  if(intvec1[vecindex1].second + k + 1 == intvec2[vecindex2].second)
					  {
						  mass = 0;
						  int kk;
						  for (kk = k; k < masses.size(); k++)
						  {
							  if(mass == secondShift + thirdShift)
								  break;
							  mass += (int) round(masses[kk] * 0.9995);
						  }
						  string addStr = s.substr(intvec1[vecindex1].second, kk + 1 - intvec1[vecindex1].second);
						  retvec.push_back(addStr);
					  }

			  }
		  }
	  }

	  // remove duplicates
	  if(!retvec.empty())
	  {
		   std::sort(retvec.begin(), retvec.end());
		   retvec.erase(unique(retvec.begin(), retvec.end()), retvec.end());
	  }
	  return retvec;
   }


  DB_fasta &DB_fasta::operator=(const DB_fasta &other)
  {
    reset();
    unsigned int sizeDB = other.size();

    desc.resize(sizeDB);
    sequences.resize(sizeDB);
    IDs.resize(sizeDB);
    masses.resize(sizeDB);
    for (unsigned int i = 0; i < sizeDB; i++)
    {
      if (other.desc[i]) {
        desc[i] = new char[strlen(other.desc[i]) + 1];
        strcpy(desc[i], other.desc[i]);
      }
      else
      {
        desc[i] = (char *) 0x0;
      }

      if (other.sequences[i]) {
        sequences[i] = new char[strlen(other.sequences[i]) + 1];
        strcpy(sequences[i], other.sequences[i]);
      }
      else
      {
        sequences[i] = (char *) 0x0;
      }

      if (other.IDs[i]) {
        IDs[i] = new char[strlen(other.IDs[i]) + 1];
        strcpy(IDs[i], other.IDs[i]);
      }
      else
      {
        IDs[i] = (char *) 0x0;
      }
      masses[i] = other.masses[i];
    }

    return *this;
  }

  void DB_fasta::replaceAA(char prevAA, char repAA)
  {
    for (unsigned int seqIdx = 0; seqIdx < sequences.size(); seqIdx++) {
      unsigned int aaIdx = 0;
      while (sequences[seqIdx][aaIdx]) {
        if (sequences[seqIdx][aaIdx] == prevAA)
          sequences[seqIdx][aaIdx] = repAA;
        aaIdx++;
      }
    }
  }

  Spectrum &DB_fasta::getMassesSpec(int index)
  {
    vector<float> tmp;
    if (masses[index].size() == 0) {
      getMassesIdx(index, tmp);
      masses[index].resize(tmp.size() + 1);
      masses[index][0].set((float) 0, (float) 0); // Zero based cumsum of amino acid masses
      for (unsigned int i = 0; i < tmp.size(); i++)
        masses[index][i + 1].set(tmp[i] + masses[index][i][0], 0);
      masses[index].parentMass = masses[index][masses[index].size() - 1][0];
      masses[index].parentCharge = 1;
    }
    return masses[index];
  }

  void DB_fasta::reset()
  {
    for (unsigned int seqIdx = 0; seqIdx < IDs.size(); seqIdx++) {
      delete[] IDs[seqIdx];
      if (desc[seqIdx] != emptyStr) {
        delete[] desc[seqIdx];
      }
      delete[] sequences[seqIdx];
    }
  }


  void DB_fasta::Add(char *ID, char *des, char *seq)
  {
    string aux(ID);
    for(int i = 0 ; i < aux.length() ; i++)
      if(aux[i] == ' ')
        aux[i] = '|';

    char * aux1 = new char[ strlen(aux.c_str()) + 1];
    char * aux2 = new char[ strlen(des) + 1];
    char * aux3 = new char[ strlen(seq) + 1];

    strcpy(aux1, aux.c_str());
    strcpy(aux2, des);
    strcpy(aux3, seq);


    IDs.push_back(aux1);
    desc.push_back(aux2);
    sequences.push_back(aux3);
    Spectrum temp;
    masses.push_back(temp);
  }



  unsigned int DB_fasta::Load(const char *filename)
  {
    BufferedLineReader blr;
    unsigned int seqIdx, lineIdx, seqSize;
    char *token;
    list<unsigned int> seqSizes;

    if (blr.Load(filename) <= 0)
      return 0;

    AAJumps jumps(1);

    // Read all the lines in the file
    seqIdx = 0;
    seqSize = 0; // tallies up total sequence size from all lines
    for (lineIdx = 0; lineIdx < blr.size(); lineIdx++)
      if (blr.getline(lineIdx)[0] == '>') {
        seqIdx++;
        if (seqIdx > 1)
          seqSizes.push_back(seqSize);
        seqSize = 0;
      }
      else if (blr.getline(lineIdx)[0] != ';')
        seqSize += strlen(blr.getline(lineIdx));
    seqSizes.push_back(seqSize);

    // Initialize data structures
    reset();
    IDs.resize(seqIdx);
    sequences.resize(seqIdx);
    desc.resize(seqIdx);
    masses.resize(seqIdx);
    for (seqIdx = 0; seqIdx < desc.size(); seqIdx++) {
      desc[seqIdx] = emptyStr;
      IDs[seqIdx] = emptyStr;
      sequences[seqIdx] = emptyStr;
    }

    // Load sequences
    seqIdx = 0;
    for (lineIdx = 0; lineIdx < blr.size() and blr.getline(lineIdx)[0] != '>'; lineIdx++)
      ;
    for (; lineIdx < blr.size(); lineIdx++) {

      // Comment line starts with '>'
      if (blr.getline(lineIdx)[0] == '>') {

        sequences[seqIdx] = new char[seqSizes.front() + 1];
        if (!sequences[seqIdx]) {
          cerr << "ERROR: Could not allocate " << seqSizes.front() + 1
              << " bytes!\n";
          return 0;
        }
        sequences[seqIdx][0] = 0;
        seqSizes.pop_front();

        token = strtok(&blr.getline(lineIdx)[1], " ");
        if (token) {
          IDs[seqIdx] = new char[strlen(token) + 1];
          if (!IDs[seqIdx]) {
            cerr << "ERROR: Could not allocate " << strlen(token) + 1
                << " bytes!\n";
            return 0;
          }
          strcpy(IDs[seqIdx], token);
        } else {
          IDs[seqIdx] = emptyStr;
        }
        token = strtok(NULL, "");
        if (token) {
          desc[seqIdx] = new char[strlen(token) + 1];
          if (!desc[seqIdx]) {
            cerr << "ERROR: Could not allocate " << strlen(token) + 1
                << " bytes!\n";
            return 0;
          }
          strcpy(desc[seqIdx], token);
        } else {
          desc[seqIdx] = emptyStr;
        }
        seqSize = 0;
        seqIdx++;
      } else if (seqIdx > 0 and blr.getline(lineIdx)[0] != 0 &&
                 blr.getline(lineIdx)[0] != ';') {
        int length = strlen(blr.getline(lineIdx));
        char * temp = new char[length + 1];
        if (temp == 0x0) {
          cerr << "ERROR: Could not allocate " << length + 1
               << " bytes for temp string!\n";
          return 0;
        }
        strcpy(temp, blr.getline(lineIdx));
        string tempChar = "X";
        size_t k = 0;
        for (size_t i = 0; i < length; i++) {
           tempChar[0] = temp[i];
           if (jumps.massLookup.find(tempChar) != jumps.massLookup.end()) {
             temp[k] = temp[i];
             k++;
           } else {
             static bool beenWarned = false;
             if (!beenWarned) {
               WARN_MSG("Removed unknown AA [" << tempChar << "] on line [" << lineIdx + 1<< "] of " << filename);
               WARN_MSG("  Further warnings of this type have been supressed");
               beenWarned = true;
             }
           }
        }
        temp[k] = '\0'; // put null on end
        strcpy(&(sequences[seqIdx - 1][seqSize]), temp);
        seqSize += strlen(temp);
	      delete[] temp;
      }
    }

    return (seqIdx);
  }



  unsigned int DB_fasta::Save(const char * filename)
  {
    FILE * pepsFID = (FILE *) fopen(filename, "wb");
    if (pepsFID == 0x0) {
      return -1;
    }
    unsigned int count = 0;
    for (unsigned int index = 0; index < IDs.size(); index++) {
      fprintf(pepsFID,
              ">%s %s\n%s\n",
              IDs[index],
              desc[index],
              sequences[index]);
      count++;
    }
    fclose(pepsFID);

    return count;
  }




  unsigned int DB_fasta::Save(const char * filename,
                              vector<bool> & saveFlags)
  {
    FILE * pepsFID = (FILE *) fopen(filename, "wb");
    if (pepsFID == 0x0) {
      return -1;
    }
    unsigned int count = 0;
    for (unsigned int index = 0; index < saveFlags.size(); index++) {
      if (saveFlags[index]) {
        fprintf(pepsFID,
                ">%s %s\n%s\n",
                IDs[index],
                desc[index],
                sequences[index]);
        count++;
      }
    }
    fclose(pepsFID);

    return count;
  }

  void DB_fasta::output(ostream &out)
  {
    for (unsigned int i = 0; i < IDs.size(); i++)
      std::cout << IDs[i] << ":" << desc[i] << "\n" << sequences[i] << "\n";
  }

  unsigned int DB_fasta::find(const char *peptide,
                              list<sps::tuple<int,float,string> > &matches)
  {
    matches.clear();
    if (peptide == (char *) 0 or strlen(peptide) == 0)
      return 0;

    sps::tuple<int,float,string> curTuple;
    curTuple.m2.resize(strlen(peptide));
    for (unsigned int aaIdx = 0; aaIdx < curTuple.m2.size(); aaIdx++)
      curTuple.m2[aaIdx] = peptide[aaIdx];

    for (curTuple.m0 = 0; curTuple.m0 < sequences.size(); curTuple.m0++) {
      if (strstr(sequences[curTuple.m0], peptide)) {
        matches.push_back(curTuple);
      }
    }
    return matches.size();
  }

  unsigned int DB_fasta::find(const string &peptide,
                              list<sps::tuple<int,float,string> > &matches)
  {
    matches.clear();

    if (peptide.compare("") == 0 )
    {
      return 0;
    }

    AAJumps jumps(1);

    for (int dbIndex = 0; dbIndex < sequences.size(); dbIndex++)
    {
      string currProtein(sequences[dbIndex]);
      size_t found = currProtein.find(peptide);
      while (found != string::npos)
      {
        string prefixPeptide = currProtein.substr(0,found);
        //DEBUG_VAR(prefixPeptide);
        float startMass = jumps.getPeptideMass(prefixPeptide);
        //DEBUG_VAR(startMass);
        sps::tuple<int,float,string> currHit;
        currHit.m0 = dbIndex;
        currHit.m1 = startMass;
        currHit.m2 = peptide;
        matches.push_back(currHit);
        found = currProtein.find(peptide, found +1);
      }
    }
    return matches.size();
  }
  // ***************************************************
  // ***  Protein alignment methods
  // ***************************************************

  bool Load_clustalw(const char *filename,
                     DB_fasta &db,
                     int protID1,
                     int protID2,
                     vector<TwoValues<int> > &matchedIndices)
  {
    BufferedLineReader blr;
    matchedIndices.resize(0);
    char *matchedP1, *matchedP2; // Alignment string for protein P1/P2
    char *token;
    unsigned int lineIdx, pivot, szP1, szP2, idxMatch, idxP1, idxP2;

    if (min(protID1, protID2) < 0 or max(protID1, protID2) >= db.size())
      return false;
    if (blr.Load(filename) <= 0) {
      cerr << "ERROR loading " << filename << "!\n";
      return false;
    }
    szP1 = strlen(db[protID1]);
    szP2 = strlen(db[protID2]);
    matchedP1 = new char[1 + szP1 + szP2];
    matchedP1[0] = 0;
    matchedP2 = new char[1 + szP1 + szP2];
    matchedP2[0] = 0;

    // Merge the alignment strings into matchedP1/matchedP2
    for (lineIdx = 3; lineIdx < blr.size(); lineIdx++) {
      token = strtok(blr.getline(lineIdx++), " ");
      token = strtok(NULL, " ");
      strcat(matchedP1, token);

      if (lineIdx >= blr.size()) {
        cerr << "ERROR: Invalid clustalw .aln format!\n";
        delete[] matchedP1;
        delete[] matchedP2;
        return false;
      }
      token = strtok(blr.getline(lineIdx++), " ");
      token = strtok(NULL, " ");
      strcat(matchedP2, token);

      lineIdx++;
    }

    // Find matching aa-breaks (b0/N-Term = index 0)
    szP1 = strlen(matchedP1);
    szP2 = strlen(matchedP2);
    if (szP1 != szP2) {
      cerr << "ERROR: Invalid clustalw .aln format!\n";
      delete[] matchedP1;
      delete[] matchedP2;
      return false;
    }
    matchedIndices.resize(1 + szP1);
    idxMatch = 0;
    idxP1 = 0;
    idxP2 = 0;
    for (pivot = 0; pivot < szP1; pivot++) {
      if (matchedP1[pivot] == '-') {
        idxP2++;
        continue;
      }
      if (matchedP2[pivot] == '-') {
        idxP1++;
        continue;
      }
      matchedIndices[idxMatch++].set((int) idxP1++, (int) idxP2++);
    }
    matchedIndices[idxMatch].set(matchedIndices[idxMatch - 1][0] + 1,
                                 matchedIndices[idxMatch - 1][1] + 1); // match smallest bn at the end of the alignment
    matchedIndices.resize(idxMatch + 1);
    delete[] matchedP1;
    delete[] matchedP2;
    return true;
  }

  //
  // Load_clustalw_multiple - Loads multiple clustalw alignments. Aligned protein indices and clustalw .aln filenames
  //      are read from the text file filename with one clustalw alignment per line:
  //
  //        protein_idx_1;protein_idx_2;clustalw_aln_filename
  //
  //   db           - database containing the sequences of the aligned proteins
  //   matchedProts - indices of the matched proteins
  //   matchedAAs   - amino acid correspondences between matched proteins
  //
  bool Load_clustalw_multiple(const char *filename,
                              DB_fasta &db,
                              vector<TwoValues<int> > &matchedProts,
                              vector<vector<TwoValues<int> > > &matchedAAs,
                              unsigned int proteinIndexOffset)
  {
    BufferedLineReader blr;
    unsigned int idxP1, idxP2, lineIdx;
    char *token;

    if (blr.Load(filename) < 0) {
      cerr << "ERROR loading " << filename << "!\n";
      return false;
    }
    matchedAAs.resize(blr.size());
    matchedProts.resize(blr.size());

    DEBUG_VAR(blr.size());

    for (lineIdx = 0; lineIdx < blr.size() and strlen(blr.getline(lineIdx)) > 0; lineIdx++) {
      idxP1 = 0;
      idxP2 = 0;
      token = strtok(blr.getline(lineIdx), ";");
      if (token)
        idxP1 = (int) atoi(token);
      idxP1 -= proteinIndexOffset; // Convert to zero-based protein indices
      if (idxP1 >= db.size()) {
        cerr << "ERROR reading " << filename << ": invalid protein index "
            << idxP1 << " in line " << lineIdx << "!\n";
        return false;
      }

      token = strtok(NULL, ";");
      if (token)
        idxP2 = (int) atoi(token);
      idxP2 -= proteinIndexOffset; // Convert to zero-based protein indices
      if (idxP2 >= db.size()) {
        cerr << "ERROR reading " << filename << ": invalid protein index "
            << idxP2 << " in line " << lineIdx << "!\n";
        return false;
      }

      token = strtok(NULL, ";");
      if (!token) {
        cerr << "ERROR reading " << filename << ": invalid filename in line "
            << lineIdx << "!\n";
        return false;
      }
      matchedProts[lineIdx].set((int) idxP1, (int) idxP2);
      if (not Load_clustalw(token, db, idxP1, idxP2, matchedAAs[lineIdx]))
        return false;
    }
    matchedAAs.resize(lineIdx);
    matchedProts.resize(lineIdx);
    return true;
  }

  // ***************************************************
  // ***  DB_index methods
  // ***************************************************

  static class Text2aa
  { // Table to convert ascii codes (0-255) to amino acid indices (0-17)
  public:
    vector<char> convTable;
    Text2aa()
    {
      convTable.resize(256);
      for (unsigned int idx = 0; idx < 256; idx++)
        convTable[idx] = (char) 0;
      convTable[65] = (char) 0;
      convTable[97] = (char) 0; // 'A', 'a'
      convTable[67] = (char) 1;
      convTable[99] = (char) 1; // 'C', 'c'
      convTable[68] = (char) 2;
      convTable[100] = (char) 2; // 'D', 'd'
      convTable[69] = (char) 3;
      convTable[101] = (char) 3; // 'E', 'e'
      convTable[70] = (char) 4;
      convTable[102] = (char) 4; // 'F', 'f'
      convTable[71] = (char) 5;
      convTable[103] = (char) 5; // 'G', 'g'
      convTable[72] = (char) 6;
      convTable[104] = (char) 6; // 'H', 'h'
      convTable[73] = (char) 7;
      convTable[105] = (char) 7; // 'I', 'i'
      convTable[75] = (char) 8;
      convTable[107] = (char) 8; // 'K', 'k'
      convTable[76] = (char) 7;
      convTable[108] = (char) 7; // 'L', 'l'
      convTable[77] = (char) 9;
      convTable[109] = (char) 9; // 'M', 'm'
      convTable[78] = (char) 10;
      convTable[110] = (char) 10; // 'N', 'n'
      convTable[80] = (char) 11;
      convTable[112] = (char) 11; // 'P', 'p'
      convTable[81] = (char) 9;
      convTable[113] = (char) 9; // 'Q', 'q'
      convTable[82] = (char) 12;
      convTable[114] = (char) 12; // 'R', 'r'
      convTable[83] = (char) 13;
      convTable[115] = (char) 13; // 'S', 's'
      convTable[84] = (char) 14;
      convTable[116] = (char) 14; // 'T', 't'
      convTable[86] = (char) 15;
      convTable[118] = (char) 15; // 'V', 'v'
      convTable[87] = (char) 16;
      convTable[119] = (char) 16; // 'W', 'w'
      convTable[89] = (char) 17;
      convTable[121] = (char) 17; // 'Y', 'y'
    }
    char &operator[](int idx)
    {
      return convTable[idx];
    }
  } text2aa;

  void DB_index::hash_init(int tagLength)
  {
    coeffs1.resize(tagLength);
    coeffs2.resize(tagLength);
    srand( time(NULL));
    for (int cIdx = 0; cIdx < tagLength; cIdx++) {
      coeffs1[cIdx] = rand() % SHRT_MAX;
      coeffs2[cIdx] = rand() % SHRT_MAX;
    }
  }

  int DB_index::hash(int hashIdx, const char *tag)
  {
    int index = 0, cIdx;
    if (hashIdx == 1) {
      for (cIdx = 0; cIdx < tagLength; cIdx++)
        index += (coeffs1[cIdx] * ((int) text2aa[tag[cIdx]])) % indexDim1;
      return (int) (index % indexDim1);
    }
    else {
      for (cIdx = 0; cIdx < tagLength; cIdx++)
        index += (coeffs2[cIdx] * ((int) text2aa[tag[cIdx]])) % indexDim2;
      return (int) (index % indexDim2);
    }
  }

  DB_index::DB_index(DB_fasta &db,
                     int newIndexDim1,
                     int newIndexDim2,
                     int newTagLength)
  {
    buildIndex(db, newIndexDim1, newIndexDim2, newTagLength);
  }

  void DB_index::buildIndex(DB_fasta &db,
                            int newIndexDim1,
                            int newIndexDim2,
                            int newTagLength)
  {
    vector < vector<list<Tag> > > tmpIndex; // Temporary index (with lists instead of vectors)
    int c1, c2; // Index coordinates 1, 2
    pair<int, int> tmpEntry; // Temporary pair to hold (proteinID,tagPos) values for insertion
    char *protSeq; // Pointer to current protein sequence
    int lastTagPos, tagIndex; // Position of last tag in current protein sequence and tag index
    list<Tag>::iterator tagIter; // Iterator of tags in a collision list
    Tag curTag; // Holder for current tag being processed

    indexDim1 = newIndexDim1;
    indexDim2 = newIndexDim2;
    tagLength = newTagLength;
    hash_init(newTagLength); //Initialize hash function coefficients for current tag length

    // Build new temporary index (with lists instead of vectors)
    tmpIndex.resize(indexDim1);
    for (c1 = 0; c1 < indexDim1; c1++) {
      tmpIndex[c1].resize(indexDim2);
      for (c2 = 0; c2 < indexDim2; c2++)
        tmpIndex[c1][c2].clear();
    }

    curTag.text = new char[sizeof(char) * (tagLength + 1)];
    curTag.text[tagLength] = char(0);
    for (int protIdx = 0; protIdx < (int) db.size(); protIdx++) {

      protSeq = db[protIdx];
      lastTagPos = strlen(protSeq) - tagLength;
      tmpEntry.first = protIdx;
      for (tagIndex = 0; tagIndex <= lastTagPos; tagIndex++) {
        strncpy(curTag.text, &protSeq[tagIndex], tagLength);
        c1 = hash(1, curTag.text);
        c2 = hash(2, curTag.text);
        tmpEntry.second = tagIndex;
        // check if current tag is already in the collision list
        bool tagProcessed = false;
        for (tagIter = tmpIndex[c1][c2].begin(); tagIter
            != tmpIndex[c1][c2].end(); tagIter++)
          if (curTag == *tagIter) {
            tagIter->insts.push_back(tmpEntry);
            tagProcessed = true;
            break;
          }
        if (!tagProcessed) {
          curTag.insts.clear();
          curTag.insts.push_back(tmpEntry);
          tmpIndex[c1][c2].push_back(curTag);
        }
      }
    }

    // Build definitive index that uses vectors instead of lists
    index.resize(indexDim1);
    for (c1 = 0; c1 < indexDim1; c1++) {
      index[c1].resize(indexDim2);
      for (c2 = 0; c2 < indexDim2; c2++) {
        if (tmpIndex[c1][c2].size() == 0)
          continue;
        index[c1][c2].resize(tmpIndex[c1][c2].size());
        tagIndex = 0;
        for (list<Tag>::iterator iter = tmpIndex[c1][c2].begin(); iter
            != tmpIndex[c1][c2].end(); iter++)
          index[c1][c2][tagIndex++] = *iter;
        tmpIndex[c1][c2].clear();
      }
    }

    for (c1 = 0; c1 < indexDim1; c1++)
      tmpIndex[c1].resize(0);
  }

  void DB_index::insert(char *tag, int proteinID, int tagPos)
  {

  }

  bool DB_index::find(const char *tag, list<pair<int, int> > **location)
  {
    int c1, c2; // Index coordinates 1, 2
    vector<Tag>::iterator tagIter; // Iterator of tags in a collision list

    c1 = hash(1, tag);
    c2 = hash(2, tag);
    for (tagIter = index[c1][c2].begin(); tagIter != index[c1][c2].end(); tagIter++)
      if (compareTags(tag, tagIter->text)) {
        (*location) = &(tagIter->insts);
        return true;
      }
    (*location) = (list<pair<int, int> > *) 0;
    return false;
  }

  bool DB_index::find(const char *tag,
                      DB_fasta &db,
                      list<sps::tuple<int,float,string> > &peptides,
                      int minMatchFlanking,
                      float flankPref,
                      float tolPref,
                      float flankSuff,
                      float tolSuff)
  {
    list < pair<int, int> > *location; // peptides.clear(); // Allow addition of peptides to an existing set
    if (flankPref < -tolPref || flankSuff < -tolSuff)
    {
      //DEBUG_TRACE;
      return false;
    }
    if (not find(tag, &location))
    {
      //DEBUG_TRACE;
      return false;
    }
    // Match flanking masses
    string curPeptide;
    curPeptide.reserve(128);
    sps::tuple<int,float,string> curTuple;
    vector<float> protMasses;
    int curMatchingFlanks;
    int aaIdx, pepStart, pepEnd;
    float cumMass;
    for (list<pair<int, int> >::iterator iter = location->begin(); iter != location->end(); iter++) {
      db.getMassesIdx(iter->first, protMasses);

      // Match prefix flanking mass
      aaIdx = iter->second;
      cumMass = 0;
      curMatchingFlanks = 0;
      while (fabs(cumMass - flankPref) > tolPref + 0.0001 and cumMass < flankPref) {
        if (aaIdx == 0)
          break;
        else
          cumMass += protMasses[--aaIdx];
      }

      if (fabs(cumMass - flankPref) <= tolPref + 0.0001) {
        curMatchingFlanks++;
        pepStart = aaIdx;
      } else {
        if (minMatchFlanking == 2)
          continue;
        else
          pepStart = aaIdx + 1;
      } // Default to modifications with positive mass offset

      // Match suffix flanking mass
      aaIdx = min(iter->second + tagLength, (int) protMasses.size());
      cumMass = 0;
      while (fabs(cumMass - flankSuff) > (tolSuff + 0.0001) and
             cumMass < flankSuff) {
        if (aaIdx == (int) protMasses.size())
          break;
        else
          cumMass += protMasses[aaIdx++];
      }

      if (fabs(cumMass - flankSuff) <= tolSuff + 0.0001) {
        curMatchingFlanks++;
        pepEnd = aaIdx;
      } // One aa past the end of the found peptide
      else {
        pepEnd = aaIdx - 1;
      }

      if (curMatchingFlanks < minMatchFlanking)
        continue;

      curPeptide.resize(pepEnd - pepStart);
      for (aaIdx = 0; aaIdx < (int) curPeptide.size(); aaIdx++)
        curPeptide[aaIdx] = db[iter->first][pepStart + aaIdx];
      curTuple.m0 = iter->first;
      curTuple.m1 = db.getMassesSpec(iter->first)[iter->second][0] - flankPref;
      curTuple.m2 = curPeptide;
      peptides.push_back(curTuple);
    }
    return peptides.size() > 0;
  }

  bool DB_index::findAll(const char *tag,
                      DB_fasta &db,
                      list<sps::tuple<int,float,string> > &peptides,
                      int minMatchFlanking,
                      float flankPref,
                      float tolPref,
                      float flankSuff,
                      float tolSuff)
  {
    list < pair<int, int> > *location; // peptides.clear(); // Allow addition of peptides to an existing set
    if (flankPref < -tolPref || flankSuff < -tolSuff)
    {
      //DEBUG_TRACE;
      return false;
    }
    if (not find(tag, &location))
    {
      //DEBUG_TRACE;
      return false;
    }
    // Match flanking masses
    string curPeptide;
    curPeptide.reserve(128);
    sps::tuple<int,float,string> curTuple;
    vector<float> protMasses;
    int curMatchingFlanks;
    int aaIdx, pepStart, pepEnd;
    float cumMass;
    
    if (DEBUG_TAG) DEBUG_VAR(tag);

    for (list<pair<int, int> >::iterator iter = location->begin(); iter != location->end(); iter++) {
      db.getMassesIdx(iter->first, protMasses);
      //DEBUG_VAR(protMasses.size());

      string dbString = db.getSequence(iter->first);
      //DEBUG_VAR(dbString.length());
      
      int prevPos = 0;
      int pos = dbString.find(tag, prevPos);

      if (DEBUG_TAG) DEBUG_VAR(iter->first);
      if (DEBUG_TAG) DEBUG_VAR(iter->second);
      if (DEBUG_TAG) DEBUG_VAR(pos);
      while (pos != string::npos) {

        // Match prefix flanking mass
        aaIdx = pos;       //iter->second
        if (DEBUG_TAG) DEBUG_VAR(aaIdx);
        cumMass = 0;
        curMatchingFlanks = 0;
        while (fabs(cumMass - flankPref) > tolPref + 0.0001 and cumMass < flankPref) {
          if (aaIdx == 0)
            break;
          else
            cumMass += protMasses[--aaIdx];
        }
        if (DEBUG_TAG) DEBUG_VAR(cumMass);

        if (fabs(cumMass - flankPref) <= tolPref + 0.0001) {
          curMatchingFlanks++;
          pepStart = aaIdx;
        } else {
          if (minMatchFlanking == 2)
            continue;
          else
            pepStart = aaIdx + 1;
        } // Default to modifications with positive mass offset
        if (DEBUG_TAG) DEBUG_VAR(pepStart);

        // Match suffix flanking mass
        aaIdx = min(pos + tagLength, (int) protMasses.size());
        cumMass = 0;
        while (fabs(cumMass - flankSuff) > (tolSuff + 0.0001) and
              cumMass < flankSuff) {
          if (aaIdx == (int) protMasses.size())
            break;
          else
            cumMass += protMasses[aaIdx++];
        }
        if (DEBUG_TAG) DEBUG_VAR(cumMass);

        if (fabs(cumMass - flankSuff) <= tolSuff + 0.0001) {
          curMatchingFlanks++;
          pepEnd = aaIdx;
        } // One aa past the end of the found peptide
        else {
          pepEnd = aaIdx - 1;
        }
        if (DEBUG_TAG) DEBUG_VAR(pepEnd);

        if (curMatchingFlanks < minMatchFlanking)
          continue;

        curPeptide.resize(pepEnd - pepStart);
        for (aaIdx = 0; aaIdx < (int) curPeptide.size(); aaIdx++)
          curPeptide[aaIdx] = db[iter->first][pepStart + aaIdx];
        curTuple.m0 = iter->first;
        curTuple.m1 = db.getMassesSpec(iter->first)[iter->second][0] - flankPref;
        curTuple.m2 = curPeptide;
        peptides.push_back(curTuple);

        prevPos = prevPos + pos + 1;
        if (DEBUG_TAG) DEBUG_VAR(prevPos);
        pos = dbString.find(tag, prevPos);
        if (DEBUG_TAG) DEBUG_VAR(pos);
        
      } // while (pos != string::npos)
      
      if (DEBUG_TAG) DEBUG_TRACE;

    } // for (list<pair<int, int> >::iterator
    
    if (DEBUG_TAG) DEBUG_TRACE;
    return peptides.size() > 0;
  }


  float MatchSpecToPeptide(Spectrum &spec,
                           const char *peptide,
                           float peakTol,
                           float ionOffset,
                           bool filterMatched,
                           bool *isReversed,
                           AAJumps *aaJumps,
                           vector<int> *pepMatches)
  {
    vector<float> pepMasses; // Peptide masses derived from peptide sequence
    Spectrum pepSpec; // Peptide spectrum derived from peptide sequence
    Spectrum specRev; // Reversed version of spec
    vector<int> idxSpec, idxSpecRev, idxPep, idxPepRev; // indices of matched peaks
    vector<int> *bestSpecIdx, *bestPepIdx; // Pointers to best-matching sets of peaks
    float score, scoreRev; // Match scores
    unsigned int peakIdx;

    spec.reverse(2 * ionOffset, &specRev);

    // Get peptide sequence as a spectrum
    string annot(peptide);
    if(aaJumps) {
//cerr<<"Using input aaJumps, first mass = "<<(*aaJumps)[0]<<endl;
    	aaJumps->getPRMMasses(annot, pepMasses);
    } else {
        AAJumps jumps(1);
        jumps.getPRMMasses(annot, pepMasses);
    }
    pepSpec.resize(pepMasses.size() + 1);
    pepSpec[0].set(0.0, 0.0);
//cerr<<"  (MatchSpecToPeptide) peptide "<<peptide<<" -> ";
    for (peakIdx = 0; peakIdx < pepMasses.size(); peakIdx++) {
      pepSpec[peakIdx + 1].set(pepMasses[peakIdx] + ionOffset, 0.0);
//cerr<<pepMasses[peakIdx]+ionOffset<<", ";
    }
//cerr<<endl;

    // Match peptide to spectrum
    score = ScoreOverlap6(spec, pepSpec, 0.0, peakTol, idxSpec, idxPep);
//cerr<<"  --> score = "<<score<<endl;
//for(unsigned int i=0; i<idxSpec.size(); i++)
//	cerr<<"  --> ["<<i<<"]: ("<<idxSpec[i]<<","<<spec[idxSpec[i]][0]<<","<<spec[idxSpec[i]][1]<<") / ("<<idxPep[i]<<","<<pepSpec[idxPep[i]][0]<<")\n";
    scoreRev = ScoreOverlap6(specRev,
                             pepSpec,
                             0.0,
                             peakTol,
                             idxSpecRev,
                             idxPepRev);
//cerr<<"  --> scoreRev = "<<scoreRev<<endl;
//for(unsigned int i=0; i<idxSpecRev.size(); i++)
//	cerr<<"  --> ["<<i<<"]: ("<<idxSpecRev[i]<<","<<specRev[idxSpecRev[i]][0]<<","<<specRev[idxSpecRev[i]][1]<<") / ("<<idxPep[i]<<","<<pepSpec[idxPepRev[i]][0]<<")\n";
    if (score >= scoreRev - 0.0001) {
      if(isReversed)
        (*isReversed) = false;
      bestSpecIdx = &idxSpec;
      bestPepIdx = &idxPep;
    }
    else {
      if (filterMatched)
        spec = specRev;
      if(isReversed)
    	  (*isReversed) = true;
      bestSpecIdx = &idxSpecRev;
      bestPepIdx = &idxPepRev;
      score = scoreRev;
    }

    if (filterMatched) {
      for (peakIdx = 0; peakIdx < bestSpecIdx->size(); peakIdx++) {
        spec[peakIdx] = spec[(*bestSpecIdx)[peakIdx]];
        (*bestSpecIdx)[peakIdx] = peakIdx;
      }
      spec.resize(bestSpecIdx->size());
    }
    if (pepMatches)
      pepMatches->assign(bestPepIdx->begin(), bestPepIdx->end());
    return score;
  }

  void unit_MatchSpecsToPeps()
  {
    SpecSet specs(3), peptides(1);

    // Set specs[0] to TIDE
    specs[0].resize(3);
    specs[0].parentMass = 477.2; // includes H2O + 1
    specs[0][0].set(102.05, 1);
    specs[0][1].set(215.13, 1);
    specs[0][2].set(330.16, 1);

    // Set specs[1] to TIDE, with a noise peak in the zero-mass range and the parent mass peak, score = 4
    specs[1].resize(5);
    specs[1].parentMass = 477.2; // includes H2O + 1
    specs[1][0].set(10.05, 1);
    specs[1][1].set(102.05, 1);
    specs[1][2].set(215.13, 1);
    specs[1][3].set(330.16, 1);
    specs[1][4].set(459.05, 1);

    // Set specs[2] to TIDE, only y-ions, score =3
    specs[2].resize(3);
    specs[2].parentMass = 477.2; // includes H2O + 1
    specs[2][0].set(148.04, 1);
    specs[2][1].set(263.07, 1);
    specs[2][2].set(376.15, 1);

    // Set peptides[0] to PEPTIDE
    peptides[0].resize(8);
    peptides[0].parentMass = 800.35; // includes H2O + 1
    peptides[0][0].set(1, 1);
    peptides[0][1].set(98.05, 1);
    peptides[0][2].set(227.1, 1);
    peptides[0][3].set(324.15, 1);
    peptides[0][4].set(425.2, 1);
    peptides[0][5].set(538.28, 1);
    peptides[0][6].set(653.31, 1);
    peptides[0][7].set(782.35, 1);

    vector < TwoValues<float> > scores(1);
    vector<float> ionOffsets(1);
    ionOffsets[0] = 0;
    scores[0].set(1, 0);

    vector < TwoValues<float> > pepMatchScores;
    vector < list<unsigned int> > pepMatchedSpecs;
    MatchSpecsToPeps(specs,
                     peptides,
                     0.5,
                     1,
                     19,
                     0,
                     scores,
                     ionOffsets,
                     pepMatchScores,
                     pepMatchedSpecs,
                     true);
    if (pepMatchScores[0][0] != 10) {
      cerr << "(Failure in unit_MatchSpecsToPeps) - expected score 10, got "
          << pepMatchScores[0][0] << endl;
    }

    // Repeat the test with no c-term H2O
    for (unsigned int i = 0; i < specs.size(); i++)
      specs[i].parentMass -= 18;
    for (unsigned int i = 0; i < specs[2].size(); i++)
      specs[2][i][0] -= 18; // Subtract H2O from y-ions
    for (unsigned int i = 0; i < peptides.size(); i++)
      peptides[i].parentMass -= 18;
    MatchSpecsToPeps(specs,
                     peptides,
                     0.5,
                     1,
                     1,
                     0,
                     scores,
                     ionOffsets,
                     pepMatchScores,
                     pepMatchedSpecs,
                     true);
    if (pepMatchScores[0][0] != 11) { // 11 because parent mass is matched twice with spectrum 1
      cerr << "(Failure in unit_MatchSpecsToPeps) - expected score 11, got "
          << pepMatchScores[0][0] << endl;
    }
  }

  //
  //  MatchSpecsToPeps - matches a set of spectra to a set of peptides (represented
  //    as a set of sets of masses in peptides) and returns a set of peptide
  //    (score,index) pairs sorted by increasing score.
  //
  //    scores     - contain 2 values per ion type: premium if present, penalty if absent (pos 0/1).
  //                   If scores is empty then MatchSpecsToPeps returns the percentage of explained intensity.
  //                   *** NOTE ***: should be SORTED by decreasing ion score (first match takes the peak assignment).
  //    ionOffsets - contains the mass offsets (from the prefix mass) where the
  //                   spectrum peaks are to be found: [0] = -18 for b-H2O, etc.
  //                   *** NOTE ***: should be in the same order as scores.
  //    ionOffset  - mass offset caused by each fragment charge (~1 for raw MS/MS, 0 for PRM spectra).
  //                   (expected on both specs and peptides)
  //    ctermMass  - value such that specs[i].parentMass-ctermMass and peptides[j].parentMass-ctermMass indicate the summed amino acid masses
  //
  //    pepMatchScores  - Summed score of peptide matches to spectra (pos.0), peptide index (pos.1); sorted by increasing match scores / peptide indices.
  //    pepMatchedSpecs - Indices of the spectra that could be matched to the corresponding peptide
  //
  bool MatchSpecsToPeps(SpecSet &specs,
                        SpecSet &peptides,
                        float peakTol,
                        float ionOffset,
                        float ctermMass,
                        float noisePenalty,
                        vector<TwoValues<float> > &scores,
                        vector<float> &ionOffsets,
                        vector<TwoValues<float> > &pepMatchScores,
                        vector<list<unsigned int> > &pepMatchedSpecs,
                        bool matchSuffix)
  {
    if (specs.size() == 0 or peptides.size() == 0) {
      pepMatchScores.resize(0);
      pepMatchedSpecs.resize(0);
      return false;
    }
    unsigned int pepIdx, specIdx, pivot, startIdx; // Index of the peptide peak where the current spectrum starts
    int endIdx, // Index of the peptide peak where the current spectrum ends
        numPepPeaks; // Number of peptide peaks between startIdx and endIdx
    Spectrum revPep; // Used to match y-ions
    float curScore, // Current score when matching a spectrum against a peptide
        bestScore, // Best score when matching a spectrum against a peptide
        specTIC, // Total spectrum intensity (Total Ion Current), used for percentage of explained intensity
        worstPenalty, // Worst possible penalty for missing a peak
        extraCTermMass, // Used to match y-ions when a spectrum is matched to a water loss from the peptide
        specPeptideMass,// Summed amino acid masses from spectrum's parent mass
        peptideMass; // Summed amino acid masses from peptide's parent mass
    vector<int> idxMatchedPep, idxMatchedSpec;
    vector<bool> matchedSpecPeaks;
    bool matchedSpec;
    //	for(specIdx=0; specIdx<specs.size(); specIdx++) specs[specIdx].parentMass -= ionOffset;
    worstPenalty = 0;
    for (pivot = 0; pivot < scores.size(); pivot++)
      worstPenalty = min(worstPenalty, scores[pivot][1]);

    //for(specIdx=0;specIdx<specs.size();specIdx++) { cerr<<"Spectrum "<<specIdx<<":\n"; specs[specIdx].output(cerr); }

    pepMatchScores.resize(peptides.size());
    pepMatchedSpecs.resize(peptides.size());
    for (pepIdx = 0; pepIdx < peptides.size(); pepIdx++) {
      peptides[pepIdx].reverse(2 * ionOffset, &revPep);
      //cerr<<"Reversed version of peptide "<<pepIdx<<":\n";   revPep.output(cerr);
      //		peptides[pepIdx].parentMass -= ionOffset;
      peptideMass = peptides[pepIdx].parentMass - ctermMass;
      pepMatchScores[pepIdx][0] = 0;
      pepMatchScores[pepIdx][1] = (float) pepIdx;

      for (specIdx = 0; specIdx < specs.size(); specIdx++) {
        // Start by assuming worst-case scenario: spectrum PM does not match peptide
        bestScore = worstPenalty * scores.size() * peptides[pepIdx].size()
            + noisePenalty * specs[specIdx].size();
        specPeptideMass = specs[specIdx].parentMass - ctermMass;
        //cerr<<"initial score peptide "<<pepIdx<<" is "<<scores[1]*2*peptides[pepIdx].size()<<" ("<<peptides[pepIdx].size()<<" peaks)\n";
        //cerr<<"initial score spectrum "<<specIdx<<" is "<<scores[2]*specs[specIdx].size()<<" ("<<specs[specIdx].size()<<" peaks)\n";
        matchedSpecPeaks.resize(specs[specIdx].size());
        specTIC = 0;
        for (pivot = 0; pivot < specs[specIdx].size(); pivot++)
          specTIC += specs[specIdx][pivot][1];

        matchedSpec = false;
        for (startIdx = 0; startIdx < peptides[pepIdx].size()
            and peptides[pepIdx][startIdx][0] - ionOffset + specPeptideMass
                <= peptideMass + peakTol; startIdx++) {
          // Find where the peptide could end
          //cerr<<"Looking for mass "<<peptides[pepIdx][startIdx][0]+specPeptideMass<<" from "<<peptides[pepIdx][startIdx][0]<<" (peak "<<startIdx<<", offset "<<specPeptideMass<<")\n";
          endIdx = peptides[pepIdx].findClosest(peptides[pepIdx][startIdx][0]
              + specPeptideMass);
          extraCTermMass = 0;
          curScore = 0;
          if (endIdx < (int) startIdx or fabs(peptides[pepIdx][endIdx][0]
              - peptides[pepIdx][startIdx][0] - specPeptideMass) > peakTol) {
            // Try to interpret spectrum as a water loss
            endIdx = peptides[pepIdx].findClosest(peptides[pepIdx][startIdx][0]
                + specPeptideMass + AAJumps::massH2O);
            extraCTermMass = -AAJumps::massH2O;
            if (endIdx < (int) startIdx or fabs(peptides[pepIdx][endIdx][0]
                - peptides[pepIdx][startIdx][0] - specPeptideMass
                + AAJumps::massH2O) > peakTol)
              continue;
          }
          //if(pepIdx==3272 or pepIdx==3492) cerr<<" -- found a matching peptide at ("<<startIdx<<","<<endIdx<<"), ("<<peptides[pepIdx][startIdx][0]<<","<<peptides[pepIdx][endIdx][0]<<"), specPeptideMass = "<<specPeptideMass<<" + extraCTermMass = "<<extraCTermMass<<"\n";
          matchedSpec = true;
          for (pivot = 0; pivot < matchedSpecPeaks.size(); pivot++)
            matchedSpecPeaks[pivot] = false;
          numPepPeaks = endIdx - startIdx - 1; // Masses zero and parent mass are always matched
          if (specs[specIdx][0][0] < 25)
            numPepPeaks++; // If there are peaks in the zero-mass range
          if (specs[specIdx][specs[specIdx].size() - 1][0] > specPeptideMass
              - 25)
            numPepPeaks++; // If there are peaks in the parent-mass range

          // Match prefix masses
          float shift = peptides[pepIdx][startIdx][0] - ionOffset;
          for (unsigned int ionIdx = 0; ionIdx < ionOffsets.size(); ionIdx++) {
            FindMatchPeaks(peptides[pepIdx], specs[specIdx], shift
                - ionOffsets[ionIdx], peakTol, idxMatchedPep, idxMatchedSpec);
            //					curScore = (scores[0]/2)*idxMatchedPep.size() + (scores[1]/2)*(numPepPeaks-idxMatchedPep.size());
            if (scores.size() > ionIdx)
              curScore += scores[ionIdx][0] * idxMatchedPep.size()
                  + scores[ionIdx][1] * (numPepPeaks - idxMatchedPep.size());
            for (pivot = 0; pivot < idxMatchedSpec.size(); pivot++)
              matchedSpecPeaks[idxMatchedSpec[pivot]] = true;
            /* if(pepIdx==3272 or pepIdx==3492) {
             cerr<<" -- -- Got "<<idxMatchedSpec.size()<<" prefix matches at offset "<<ionOffsets[ionIdx]<<", curScore = "<<curScore;
             if(idxMatchedSpec.size()>0) cerr<<": ";
             for(unsigned int i=0; i<idxMatchedPep.size();i++) cerr<<"("<<peptides[pepIdx][idxMatchedPep[i]][0]<<","<<specs[specIdx][idxMatchedSpec[i]][0]<<")";
             cerr<<"\n";
             }*/
          }

          if (matchSuffix) {
            // Match suffix masses
            //					shift = (peptides[pepIdx].parentMass-(peptides[pepIdx][startIdx][0]-ionOffset)-specs[specIdx].parentMass)-ctermMass-extraCTermMass;
            shift = peptideMass - (peptides[pepIdx][endIdx][0] - ionOffset)
                - extraCTermMass; // shift to apply to the spectrum masses
            for (unsigned int ionIdx = 0; ionIdx < ionOffsets.size(); ionIdx++) {
              FindMatchPeaks(revPep,
                             specs[specIdx],
                             shift - ionOffsets[ionIdx],
                             peakTol,
                             idxMatchedPep,
                             idxMatchedSpec);
              //						curScore += (scores[0]/2)*idxMatchedPep.size() + (scores[1]/2)*(numPepPeaks-idxMatchedPep.size());
              if (scores.size() > ionIdx)
                curScore += scores[ionIdx][0] * idxMatchedPep.size()
                    + scores[ionIdx][1] * (numPepPeaks - idxMatchedPep.size());
              for (pivot = 0; pivot < idxMatchedSpec.size(); pivot++)
                matchedSpecPeaks[idxMatchedSpec[pivot]] = true;
              /* if(pepIdx==3272 or pepIdx==3492) {
               cerr<<" -- -- Got "<<idxMatchedSpec.size()<<" suffix matches at offset "<<ionOffsets[ionIdx]<<", shift = "<<shift<<", curScore = "<<curScore;
               if(idxMatchedSpec.size()>0) cerr<<": ";
               for(unsigned int i=0; i<idxMatchedPep.size();i++) cerr<<"("<<revPep[idxMatchedPep[i]][0]<<","<<specs[specIdx][idxMatchedSpec[i]][0]<<")";
               cerr<<"\n";
               }*/
            }
          }

          //				for(pivot=0; pivot<matchedSpecPeaks.size(); pivot++)
          //					curScore+=(matchedSpecPeaks[pivot]?scores[0]:scores[2]);

          if (scores.size() == 0) { // Explained intensity in matched peaks
            curScore = 0;
            for (pivot = 0; pivot < matchedSpecPeaks.size(); pivot++)
              if (matchedSpecPeaks[pivot])
                curScore += specs[specIdx][pivot][1];
            curScore /= specTIC;
          }

          // Keep track of best match
          //cerr<<" -- curScore = "<<curScore<<", bestScore = "<<bestScore<<endl;
          if (curScore > bestScore)
            bestScore = curScore;
        }
        pepMatchScores[pepIdx][0] += bestScore;
        if (matchedSpec)
          pepMatchedSpecs[pepIdx].push_back(specIdx);
      }
    }
    sort(pepMatchScores.begin(), pepMatchScores.end());

    return true;
  }

  //
  //  MatchSpecsToPeps_old - matches a set of spectra to a set of peptides (represented
  //    as a set of sets of masses in peptides) and returns a set of peptide
  //    (score,index) pairs sorted by increasing score. The vector scores is supposed
  //    to contain 4 values for the following possibilities:
  //       Peak in peptide      Peak in spectrum     scores entry
  //            yes                   yes             scores[0]
  //            yes                   no              scores[1]
  //            no                    yes             scores[2]
  //    ionOffsets - contains the mass offsets (from the prefix mass) where the
  //       spectrum peaks are to be found: [0] = -18 for b-H2O, etc.
  //
  //    offsetScores[i][j] - ion mass offset, premium if present, penalty if absent (j=0/1/2)
  //
  bool MatchSpecsToPeps_old(SpecSet &specs,
                            SpecSet &peptides,
                            float peakTol,
                            float ionOffset,
                            float ctermMass,
                            vector<float> &scores,
                            vector<float> &ionOffsets,
                            vector<TwoValues<float> > &pepMatchScores,
                            vector<list<unsigned int> > &pepMatchedSpecs,
                            bool matchSuffix)
  {
    if (specs.size() == 0 or peptides.size() == 0 or scores.size() != 3) {
      pepMatchScores.resize(0);
      pepMatchedSpecs.resize(0);
      return false;
    }
    unsigned int pepIdx, specIdx, pivot, startIdx; // Index of the peptide peak where the current spectrum starts
    int endIdx, // Index of the peptide peak where the current spectrum ends
        numPepPeaks; // Number of peptide peaks between startIdx and endIdx
    Spectrum revPep; // Used to match y-ions
    float curScore, // Current score when matching a spectrum against a peptide
        bestScore, // Best score when matching a spectrum against a peptide
        extraCTermMass; // Used to match y-ions when a spectrum is matched to a water loss from the peptide
    vector<int> idxMatchedPep, idxMatchedSpec;
    vector<bool> matchedSpecPeaks;
    bool matchedSpec;
    for (specIdx = 0; specIdx < specs.size(); specIdx++)
      specs[specIdx].parentMass -= ionOffset;

    pepMatchScores.resize(peptides.size());
    pepMatchedSpecs.resize(peptides.size());
    for (pepIdx = 0; pepIdx < peptides.size(); pepIdx++) {
      peptides[pepIdx].reverse(2 * ionOffset, &revPep);
      //cerr<<"Reversed version of peptide "<<pepIdx<<":\n";   revPep.output(cerr);
      peptides[pepIdx].parentMass -= ionOffset;
      pepMatchScores[pepIdx][0] = 0;
      pepMatchScores[pepIdx][1] = (float) pepIdx;

      for (specIdx = 0; specIdx < specs.size(); specIdx++) {
        // Start by assuming worst-case scenario: spectrum PM does not match peptide
        bestScore = scores[1] * 2 * peptides[pepIdx].size() + scores[2]
            * specs[specIdx].size();
        //cerr<<"initial score peptide "<<pepIdx<<" is "<<scores[1]*2*peptides[pepIdx].size()<<" ("<<peptides[pepIdx].size()<<" peaks)\n";
        //cerr<<"initial score spectrum "<<specIdx<<" is "<<scores[2]*specs[specIdx].size()<<" ("<<specs[specIdx].size()<<" peaks)\n";
        matchedSpecPeaks.resize(specs[specIdx].size());

        matchedSpec = false;
        for (startIdx = 0; startIdx < peptides[pepIdx].size()
            and peptides[pepIdx][startIdx][0] - ionOffset
                + specs[specIdx].parentMass <= peptides[pepIdx].parentMass
                + peakTol; startIdx++) {
          // Find where the peptide could end
          //cerr<<"Looking for mass "<<peptides[pepIdx][startIdx][0]+specs[specIdx].parentMass-ctermMass<<" from "<<peptides[pepIdx][startIdx][0]<<" (peak "<<startIdx<<", offset "<<specs[specIdx].parentMass-ctermMass<<")\n";
          endIdx = peptides[pepIdx].findClosest(peptides[pepIdx][startIdx][0]
              + specs[specIdx].parentMass - ctermMass);
          extraCTermMass = 0;
          if (endIdx < (int) startIdx or fabs(peptides[pepIdx][endIdx][0]
              - peptides[pepIdx][startIdx][0] - specs[specIdx].parentMass)
              > peakTol) {
            // Try to interpret spectrum as a water loss
            endIdx = peptides[pepIdx].findClosest(peptides[pepIdx][startIdx][0]
                + specs[specIdx].parentMass - ctermMass + 18.010564686);
            extraCTermMass = -18.010564686;
            if (endIdx < (int) startIdx or fabs(peptides[pepIdx][endIdx][0]
                - peptides[pepIdx][startIdx][0] - specs[specIdx].parentMass
                + AAJumps::massH2O) > peakTol)
              continue;
          }
          matchedSpec = true;
          for (pivot = 0; pivot < matchedSpecPeaks.size(); pivot++)
            matchedSpecPeaks[pivot] = false;
          numPepPeaks = endIdx - startIdx + 1;

          // Match prefix masses
          float shift = peptides[pepIdx][startIdx][0] - ionOffset;
          FindMatchPeaks(peptides[pepIdx],
                         specs[specIdx],
                         shift,
                         peakTol,
                         idxMatchedPep,
                         idxMatchedSpec);
          curScore = scores[0] * idxMatchedPep.size() + scores[1]
              * (numPepPeaks - idxMatchedPep.size());
          for (pivot = 0; pivot < idxMatchedSpec.size(); pivot++)
            matchedSpecPeaks[idxMatchedSpec[pivot]] = true;

          // Match prefix neutral losses
          for (unsigned int ionIdx = 0; ionIdx < ionOffsets.size(); ionIdx++) {
            FindMatchPeaks(peptides[pepIdx], specs[specIdx], shift
                + ionOffsets[ionIdx], peakTol, idxMatchedPep, idxMatchedSpec);
            curScore = (scores[0] / 2) * idxMatchedPep.size() + (scores[1] / 2)
                * (numPepPeaks - idxMatchedPep.size());
            for (pivot = 0; pivot < idxMatchedSpec.size(); pivot++)
              matchedSpecPeaks[idxMatchedSpec[pivot]] = true;
          }

          if (matchSuffix) {
            // Match suffix masses
            shift = (peptides[pepIdx].parentMass
                - (peptides[pepIdx][startIdx][0] - ionOffset)
                - specs[specIdx].parentMass) - ctermMass - extraCTermMass;
            FindMatchPeaks(revPep,
                           specs[specIdx],
                           shift,
                           peakTol,
                           idxMatchedPep,
                           idxMatchedSpec);
            curScore += scores[0] * idxMatchedPep.size() + scores[1]
                * (numPepPeaks - idxMatchedPep.size());
            for (pivot = 0; pivot < idxMatchedSpec.size(); pivot++)
              matchedSpecPeaks[idxMatchedSpec[pivot]] = true;

            // Match suffix neutral losses
            for (unsigned int ionIdx = 0; ionIdx < ionOffsets.size(); ionIdx++) {
              FindMatchPeaks(revPep,
                             specs[specIdx],
                             shift + ionOffsets[ionIdx],
                             peakTol,
                             idxMatchedPep,
                             idxMatchedSpec);
              curScore += (scores[0] / 2) * idxMatchedPep.size() + (scores[1]
                  / 2) * (numPepPeaks - idxMatchedPep.size());
              for (pivot = 0; pivot < idxMatchedSpec.size(); pivot++)
                matchedSpecPeaks[idxMatchedSpec[pivot]] = true;
            }
          }

          for (pivot = 0; pivot < matchedSpecPeaks.size(); pivot++)
            curScore += (matchedSpecPeaks[pivot] ? scores[0] : scores[2]);

          // Keep track of best match
          if (curScore > bestScore)
            bestScore = curScore;
        }
        pepMatchScores[pepIdx][0] += bestScore;
        if (matchedSpec)
          pepMatchedSpecs[pepIdx].push_back(specIdx);
      }
    }
    sort(pepMatchScores.begin(), pepMatchScores.end());

    return true;
  }
}

