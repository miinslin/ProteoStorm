#ifndef DB_FASTA_H

/** */
#define DB_FASTA_H

#include "aminoacid.h"
#include "spectrum.h"
#include "SpecSet.h"
#include "tuple.h"
#include <vector>
#include <cstdio>
#include <iostream>

namespace specnets
{
	using namespace std;

	/**
	 * Database of protein sequences with load/save in FASTA format
	 */
	class DB_fasta {
	protected:

		/**
		 * TODO: add description
		 *
		 */
		static char emptyStr[];

		/**
		 * Protein descriptions stored as \0-terminated strings
		 *
		 */
		vector<char *> desc;

		/**
		 * Protein sequences stored as \0-terminated strings
		 *
		 */
		vector<char *> sequences;

		/*
		 * old index: not currently used
		 *
		 */
		vector < map <int, vector< pair<int, int> > > > massOffsets;
		map<int, vector<pair<int, int> > >::iterator it;


		 // Mass index
		 vector < vector< pair<int, int> > > massesIndex;



		 vector< pair< int, int > > inline intersectSets(vector< pair< int, int> > setA, vector< pair< int, int> > setB )
		 {
			/* DEBUG_MSG("A:");
			 for(int a = 0; a < setA.size(); a++)
			 {
				 DEBUG_MSG(setA[a].first << " " <<  setA[a].second);

			 }

			 DEBUG_MSG("B:");
			 for(int a = 0; a < setB.size(); a++)
			 {
				 DEBUG_MSG(setB[a].first << " " <<  setB[a].second);

			 }*/
			 vector< pair< int, int > > retvec;

			 for(int a = 0; a < setA.size(); a++)
			 {

				 for(int b = 0; b < setB.size(); b++)
				 {
					 if(setA[a].first == setB[b].first && setA[a].second == setB[b].second )
					 {
						 retvec.push_back(setA[a]);

					 }
				 }
			 }
			 return retvec;

		 }

		 bool inline checkSingleton(int m1, int m2, int m3)
		 {

			 if(m1 <= 0 && m2 <= 0)
				 return true;
			 else if(m1 <= 0 && m3 <= 0)
				 return true;
			 else if(m2 <= 0 && m3 <= 0)
				 return true;
			 else
				 return false;

		 }

		 int inline returnSingleton(int m1, int m2, int m3)
		 {

			 if(m1 <= 0 && m2 <= 0)
				 return m3;
			 else if(m1 <= 0 && m3 <= 0)
				 return m2;
			 else if(m2 <= 0 && m3 <= 0)
				 return m1;
			 else
				 return false;

		 }


	public:

		 void addDecoyReversed();
		 void addDecoyShuffled();
     void replaceDecoyShuffled();
     void shuffleSequence(char * original, char * shuffled);

		/**
		 * Protein ID descriptors stored as \0-terminated strings
		 *
		 */
		vector<char *> IDs;

		/**
		 * Contains each protein's cumulative masses as a Spectrum.
		 */
		vector<Spectrum> masses;

		/**
		 * TODO: add description
		 *
		 */
		~DB_fasta();

		/**
		 * TODO: add description
		 *
		 *@param id
		 *@return
		 */
		char *operator[](char *id) const;

		/**
		 * TODO: add description
		 *
		 *@param index
		 *@return
		 */
		char *operator[](int index) const {
			return sequences[index];
		}

		/**
		 * Copies content from another DB_fasta object
		 *
		 *@param other Contents to be copied into this object
		 *@return
		 */
		DB_fasta &operator=(const DB_fasta &other);

		/**
		 * TODO: add description
		 *
		 *@param index
		 *@param putHere
		 */
		void getMassesIdx(int index, vector<float> &putHere) {
			getMasses(sequences[index], putHere);
		}

		/**
		 * TODO: add description
		 *
		 *@param index
		 *@return
		 */
		Spectrum &getMassesSpec(int index);

		/**
		 * Replaces all amino acids prevAA with repAA (e.g., all 'I' become 'L')
		 * on all protein sequences.
		 *
		 *@param prevAA Previous amino acid character
		 *@param repAA Replacement amino acid character
		 */
		void replaceAA(char prevAA, char repAA);


		/**
		 * Christina Boucher: returns true if there the integral pair
		 * is found in the massOffsets data structure; otherwise, false
		 * is returned.
		 */
		 vector<string> findPeptideMultiDimIndex(int, int, int);


		 /**
		  * Christina Boucher: populates massOffsets data structure.
		  *
		  *@param
		  *@param
		  */
		 void populateIndex(bool);

		 vector<int> peptideIndices;

		 /**
		  * Christina Boucher: returns true if there the integral pair
		  * is found in the massOffsets data structure; otherwise, false
		  * is returned.
		  */
		 vector<string> findPeptides(int, int, int);

		 vector<string> findPeptidesForTriple(int, int, int, int, int, int );

		 vector<string> findPeptideForSingleMass(int);	 // need to implement this


		/**
		 * Get protein ID text from protein index.
		 *
		 *@param index Protein index
		 *@return Protein ID string (char *)
		 */
		char *getID(unsigned int index) const {
			if (index<0 or index>=IDs.size()) return (char *)0;
			else return IDs[index];
		}

		/**
		 * Provide access to protein sequences
		 *
		 *@param index - Protein index
		 *@return (char *) to \0-terminated string (protein sequence)
		 */
		char *getSequence(unsigned int index) const {
			if (index<0 or index>=sequences.size()) return (char *)0;
			else return sequences[index];
		}

		/**
		 * TODO: add description
		 *
		 *@param index
		 *@return
		 */
		char *getDesc(unsigned int index) const {
			if (index<0 or index>=desc.size()) return (char *)0;
			else return desc[index];
		}

		bool checkIfDecoy(int index)
		{

			if(index < 0 || index >= sequences.size())
			{
				DEBUG_MSG("index out of range. index: " << index);
				exit(1);
			}
			else{

			    string c = desc[index];

				if(c[0] == 'X' && c[1] == 'X' && c[2] == 'X')
				{
					return true;
				}
				return false;

			}
		}

		/**
		 * TODO: add description
		 *
		 *@return
		 */
		unsigned int size() const {
			return IDs.size();
		}

		/**
		 * TODO: add description
		 *
		 */
		void reset();

		/**
		 * TODO: add description
		 *
		 *@param filename
		 *@return
		 */
		void Add(char *ID, char *des, char *seq);

		/**
		 * TODO: add description
		 *
		 *@param filename
		 *@return
		 */
		unsigned int Load(const char *filename);

		/**
		 * TODO: add description
		 *
		 *@param filename
		 *@return
		 */
		unsigned int Save(const char *filename);

		/**
		 * TODO: add description
		 *
		 *@param filename
		 *@return
		 */
    unsigned int Save(const char * filename, vector<bool> & saveFlags);

		/**
		 * TODO: add description
		 *
		 *@param out
		 */
		void output(ostream &out);

		/**
		 * Searches all protein sequences for the string peptide.
		 *
		 *@param peptide
		 *@param matches List of (protein index, start index, peptide) for matched strings
		 *@return Number of matches
		 */
		unsigned int find(const char *peptide, list<sps::tuple<int,float,string> > &matches);

    /**
     * Searches all protein sequences for the string peptide.
     *
     *@param peptide
     *@param matches List of (protein index, start mass, peptide) for matched strings
     *@return Number of matches
     */
    unsigned int find(const string &peptide, list<sps::tuple<int,float,string> > &matches);

	};

	/**
	 * TODO: add description
	 *
	 *@param filename
	 *@param db
	 *@param protID1
	 *@param protID2
	 *@param matchedIndices
	 *@return
	 */
	bool Load_clustalw(const char *filename, DB_fasta &db, int protID1, int protID2,
			vector<TwoValues<int> > &matchedIndices);


	/**
	 * TODO: add description
	 *
	 *@param filename
	 *@param db
	 *@param matchedProts
	 *@param matchedAAs
	 *@param proteinIndexOffset Number to subtract from indices in index file to obtain zero-based indices
	 *@return
	 */
	bool Load_clustalw_multiple(const char *filename, DB_fasta &db,
			vector<TwoValues<int> > &matchedProts,
			vector<vector<TwoValues<int> > > &matchedAAs,
			unsigned int proteinIndexOffset = 0);

	/**
	 * TODO: add description
	 *
	 */
	class DB_index {


		/**
		 * TODO: add description
		 *
		 */
		class Tag {
		public:


			/**
			 * Tag text.
			 */
			char *text;

			/**
			 * List of tag instances as (protein ID, starting amino acid position).
			 */
			list<pair<int, int> > insts;

			/**
			 * TODO: add description
			 *
			 */
			Tag() {
				text = (char *) 0;
			}

			/**
			 * TODO: add description
			 *
			 *@param other
			 */
			Tag(const Tag &other) {
				if (other.text == (char *) 0) {
					text = (char *) 0;
					insts.clear();
					return;
				}
				text = new char[strlen(other.text) + 1];
				strcpy(text, other.text);
				insts = other.insts;
			}

			/**
			 * Destructor.
			 */
			~Tag() {
				if (text != (char *) 0)
					delete[] text;
			}

			/**
			 * TODO: add description
			 *
			 *@param newText
			 *@param textLen
			 */
			void setTag(char *newText, int textLen) {
				if (text != (char *) 0)
					delete[] text;
				text = new char[textLen + 1];
				strcpy(text, newText);
			}

			/**
			 * TODO: add description
			 *
			 *@param protId
			 *@param tagPos
			 */
			void addInst(int protId, int tagPos) {
				pair<int, int> tmp(protId, tagPos);
				insts.push_back(tmp);
			}

			/**
			 * TODO: add description
			 *
			 *@param newInstance
			 */
			void addInst(pair<int, int> &newInstance) {
				insts.push_back(newInstance);
			}

			/**
			 * TODO: add description
			 *
			 *@param other
			 *@return
			 */
			bool operator==(const Tag &other) {
				return strcmp(text, other.text) == 0;
			}

			/**
			 * TODO: add description
			 *
			 *@param other
			 *@return
			 */
			Tag &operator=(const Tag &other) {
				if (text != (char *) 0)
					delete[] text;
				text = new char[strlen(other.text) + 1];
				strcpy(text, other.text);
				insts = other.insts;
				return *this;
			}
		};

		/**
		 * TODO: add description
		 *
		 */
		vector<vector<vector<Tag> > > index;

		/**
		 * Coefficients for hash1.
		 */
		vector<int> coeffs1;

		/**
		 *  Coefficients for hash2.
		 */
		vector<int> coeffs2;


		/**
		 * TODO: add description
		 *
		 */
		int tagLength;

		/**
		 * TODO: add description
		 *
		 */
		int indexDim1;

		/**
		 * TODO: add description
		 *
		 */
		int indexDim2;

		/**
		 * TODO: add description
		 *
		 *@param tagLength
		 */
		void hash_init(int tagLength);

		/**
		 * TODO: add description
		 *
		 *@param hashIdx
		 *@param tag
		 *@return
		 */
		int hash(int hashIdx, const char *tag);

		/**
		 * TODO: add description
		 *
		 *@param tag1
		 *@param tag2
		 *@return
		 */
		bool compareTags(const char *tag1, const char *tag2) {
			for (int tagIndex = 0; tagIndex < tagLength; tagIndex++)
				if (tag1[tagIndex] != tag2[tagIndex])
					return false;
			return true;
		}

		/**
		 * TODO: add description
		 *
		 *@param t
		 */
		void doNothing(char *t) {}

	public:

		/**
		 * TODO: add description
		 *
		 *@param db
		 *@param newIndexDim1
		 *@param newIndexDim2
		 *@param newTagLength
		 */
		DB_index(DB_fasta &db, int newIndexDim1, int newIndexDim2,
				int newTagLength);


		/**
		 * TODO: add description
		 *
		 *@param db
		 *@param newIndexDim1
		 *@param newIndexDim2
		 *@param newTagLength
		 */
		void buildIndex(DB_fasta &db, int newIndexDim1, int newIndexDim2,
				int newTagLength);

		/**
		 * TODO: add description
		 *
		 *@param tag
		 *@param proteinID
		 *@param tagPos
		 */
		void insert(char *tag, int proteinID, int tagPos);

		/**
		 * TODO: add description
		 *
		 *@param tag
		 *@param location
		 *@return
		 */
		bool find(const char *tag, list<pair<int, int> > **location);

		/**
		 * TODO: add description
		 *
		 *@param tag
		 *@param db
		 *@param peptides List of (protein index, peptide) for matched tags
		 *@param minMatchFlanking
		 *@param flankPref
		 *@param tolPref
		 *@param flankSuff
		 *@param tolSuff
		 *@return
		 */
		bool find(const char *tag,
		          DB_fasta &db, 
		          list<sps::tuple<int, float, string> > &peptides,
                        int minMatchFlanking, 
                        float flankPref, 
                        float tolPref,
                        float flankSuff, 
                        float tolSuff);

              bool findAll(const char *tag,
		          DB_fasta &db, 
		          list<sps::tuple<int, float, string> > &peptides,
                        int minMatchFlanking, 
                        float flankPref, 
                        float tolPref,
                        float flankSuff, 
                        float tolSuff);
	};


	/**
	 * MatchSpecToPeptide - scores the direct/reversed match of the spectrum spec
	 * against a given peptide sequence; the match scores are defined as the summed
	 * intensity of matched peaks in spec.
	 *
	 *@param spec PRM or MS/MS spectrum
	 *@param peptide Peptide string
	 *@param peakTol Peak mass tolerance
	 *@param ionOffset Ion offset added to theoretical peptide masses to match masses in spec
	 *@param filterMatched Select matched peaks in spec
	 *@param isReversed Indicates whether the maximum score was with a match to the reversed spectrum
	 *@param aaJumps Use these amino acid masses instead of the global masses
	 *@param pepMatches Indices of matched peptide masses (1-to-1 correspondence to spec if filterMatched set to true)
	 *@return Summed score of matched peaks (in spec)
	 */
	float MatchSpecToPeptide(Spectrum &spec, const char *peptide, float peakTol,
			float ionOffset, bool filterMatched, bool *isReversed = 0,
			AAJumps *aaJumps = 0, vector<int> *pepMatches = 0);

	/**
	 * MatchSpecsToPeps - scores the match of each spectrum in specs against all entries
	 *  in peptides (peptide substrings are also matched); each spectrum is assigned to a
	 *  peptide with maximal match score.
	 *
	 *@param specs
	 *@param peptides
	 *@param peakTol
	 *@param ionOffset
	 *@param ctermMass
	 *@param noisePenalty
	 *@param scores
	 *@param ionOffsets
	 *@param pepMatchScores
	 *@param pepMatchedSpecs
	 *@param matchSuffix
	 *@return
	 */
	bool MatchSpecsToPeps(SpecSet &specs, SpecSet &peptides, float peakTol,
			float ionOffset, float ctermMass, float noisePenalty,
			vector<TwoValues<float> > &scores, vector<float> &ionOffsets,
			vector<TwoValues<float> > &pepMatchScores,
			vector<list<unsigned int> > &pepMatchedSpecs, bool matchSuffix = true);

	/**
	 * TODO: add description
	 *
	 *@param specs
	 *@param peptides
	 *@param peakTol
	 *@param ionOffset
	 *@param ctermMass
	 *@param scores
	 *@param ionOffsets
	 *@param pepMatchScores
	 *@param pepMatchedSpecs
	 *@param matchSuffix
	 *@
	 *@return
	 */
	bool MatchSpecsToPeps_old(SpecSet &specs, SpecSet &peptides, float peakTol,
			float ionOffset, float ctermMass, vector<float> &scores,
			vector<float> &ionOffsets, vector<TwoValues<float> > &pepMatchScores,
			vector<list<unsigned int> > &pepMatchedSpecs, bool matchSuffix = true);

	/**
	 * TODO: add description
	 *
	 */
	void unit_MatchSpecsToPeps();

} // namespace specnets

#endif
