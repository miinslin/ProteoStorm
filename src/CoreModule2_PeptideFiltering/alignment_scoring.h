#ifndef ALIGNMENT_SCORING
/** */
#define ALIGNMENT_SCORING

#include <map>
#include "aminoacid.h"
#include "SpecSet.h"
#include "spectrum.h"

namespace specnets
{
	using namespace std;
	/**
	 * Finds all _pairs_ of matching peaks (i.e. peak masses within tolerance of
	 * one another).
	 *
	 *@param spec1
	 *@param spec2
	 *@param shift
	 *@param tolerance
	 *@param idx1
	 *@param idx2
	 *@return
	 */
	int FindMatchPeaksAll(const Spectrum &spec1, const Spectrum &spec2, float shift,
			float tolerance, vector<int> &idx1, vector<int> &idx2);

	/**
	 * TODO: add description
	 *
	 *@param spec1
	 *@param spec2
	 *@param shift
	 *@param tolerance
	 *@param idx1
	 *@param idx2
	 *@return
	 */
	int FindMatchPeaksAll2(const Spectrum &spec1, const Spectrum &spec2, float shift,
			float tolerance, vector<int> &idx1, vector<int> &idx2);

	
	int FindMatchPeaksAll2_efficient(const Spectrum &spec1, const Spectrum &spec2, float shift,
			float tolerance, vector<int> &idx1, vector<int> &idx2);

	/**
	 * Finds all _peaks_ that match at least one peak on the other spectrum
	 * (i.e. peak masses within tolerance of one another)
	 *
	 * **Note: that the output of this function does not have a one-to-one
	 * correspondence like FindMatchPeaksAll.
	 *
	 *@param spec1
	 *@param spec2
	 *@param shift
	 *@param tolerance
	 *@param idx1
	 *@param idx2
	 */
	void FindMatchPeaks(Spectrum &spec1, Spectrum &spec2, float shift,
			float tolerance, vector<int> &idx1, vector<int> &idx2);

	/**
	 * TODO: add description
	 *
	 *@param spec1
	 *@param spec2
	 *@param minShift
	 *@param minAAmass
	 *@param ionOffset
	 *@param peakTol
	 *@param resolution
	 *@param ctermMass
	 *@param scoredShifts
	 */
	void CyclicAlign(Spectrum &spec1, Spectrum &spec2, float minShift,
			float minAAmass, float ionOffset, float peakTol, float resolution,
			float ctermMass, vector<TwoValues<float> > &scoredShifts);


	/**
	 * Temporarily here for testing - this is supposed to be a function internal
	 * to this module.
	 *
	 *@param spec1
	 *@param spec2
	 *@param peakTol
	 *@param pmTol
	 *@param minRatio
	 *@param minPeakAreaOvlp
	 *@param minNumMatchedPeaks
	 *@param validJumps
	 *@param shiftScores
	 *@param shiftPairs
	 *@param shiftMatchedPeaks
	 *@param bestCandidateScores
	 *@param bestCandidateMP
	 *@param minAbsShift
	 */

	TwoValues<int> computeShifts(Spectrum &spec1, Spectrum &spec2, float peakTol, float pmTol,
								 float minOvlpMassSize, float minScore1, float minScore2, int minNumMatchedPeaks,
								 AAJumps &validJumps, float resolution,
								 list<float> &shiftScores,
								 list<TwoValues<unsigned int> > &shiftPairs,
								 vector<list<TwoValues<int> > > &shiftMatchedPeaks,
								 TwoValues<float> &bestCandidateScores, TwoValues<int> &bestCandidateMP,
								 float minAbsShift = 0, bool addSymmetric=true);

	TwoValues<int> computeShiftsWithFiltered(Spectrum &spec1, Spectrum &spec2, float peakTol, float pmTol,
											 float minOvlpMassSize, float minScore1, float minScore2, int minNumMatchedPeaks,
											 AAJumps &validJumps, float resolution,
											 set<int> &possibleDeltas, list<float> &shiftScores,
											 list<TwoValues<unsigned int> > &shiftPairs,
											 vector<list<TwoValues<int> > > &shiftMatchedPeaks,
											 TwoValues<float> &bestCandidateScores, TwoValues<int> &bestCandidateMP,
											 float minAbsShift = 0, bool addSymmetric=true);

	float computeBestShift(Spectrum& spec1, Spectrum& spec2, float peakTol,
	            float resolution, float shiftOffset, int minMatchedPeaks, list<TwoValues<int> >& matched);

	void computeShiftsRaw(Spectrum& spec1, Spectrum& spec2, float peakTol,
	            float resolution, float shiftOffset, map<int, list<TwoValues<int> > >& shifts, unsigned int minNumMatchedPeaks, bool usePPM = false);

	/**
	 * Cleans up the static variables used by computeShifts().
	 */
	void computeShifts_cleanup();


	/**
	 * TODO: add description
	 *
	 *@param spec1
	 *@param spec2
	 *@param shiftsList
	 *@param resolution
	 */
	void computeShifts2(Spectrum &spec1, Spectrum &spec2,
			vector<float> &shiftsList, float resolution = 0.1);

	/**
	 * TODO: add description
	 *
	 *@param spec1
	 *@param spec2
	 *@param shift
	 *@param tolerance
	 *@param idxMatched1
	 *@param idxMatched2
	 *@param minIPdist
	 *@param offsetPenalty
	 *@return
	 */
	float ScoreOverlap6(const Spectrum &spec1, const Spectrum &spec2, float shift,
			float tolerance, vector<int> &idxMatched1, vector<int> &idxMatched2,
			float minIPdist = AAJumps::minAAmass, float *offsetPenalty = 0);

	/**
	 * TODO: add description
	 *
	 *@param spec1
	 *@param idx1all
	 *@param spec2
	 *@param idx2all
	 *@param shift
	 *@param tolerance
	 *@param idxMatched1
	 *@param idxMatched2
	 *@param minIPdist
	 *@param offsetPenalty
	 *@return
	 */
	float ScoreOverlap6mp(const Spectrum &spec1, vector<int> idx1all, const Spectrum &spec2,
			vector<int> idx2all, float shift, float tolerance,
			vector<int> &idxMatched1, vector<int> &idxMatched2, float minIPdist =
					AAJumps::minAAmass, float *offsetPenalty = 0);

	/**
	 * TODO: add description
	 *
	 *@param spec1
	 *@param idx1all
	 *@param spec2
	 *@param idx2all
	 *@param shift
	 *@param tolerance
	 *@param idxMatched1
	 *@param idxMatched2
	 *@param symmetryOffset
	 *@return
	 */
	float ScoreOverlap7(Spectrum &spec1, vector<int> idx1all, Spectrum &spec2,
			vector<int> idx2all, float shift, float tolerance,
			vector<int> &idxMatched1, vector<int> &idxMatched2,
			float symmetryOffset = 0);

	/**
	 * TODO: add description
	 *
	 *@param spec
	 *@param tolerance
	 *@param pmOffset
	 *@param idxMatched
	 *@param includeSymmetric
	 *@return
	 */
	float getMaxSparseSet(Spectrum &spec, float tolerance, float pmOffset, vector<
			int> &idxMatched, bool includeSymmetric = false);


	void getALGFProbDensity(vector<float>& probDensity,
			                Spectrum &spec,
	                        int scoreMax,
	                        float peakTol,
	                        float minDistance = AAJumps::minAAmass);

	double*** getALGFProbDensity(Spectrum &spec, int scoreMax, float peakTol, float minDistance = AAJumps::minAAmass);
	void freeALGFTable(double ***dpTable, Spectrum &spec);
	double getALGFPValue(double ***dpTable, Spectrum &spec, int fwScore, int rvScore, double mass);
	double getForwardALGFPValue(double ***dpTable, Spectrum &spec, int score, double mass);
	double getReverseALGFPValue(double ***dpTable, Spectrum &spec, int score, double mass);
	double getALGFPValueForInternalRegion(Spectrum &spec, int scoreMax, float peakTol,
			float startMass, float endMass, int score, float minDistance = AAJumps::minAAmass);

	// #define DEBUG 1

	#ifdef DEBUG

	/**
	 * TODO: add description
	 *
	 *@param spec1
	 *@param spec2
	 *@param peakTol
	 *@param matched1
	 *@param matched2
	 *@param maxAAJump
	 *@param sameVertexPenalty
	 *@param ptmPenalty
	 *@param forceSymmetry
	 *@param adZPMmatches
	 *@param debug
	 *@return
	 */
	float SpectrumAlignment(Spectrum *spec1, Spectrum *spec2, float peakTol,
			Spectrum *matched1, Spectrum *matched2, int maxAAJump,
			float sameVertexPenalty=-5, float ptmPenalty=-5, bool forceSymmetry=true,
			bool addZPMmatches=false, ostream &debug=cerr);
	#else

	/**
	 * TODO: add description
	 *
	 *@param spec1
	 *@param spec2
	 *@param peakTol
	 *@param matched1
	 *@param matched2
	 *@param maxAAJump
	 *@param sameVertexPenalty
	 *@param ptmPenalty
	 *@param forceSymmetry
	 *@param addZPMmatches
	 *@return
	 */
	float SpectrumAlignment(Spectrum *spec1, Spectrum *spec2, float peakTol,
			Spectrum *matched1, Spectrum *matched2, int maxAAJump,
			float sameVertexPenalty = -5, float ptmPenalty = -5,
			bool forceSymmetry = true, bool addZPMmatches = false);
	#endif
	
        // Helper class for getPairCosines
        class GPCAux {
        public:
            float score;  // Product of peak intensities
            int i, j;     // Peak indices

            bool operator<(GPCAux & other);
        };
        bool GPCAux_cmp(const GPCAux &a, const GPCAux &b);

}

#endif
