#include "Logger.h"
#include "alignment_scoring.h"
#include "inputParams.h"
#include "aminoacid.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <algorithm>
#include <cmath>

namespace DEKEL
{

  const int None = -1;
  typedef pair<vector<int>, vector<int> > pair2;

  pair<double, pair2> align_simple(double parent_mass1,
                                   double parent_mass2,
                                   vector<double> peaks,
                                   vector<double> peaks2,
                                   vector<int> & common,
                                   vector<int> & common2,
                                   vector<double> & common_scores,
                                   vector<double> & common2_scores,
                                   vector<int> & prev,
                                   vector<int> & next,
                                   vector<vector<int> > & left_jumps,
                                   vector<vector<int> > & right_jumps,
                                   double ptm_penalty,
                                   double ambiguous_penalty);

  pair<double, pair2> align(double parent_mass1,
                            double parent_mass2,
                            vector<double> peaks,
                            vector<double> peaks2,
                            vector<int> & common,
                            vector<int> & common2,
                            vector<double> & common_scores,
                            vector<double> & common2_scores,
                            vector<int> & prev,
                            vector<int> & next,
                            vector<int> & prev2,
                            vector<int> & next2,
                            vector<vector<int> > & left_jumps,
                            vector<vector<int> > & right_jumps,
                            vector<vector<int> > & left_jumps2,
                            vector<vector<int> > & right_jumps2,
                            vector<vector<int> > & left_neighbors,
                            vector<vector<int> > & right_neighbors,
                            double same_vertex_penalty,
                            double ptm_penalty);
}

namespace specnets
{
  //
  //  FindMatchPeaksAll - like findMatchPeaksAll.m. Indices of the matched PRMs are returned in idx1 and idx2.
  //
  int FindMatchPeaksAll(const Spectrum &spec1,
                        const Spectrum &spec2,
                        float shift,
                        float tolerance,
                        vector<int> &idx1,
                        vector<int> &idx2)
  {
    int i, j;            // Iterators over the peaks indices
    int low = 0, high = 0; // Index bounds of the peaks in spec2 the lie within tolerance of current peak in spec1

    idx1.resize(0);
    idx2.resize(0);
    for (i = 0; i < (int)spec1.size(); i++)
    {

      float lowMass = spec1[i][0] - tolerance - 0.000001;
      float highMass = spec1[i][0] + tolerance + 0.000001;

      while ((low < (int)spec2.size()) && lowMass > (spec2[low][0] + shift))
        low++;
      if (low > 0) low--;
      while ((high < (int)spec2.size()) && highMass >= (spec2[high][0] + shift))
        high++;  // high is index of first unreachable peak
      if (high > 0) high--;


      for (j = low; j <= high; j++)
      {
        if (spec2[j][0] + shift > lowMass && spec2[j][0] + shift <= highMass)
        {
          idx1.push_back(i);
          idx2.push_back(j);
        }
      }
    }
    return idx1.size();
  }

  //
  //  FindMatchPeaksAll2 - like findMatchPeaksAll.m. Indices of the matched PRMs are returned in idx1 and idx2.
  //  Shift is assumed to be spec1 - spec2
  //
  int FindMatchPeaksAll2_efficient(const Spectrum &spec1,
                         const Spectrum &spec2,
                         float shift,
                         float tolerance,
                         vector<int> &idx1,
                         vector<int> &idx2)
  {

    float adj_tolerance = tolerance + 0.000001;
    vector<float> peak2_masses;
    for(int i = 0; i < (int)spec2.size(); i++){
        peak2_masses.push_back(spec2[i][0]);
    }

    int i, j;            // Iterators over the peaks indices
    int low = 0, high = 0; // Index bounds of the peaks in spec2 the lie within tolerance of current peak in spec1

    idx1.resize(0);
    idx2.resize(0);
    for (i = 0; i < (int)spec1.size(); i++)
    {
      float left_mz_bound = spec1[i][0] - shift - adj_tolerance;
      float right_mz_bound = spec1[i][0] - shift + adj_tolerance;

      std::vector<float>::iterator low_it,up_it;
      low_it = std::lower_bound (peak2_masses.begin(), peak2_masses.end(), left_mz_bound);
      up_it = std::upper_bound (peak2_masses.begin(), peak2_masses.end(), right_mz_bound);

      low = low_it - peak2_masses.begin();
      high = up_it - peak2_masses.begin();
      
      for (j = low; j < high; j++)
      {
        idx1.push_back(i);
        idx2.push_back(j);
      }
    }
    return idx1.size();
  }

  //
  //  FindMatchPeaksAll2 - like findMatchPeaksAll.m. Indices of the matched PRMs are returned in idx1 and idx2.
  //
  int FindMatchPeaksAll2(const Spectrum &spec1,
                         const Spectrum &spec2,
                         float shift,
                         float tolerance,
                         vector<int> &idx1,
                         vector<int> &idx2)
  {
    int i, j;            // Iterators over the peaks indices
    int low = 0, high = 0; // Index bounds of the peaks in spec2 the lie within tolerance of current peak in spec1

    idx1.resize(0);
    idx2.resize(0);
    for (i = 0; i < (int)spec1.size(); i++)
    {
      while (low > 0
          && (spec1[i][0] - tolerance - 0.000001) < (spec2[low][0] + shift))
        low--;
      while ((low < (int)spec2.size())
          && (spec1[i][0] - tolerance - 0.000001) > (spec2[low][0] + shift))
        low++;
      while ((high < (int)spec2.size())
          && (spec1[i][0] + tolerance + 0.000001) >= (spec2[high][0] + shift))
        high++;  // high is index of first unreachable peak
      for (j = low; j < high; j++)
      {
        idx1.push_back(i);
        idx2.push_back(j);
      }
    }
    return idx1.size();
  }

  //
  //  FindMatchPeaks - Finds all _peaks_ that match at least one peak on the
  //        other spectrum (i.e. peak masses within tolerance of one another).
  //        Note that the output of this function does not have a one-to-one
  //        correspondence like FindMatchPeaksAll.
  //        Indices of the matched PRMs are returned in idx1 and idx2.
  //
  //  NOTE: idx1 and idx2 should be pre-allocated to have enough storage and avoid resizing operations in the middle
  //
  void FindMatchPeaks(Spectrum &spec1,
                      Spectrum &spec2,
                      float shift,
                      float tolerance,
                      vector<int> &idx1,
                      vector<int> &idx2)
  {
    unsigned int i, j;            // Iterators over the peaks indices
    unsigned int low = 0, high = 0; // Index bounds of the peaks in spec2 the lie within tolerance of current peak in spec1
    vector<bool> match1(spec1.size()), match2(spec2.size());
    for (i = 0; i < spec1.size(); i++)
      match1[i] = false;
    for (i = 0; i < spec2.size(); i++)
      match2[i] = false;

    for (i = 0; i < spec1.size(); i++)
    {
      while ((low < spec2.size())
          && (spec1[i][0] - tolerance - 0.000001) > (spec2[low][0] + shift))
        low++;
      while ((high < spec2.size())
          && (spec1[i][0] + tolerance + 0.000001) >= (spec2[high][0] + shift))
        high++;  // high is index of first unreachable peak
      for (j = low; j < high; j++)
      {
        match1[i] = true;
        match2[j] = true;
      }
    }

    idx1.resize(0);
    idx2.resize(0);
    for (i = 0; i < spec1.size(); i++)
      if (match1[i])
        idx1.push_back(i);
    for (i = 0; i < spec2.size(); i++)
      if (match2[i])
        idx2.push_back(i);
  }

  void CyclicAlign(Spectrum &spec1,
                   Spectrum &spec2,
                   float minShift,
                   float minAAmass,
                   float ionOffset,
                   float peakTol,
                   float resolution,
                   float ctermMass,
                   vector<TwoValues<float> > &scoredShifts)
  {
    unsigned int pivot, peakIdx;
    Spectrum specDouble;
    specDouble = spec1;
    unsigned int numPeaks = spec1.size(), numPeaks2 = 2 * numPeaks;
    float massOffset = spec1.parentMass - ctermMass - AAJumps::massHion
        + ionOffset; // Offset that guarantees that bn overlaps with b0

    specDouble.resize(numPeaks2);
    specDouble.parentMass = specDouble.parentMass * 2 - ctermMass
        - AAJumps::massHion;
    for (peakIdx = numPeaks; peakIdx < numPeaks2; peakIdx++)
    {
      specDouble[peakIdx] = spec1[peakIdx - numPeaks];
      specDouble[peakIdx][0] += massOffset;
    }
    /*	specDouble.resize(3*spec1.size());   specDouble.parentMass=specDouble.parentMass*3-2*(ctermMass+AAJumps::massHion);
     for(peakIdx=0; peakIdx<numPeaks; peakIdx++) { specDouble[peakIdx]=spec1[peakIdx]; specDouble[peakIdx][0]-=massOffset; }
     for(peakIdx=0; peakIdx<numPeaks; peakIdx++) specDouble[numPeaks+peakIdx]=spec1[peakIdx];
     for(peakIdx=0; peakIdx<numPeaks; peakIdx++) { specDouble[numPeaks2+peakIdx]=spec1[peakIdx]; specDouble[numPeaks2+peakIdx][0]+=massOffset; }
     */
    vector<float> shiftsList;
    computeShifts2(spec1, spec2, shiftsList, resolution);
    scoredShifts.resize(shiftsList.size());

    vector<int> idx1, idx2, idxMatched1, idxMatched2;
    unsigned int shiftsIdx = 0;
    for (pivot = 0; pivot < shiftsList.size(); pivot++)
    {
      if (shiftsList[pivot] < 0)
        continue;
      scoredShifts[shiftsIdx].set(0, shiftsList[pivot]);
      FindMatchPeaksAll2(specDouble,
                         spec2,
                         shiftsList[pivot],
                         peakTol,
                         idx1,
                         idx2);
      for (peakIdx = 0; peakIdx < idx1.size(); peakIdx++)
        scoredShifts[shiftsIdx][0] += specDouble[idx1[peakIdx]][1]
            * spec2[idx2[peakIdx]][1];
      //			scoredShifts[shiftsIdx][0]+=specDouble[idx1[peakIdx]][1]+spec2[idx2[peakIdx]][1];
      //		ScoreOverlap6(specDouble, idx1, spec2, idx2, shiftsList[pivot], peakTol, idxMatched1, idxMatched2, minAAmass);
      //		for(peakIdx=0; peakIdx<idxMatched1.size(); peakIdx++)
      //			scoredShifts[shiftsIdx][0]+=specDouble[idxMatched1[peakIdx]][1]*spec2[idxMatched2[peakIdx]][1];
      //			scoredShifts[shiftsIdx][0]+=specDouble[idxMatched1[peakIdx]][1]+spec2[idxMatched2[peakIdx]][1];
      shiftsIdx++;
    }
  }

  // Vectors are too slow for this computShifts - this is the "gatekeeper" so it must be as fast as possible
  /* static bool computeShiftsFirstUse=true;
   static vector<TwoValues<float> > shiftScoresV;   // TwoValues to hold match scores per spectrum (spec1=[0],spec2=[1])
   static vector<TwoValues<int> > shiftMPcount;     // Counts the number of matched peaks per shift (spec1=[0],spec2=[1])
   static vector<unsigned int> shiftPairsV;         // index of the highest scoring symmetric shift
   //static vector<bool> validShift;                  // vector<bool> is way too slow to access
   static vector<char> validShift;                  // Used to signal whether a shift is valid (in terms of start/end-point AA offsets)
   //   and whether a shift's maximum match score was above the threshold (thus keeping matched peak lists in shiftMatchedPeaks.
   */

  //
  //  computeShifts - Computes all eligible shifts of spec2 in relation to spec1.
  //
  //  peakTol, pmTol - peak and parent mass tolerances (typically 0.5 and 1Da)
  //  minRatio - minimum acceptable ratio of matched peaks score to total spectrum score (typically 0.3-0.4)
  //  minPeakAreaOvlp  - If >0 require the overlaps to span at least this much area (peak area = largest peak mass-smallest peak mass)
  //  validJumps       - If != NULL then enforce that shifts and parent mass differences must match some value in validJumps
  //
  //  shiftScores - list of shift pair scores, sorted by decreasing shift score
  //  shiftPairs  - list of shift pairs (shift,shiftSym), both including shiftOffset, in the same order as shiftScores
  //  shiftMatchedPeaks - list of matched peaks in (spec1,spec2) for the every shift in any shiftPairs entry. Directly addressable by round(shift mass offset/InputParams::Resolution)
  //  bestCandidateScores - maximal pair of matched peak scores. Matched score in spec1/pos[0], spec2/pos[1]
  //  bestCandidateMP - maximal number of matched peaks. Num matched peaks in spec1/pos[0], spec2/pos[1], including best symmetric shift
  //
  // Returns a pair of values:
  //   pos[0]: integer shifts offset from the minimum shift considered
  //   pos[1]: integer index of the maximum eligible symmetric shift
  //
  TwoValues<int> computeShifts(Spectrum &spec1,
                               Spectrum &spec2,
                               float peakTol,
                               float pmTol,
                               float minOvlpMassSize,
                               float minScore1,
                               float minScore2,
                               int minNumMatchedPeaks,
                               AAJumps &validJumps,
                               float resolution,
                               list<float> &shiftScores,
                               list<TwoValues<unsigned int> > &shiftPairs,
                               //vector<vector<TwoValues<int> > > &shiftMatchedPeaks,
                               vector<list<TwoValues<int> > > &shiftMatchedPeaks,
                               TwoValues<float> &bestCandidateScores,
                               TwoValues<int> &bestCandidateMP,
                               float minAbsShift,
                               bool addSymmetric)
  {

    // Helper global variables for computeShifts (to avoid large resizes every time the function is called)

    unsigned int szSpec1 = spec1.size(), szSpec2 = spec2.size(), szVecs = 0; // Size of the vectors used in computeShifts. Determined by minimum/maximum shifts to consider

/*    float maxPM = max(spec1.parentMass, spec2.parentMass);
    if (minPeakAreaOvlp > 0)
      minPeakAreaOvlp = minPeakAreaOvlp * maxPM;
    if (spec2.parentMass < minPeakAreaOvlp
        || spec1.parentMass < minPeakAreaOvlp)
      return (TwoValues<int>(0, 0)); // Test if required overlap areas are feasible //*/

    int shiftsOffset = (int)ceil((spec2.parentMass + peakTol)
        / resolution);
    szVecs = 1 + shiftsOffset
        + (int)ceil((spec1.parentMass + peakTol) / resolution);

/*    unsigned int csMaxVecSize = szVecs;   // Allocated size for the arrays below
    TwoValues<float> shiftScoresV[szVecs]; // TwoValues to hold match scores per spectrum (spec1=[0],spec2=[1])
    TwoValues<int> shiftMPcount[szVecs]; // Counts the number of matched peaks per shift (spec1=[0],spec2=[1])
    unsigned int shiftPairsV[szVecs]; // index of the highest scoring symmetric shift
    char validShift[szVecs]; // Used to signal whether a shift is valid (in terms of start/end-point AA offsets)//*/

    unsigned int csMaxVecSize = szVecs;   // Allocated size for the arrays below
	vector<TwoValues<float> > shiftScoresV(szVecs); // TwoValues to hold match scores per spectrum (spec1=[0],spec2=[1])
	vector<TwoValues<int> > shiftMPcount(szVecs); // Counts the number of matched peaks per shift (spec1=[0],spec2=[1])
	vector<unsigned int> shiftPairsV(szVecs); // index of the highest scoring symmetric shift
	vector<char> validShift(szVecs); // Used to signal whether a shift is valid (in terms of start/end-point AA offsets)

    shiftMatchedPeaks.resize(szVecs);
    unsigned int leftShiftIdx = (unsigned int)(
    		minOvlpMassSize <= 2 * peakTol ? 0 :
            floor((minOvlpMassSize - 2 * peakTol) / resolution));
    unsigned int rightShiftIdx = szVecs - 1 - leftShiftIdx;

    // Initializations
    for (unsigned int i = 0; i < szVecs; i++)
    {
      shiftScoresV[i].set(0, 0);
      shiftMPcount[i].set(0, 0);
    }

/*    float totalScore1 = 0;
    for (unsigned int i = 0; i < szSpec1; i++)
      totalScore1 += spec1[i][1];
    float totalScore2 = 0;
    for (unsigned int i = 0; i < szSpec2; i++)
      totalScore2 += spec2[i][1];
    float minScore1 = minRatio * totalScore1, minScore2 = minRatio
        * totalScore2;//*/

    shiftScores.clear();
    shiftPairs.clear();
    for (unsigned int i = 0; i < szVecs; i++)
      shiftMatchedPeaks[i].clear(); //.resize(0);
    int shiftIndex, intPeakTol = (int)round(peakTol / resolution),
        intPMTol = (int)round(pmTol / resolution);

    // Populate shiftScoresV and shiftMatchedPeaks - spectra should have peaks at 0/19/PM-19/PM to match endpoints to internal peaks
    TwoValues<int> peakPair;
    for (unsigned int idxSpec1 = 0; idxSpec1 < szSpec1; idxSpec1++)
      for (unsigned int idxSpec2 = 0; idxSpec2 < szSpec2; idxSpec2++)
      {
        shiftIndex = shiftsOffset
            + (int)round((spec1[idxSpec1][0] - spec2[idxSpec2][0])
                / resolution);
        peakPair.set((int)idxSpec1, (int)idxSpec2);

        for (int tolIdx = -intPeakTol; tolIdx <= intPeakTol; tolIdx++)
        {
          int shiftIndexTol = shiftIndex + tolIdx; // float curPenalty=abs(tolIdx/10.0);
          if (shiftIndexTol >= 0 && shiftIndexTol < szVecs)
            shiftMatchedPeaks[shiftIndexTol].push_back(peakPair);
        }
      }

    // Add peak scores to all reachable shift positions
    int lastPeak;                // Index of last peak in spectrum 1 (per shift)
    char *match2 = new char[szSpec2]; // Matched peaks in spectrum 2 (per shift)
    for (unsigned int shiftIndex = leftShiftIdx; shiftIndex < rightShiftIdx;
        shiftIndex++)
    {  // Add the scores of all matched peaks without double-counting errors
      if (shiftMatchedPeaks[shiftIndex].size() > 0)
      {
        lastPeak = -1;
        for (unsigned int j = 0; j < szSpec2; j++)
          match2[j] = 0;
        list<TwoValues<int> >::iterator matchStart =
            shiftMatchedPeaks[shiftIndex].begin(), matchEnd =
            shiftMatchedPeaks[shiftIndex].end();
        for (list<TwoValues<int> >::iterator matchIter = matchStart;
            matchIter != matchEnd; matchIter++)
        {
          if ((*matchIter)[0] > lastPeak)
          {
            shiftScoresV[shiftIndex][0] += spec1[(*matchIter)[0]][1];
            shiftMPcount[shiftIndex][0]++;
            lastPeak = (*matchIter)[0];
          }
          if (!match2[(*matchIter)[1]])
          {
            shiftScoresV[shiftIndex][1] += spec2[(*matchIter)[1]][1];
            shiftMPcount[shiftIndex][1]++;
            match2[(*matchIter)[1]] = 1;
          }
        }
      }
    }
    delete[] match2;

    // Find and mark valid shifts
    float shiftMass = 0, maxMass = validJumps.masses[validJumps.size() - 1];
    //  	for(shiftIndex=0; shiftIndex<(int)validShift.size(); shiftIndex++) validShift[shiftIndex]=true;
    for (shiftIndex = 0; shiftIndex < leftShiftIdx; shiftIndex++)
      validShift[shiftIndex] = 0;
    for (shiftIndex = rightShiftIdx + 1; shiftIndex < szVecs; shiftIndex++)
      validShift[shiftIndex] = 0;

    for (shiftIndex = leftShiftIdx; shiftIndex <= (int)rightShiftIdx;
        shiftIndex++)
    {
      validShift[shiftIndex] = 1;

/*      validShift[shiftIndex] = 0;
      shiftMass = fabs((shiftIndex - shiftsOffset) * resolution);
      if (shiftMass >= minAbsShift - peakTol
          && (shiftMass <= peakTol + 0.000001
              || shiftMass > maxMass + peakTol - 0.000001
              || validJumps.isValid(shiftMass, peakTol)))
      {
        validShift[shiftIndex] = 1;
      }//*/

    }

    // Create shiftsTmp with entries observing the minRatio and minPeakAreaOvlp restrictions.
    int middleShift = (int)round((spec1.parentMass - spec2.parentMass)
        / (2 * resolution));
    int middleTimesTwo = (int)round((spec1.parentMass - spec2.parentMass)
        / resolution);
    int upperShiftLimit = shiftsOffset + middleShift;
    int shiftSym;
    list<TwoValues<float> > shiftsTmp; // List of all curShift pairs (see line below)
    TwoValues<float> curShift; // Maximum score of the shift pair (pos[0]), index of the base shift (pos[1])

    bestCandidateScores.set(0, 0);
    bestCandidateMP.set(0, 0);
    for (shiftIndex = leftShiftIdx; shiftIndex <= upperShiftLimit; shiftIndex++)
    {
      if (!validShift[shiftIndex])
        continue;
      shiftSym = shiftsOffset + (middleTimesTwo - (shiftIndex - shiftsOffset));

      // Choose highest scoring symmetric shift within tolerance
      float maxScore = 0, score1, score2, curPenalty;
      int maxScoreIdx = -1;
      int numPeaks1, numPeaks2;
      if (addSymmetric)
      {
        for (int symIdx = max(0, shiftSym - intPMTol);
            symIdx <= shiftSym + intPMTol && symIdx < (int)szVecs; symIdx++)
        {
          if (!validShift[symIdx])
            continue;
          score1 = shiftScoresV[shiftIndex][0] + shiftScoresV[symIdx][0]; // +shiftScoresV[shiftSym][0];  Bug fixed 2006/06/19
          score2 = shiftScoresV[shiftIndex][1] + shiftScoresV[symIdx][1]; // +shiftScoresV[shiftSym][1];
          numPeaks1 = shiftMPcount[shiftIndex][0] + shiftMPcount[symIdx][0];
          numPeaks2 = shiftMPcount[shiftIndex][1] + shiftMPcount[symIdx][1];
          if (score1 >= minScore1 && score2 >= minScore2
              && numPeaks1 >= minNumMatchedPeaks
              && numPeaks2 >= minNumMatchedPeaks)
          {
            curPenalty = round(abs(shiftSym - symIdx)
                * resolution);
            if (score1 + score2 - curPenalty > maxScore)
            {
              maxScore = score1 + score2 - curPenalty;
              maxScoreIdx = symIdx;
            }

            if ((score1 + score2)
                > (bestCandidateScores[0] + bestCandidateScores[1]))
              bestCandidateScores.set(score1, score2);
            if ((numPeaks1 + numPeaks2)
                > (bestCandidateMP[0] + bestCandidateMP[1]))
              bestCandidateMP.set(numPeaks1, numPeaks2);
          }
        }

      }
      else
      {
        score1 = shiftScoresV[shiftIndex][0];
        score2 = shiftScoresV[shiftIndex][1];
        numPeaks1 = shiftMPcount[shiftIndex][0];
        numPeaks2 = shiftMPcount[shiftIndex][1];
        if (score1 >= minScore1 && score2 >= minScore2
            && numPeaks1 >= minNumMatchedPeaks
            && numPeaks2 >= minNumMatchedPeaks)
        {
          if ((score1 + score2)
              > (bestCandidateScores[0] + bestCandidateScores[1]))
            bestCandidateScores.set(score1, score2);
          if ((numPeaks1 + numPeaks2)
              > (bestCandidateMP[0] + bestCandidateMP[1]))
            bestCandidateMP.set(numPeaks1, numPeaks2);
          maxScore = score1 + score2;
          maxScoreIdx = 0;
        }
      }
      if (maxScoreIdx == -1)
        continue;

      curShift.set(maxScore, shiftIndex);
      shiftsTmp.push_back(curShift);
      //		shiftPairsV[shiftIndex] = maxScoreIdx;  // This function should not decide the best symmetric shift - that's for the DP function to decide
      shiftPairsV[shiftIndex] = shiftSym;
    }

    // Sort shiftsTmp, create output structures
    shiftsTmp.sort();
    int maxUsedShiftIndex = 0;
    for (unsigned int i = 0; i < szVecs; i++)
      validShift[i] = 0;
    list<TwoValues<float> >::reverse_iterator iter = shiftsTmp.rbegin();
    for (; iter != shiftsTmp.rend(); iter++)
    {
      shiftIndex = (int)round((*iter)[1]);
      shiftSym = shiftPairsV[shiftIndex];
      if (maxUsedShiftIndex < shiftSym)
        maxUsedShiftIndex = shiftSym;
      shiftScores.push_back((*iter)[0]);
      shiftPairs.push_back(TwoValues<unsigned int>(shiftIndex, shiftSym));
      validShift[shiftIndex] = 1;
      for (int shiftSymTol = max(0, shiftSym - intPMTol);
          shiftSymTol <= shiftSym + intPMTol && shiftSymTol < (int)szVecs;
          shiftSymTol++)
        validShift[shiftSymTol] = 1;
    }

    // Clear list of matched peaks for unused shifts
    for (shiftIndex = 0; shiftIndex < szVecs; shiftIndex++)
      if (!validShift[shiftIndex])
        shiftMatchedPeaks[shiftIndex].clear(); //.resize(0);

    return TwoValues<int>(shiftsOffset, maxUsedShiftIndex);
  }

  TwoValues<int> computeShiftsWithFiltered(Spectrum &spec1,
										   Spectrum &spec2,
										   float peakTol,
										   float pmTol,
										   float minOvlpMassSize,
										   float minScore1,
										   float minScore2,
										   int minNumMatchedPeaks,
										   AAJumps &validJumps,
										   float resolution,
										   set<int> &possibleDeltas, //selected deltas by tag-alignment
										   list<float> &shiftScores,
										   list<TwoValues<unsigned int> > &shiftPairs,
										   vector<list<TwoValues<int> > > &shiftMatchedPeaks,
										   TwoValues<float> &bestCandidateScores,
										   TwoValues<int> &bestCandidateMP,
										   float minAbsShift,
										   bool addSymmetric)
  {

    // Helper global variables for computeShifts (to avoid large resizes every time the function is called)

    unsigned int szSpec1 = spec1.size(), szSpec2 = spec2.size(); // Size of the vectors used in computeShifts. Determined by minimum/maximum shifts to consider
    int shiftsOffset = (int)ceil((spec2.parentMass + peakTol) / resolution);
    int szVecs = 1 + shiftsOffset + (int)ceil((spec1.parentMass + peakTol) / resolution);
    unsigned int leftShiftIdx = (unsigned int)(minOvlpMassSize <= 2*peakTol ? 0 : floor((minOvlpMassSize - 2*peakTol) / resolution));
    unsigned int rightShiftIdx = szVecs - 1 - leftShiftIdx;
    vector<char> validShift(szVecs); // Used to signal whether a shift is valid (in terms of start/end-point AA offsets)

    int shiftIndex;//common var;
    // Find and mark valid shifts from deltas by tag-alignment
    for (shiftIndex = 0; shiftIndex < szVecs; shiftIndex++) validShift[shiftIndex] = 0;

    int iso_err_buf = 3;
    int valid_dt_count=0, max_delta=-szVecs, min_delta=szVecs;
    for (set<int>::iterator itr=possibleDeltas.begin(); itr!=possibleDeltas.end(); itr++)
    {
    	int sci = (*itr-iso_err_buf)/resolution + shiftsOffset;
    	int eci = (*itr+iso_err_buf)/resolution + shiftsOffset;
    	if( eci < leftShiftIdx || rightShiftIdx < sci ) continue;

    	if( sci < leftShiftIdx ) sci = leftShiftIdx;
    	if( rightShiftIdx < eci ) eci = rightShiftIdx;

    	for( shiftIndex=sci; shiftIndex<=eci; shiftIndex++ ){
    		validShift[shiftIndex] = 1;
    	}
    	valid_dt_count++;
    	if( max_delta < *itr ) max_delta = *itr;
    	if( *itr < min_delta ) min_delta = *itr;
    }
    if( valid_dt_count == 0 ) return (TwoValues<int>(0, 0));
    max_delta += iso_err_buf;
    min_delta -= iso_err_buf;

   	vector<TwoValues<float> > shiftScoresV(szVecs); // TwoValues to hold match scores per spectrum (spec1=[0],spec2=[1])
   	vector<TwoValues<int> > shiftMPcount(szVecs); // Counts the number of matched peaks per shift (spec1=[0],spec2=[1])
   	vector<unsigned int> shiftPairsV(szVecs); // index of the highest scoring symmetric shift
    shiftMatchedPeaks.resize(szVecs);

    // Initializations
	for (unsigned int i = 0; i < szVecs; i++)
	{
	  shiftScoresV[i].set(0, 0);
	  shiftMPcount[i].set(0, 0);
	}

	shiftScores.clear();
	shiftPairs.clear();
	for (unsigned int i = 0; i < szVecs; i++)
	  shiftMatchedPeaks[i].clear(); //.resize(0);
	int intPeakTol = (int)round(peakTol / resolution), intPMTol = (int)round(pmTol / resolution);

    // Populate shiftScoresV and shiftMatchedPeaks - spectra should have peaks at 0/19/PM-19/PM to match endpoints to internal peaks
    TwoValues<int> peakPair;
    for (unsigned int idxSpec1 = 0; idxSpec1 < szSpec1; idxSpec1++)
    {
      for (unsigned int idxSpec2 = 0; idxSpec2 < szSpec2; idxSpec2++)
      {
    	float delta = spec1[idxSpec1][0] - spec2[idxSpec2][0];
    	if( delta < min_delta ) break;
    	if( max_delta < delta ) continue;

        shiftIndex = shiftsOffset + (int)round(delta/resolution);

        if (!validShift[shiftIndex]) continue;

        peakPair.set((int)idxSpec1, (int)idxSpec2);
        for (int tolIdx = -intPeakTol; tolIdx <= intPeakTol; tolIdx++)
        {
          int shiftIndexTol = shiftIndex + tolIdx; // float curPenalty=abs(tolIdx/10.0);
          if (shiftIndexTol >= 0 && shiftIndexTol < szVecs)
            shiftMatchedPeaks[shiftIndexTol].push_back(peakPair);
        }
      }
    }

    // Add peak scores to all reachable shift positions
    int lastPeak;                // Index of last peak in spectrum 1 (per shift)
    char *match2 = new char[szSpec2]; // Matched peaks in spectrum 2 (per shift)
    for (unsigned int shiftIndex = leftShiftIdx; shiftIndex < rightShiftIdx; shiftIndex++)
    {  // Add the scores of all matched peaks without double-counting errors
      if (shiftMatchedPeaks[shiftIndex].size() > 0)
      {
        lastPeak = -1;
        for (unsigned int j = 0; j < szSpec2; j++)
          match2[j] = 0;
        list<TwoValues<int> >::iterator matchStart =
            shiftMatchedPeaks[shiftIndex].begin(), matchEnd =
            shiftMatchedPeaks[shiftIndex].end();
        for (list<TwoValues<int> >::iterator matchIter = matchStart;
            matchIter != matchEnd; matchIter++)
        {
          if ((*matchIter)[0] > lastPeak)
          {
            shiftScoresV[shiftIndex][0] += spec1[(*matchIter)[0]][1];
            shiftMPcount[shiftIndex][0]++;
            lastPeak = (*matchIter)[0];
          }
          if (!match2[(*matchIter)[1]])
          {
            shiftScoresV[shiftIndex][1] += spec2[(*matchIter)[1]][1];
            shiftMPcount[shiftIndex][1]++;
            match2[(*matchIter)[1]] = 1;
          }
        }
      }
    }
    delete[] match2;

    // Create shiftsTmp with entries observing the minRatio and minPeakAreaOvlp restrictions.
    int middleShift = (int)round((spec1.parentMass - spec2.parentMass)
        / (2 * resolution));
    int middleTimesTwo = (int)round((spec1.parentMass - spec2.parentMass)
        / resolution);
    int upperShiftLimit = shiftsOffset + middleShift;
    int shiftSym;
    list<TwoValues<float> > shiftsTmp; // List of all curShift pairs (see line below)
    TwoValues<float> curShift; // Maximum score of the shift pair (pos[0]), index of the base shift (pos[1])

    bestCandidateScores.set(0, 0);
    bestCandidateMP.set(0, 0);
    for (shiftIndex = leftShiftIdx; shiftIndex <= upperShiftLimit; shiftIndex++)
    {
      if (!validShift[shiftIndex]) continue;

      shiftSym = shiftsOffset + (middleTimesTwo - (shiftIndex - shiftsOffset));

      // Choose highest scoring symmetric shift within tolerance
      float maxScore = 0, score1, score2, curPenalty;
      int maxScoreIdx = -1;
      int numPeaks1, numPeaks2;
      if (addSymmetric)
      {
        for (int symIdx = max(0, shiftSym - intPMTol);
            symIdx <= shiftSym + intPMTol && symIdx < (int)szVecs; symIdx++)
        {
          if (!validShift[symIdx]) continue;

          score1 = shiftScoresV[shiftIndex][0] + shiftScoresV[symIdx][0]; // +shiftScoresV[shiftSym][0];  Bug fixed 2006/06/19
          score2 = shiftScoresV[shiftIndex][1] + shiftScoresV[symIdx][1]; // +shiftScoresV[shiftSym][1];
          numPeaks1 = shiftMPcount[shiftIndex][0] + shiftMPcount[symIdx][0];
          numPeaks2 = shiftMPcount[shiftIndex][1] + shiftMPcount[symIdx][1];
          if (score1 >= minScore1 && score2 >= minScore2
              && numPeaks1 >= minNumMatchedPeaks
              && numPeaks2 >= minNumMatchedPeaks)
          {
            curPenalty = round(abs(shiftSym - symIdx)
                * resolution);
            if (score1 + score2 - curPenalty > maxScore)
            {
              maxScore = score1 + score2 - curPenalty;
              maxScoreIdx = symIdx;
            }

            if ((score1 + score2)
                > (bestCandidateScores[0] + bestCandidateScores[1]))
              bestCandidateScores.set(score1, score2);
            if ((numPeaks1 + numPeaks2)
                > (bestCandidateMP[0] + bestCandidateMP[1]))
              bestCandidateMP.set(numPeaks1, numPeaks2);
          }
        }

      }
      else
      {
        score1 = shiftScoresV[shiftIndex][0];
        score2 = shiftScoresV[shiftIndex][1];
        numPeaks1 = shiftMPcount[shiftIndex][0];
        numPeaks2 = shiftMPcount[shiftIndex][1];
        if (score1 >= minScore1 && score2 >= minScore2
            && numPeaks1 >= minNumMatchedPeaks
            && numPeaks2 >= minNumMatchedPeaks)
        {
          if ((score1 + score2)
              > (bestCandidateScores[0] + bestCandidateScores[1]))
            bestCandidateScores.set(score1, score2);
          if ((numPeaks1 + numPeaks2)
              > (bestCandidateMP[0] + bestCandidateMP[1]))
            bestCandidateMP.set(numPeaks1, numPeaks2);
          maxScore = score1 + score2;
          maxScoreIdx = 0;
        }
      }
      if (maxScoreIdx == -1)
        continue;

      curShift.set(maxScore, shiftIndex);
      shiftsTmp.push_back(curShift);
      //		shiftPairsV[shiftIndex] = maxScoreIdx;  // This function should not decide the best symmetric shift - that's for the DP function to decide
      shiftPairsV[shiftIndex] = shiftSym;
    }

    // Sort shiftsTmp, create output structures
    shiftsTmp.sort();
    int maxUsedShiftIndex = 0;
    for (unsigned int i = 0; i < szVecs; i++)
      validShift[i] = 0;

    list<TwoValues<float> >::reverse_iterator iter = shiftsTmp.rbegin();
    for (; iter != shiftsTmp.rend(); iter++)
    {
      shiftIndex = (int)round((*iter)[1]);
      shiftSym = shiftPairsV[shiftIndex];
      if (maxUsedShiftIndex < shiftSym)
        maxUsedShiftIndex = shiftSym;
      shiftScores.push_back((*iter)[0]);
      shiftPairs.push_back(TwoValues<unsigned int>(shiftIndex, shiftSym));
      validShift[shiftIndex] = 1;
      for (int shiftSymTol = max(0, shiftSym - intPMTol);
          shiftSymTol <= shiftSym + intPMTol && shiftSymTol < (int)szVecs;
          shiftSymTol++)
        validShift[shiftSymTol] = 1;
    }

    // Clear list of matched peaks for unused shifts
    for (shiftIndex = 0; shiftIndex < szVecs; shiftIndex++)
      if (!validShift[shiftIndex])
        shiftMatchedPeaks[shiftIndex].clear(); //.resize(0);

    return TwoValues<int>(shiftsOffset, maxUsedShiftIndex);
  }


  float computeBestShift(Spectrum& spec1,
                         Spectrum& spec2,
                         float peakTol,
                         float resolution,
                         float shiftOffset,
                         int minMatchedPeaks,
                         list<TwoValues<int> >& matched)
  {
    matched.clear();

    map<int, list<TwoValues<int> > > shifts;
    list<TwoValues<int> > shiftMP;
    TwoValues<int> MP;
    map<int, float> shiftScore;
    int intPkTol;
    float bestShift;
    float bestScore = 0;
    intPkTol = 2
        * (int)(round(peakTol / resolution) + 0.01);

    for (int i = 0; i < spec1.size(); i++)
    {
      for (int j = 0; j < spec2.size(); j++)
      {
        float shift = spec1[i][0] - spec2[j][0] + shiftOffset;
        int intShift = (int)(round(shift / resolution) + 0.01);
        MP[0] = i;
        MP[1] = j;

        for (int sShift = intShift - intPkTol; sShift <= intShift + intPkTol;
            sShift++)
        {
          float score = spec1[i][1] + spec2[j][1]
              - (((float)abs(sShift - intShift)) * 0.001);
          if (shiftScore.count(sShift) > 0)
          {
            shiftScore[sShift] += score;
            shifts[sShift].push_back(MP);
          }
          else
          {
            shiftScore[sShift] = score;
            shiftMP.clear();
            shiftMP.push_back(MP);
            shifts[sShift] = shiftMP;
          }

          if (shiftScore[sShift] > bestScore
              && shifts[sShift].size() >= minMatchedPeaks)
          {
            bestScore = shiftScore[sShift];
            bestShift = ((float)sShift) * resolution;
            matched = shifts[sShift];
          }
        }
      }
    }
    return bestShift;
  }

  void computeShiftsRaw(Spectrum& spec1,
                        Spectrum& spec2,
                        float peakTol,
                        float resolution,
                        float shiftOffset,
                        map<int, list<TwoValues<int> > >& bestShifts,
                        unsigned int minNumMatchedPeaks,
                        bool usePPM)
  {
    bestShifts.clear();
    map<int, list<TwoValues<int> > > shifts;
    list<TwoValues<int> > shiftMP;
    TwoValues<int> MP;
    map<int, float> shiftScore;
    int intPkTol;
    if (!usePPM)
    {
      intPkTol = 2
          * (int)(round(peakTol / resolution) + 0.01);
    }

    //cout << "finding shifts ... "; cout.flush();
    for (int i = 0; i < spec1.size(); i++)
    {
      for (int j = 0; j < spec2.size(); j++)
      {
        float shift = spec1[i][0] - spec2[j][0] + shiftOffset;
        int intShift = (int)(round(shift / resolution) + 0.01);
        MP[0] = i;
        MP[1] = j;

        if (usePPM)
        {
          float shiftol = (InputParams::PPM * shift) / 1000000.0;
          float resol = getResolution(shiftol);
          intPkTol = (int)(round(shiftol / resol) + 0.01);
        }

        for (int sShift = intShift - intPkTol; sShift <= intShift + intPkTol;
            sShift++)
        {
          float score = spec1[i][1] + spec2[j][1]
              - (((float)abs(sShift - intShift)) * 0.001);
          if (shiftScore.count(sShift) > 0)
          {
            shiftScore[sShift] += score;
            shifts[sShift].push_back(MP);
          }
          else
          {
            shiftScore[sShift] = score;
            shiftMP.clear();
            shiftMP.push_back(MP);
            shifts[sShift] = shiftMP;
          }
        }
      }
    }

    /*
     cout << "\n\nShifts\n";
     int mx = 0;
     for (map<int, float>::iterator shiftScoreIt = shiftScore.begin(); shiftScoreIt != shiftScore.end(); shiftScoreIt ++) {
     int shift = shiftScoreIt->first;
     float score = shiftScoreIt->second;
     cout << shift << " - " << score << " : ";
     if (shifts[shift].size() > mx) mx = shifts[shift].size();
     for (list<TwoValues<int> >::iterator it = shifts[shift].begin(); it != shifts[shift].end(); it++) {
     cout << spec1[(*it)[0]][0] << "," << spec2[(*it)[1]][0] << "(" << (*it)[0] << "," << (*it)[1] << "); ";
     }
     cout << "\n";
     }

     cout << "\nMAX = " << mx << "\n\n";
     */

    // only keep shifts that are at the center of their resolution distribution
    float lastScore = -1.0;
    int lastShift;
    bool increasing = true, decreasing = false, peaking = false;

    list<int> prevShifts;
    //cout << "finished, picking top shifts ... "; cout.flush();
    for (map<int, float>::iterator shiftScoreIt = shiftScore.begin();
        shiftScoreIt != shiftScore.end(); shiftScoreIt++)
    {
      int shift = shiftScoreIt->first;
      float score = shiftScoreIt->second;

      if (shifts[shift].size() < minNumMatchedPeaks)
      {
        increasing = true;
        decreasing = false;
        peaking = false;
        continue;
      }

      if (lastScore > 0 && score > lastScore && lastShift + 1 == shift)
      {
        increasing = true;
        decreasing = false;
        peaking = false;
      }
      else if (lastScore > 0 && score < lastScore && lastShift + 1 == shift)
      {
        if (increasing || peaking)
        {
          prevShifts.push_back(lastShift);
        }
        decreasing = true;
        increasing = false;
        peaking = false;
      }
      else if (lastScore > 0 && score == lastScore && lastShift + 1 == shift)
      {
        increasing = false;
        decreasing = false;
        peaking = true;
        prevShifts.push_back(lastShift);
      }
      else
      {
        if (increasing || peaking)
        {
          prevShifts.push_back(lastShift);
        }
        increasing = true;
        decreasing = false;
        peaking = false;
      }

      if (!peaking && prevShifts.size() > 0)
      {
        float bestShift = 0;
        for (list<int>::iterator lit = prevShifts.begin();
            lit != prevShifts.end(); lit++)
        {
          bestShift += *lit;
        }
        bestShift /= (float)prevShifts.size();
        int bestShiftInt = (int)(round(bestShift) + 0.01);
        bestShifts[bestShiftInt] = shifts[bestShiftInt];
        prevShifts.clear();
      }

      lastScore = score;
      lastShift = shift;
    }

    /*
     cout << "\n\nBest Shifts\n";
     mx = 0;
     for (map<int, list<TwoValues<int> > >::iterator shiftScoreIt = bestShifts.begin(); shiftScoreIt != bestShifts.end(); shiftScoreIt ++) {
     int shift = shiftScoreIt->first;
     float score = shiftScore[shift];
     cout << shift << " - " << score << " : ";
     if (bestShifts[shift].size() > mx) mx = bestShifts[shift].size();
     for (list<TwoValues<int> >::iterator it = bestShifts[shift].begin(); it != bestShifts[shift].end(); it++) {
     cout << spec1[(*it)[0]][0] << "," << spec2[(*it)[1]][0] << "(" << (*it)[0] << "," << (*it)[1] << "); ";
     }
     cout << "\n";
     }

     cout << "\nMAX = " << mx << "\n\n";
     cout.flush();
     */
    //cout << "finished, "; cout.flush();
  }

  // Helper global variables for computeShifts2
  static bool computeShifts2firstUse = true;
  void computeShifts2(Spectrum &spec1,
                      Spectrum &spec2,
                      vector<float> &shiftsList,
                      float resolution)
  {
    if (computeShifts2firstUse)
    {
      computeShifts2firstUse = false;
    }
    shiftsList.resize(spec1.size() * spec2.size() + 1);
    unsigned int shiftsIdx = 0, spec1idx, spec2idx;
    shiftsList[shiftsIdx++] = 0;
    //	for(spec1idx=0; spec1idx<spec1.size(); spec1idx++) shiftsList[shiftsIdx++]=spec1[spec1idx][0]; // One peak is always matched
    //	for(spec1idx=0; spec1idx<spec1.size(); spec1idx++) shiftsList[shiftsIdx++]=spec1[spec1idx][0]-spec2.parentMass; // One peak is always matched
    for (spec1idx = 0; spec1idx < spec1.size(); spec1idx++)
      for (spec2idx = 0; spec2idx < spec2.size(); spec2idx++)
        shiftsList[shiftsIdx++] = (spec1[spec1idx][0] - spec2[spec2idx][0]);

    // Remove duplicate shift entries
    sort(shiftsList.begin(), shiftsList.end());
    unsigned int shiftsIdxUnique = 0;
    for (shiftsIdx = 1; shiftsIdx < shiftsList.size(); shiftsIdx++)
      if (fabs(shiftsList[shiftsIdx] - shiftsList[shiftsIdxUnique])
          > resolution)
      {
        shiftsIdxUnique++;
        if (shiftsIdx > shiftsIdxUnique)
          shiftsList[shiftsIdxUnique] = shiftsList[shiftsIdx];
      }
    shiftsList.resize(shiftsIdxUnique + 1);
  }

  //
  //  ScoreOverlap6 - Front-end function for ScoreOverlap6
  //
  //  minIPdist - minimum allowed mass distance between 2 consecutive peaks (defaults to AAJumps::minAAmass)
  //  idxMatched       - Indices of the matched PRMs (sparse set), col 0 for spec1 and col 1 for spec2
  //
  float ScoreOverlap6(const Spectrum &spec1,
                      const Spectrum &spec2,
                      float shift,
                      float tolerance,
                      vector<int> &idxMatched1,
                      vector<int> &idxMatched2,
                      float minIPdist,
                      float *offsetPenalty)
  {
    vector<int> idx1all, idx2all;
    FindMatchPeaksAll(spec1, spec2, shift, tolerance, idx1all, idx2all);

    return ScoreOverlap6mp(spec1,
                           idx1all,
                           spec2,
                           idx2all,
                           shift,
                           tolerance,
                           idxMatched1,
                           idxMatched2,
                           minIPdist,
                           offsetPenalty);
  }

  //
  //  ScoreOverlap6mp - like ScoreOverlapE.m with distance function 6. Min distance between PRMs is 57 - 2*tolerance
  //
  //  minIPdist - minimum allowed mass distance between 2 consecutive peaks (defaults to AAJumps::minAAmass)
  //  idx1all, idx2all - as returned by FindMatchPeaksAll
  //  idxMatched       - Indices of the matched PRMs (sparse set), col 0 for spec1 and col 1 for spec2
  //
  float ScoreOverlap6mp(const Spectrum &spec1,
                        vector<int> idx1all,
                        const Spectrum &spec2,
                        vector<int> idx2all,
                        float shift,
                        float tolerance,
                        vector<int> &idxMatched1,
                        vector<int> &idxMatched2,
                        float minIPdist,
                        float *offsetPenalty)
  {
    vector<TwoValues<float> > values(min(idx1all.size(), idx2all.size()) + 1); // Keeps the values for the dynamic programming recursion: predecessor (col 0) and predecessor score (col 1)
    // +1 because first line is (0,0) for DP initialization
    int forbiddenIdx; // Index of the first forbidden PRM (because it is too close to the current PRM
    int maxIdx;     // Index of the best PRM that the current PRM can connect to
    float bestMatch;   // Best path score so far
    int bestMatchIdx;  // Index of the last PRM in the best path so far
    int i, j;         // Iterator vars

    idxMatched1.resize(0);
    idxMatched2.resize(0);

    values[0][0] = 0;
    values[0][1] = 0;
    maxIdx = 0;
    forbiddenIdx = 1;
    bestMatch = 0;
    bestMatchIdx = 0;
    for (i = 1; i < (int)values.size(); i++)
    {
      while (forbiddenIdx < i
          && spec1[idx1all[forbiddenIdx - 1]][0]
              <= (spec1[idx1all[i - 1]][0] - minIPdist + 2 * tolerance)
          && spec2[idx2all[forbiddenIdx - 1]][0]
              <= (spec2[idx2all[i - 1]][0] - minIPdist + 2 * tolerance))
      {
        // This is executed only when forbidden is advanced
        if (values[forbiddenIdx][1] > values[maxIdx][1])
        {
          maxIdx = forbiddenIdx;
        }
        forbiddenIdx++;
      }
      values[i][0] = maxIdx;
      values[i][1] = values[maxIdx][1] + spec1[idx1all[i - 1]][1]
          + spec2[idx2all[i - 1]][1];
      if (offsetPenalty)
        values[i][1] -= abs(spec1[idx1all[i - 1]][0]
            - (spec2[idx2all[i - 1]][0] + shift));
      if (values[i][1] > bestMatch)
      {
        bestMatch = values[i][1];
        bestMatchIdx = i;
      }  // Keep track of where the best path ends
    }

    list<int> bestPath;
    while (bestMatchIdx > 0)
    {
      if (offsetPenalty)
        (*offsetPenalty) += abs(spec1[idx1all[bestMatchIdx - 1]][0]
            - (spec2[idx2all[bestMatchIdx - 1]][0] + shift));
      bestPath.push_back(bestMatchIdx - 1);
      bestMatchIdx = (int)values[bestMatchIdx][0];
    }

    // ******************************************
    //    Populate final lists of matched PRMs
    // ******************************************
    unsigned int bestPathLength = bestPath.size();
    idxMatched1.resize(bestPathLength);
    idxMatched2.resize(bestPathLength);

    for (unsigned int idxPath = 0; idxPath < bestPathLength; idxPath++)
    {
      idxMatched1[idxPath] = idx1all[bestPath.back()];
      idxMatched2[idxPath] = idx2all[bestPath.back()];
      bestPath.pop_back();
    }

    return bestMatch;
  }

  //
  //  ScoreOverlap7 - like ScoreOverlapE.m with distance function 7 but gets matching PRMs from findMatchPeaksAll instead of findMatchPeaks2a. Min distance between PRMs is 57 - 2*tolerance
  //
  //  idx1all, idx2all - as returned by FindMatchPeaksAll
  //  idxMatched1,2    - Indices of the matched PRMs (sparse set) for spec1 and for spec2
  //  symmetryOffset   - Masses of symmetric peaks should add up to sum(peptide masses)+18+symmetryOffset
  //
  float ScoreOverlap7(Spectrum &spec1,
                      vector<int> idx1all,
                      Spectrum &spec2,
                      vector<int> idx2all,
                      float shift,
                      float tolerance,
                      vector<int> &idxMatched1,
                      vector<int> &idxMatched2,
                      float symmetryOffset)
  {
    Spectrum tmpSpec;
    vector<int> idxMatched;
    float score;

    // Use FindMatchPeaksAll results to construct a single spectrum with all possible peak matches
    tmpSpec.parentMass = (spec1.parentMass + spec2.parentMass) / 2;
    tmpSpec.resize(idx1all.size());
    tmpSpec.idDist = spec1.idDist;

    //cerr<<"tmpSpec: (parent mass = "<<tmpSpec.parentMass<<")\n";
    for (unsigned int i = 0; i < tmpSpec.size(); i++)
    {
   // 	tmpSpec[i].set((spec1[idx1all[i]][0] + spec2[idx2all[i]][0] + shift) / 2,
    	tmpSpec[i].set((spec1[idx1all[i]][0] + spec2[idx2all[i]][0]) / 2,
                     spec1[idx1all[i]][1] + spec2[idx2all[i]][1]);
      //cerr<<i<<": "<<spec1[idx1all[i]][0]<<"+"<<spec2[idx2all[i]][0]+shift<<" -> ["<<tmpSpec[i][0]<<","<<tmpSpec[i][1]<<"]\n";
    }

    // Use getMaxSparseSet to select which peak matches to keep and determine match score
    //	score = getMaxSparseSet(tmpSpec, tolerance, 0, idxMatched);
    score = getMaxSparseSet(tmpSpec,
                            tolerance,
                            symmetryOffset,
                            idxMatched,
                            true);

    /*cerr << "tmpSpec's matched peaks:\n";
     for(int i=0;i<idxMatched.size();i++)
     cerr <<idxMatched[i]<<"\t"<<tmpSpec[idxMatched[i]][0]<<"\t"<<tmpSpec[idxMatched[i]][1]<<"\n";
     */
    idxMatched1.resize(idxMatched.size());
    idxMatched2.resize(idxMatched.size());
    for (unsigned int i = 0; i < idxMatched.size(); i++)
    {
      idxMatched1[i] = idx1all[idxMatched[i]];
      idxMatched2[i] = idx2all[idxMatched[i]];
    }

    return score;
  }

  //
  //  getMaxSparseSet - like ScoreOverlap7 but takes only one spectrum as input
  //
  //  pmOffset=0 for PRM spectra and pmOffset=2 for MS/MS spectra.
  //  idxMatched - Indices of the matched PRMs (sparse set)
  //  includeSymmetric - set to true if the sparse set should also include the symmetric
  //                      peaks of the matched peaks (e.g. ScoreOverlap7)
  //
  //  NOTE: Make sure that spec.idDist is adequately set.
  //
  float getMaxSparseSet(Spectrum &spec,
                        float tolerance,
                        float pmOffset,
                        vector<int> &idxMatched,
                        bool includeSymmetric)
  {
    vector<vector<TwoValues<int> > > prevPair, // Previous pair for every possible pair
        bestPair;     // Best pair so far
    vector<vector<float> > scores,         // Scores for every possible path
        bestScores;     // Best path score so far
    vector<float> peakScores; // Adjusted peak scores (including the score of the symmetric peaks)
    vector<int> peakPairs; // Indices of the other paired peak (if any) used in peakScores
    vector<TwoValues<int> > prefixPRMs; // PRMs with mass<=aaMass/2 and sorted by increasing mass
                                        //  Col 1 is PRM index, col 2 is index of closest PRM >=57-2*tolerance Da away
    vector<TwoValues<int> > suffixPRMs; // PRMs with mass>aaMass/2 and sorted by increasing distance to aaMass
    TwoValues<int> globalBestPair(0, 0);
    float aaMass = spec.parentMass
        + (pmOffset - spec.idDist) * AAJumps::massHion, minInterPeakDist = 57
        * spec.idDist - 2 * tolerance, globalBestScore = 0;
    int i, j, k, p, idxPref, idxSuff;

    if (spec.size() == 0)
    {
      idxMatched.resize(0);
      return 0;
    }

    // Make sure that the spectrum has PRMs at zero and at aaMass
    short addZero = 0, addPM = 0;
    Spectrum oldSpec;    // oldSpec is used to keep a copy of the input spectrum
                         //   whenever addZero or addPM are >0
    if (spec[0][0] > tolerance)
      addZero++;
    if (spec[spec.size() - 1][0]
        < aaMass - AAJumps::massH2O * spec.idDist - tolerance)
      addPM++;
    if (addZero + addPM > 0)
    {
      oldSpec = spec;
      spec.resize(spec.size() + addZero + addPM);
      if (addZero)
      {
        for (i = spec.size() - 1 - addPM; i >= 1; i--)
          spec[i] = spec[i - 1];
        spec[0].set(0, 0);
      }
      if (addPM)
        spec[spec.size() - 1].set(aaMass - AAJumps::massH2O * spec.idDist, 0);
    }

    //for(unsigned int peakIdx=0; peakIdx<spec.size(); peakIdx++)
    //	cerr<<peakIdx<<": "<<spec[peakIdx][0]<<", "<<spec[peakIdx][1]<<endl;

    // Populate prefixPRMs and suffixPRMs - peak index (col.1), predecessor index (col.2)
    prefixPRMs.resize(spec.size());
    for (i = 0; i < (int)spec.size() && spec[i][0] <= aaMass / 2; i++)
    {
      for (p = i - 1; p >= 0 && spec[i][0] - spec[p][0] < minInterPeakDist; p--)
        ;
      prefixPRMs[i].set(i, max((int)p, 0));
    }
    prefixPRMs.resize(i);
    suffixPRMs.resize(spec.size() - prefixPRMs.size());
    suffixPRMs[suffixPRMs.size() - 1].set(spec.size() - 1, 0);
    for (j = spec.size() - 1, k = 0; j >= 0 && spec[j][0] > aaMass / 2; j--)
    {
      suffixPRMs[k][0] = j;
      for (p = 0;
          suffixPRMs[p][0] > j
              && spec[suffixPRMs[p][0]][0] - spec[j][0] > minInterPeakDist; p++)
        ;
      suffixPRMs[k++].set(j, max(p - 1, 0));
    }

    // Resize and initialize all the DP variables
    prevPair.resize(prefixPRMs.size());
    bestPair.resize(prefixPRMs.size());
    scores.resize(prefixPRMs.size());
    bestScores.resize(prefixPRMs.size());
    for (int i = 0; i < (int)prefixPRMs.size(); i++)
    {
      prevPair[i].resize(suffixPRMs.size());
      bestPair[i].resize(suffixPRMs.size());
      scores[i].resize(suffixPRMs.size());
      bestScores[i].resize(suffixPRMs.size());
      for (j = 0; j < suffixPRMs.size(); j++)
      {
        scores[i][j] = 0;
        bestScores[i][j] = 0;
      }
    }

    // Consolidate the scores of symmetric peaks
    peakScores.resize(spec.size());
    peakPairs.resize(spec.size());
    for (k = 0; k < spec.size(); k++)
    {
      peakScores[k] = spec[k][1];
      peakPairs[k] = -1;
    }
    vector<vector<float> > pairs;
    vector<vector<int> > pairsIdx;
    spec.getSymmetricPeakPairs(pmOffset, tolerance, pairs, pairsIdx);
    for (k = 0; k < pairsIdx.size(); k++)
    {
      if (pairsIdx[k][0] < 0 || pairsIdx[k][1] < 0)
        continue;
      if (peakScores[pairsIdx[k][0]] < pairs[k][2])
      {
        peakScores[pairsIdx[k][0]] = pairs[k][2];
        peakPairs[pairsIdx[k][0]] = pairsIdx[k][1];
      }
      if (peakScores[pairsIdx[k][1]] < pairs[k][2])
      {
        peakScores[pairsIdx[k][1]] = pairs[k][2];
        peakPairs[pairsIdx[k][1]] = pairsIdx[k][0];
      }
    }

    //
    // DP to compute maximum sparse set
    //
    scores[0][0] = peakScores[prefixPRMs[0][0]] + peakScores[suffixPRMs[0][0]];
    bestScores[0][0] = scores[0][0];
    bestPair[0][0].set(0, 0);
    globalBestScore = scores[0][0];
    globalBestPair.set(0, 0);
    for (i = 0; i < prefixPRMs.size(); i++)
      for (j = 0; j < suffixPRMs.size(); j++)
      {
        if ((i == 0 && j == 0)
            || spec[suffixPRMs[j][0]][0] - spec[prefixPRMs[i][0]][0]
                < minInterPeakDist)
          continue;

        // Set default values of best scores/pairs for position [i][j]
        if (i > 0)
        {
          bestScores[i][j] = bestScores[i - 1][j];
          bestPair[i][j] = bestPair[i - 1][j];
        }
        if (j > 0 && bestScores[i][j - 1] > bestScores[i][j])
        {
          bestScores[i][j] = bestScores[i][j - 1];
          bestPair[i][j] = bestPair[i][j - 1];
        }

        idxPref = prefixPRMs[i][0];
        idxSuff = suffixPRMs[j][0];
        // j ranges over suffixes whose masses differ from spec[i][0]
        if (fabs(spec[idxPref][0] + spec[idxSuff][0] - aaMass) > 2 * tolerance)
        {
          //cerr<<" --- "<<spec[idxPref][0]+spec[idxSuff][0]<<", "<<aaMass<<endl;
          if (spec[idxPref][0] > aaMass - spec[idxSuff][0])
          {  // last jump was on the prefix side
            scores[i][j] = peakScores[idxPref]
                + bestScores[prefixPRMs[i][1]][j];
            prevPair[i][j] = bestPair[prefixPRMs[i][1]][j];
            //cerr<<"["<<i<<","<<j<<"] jumping from ["<<prefixPRMs[i][1]<<","<<j<<"], score "<<scores[i][j]<<"\n";
          }
          else
          {    // last jump was on the suffix side
            scores[i][j] = peakScores[idxSuff]
                + bestScores[i][suffixPRMs[j][1]];
            prevPair[i][j] = bestPair[i][suffixPRMs[j][1]];
            //cerr<<"["<<i<<","<<j<<"] jumping from ["<<i<<","<<suffixPRMs[j][1]<<"], score "<<scores[i][j]<<"\n";
          }
        }
        else
        {  // still consider these pairs but don't increase the score
          scores[i][j] = bestScores[prefixPRMs[i][1]][j];
          prevPair[i][j] = bestPair[prefixPRMs[i][1]][j];
          if (scores[i][j] < bestScores[i][suffixPRMs[j][1]])
          {
            scores[i][j] = bestScores[i][suffixPRMs[j][1]];
            prevPair[i][j] = bestPair[i][suffixPRMs[j][1]];
          }
          //cerr<<"["<<i<<","<<j<<"] (same!) jumping from ["<<i<<","<<suffixPRMs[j][1]<<"] or ["<<prefixPRMs[i][1]<<","<<j<<"], score "<<scores[i][j]<<"\n";
        }

        if (scores[i][j] > bestScores[i][j])
        {
          bestScores[i][j] = scores[i][j];
          bestPair[i][j].set(i, j);
        }
        if (bestScores[i][j] > globalBestScore)
        {
          globalBestScore = bestScores[i][j];
          globalBestPair.set(i, j);
        }
      }

    //cerr<<"Global best score is "<<globalBestScore<<", pair ["<<globalBestPair[0]<<","<<globalBestPair[1]<<"]\n";

    // Construct idxMatched
    TwoValues<int> tmpPair;
    int curPeakIdx;
    vector<bool> idxMatchedBool(spec.size()); // Boolean vector used to mark matched peaks
    for (int i = 0; i < spec.size(); i++)
      idxMatchedBool[i] = false;
    if (globalBestPair[0] > 0)
    {
      idxMatchedBool[prefixPRMs[globalBestPair[0]][0] - addZero] = true;
      if (includeSymmetric && peakPairs[prefixPRMs[globalBestPair[0]][0]] > 0)
        idxMatchedBool[peakPairs[prefixPRMs[globalBestPair[0]][0]] - addZero] =
            true;
    }
    if (globalBestPair[1] > 0)
    {
      idxMatchedBool[suffixPRMs[globalBestPair[1]][0] - addZero] = true;
      if (includeSymmetric && peakPairs[suffixPRMs[globalBestPair[1]][0]] > 0)
        idxMatchedBool[peakPairs[suffixPRMs[globalBestPair[1]][0]] - addZero] =
            true;
    }
    while (globalBestPair[0] > 0 || globalBestPair[1] > 0)
    {
      tmpPair = prevPair[globalBestPair[0]][globalBestPair[1]];
      curPeakIdx = -1;
      if (tmpPair[0] != globalBestPair[0])
      {
        if (tmpPair[0] > 0)
          curPeakIdx = prefixPRMs[tmpPair[0]][0];
      }
      else if (tmpPair[1] > 0)
        curPeakIdx = suffixPRMs[tmpPair[1]][0];
      if (curPeakIdx > 0)
      {
        idxMatchedBool[curPeakIdx - addZero] = true;
        if (includeSymmetric && peakPairs[curPeakIdx] > 0)
          idxMatchedBool[peakPairs[curPeakIdx] - addZero] = true;
        //cerr << "Marked ["<<curPeakIdx-addZero<<","<<peakPairs[curPeakIdx]-addZero<<"]\n";
      }
      globalBestPair = tmpPair;
    }
    if (addZero == 0)
      idxMatchedBool[prefixPRMs[0][0]] = true; // If these were already in the spectrum
    if (addPM == 0)
      idxMatchedBool[suffixPRMs[0][0] - addZero] = true; //  then their scores are part of the match

    // Copy matches from idxMatchedBool to idxMatched
    idxMatched.resize(spec.size());
    k = 0;
/*    for (int i = 0; i < spec.size(); i++)
      if (idxMatchedBool[i])
        idxMatched[k++] = i;//*/
    for (int i = addZero; i < spec.size()-addPM; i++){
	  if ( idxMatchedBool[i] ){
		idxMatched[k++] = i-addZero;
	  }
	}
    idxMatched.resize(k);

    if (addZero + addPM > 0)
      spec = oldSpec; // Reverse the addition of peaks at masses zero and parentMass

    return globalBestScore;
  }

  static double** forwardAlignGF(Spectrum &spec,
                                 int scoreMax,
                                 float minDistance)
  {

    static double ProbOfMatchingRandomPeak = 0.05;

    int rowMax = spec.size();
    int colMax = scoreMax + 1;

    double **dpTable = (double**)malloc(sizeof(double*) * rowMax);
    for (int i = 0; i < rowMax; i++)
    {
      dpTable[i] = (double*)malloc(sizeof(double) * colMax);
    }

    for (int i = 0; i < rowMax; i++)
    {
      for (int j = 0; j < colMax; j++)
      {
        dpTable[i][j] = 0;
      }
    }
//	  dpTable[0][0] = 1;
    dpTable[0][(int)spec[0][1]] = 1;

    for (int peak = 1; peak < rowMax; peak++)
    {
      for (int score = 0; score < colMax; score++)
      {
        dpTable[peak][score] = dpTable[peak - 1][score]
            * (1 - ProbOfMatchingRandomPeak);
        int preScore = score - spec[peak][1];
        if (preScore < 0)
          continue;

        int prePeak = 0;
        for (int i = peak - 1; i > -1; i--)
        {
          if (spec[peak][0] - spec[i][0] >= minDistance)
          {
            prePeak = i;
            break;
          }
        }
        dpTable[peak][score] += ProbOfMatchingRandomPeak
            * pow(1 - ProbOfMatchingRandomPeak, peak - prePeak - 1)
            * dpTable[prePeak][preScore];
      }
    }
    /*	  for(int score=0; score<colMax; score++){
     dpTable[rowMax-1][score] = dpTable[rowMax-2][score];
     }//*/

    for (int peak = 0; peak < rowMax; peak++)
    {
      double pSum = 0;
      for (int score = colMax - 1; score > -1; score--)
      {
        pSum += dpTable[peak][score];
        dpTable[peak][score] = pSum;
      }
      double norm = 1. / pSum;
      for (int score = 0; score < colMax; score++)
      {
        dpTable[peak][score] *= norm;
      }
    }
    return dpTable;
  }

  static double** reverseAlignGF(Spectrum &spec,
                                 int scoreMax,
                                 float minDistance)
  {

    static double ProbOfMatchingRandomPeak = 0.05;

    int rowMax = spec.size();
    int colMax = scoreMax + 1;

    double **dpTable = (double**)malloc(sizeof(double*) * rowMax);
    for (int i = 0; i < rowMax; i++)
    {
      dpTable[i] = (double*)malloc(sizeof(double) * colMax);
    }

    for (int i = 0; i < rowMax; i++)
    {
      for (int j = 0; j < colMax; j++)
      {
        dpTable[i][j] = 0;
      }
    }
    //  dpTable[rowMax-1][0] = 1;
    dpTable[rowMax - 1][(int)spec[rowMax - 1][1]] = 1;

    for (int peak = rowMax - 2; peak > -1; peak--)
    {
      for (int score = 0; score < colMax; score++)
      {
        dpTable[peak][score] = dpTable[peak + 1][score]
            * (1 - ProbOfMatchingRandomPeak);
        int preScore = score - spec[peak][1];
        if (preScore < 0)
          continue;

        int prePeak = rowMax - 1;
        for (int i = peak + 1; i < rowMax; i++)
        {
          if (spec[i][0] - spec[peak][0] >= minDistance)
          {
            prePeak = i;
            break;
          }
        }
        dpTable[peak][score] += ProbOfMatchingRandomPeak
            * pow(1 - ProbOfMatchingRandomPeak, prePeak - peak - 1)
            * dpTable[prePeak][preScore];
      }
    }
    /*  for(int score=0; score<colMax; score++){
     dpTable[0][score] = dpTable[1][score];
     }//*/

    for (int peak = 0; peak < rowMax; peak++)
    {
      double pSum = 0;
      for (int score = colMax - 1; score > -1; score--)
      {
        pSum += dpTable[peak][score];
        dpTable[peak][score] = pSum;
      }
      double norm = 1. / pSum;
      for (int score = 0; score < colMax; score++)
      {
        dpTable[peak][score] *= norm;
      }
    }
    return dpTable;
  }

  void getALGFProbDensity(vector<float>& probDensity,
		                  Spectrum &spec,
                          int scoreMax,
                          float peakTol,
                          float minDistance)
  {
    double **dpTable = forwardAlignGF(spec, scoreMax, minDistance - 2 * peakTol);

    int size = scoreMax+1, index = spec.size()-1;
    probDensity.resize(size);
    for(int i=0; i<size; i++){
    	probDensity[i] = dpTable[index][i];
    }

    for (int i = 0; i < spec.size(); i++)
    	free(dpTable[i]);
    free(dpTable);
  }


  double*** getALGFProbDensity(Spectrum &spec,
                               int scoreMax,
                               float peakTol,
                               float minDistance)
  {
    double ***dpTable = (double***)malloc(sizeof(double**) * 2);
    dpTable[0] = forwardAlignGF(spec, scoreMax, minDistance - 2 * peakTol);
    dpTable[1] = reverseAlignGF(spec, scoreMax, minDistance - 2 * peakTol);
    return dpTable;
  }
  void freeALGFTable(double ***dpTable, Spectrum &spec)
  {
    for (int i = 0; i < spec.size(); i++)
    {
      free(dpTable[0][i]);
      free(dpTable[1][i]);
    }
    free(dpTable[0]);
    free(dpTable[1]);

    free(dpTable);
  }
  double getALGFPValue(double ***dpTable,
                       Spectrum &spec,
                       int fwScore,
                       int rvScore,
                       double mass)
  {
    int index = 0;
    double pvalue = 0;
    for (int i = spec.size() - 1; i > -1; i--)
    {
      if (spec[i][0] < mass)
      {
        index = i;
        break;
      }
    }
    pvalue = dpTable[0][index][fwScore];

    //reverse score
    double deltaMass = spec[spec.size() - 1][0] - mass;
    index = 0;
    for (int i = 0; i < spec.size(); i++)
    {
      if (spec[i][0] > deltaMass)
      {
        index = i;
        break;
      }
    }
    if (dpTable[1][index][rvScore] < pvalue)
      pvalue = dpTable[1][index][rvScore];  //*/
    return pvalue;
  }
  double getForwardALGFPValue(double ***dpTable,
                              Spectrum &spec,
                              int score,
                              double mass)
  {
    int index = 0;
    for (int i = spec.size() - 1; i > -1; i--)
    {
      if (spec[i][0] < mass)
      {
        index = i;
        break;
      }
    }
    return dpTable[0][index][score];
  }
  double getReverseALGFPValue(double ***dpTable,
                              Spectrum &spec,
                              int score,
                              double mass)
  {
    //reverse score
    double deltaMass = spec[spec.size() - 1][0] - mass;
    int index = 0;
    for (int i = 0; i < spec.size(); i++)
    {
      if (spec[i][0] > deltaMass)
      {
        index = i;
        break;
      }
    }
    return dpTable[1][index][score];
  }

  double getALGFPValueForInternalRegion(Spectrum &spec,
                                        int scoreMax,
                                        float peakTol,
                                        float startMass,
                                        float endMass,
                                        int threscore,
                                        float minDistance)
  {

    static double ProbOfMatchingRandomPeak = 0.05;
    float minJump = minDistance - peakTol;

    int rowMax = spec.size();
    for (int i = spec.size() - 1; i > -1; i--)
    {
      if (fabs(spec[i][0] - endMass) < peakTol || spec[i][0] < endMass)
      {
        rowMax = i + 1;
        break;
      }
    }

    int colMax = scoreMax + 1;

    double **dpTable = (double**)malloc(sizeof(double*) * rowMax);
    for (int i = 0; i < rowMax; i++)
    {
      dpTable[i] = (double*)malloc(sizeof(double) * colMax);
    }
    for (int i = 0; i < rowMax; i++)
    {
      for (int j = 0; j < colMax; j++)
      {
        dpTable[i][j] = 0;
      }
    }

    int startPoint = 0, baseScore = 0;
    for (int i = 0; i < spec.size(); i++)
    {
      if (fabs(spec[i][0] - startMass) < peakTol)
      {
        startPoint = i;
        baseScore = spec[i][1];
        break;
      }
      else if (spec[i][0] > startMass)
      {
        startPoint = i - 1;
        break;
      }
    }

    dpTable[startPoint][baseScore] = 1;

    for (int peak = startPoint + 1; peak < rowMax; peak++)
    {
      for (int score = 0; score < colMax; score++)
      {
        dpTable[peak][score] = dpTable[peak - 1][score]
            * (1 - ProbOfMatchingRandomPeak);
        int preScore = score - spec[peak][1];
        if (preScore < 0)
          continue;

        int prePeak = startPoint;
        for (int i = peak - 1; i >= startPoint; i--)
        {
          if (spec[peak][0] - spec[i][0] >= minJump)
          {
            prePeak = i;
            break;
          }
        }
        dpTable[peak][score] += ProbOfMatchingRandomPeak
            * pow(1 - ProbOfMatchingRandomPeak, peak - prePeak - 1)
            * dpTable[prePeak][preScore];
      }
    }

    double pSum = 0;
    for (int score = colMax - 1; score > -1; score--)
    {
      pSum += dpTable[rowMax - 1][score];
      dpTable[rowMax - 1][score] = pSum;
    }
    double norm = 1. / pSum;
    for (int score = 0; score < colMax; score++)
    {
      dpTable[rowMax - 1][score] *= norm;
    }
    double pvalue = dpTable[rowMax - 1][threscore];

    for (int i = 0; i < rowMax; i++)
    {
      free(dpTable[i]);
    }
    free(dpTable);
    return pvalue;
  }

  //
  //  findMatchingJumps - Finds all peaks to the left and to the right of peaks[peakIdx]
  //    whose peak masses correspond to a valid jump from peaks[peakIdx]. The indices
  //    of these peaks are returned in leftMatches and rightMatches.
  //
  void findMatchingJumps(int peakIdx,
                         vector<double> &peaks,
                         vector<float> &jumps,
                         float peakTol,
                         vector<int> &leftMatches,
                         vector<int> &rightMatches);

#ifdef DEBUG
  float SpectrumAlignment(Spectrum *spec1, Spectrum *spec2, float peakTol,
      Spectrum *matched1, Spectrum *matched2, int maxAAJump,
      float sameVertexPenalty, float ptmPenalty, bool forceSymmetry,
      bool addZPMmatches, ostream &debug)
  {
#else
  float SpectrumAlignment(Spectrum *spec1,
                          Spectrum *spec2,
                          float peakTol,
                          Spectrum *matched1,
                          Spectrum *matched2,
                          int maxAAJump,
                          float sameVertexPenalty,
                          float ptmPenalty,
                          bool forceSymmetry,
                          bool addZPMmatches)
  {
#endif

    const float MIN_AA_MASS = 57.0214637230;
    //	const float MIN_AA_MASS = 1;

    if (spec1->size() == 0 || spec2->size() == 0)
    {
      if (matched1)
        matched1->resize(0);
      if (matched2)
        matched1->resize(0);
      return 0;
    }

    AAJumps jumps(0);
    float jumpsSupremum = 0;
    if (maxAAJump > 0)
    {
      jumps.getjumps(maxAAJump);
      jumpsSupremum = jumps.masses[jumps.masses.size() - 1] + 2 * peakTol
          + .00001;
    }

    //#ifdef DEBUG
    //	ofstream debug("dekel_align_debug.txt");
    //#endif

    if (spec1->parentMass > spec2->parentMass)
    {
      Spectrum *tmp = spec1;
      spec1 = spec2;
      spec2 = tmp;
      tmp = matched1;
      matched1 = matched2;
      matched2 = tmp;
    }

    vector<double> peaks;
    peaks.reserve(spec1->size() + 1);
    vector<double> peaks2(spec2->size());
    for (int i = 0; i < spec2->size(); i++)
      peaks2[i] = (*spec2)[i][0];
    vector<int> common;
    common.reserve(spec1->size() + 1);
    vector<int> common2;
    common2.reserve(spec1->size() + 1);
    vector<double> common_scores;
    common_scores.reserve(spec1->size() + 1);
    vector<double> common2_scores;
    common2_scores.reserve(spec1->size() + 1);
    vector<int> prev;
    prev.reserve(spec1->size());
    vector<int> next;
    next.reserve(spec1->size());
    vector<int> prev2;
    prev2.reserve(spec1->size());
    vector<int> next2;
    next2.reserve(spec1->size());
    vector<vector<int> > left_neighbors;
    left_neighbors.reserve(spec1->size());
    vector<vector<int> > right_neighbors;
    right_neighbors.reserve(spec1->size());

    Spectrum spec1sym;   // Symmetric version of spectrum1

    if (forceSymmetry)
    {
      spec1sym = *spec1;
      spec1sym.makeSymmetric(0, peakTol);
      //MakeSymmetric(spec1, spec1->parentMass, peakTol, &spec1sym);
    }
    else
    {
      spec1sym.resize(spec1->size());
      for (unsigned int i = 0; i < spec1sym.size(); i++)
        spec1sym[i] = (*spec1)[i];
    }

    // Fill in peaks, common, common2, common_scores, common2_scores
    int idx1,       // Index in spec1sym
        idx2 = 0,     // Index in spectrum 2
        idxPeaks = 0; // Index in peaks
    float pmDelta = spec2->parentMass - spec1->parentMass;
    for (idx1 = 0; idx1 < (int)spec1sym.size(); idx1++)
    {
      peaks.resize(idxPeaks + 1);
      common.resize(idxPeaks + 1);
      common2.resize(idxPeaks + 1);
      common_scores.resize(idxPeaks + 1);
      common2_scores.resize(idxPeaks + 1);

      // Fill in common
      if (idx2 >= spec2->size())
        idx2 = (int)spec2->size() - 1;
      for (; idx2 >= 0 && (*spec2)[idx2][0] >= spec1sym[idx1][0] - peakTol;
          idx2--)
        ;  // Move below potential matches
      for (idx2 = max(0, idx2);
          idx2 < spec2->size()
              && (*spec2)[idx2][0] < spec1sym[idx1][0] - peakTol; idx2++)
        ;  // Find first potential match
      common_scores[idxPeaks] = 0;
      common[idxPeaks] = DEKEL::None;
      for (;
          idx2 < (int)spec2->size()
              && (*spec2)[idx2][0] <= spec1sym[idx1][0] + peakTol; idx2++)
        if ((*spec2)[idx2][1] > common_scores[idxPeaks])
        {
          common_scores[idxPeaks] = (*spec2)[idx2][1];
          common[idxPeaks] = idx2;
        }
      if (common_scores[idxPeaks] > 0)
        common_scores[idxPeaks] += spec1sym[idx1][1];

      //Fill in common2
      if (idx2 >= (int)spec2->size())
        idx2 = (int)spec2->size() - 1;
      for (;
          idx2 > 0 && (*spec2)[idx2][0] >= spec1sym[idx1][0] + pmDelta - peakTol;
          idx2--)
        ;
      for (;
          idx2 < (int)spec2->size()
              && (*spec2)[idx2][0] < spec1sym[idx1][0] + pmDelta - peakTol;
          idx2++)
        ;
      common2_scores[idxPeaks] = 0;
      common2[idxPeaks] = DEKEL::None;
      for (;
          idx2 < (int)spec2->size()
              && (*spec2)[idx2][0] <= spec1sym[idx1][0] + pmDelta + peakTol;
          idx2++)
        if ((*spec2)[idx2][1] > common2_scores[idxPeaks])
        {
          common2_scores[idxPeaks] = (*spec2)[idx2][1];
          common2[idxPeaks] = idx2;
        }
      if (common2_scores[idxPeaks] > 0)
        common2_scores[idxPeaks] += spec1sym[idx1][1];

#ifdef DEBUG
      debug<<"spec1["<<idx1<<"], common_scores = "<<common_scores[idxPeaks]<<", common2_scores = "<<common2_scores[idxPeaks]<<endl;
#endif
      if (common_scores[idxPeaks] > 0 || common2_scores[idxPeaks] > 0)
      {
        peaks[idxPeaks++] = spec1sym[idx1][0];
#ifdef DEBUG
        debug<<"--- added peaks["<<idxPeaks-1<<"] = "<<peaks[idxPeaks-1]<<endl;
#endif
      }
    }
    peaks.resize(idxPeaks);

    // Make peaks symmetric and change the common/common2 structures accordingly
    Spectrum peaksTmp;
    peaksTmp.copyNP(spec1sym);
    for (int i = 0; i < peaks.size(); i++)
      peaksTmp.insertPeak(peaks[i], common_scores[i], 0);
    vector<int> indices;
    if (forceSymmetry)
    {
      peaksTmp.makeSymmetric(0, peakTol, &indices);
      //MakeSymmetric(&peaksTmp, spec1->parentMass, peakTol, NULL, &indices);
    }
    else
    {
      indices.resize(peaksTmp.size());
      for (unsigned int i = 0; i < indices.size(); i++)
        indices[i] = i;
    }

#ifdef DEBUG
    debug << "After makeSymmetric: \n";
    for(int i=0; i<peaksTmp.size(); i++)
    { debug<<"["<<peaksTmp[i][0]<<","<<peaksTmp[i][1]<<", idx = "<<indices[i]<<"]\n";}
#endif

    peaks.resize(peaksTmp.size());
    common.resize(peaksTmp.size());
    common_scores.resize(peaksTmp.size());
    common2.resize(peaksTmp.size());
    common2_scores.resize(peaksTmp.size());
    idx1 = spec1->size() - 1;    // Keeps track of the nearby peaks in spec1
    for (int i = peaksTmp.size() - 1; i >= 0; i--)
    {  // iterate backwards to avoid overwriting entries in common*
      peaks[i] = peaksTmp[i][0];
      if (indices[i] >= 0)
      {
        common[i] = common[indices[i]];
        common_scores[i] = common_scores[indices[i]]; // Note that indices[i]<=i
        common2[i] = common2[indices[i]];
        common2_scores[i] = common2_scores[indices[i]];
      }
      else
      {
        // Look in spec1 for closest peak with highest score
        int j;
        float bestScore = 0;
        for (j = idx1;
            j < spec1->size() and peaks[i] + peakTol > (*spec1)[j][0]; j++)
          ; // Find peaks after peaks[i]
        if (j == spec1->size())
          j--;
        for (; j >= 0 && peaks[i] + peakTol < (*spec1)[j][0]; j--)
          ; // Find first peak within tolerance of peaks[i]
        idx1 = j;
        for (; j >= 0 && abs(peaks[i] - (*spec1)[j][0]) <= peakTol; j--)
          bestScore = max(bestScore, (*spec1)[j][1]); // Get score of best peak within tolerance
        common[i] = DEKEL::None;
        common_scores[i] = bestScore;
        common2[i] = DEKEL::None;
        common2_scores[i] = bestScore;
      }
    }

#ifdef DEBUG
    debug<<"i: [peaks[i]] [common[i],common_scores[i]] [common2[i],common2_scores[i]]\n";
    for(int i=0; i<peaks.size(); i++)
    { debug<<i<<":\t["<<peaksTmp[i][0]<<"]\t["<<common[i]<<","<<common_scores[i]<<"]\t["<<common2[i]<<","<<common2_scores[i]<<"]\n";}
#endif

    // Fill in prev, prev2, next, next2
    prev.resize(peaks.size());
    prev2.resize(peaks.size());
    next.resize(peaks.size());
    next2.resize(peaks.size());
    left_neighbors.resize(peaks.size());
    right_neighbors.resize(peaks.size());
    vector<vector<int> > left_jumps(peaks.size()), left_jumps2(peaks.size()),
        right_jumps(peaks.size()), right_jumps2(peaks.size());
    int idxPrev = 0, // Keeps track of the righmost peak that trails the current by >=57 Da
        idxPrevDelta = 0, // Keeps track of the righmost peak that trails the current by >=delta Da
        idxNext = 0, // Keeps track of the leftmost peak ahead of current by >=57 Da
        idxNextDelta = 0; // Keeps track of the leftmost peak ahead of current by >=delta Da
    // Sep.8,05	float delta = max(MIN_AA_MASS,pmDelta)-peakTol;
    float delta = min(max(MIN_AA_MASS - peakTol, pmDelta + peakTol),
                      2 * MIN_AA_MASS - peakTol);
    float delta2 = jumpsSupremum - peakTol;
    for (idxPeaks = 0; idxPeaks < (int)peaks.size(); idxPeaks++)
    {
      if (idxPrev < 0)
        idxPrev = 0;
      for (;
          idxPrev < peaks.size()
              && peaks[idxPrev]
                  < peaks[idxPeaks] - max(MIN_AA_MASS - peakTol, delta2);
          idxPrev++)
        ;
      idxPrev--;
      if (idxPrev >= 0)
        prev[idxPeaks] = idxPrev;
      else
        prev[idxPeaks] = DEKEL::None;

      for (;
          idxPrevDelta < peaks.size()
              && peaks[idxPrevDelta] < peaks[idxPeaks] - max(delta, delta2);
          idxPrevDelta++)
        ;
      idxPrevDelta--;
      if (idxPrevDelta >= 0)
        prev2[idxPeaks] = idxPrevDelta;
      else
      {
        prev2[idxPeaks] = DEKEL::None;
        idxPrevDelta = 0;
      }

      // Similar for next/next2
      for (idxNext = idxPeaks;
          idxNext < peaks.size()
              && peaks[idxNext]
                  <= peaks[idxPeaks] + max(MIN_AA_MASS - peakTol, delta2);
          idxNext++)
        ;
      if (idxNext == peaks.size())
        next[idxPeaks] = DEKEL::None;
      else
        next[idxPeaks] = idxNext;

      for (idxNextDelta = idxNext;
          idxNextDelta < peaks.size()
              && peaks[idxNextDelta] <= peaks[idxPeaks] + max(delta, delta2);
          idxNextDelta++)
        ;
      if (idxNextDelta == peaks.size())
        next2[idxPeaks] = DEKEL::None;
      else
        next2[idxPeaks] = idxNextDelta;

      if (maxAAJump <= 0)
      {
        if (pmDelta <= MIN_AA_MASS - 2 * peakTol)
        {
          left_neighbors[idxPeaks].resize(0);
          right_neighbors[idxPeaks].resize(0);
        }
        else
        {
          unsigned int i;  // First valid peak in the target interval
          if (prev[idxPeaks] != DEKEL::None)
          {
            if (prev2[idxPeaks] != DEKEL::None)
              left_neighbors[idxPeaks].resize(prev[idxPeaks] - prev2[idxPeaks]);
            else
              left_neighbors[idxPeaks].resize(prev[idxPeaks] + 1);
            for (i = 0;
                prev[idxPeaks] >= i
                    && peaks[prev[idxPeaks] - i] >= peaks[idxPeaks] - delta;
                i++)
              left_neighbors[idxPeaks][i] = prev[idxPeaks] - i;
            left_neighbors[idxPeaks].resize(i);
          }

          // Similar for right_neighbors
          if (next[idxPeaks] != DEKEL::None)
          {
            if (next2[idxPeaks] != DEKEL::None)
              right_neighbors[idxPeaks].resize(next2[idxPeaks]
                  - next[idxPeaks]);
            else
              right_neighbors[idxPeaks].resize(peaks.size() - next[idxPeaks]);
            for (i = 0;
                next[idxPeaks] + i < peaks.size()
                    && peaks[next[idxPeaks] + i] <= peaks[idxPeaks] + delta;
                i++)
              right_neighbors[idxPeaks][i] = next[idxPeaks] + i;
            right_neighbors[idxPeaks].resize(i);
          }
        }
        left_jumps[idxPeaks].resize(0);
        left_jumps2[idxPeaks].resize(0);
        right_jumps[idxPeaks].resize(0);
        right_jumps2[idxPeaks].resize(0);
      }
      else
      {
        findMatchingJumps(idxPeaks,
                          peaks,
                          jumps.masses,
                          peakTol,
                          left_jumps[idxPeaks],
                          right_jumps[idxPeaks]);
        // left_neighbors[idxPeaks], right_neighbors[idxPeaks]
        unsigned int neighCount = 0, jumpsCount = 0;
        left_neighbors[idxPeaks].resize(left_jumps[idxPeaks].size());
        left_jumps2[idxPeaks].resize(left_jumps[idxPeaks].size());
        for (unsigned int i = 0; i < left_jumps[idxPeaks].size(); i++)
        {
          if (peaks[idxPeaks] - peaks[left_jumps[idxPeaks][i]] <= delta)
            left_neighbors[idxPeaks][neighCount++] = left_jumps[idxPeaks][i];
          if (peaks[idxPeaks] - peaks[left_jumps[idxPeaks][i]] >= delta)
            left_jumps2[idxPeaks][jumpsCount++] = left_jumps[idxPeaks][i];
        }
        left_neighbors[idxPeaks].resize(neighCount);
        left_jumps2[idxPeaks].resize(jumpsCount);

        neighCount = 0, jumpsCount = 0;
        right_neighbors[idxPeaks].resize(right_jumps[idxPeaks].size());
        right_jumps2[idxPeaks].resize(right_jumps[idxPeaks].size());
        for (unsigned int i = 0; i < right_jumps[idxPeaks].size(); i++)
        {
          if (peaks[right_jumps[idxPeaks][i]] - peaks[idxPeaks] <= delta)
            right_neighbors[idxPeaks][neighCount++] = right_jumps[idxPeaks][i];
          if (peaks[right_jumps[idxPeaks][i]] - peaks[idxPeaks] >= delta)
            right_jumps2[idxPeaks][jumpsCount++] = right_jumps[idxPeaks][i];
        }
        right_neighbors[idxPeaks].resize(neighCount);
        right_jumps2[idxPeaks].resize(jumpsCount);
      }
      /*
       #ifdef DEBUG
       debug <<"peaks["<<idxPeaks<<"]= "<<peaks[idxPeaks]<<", prev = "<<prev[idxPeaks]<<", prev2 = "<<prev2[idxPeaks]<<", next = "<<next[idxPeaks]<<", next2 = "<<next2[idxPeaks]<<endl;
       debug <<"--- left_jumps("<<left_jumps[idxPeaks].size()<<"): "; for(int i=0; i<left_jumps[idxPeaks].size(); i++) debug << left_jumps[idxPeaks][i] <<" "; debug << endl;
       debug <<"--- right_jumps("<<right_jumps[idxPeaks].size()<<"): "; for(int i=0; i<right_jumps[idxPeaks].size(); i++) debug << right_jumps[idxPeaks][i] <<" "; debug << endl;
       debug <<"--- left_jumps2("<<left_jumps2[idxPeaks].size()<<"): "; for(int i=0; i<left_jumps2[idxPeaks].size(); i++) debug << left_jumps2[idxPeaks][i] <<" "; debug << endl;
       debug <<"--- right_jumps2("<<right_jumps2[idxPeaks].size()<<"): "; for(int i=0; i<right_jumps2[idxPeaks].size(); i++) debug << right_jumps2[idxPeaks][i] <<" "; debug << endl;
       debug <<"--- left_neighbors("<<left_neighbors[idxPeaks].size()<<"): "; for(int i=0; i<left_neighbors[idxPeaks].size(); i++) debug << left_neighbors[idxPeaks][i] <<" "; debug << endl;
       debug <<"--- right_neighbors("<<right_neighbors[idxPeaks].size()<<"): "; for(int i=0; i<right_neighbors[idxPeaks].size(); i++) debug << right_neighbors[idxPeaks][i] <<" "; debug << endl;
       #endif
       */
      /*
       cerr <<"peaks["<<idxPeaks<<"]= "<<peaks[idxPeaks]<<", prev = "<<prev[idxPeaks]<<", prev2 = "<<prev2[idxPeaks]<<", next = "<<next[idxPeaks]<<", next2 = "<<next2[idxPeaks]<<endl;
       cerr <<"--- left_jumps("<<left_jumps[idxPeaks].size()<<"): "; for(int i=0; i<left_jumps[idxPeaks].size(); i++) cerr << left_jumps[idxPeaks][i] <<" "; cerr << endl;
       cerr <<"--- right_jumps("<<right_jumps[idxPeaks].size()<<"): "; for(int i=0; i<right_jumps[idxPeaks].size(); i++) cerr << right_jumps[idxPeaks][i] <<" "; cerr << endl;
       cerr <<"--- left_jumps2("<<left_jumps2[idxPeaks].size()<<"): "; for(int i=0; i<left_jumps2[idxPeaks].size(); i++) cerr << left_jumps2[idxPeaks][i] <<" "; cerr << endl;
       cerr <<"--- right_jumps2("<<right_jumps2[idxPeaks].size()<<"): "; for(int i=0; i<right_jumps2[idxPeaks].size(); i++) cerr << right_jumps2[idxPeaks][i] <<" "; cerr << endl;
       cerr <<"--- left_neighbors("<<left_neighbors[idxPeaks].size()<<"): "; for(int i=0; i<left_neighbors[idxPeaks].size(); i++) cerr << left_neighbors[idxPeaks][i] <<" "; cerr << endl;
       cerr <<"--- right_neighbors("<<right_neighbors[idxPeaks].size()<<"): "; for(int i=0; i<right_neighbors[idxPeaks].size(); i++) cerr << right_neighbors[idxPeaks][i] <<" "; cerr << endl;
       */
    }

    pair<double, DEKEL::pair2> res;

    if (forceSymmetry)
      res = DEKEL::align(spec1->parentMass - AAJumps::massMH,
                         spec2->parentMass - AAJumps::massMH,
                         peaks,
                         peaks2,
                         common,
                         common2,
                         common_scores,
                         common2_scores,
                         prev,
                         next,
                         prev2,
                         next2,
                         left_jumps,
                         right_jumps,
                         left_jumps2,
                         right_jumps2,
                         left_neighbors,
                         right_neighbors,
                         sameVertexPenalty,
                         ptmPenalty);
    else
      res = DEKEL::align_simple(spec1->parentMass - AAJumps::massMH,
                                spec2->parentMass - AAJumps::massMH,
                                peaks,
                                peaks2,
                                common,
                                common2,
                                common_scores,
                                common2_scores,
                                prev,
                                next,
                                left_jumps,
                                right_jumps,
                                ptmPenalty,
                                0);

    //cerr << "Alignment score: " << res.first << endl; cerr.flush();

    int idxMatch = 0;
    float modPos;  // Mass value of where the mod was placed
    matched1->copyNP(*spec1);
    matched2->copyNP(*spec2);
    matched1->resize(res.second.first.size() + res.second.second.size());
    matched2->resize(res.second.first.size() + res.second.second.size());
#ifdef DEBUG
    debug << "Peaks matched before the modification:\n";
#endif

    for (unsigned int i = 0; i < res.second.first.size(); i++)
    {
      idxPeaks = res.second.first[i];
      (*matched1)[idxMatch][0] = peaks[idxPeaks];
      if (common[idxPeaks] >= 0)
        (*matched2)[idxMatch].set(peaks2[common[idxPeaks]],
                                  (*spec2)[common[idxPeaks]][1]);
      else
        (*matched2)[idxMatch].set((*matched1)[idxMatch][0], 0);
      (*matched1)[idxMatch][1] = common_scores[idxPeaks]
          - (*matched2)[idxMatch][1];
#ifdef DEBUG
      debug << " -1- "<<res.second.first[i]<<"\t["<<(*matched1)[idxMatch][0]<<","<<(*matched1)[idxMatch][1]<<"] \t["<<(*matched2)[idxMatch][0]<<","<<(*matched2)[idxMatch][1]<<"]\n";
#endif
      idxMatch++;
    }

#ifdef DEBUG
    debug << "Peaks matched after the modification:\n";
#endif
    for (unsigned int i = 0; i < res.second.second.size(); i++)
    {
      idxPeaks = res.second.second[i];
      (*matched1)[idxMatch][0] = peaks[idxPeaks];
      if (common2[idxPeaks] >= 0)
        (*matched2)[idxMatch].set(peaks2[common2[idxPeaks]],
                                  (*spec2)[common2[idxPeaks]][1]);
      else
        (*matched2)[idxMatch].set((*matched1)[idxMatch][0] + pmDelta, 0);
      (*matched1)[idxMatch][1] = common2_scores[idxPeaks]
          - (*matched2)[idxMatch][1];
#ifdef DEBUG
      debug << " -2- "<<res.second.second[i]<<"\t["<<(*matched1)[idxMatch][0]<<","<<(*matched1)[idxMatch][1]<<"] \t["<<(*matched2)[idxMatch][0]<<","<<(*matched2)[idxMatch][1]<<"]\n";
#endif
      idxMatch++;
    }

    if (addZPMmatches && matched1->size() > 0 && matched2->size() > 0)
    {
      if (spec1->size() > 1 && spec2->size() > 1
          && ((*matched1)[0] == (*spec1)[0] || (*matched1)[0] == (*spec1)[1])
          && ((*matched2)[0] == (*spec2)[0] || (*matched2)[0] == (*spec2)[1]))
      {
        matched1->resize(matched1->size() + 1);
        for (unsigned int i = matched1->size() - 1; i > 0; i--)
          (*matched1)[i] = (*matched1)[i - 1];
        (*matched1)[0] = (*spec1)[0];
        (*matched1)[1] = (*spec1)[1];
        matched2->resize(matched2->size() + 1);
        for (unsigned int i = matched2->size() - 1; i > 0; i--)
          (*matched2)[i] = (*matched2)[i - 1];
        (*matched2)[0] = (*spec2)[0];
        (*matched2)[1] = (*spec2)[1];
#ifdef DEBUG
        debug << " -x- Extra match at the start: ["<<(*matched1)[0][0]<<","<<(*matched2)[0][0]<<"],["<<(*matched1)[1][0]<<","<<(*matched2)[1][0]<<"]\n";
#endif
      }
      if (spec1->size() > 1 && spec2->size() > 1
          && ((*matched1)[matched1->size() - 1] == (*spec1)[spec1->size() - 1]
              || (*matched1)[matched1->size() - 1]
                  == (*spec1)[spec1->size() - 2])
          && ((*matched2)[matched2->size() - 1] == (*spec2)[spec2->size() - 1]
              | (*matched2)[matched2->size() - 1] == (*spec2)[spec2->size() - 2]))
      {
        matched1->resize(matched1->size() + 1);
        (*matched1)[matched1->size() - 1] = (*spec1)[spec1->size() - 1];
        (*matched1)[matched1->size() - 2] = (*spec1)[spec1->size() - 2];
        matched2->resize(matched2->size() + 1);
        (*matched2)[matched2->size() - 1] = (*spec2)[spec2->size() - 1];
        (*matched2)[matched2->size() - 2] = (*spec2)[spec2->size() - 2];
#ifdef DEBUG
        debug << " -x- Extra match at the end: ["<<(*matched1)[matched1->size()-2][0]<<","<<(*matched2)[matched2->size()-2][0]<<"],["<<(*matched1)[matched1->size()-1][0]<<","<<(*matched2)[matched2->size()-1][0]<<"]\n";
#endif
      }
    }

    if (res.second.first.size() == 0)
      modPos = 0;  // No peaks in idx1 => mod at the start
    else
    {
      if (res.second.second.size() == 0)
        modPos = spec1->parentMass;  // No peaks in idx2 => mod at the end
      else
      {
        modPos = peaks[res.second.second[0]]; // Otherwise mod was placed at the first mass of spec1 in idx2
      }
    }
#ifdef DEBUG
    debug << "modPos = " << modPos<<", ptmPenalty = "<<ptmPenalty<<", index sizes = ["<<res.second.first.size()<<", "<<res.second.second.size()<<"]\n";
#endif
    return modPos;
  }

  void findMatchingJumps(int peakIdx,
                         vector<double> &peaks,
                         vector<float> &jumps,
                         float peakTol,
                         vector<int> &leftMatches,
                         vector<int> &rightMatches)
  {
    int peaksPivot, jPivot, matchCount;
    float maxJump = jumps[jumps.size() - 1] + peakTol;
    float peakDiff;

    if (peakIdx > 0)
    {
      leftMatches.resize(peakIdx);
      matchCount = 0;
      for (peaksPivot = peakIdx - 1, jPivot = 0;
          peaksPivot >= 0 && peaks[peakIdx] - peaks[peaksPivot] <= maxJump;
          peaksPivot--)
      {
        peakDiff = (float)peaks[peakIdx] - peaks[peaksPivot];
        // Find a jump with a mass lower than the peak difference
        for (jPivot = min(jPivot, (int)jumps.size() - 1);
            jPivot >= 0 && jumps[jPivot] >= peakDiff - peakTol; jPivot--)
          ;
        jPivot = max(0, jPivot);
        // Find first jump with a mass not smaller than the peak difference - peakTol
        for (; jPivot < jumps.size() && jumps[jPivot] < peakDiff - peakTol;
            jPivot++)
          ;
        if (jPivot < jumps.size()
            && abs(jumps[jPivot] - peakDiff) <= peakTol + .00001)
          leftMatches[matchCount++] = peaksPivot;
      }
      leftMatches.resize(matchCount);
    }

    rightMatches.resize(peaks.size() - peakIdx);
    matchCount = 0;
    for (peaksPivot = peakIdx + 1, jPivot = 0;
        peaksPivot < peaks.size()
            && peaks[peaksPivot] - peaks[peakIdx] <= maxJump; peaksPivot++)
    {
      peakDiff = (float)peaks[peaksPivot] - peaks[peakIdx];
      // Find a jump with a mass lower than the peak difference
      for (jPivot = min(jPivot, (int)jumps.size() - 1);
          jPivot >= 0 && jumps[jPivot] >= peakDiff - peakTol; jPivot--)
        ;
      jPivot = max(0, jPivot);
      // Find first jump with a mass not smaller than the peak difference - peakTol
      for (; jPivot < jumps.size() && jumps[jPivot] < peakDiff - peakTol;
          jPivot++)
        ;
      if (jPivot < jumps.size()
          && abs(jumps[jPivot] - peakDiff) <= peakTol + .00001)
        rightMatches[matchCount++] = peaksPivot;
    }
    rightMatches.resize(matchCount);
  }

  bool GPCAux::operator<(GPCAux & other)
  {
    return score < other.score;
  }
  bool GPCAux_cmp(const GPCAux &a, const GPCAux &b)
  {
    return a.score < b.score;
  }

} // namespace specnets

namespace DEKEL
{
  const double infinity = 1e10;

  const double mass_tolerance = 0.5;

  using std::vector;
  using std::pair;
  using std::max;
  using std::min;
  using std::max_element;

  enum celltype
  {
    cell_prefix1,
    cell_suffix1,
    cell_prefix2,
    cell_suffix2,
    cell_prefix1_L,
    cell_suffix1_R,
    cell_prefix2_L,
    cell_suffix2_R,
    cell_D1,
    cell_D2,
    cell_D3,
    cell_M1_L,
    cell_M1_R,
    cell_M2_L,
    cell_M2_R,
    cell_M3_L,
    cell_M3_R,
    INVALID
  };

  typedef vector<vector<double> > vector2;
  typedef vector<vector<vector<double> > > vector3;

  //////////////////////////////////////////////////////////////////////////////
  pair<double, pair2> align_simple(double parent_mass1,
                                   double parent_mass2,
                                   vector<double> peaks,
                                   vector<double> peaks2,
                                   vector<int> & common,
                                   vector<int> & common2,
                                   vector<double> & common_scores,
                                   vector<double> & common2_scores,
                                   vector<int> & prev,
                                   vector<int> & next,
                                   vector<vector<int> > & left_jumps,
                                   vector<vector<int> > & right_jumps,
                                   double ptm_penalty,
                                   double ambiguous_penalty)
  // parent_mass1 : parent mass of first spectrum
  // parent mass2 : parent mass of 2nd spectrum.
  // Note: parent_mass2 can be either bigger or smaller than parent_mass1
  //
  // peaks        : All the peaks m in the 1st spectrum such that
  //                either the peak m or the peak m+(parent_mass2-parent_mass1)
  // peaks2       : All the peaks of the 2st spectrum
  //
  // common       : common[i] is the index j of the peak in spectrum2 whose mass
  //                is equal to peaks[i].
  //                If no such peak exists, common[i] = None
  // common2      : common[i] is the index j of the peak in spectrum2 whose mass
  //                is equal to peaks[i]+(parent_mass2-parent_mass1).
  // common_scores: common_scores[i] is the score of the vertex
  //                (peaks[i],peaks[i]) in the alignment graph.
  //                If the vertex doesn't exist, then common_scores[i]=-infinity
  // common2_scores: common_scores2[i] is the score of the vertex
  //                (peaks[i],peaks[i]+(parent_mass2-parent_mass1).
  //
  // Let  T = the mass tolerance (e.g. 0.5)
  //
  // Let J be a set of masses (possibly empty), and M be a mass.
  // We say that a mass jump m is valid if either |m-m'| < T for some m' from J,
  // or m >= M.
  // Assume that if J is not empty, then min(J) >= 57 and max(J) < M
  //
  // Define
  //   delta2 = M-T
  //
  // prev         : prev[i] is the maximum index j such that
  //                peaks[j] < peaks[i]-delta2
  //                If no such index exists, prev[i] = None
  // next         : next[i] is the minimum index j such that
  //                peaks[j] > peaks[i]+delta2
  //
  // left_jumps   : left_jumps[i] is a list of all the indices j such that
  //   (1) peaks[i]-delta2 <= peaks[j] < peaks[i]-(57-T)
  //   (2) peaks[i]-peaks[j] is a valid mass jump
  //   If J is empty, then left_jumps[i] is an empty list
  // right_jumps   : right_jumps[i] is a list of all the indices j such that
  //   (1) peaks[i]+57-T < peaks[j] <= peaks[i]+delta2
  //   (2) peaks[j]-peaks[i] is a valid mass jump
  //   If J is empty, then right_jumps[i] is an empty list
  //
  // ptm_penalty:   penalty for having a PTM (negative number)
  // ambiguous_penalty: penalty if the first peak i before (or after) the
  //                    modification satisfies
  //                        common[i] != None and common2[i] != None
  {
    int n = common.size();

    vector<double> prefix(n, -infinity);
    vector<double> suffix(n, -infinity);
    vector<double> prefix_L(n, -infinity);
    vector<double> suffix_R(n, -infinity);

    // compute prefix
    for (int i = 0; i < n; ++i)
    {
      if (common[i] != None)
      {
        double prev_score = 0.0;
        for (vector<int>::const_iterator it = left_jumps[i].begin();
            it != left_jumps[i].end(); ++it)
        {
          prev_score = max(prev_score, prefix[*it]);
        }
        int i2 = prev[i];
        if (i2 != None)
        {
          prev_score = max(prev_score, prefix_L[i2]);
        }
        prefix[i] = common_scores[i] + prev_score;
      }

      if (i != 0)
        prefix_L[i] = max(prefix_L[i - 1], prefix[i]);
      else
        prefix_L[i] = prefix[i];
    }

    // compute suffix
    for (int j = n - 1; j >= 0; --j)
    {
      if (common2[j] != None)
      {
        double next_score = 0.0;
        for (vector<int>::const_iterator it = right_jumps[j].begin();
            it != right_jumps[j].end(); ++it)
        {
          next_score = max(next_score, suffix[*it]);
        }
        int j2 = next[j];
        if (j2 != None)
        {
          next_score = max(next_score, suffix_R[j2]);
        }
        suffix[j] = common2_scores[j] + next_score;
      }
      if (j != n - 1)
        suffix_R[j] = max(suffix_R[j + 1], suffix[j]);
      else
        suffix_R[j] = suffix[j];
    }

    double best_score = 0.0;
    int best_i = None;
    int best_j = None;
    celltype best_t = INVALID;
    for (int i = 0; i < n; ++i)
    {
      double penalty = ptm_penalty;
      if (common2[i] != None)
        penalty += ambiguous_penalty;
      for (vector<int>::const_iterator it = right_jumps[i].begin();
          it != right_jumps[i].end(); ++it)
      {
        int j = *it;
        double penalty2 = penalty;
        if (common[j] != None)
          penalty2 += ambiguous_penalty;
        if (best_score < prefix[i] + suffix[j] + penalty2)
        {
          best_i = i;
          best_j = j;
          best_t = cell_suffix2;
          best_score = prefix[i] + suffix[j] + penalty2;
        }
      }
      int j = next[i];
      if (j != None)
      {
        if (best_score < prefix[i] + suffix_R[j] + penalty)
        {
          best_i = i;
          best_j = j;
          best_t = cell_suffix2_R;
          best_score = prefix[i] + suffix_R[j] + penalty;
        }
      }
    }

    for (int i = 0; i < n; ++i)
    {
      double penalty = 0;
      if (common2[i] != None)
        penalty += ambiguous_penalty;
      if (best_score < prefix[i] + penalty)
      {
        best_i = i;
        best_j = None;
        best_t = cell_suffix2;
        best_score = prefix[i] + penalty;
      }
    }
    for (int j = 0; j < n; ++j)
    {
      double penalty = 0;
      if (common[j] != None)
        penalty += ambiguous_penalty;
      if (best_score < suffix[j] + penalty)
      {
        best_i = None;
        best_j = j;
        best_t = cell_suffix2;
        best_score = suffix[j] + penalty;
      }
    }

    vector<int> path1;
    vector<int> path2;

    if (best_t == INVALID)
      return pair<double, pair2>(best_score / parent_mass1, pair2(path1, path2));

    int i = best_i;
    while (i != None)
    {
      path1.insert(path1.begin(), i);
      double prev_score = 0.0;
      int index = 0;
      int next_i = 0;
      for (vector<int>::const_iterator it = left_jumps[i].begin();
          it != left_jumps[i].end(); ++it)
      {
        int i2 = *it;
        if (prefix[i2] > prev_score)
        {
          prev_score = prefix[i2];
          index = 1;
          next_i = i2;
        }
      }
      int i2 = prev[i];
      if (i2 != None && prefix_L[i2] > prev_score)
      {
        prev_score = prefix_L[i2];
        index = 2;
        next_i = i2;
      }

      i = next_i;
      if (index == 0)
      {
        break;
      }
      else if (index == 1)
      {
        // do nothing
      }
      else
      {
        while (prefix_L[i] != prefix[i])
          i -= 1;
      }
    }

    int j = best_j;
    if (j != None && best_t == cell_suffix2_R)
    {
      while (suffix_R[j] != suffix[j])
      {
        j += 1;
      }
    }

    while (j != None)
    {
      path2.push_back(j);
      double next_score = 0.0;
      int index = 0;
      int next_j = 0;
      for (vector<int>::const_iterator it = right_jumps[j].begin();
          it != right_jumps[j].end(); ++it)
      {
        int j2 = *it;
        if (suffix[j2] > next_score)
        {
          next_score = suffix[j2];
          index = 1;
          next_j = j2;
        }
      }
      int j2 = next[j];
      if (j2 != None && suffix_R[j2] > next_score)
      {
        next_score = suffix_R[j2];
        index = 2;
        next_j = j2;
      }

      j = next_j;
      if (index == 0)
      {
        break;
      }
      else if (index == 1)
      {
        // do nothing
      }
      else
      {
        while (suffix_R[j] != suffix[j])
        {
          j += 1;
        }
      }
    }

    return pair<double, pair2>(best_score / parent_mass1, pair2(path1, path2));

  }

  //////////////////////////////////////////////////////////////////////////////

  pair<double, pair2> find_paths(double parent_mass1,
                                 double parent_mass2,
                                 vector<double> & common_scores,
                                 vector<double> & common2_scores,
                                 double ptm_penalty,
                                 vector<int> & prev,
                                 vector<int> & next,
                                 vector<int> & prev2,
                                 vector<int> & next2,
                                 vector<vector<int> > & left_jumps,
                                 vector<vector<int> > & right_jumps,
                                 vector<vector<int> > & left_jumps2,
                                 vector<vector<int> > & right_jumps2,
                                 vector<vector<int> > & left_neighbors,
                                 vector<vector<int> > & right_neighbors,
                                 int best_i,
                                 int best_j,
                                 int best_s,
                                 celltype best_t,
                                 double best_score,
                                 vector2 & D1,
                                 vector3 & D2,
                                 vector3 & D3,
                                 vector2 & D2max,
                                 vector2 & D3max,
                                 vector2 & M1_L,
                                 vector2 & M1_R,
                                 vector2 & M2_L,
                                 vector3 & M2_R,
                                 vector3 & M3_L,
                                 vector2 & M3_R,
                                 vector<double> & prefix1,
                                 vector2 & suffix1,
                                 vector2 & prefix2,
                                 vector<double> & suffix2,
                                 vector<double> & prefix1_L,
                                 vector<double> & suffix1max,
                                 vector<double> & suffix1_R,
                                 vector<double> & prefix2max,
                                 vector<double> & prefix2_L,
                                 vector<double> & suffix2_R);

  //////////////////////////////////////////////////////////////////////////////

  int max_ind(vector<double> & v)
  {
    return max_element(v.begin(), v.end()) - v.begin();
  }

  //////////////////////////////////////////////////////////////////////////////

  pair<double, pair2> align(double parent_mass1,
                            double parent_mass2,
                            vector<double> peaks,
                            vector<double> peaks2,
                            vector<int> & common,
                            vector<int> & common2,
                            vector<double> & common_scores,
                            vector<double> & common2_scores,
                            vector<int> & prev,
                            vector<int> & next,
                            vector<int> & prev2,
                            vector<int> & next2,
                            vector<vector<int> > & left_jumps,
                            vector<vector<int> > & right_jumps,
                            vector<vector<int> > & left_jumps2,
                            vector<vector<int> > & right_jumps2,
                            vector<vector<int> > & left_neighbors,
                            vector<vector<int> > & right_neighbors,
                            double same_vertex_penalty,
                            double ptm_penalty)

  // parent_mass1 : parent mass of first spectrum
  // parent mass2 : parent mass of 2nd spectrum. We assume that
  //                parent_mass2 >= parent_mass1
  // peaks        : All the peaks m in the 1st spectrum such that
  //                either the peak m or the peak m+(parent_mass2-parent_mass1)
  // Note: peaks must be symmetric, namely, peaks[i]+peaks[n-1-i]=parent_mass1+18
  //       for all i, where n is the length of peaks
  //
  // peaks2       : All the peaks of the 2st spectrum
  // common       : common[i] is the index j of the peak in spectrum2 whose mass
  //                is equal to peaks[i].
  //                If no such peak exists, common[i] = None
  // common2      : common[i] is the index j of the peak in spectrum2 whose mass
  //                is equal to peaks[i]+(parent_mass2-parent_mass1).
  // common_scores: common_scores[i] is the score of the vertex
  //                (peaks[i],peaks[i]) in the alignment graph.
  //                If the vertex doesn't exist, then common_scores[i]=-infinity
  // common2_scores: common_scores2[i] is the score of the vertex
  //                (peaks[i],peaks[i]+(parent_mass2-parent_mass1).
  //
  //
  // Let  T = the mass tolerance (e.g. 0.5)
  //
  // Let J be a set of masses (possibly empty), and M be a mass.
  // We say that a mass jump m is valid if either |m-m'| < T for some m' from J,
  // or m >= M.
  // Assume that if J is not empty, then min(J) >= 57 and max(J) < M
  //
  // Define
  //   delta0 = parent_mass2-parent_mass1
  //   delta = min(max(delta0+T, 57-T),114-T)
  //   delta2 = M-T
  //
  // prev         : prev[i] is the maximum index j such that
  //                peaks[j] < peaks[i]-delta2
  //                If no such index exists, prev[i] = None
  // next         : next[i] is the minimum index j such that
  //                peaks[j] > peaks[i]+delta2
  // prev2        : prev2[i] is the maximum index j such that
  //                peaks[j] < peaks[i]-max(delta, delta2).
  // next2        : next[i] is the minimum index j such that
  //                peaks[j] > peaks[i]+max(delta, delta2).
  //
  // left_jumps   : left_jumps[i] is a list of all the indices j such that
  //   (1) peaks[i]-delta2 <= peaks[j] < peaks[i]-(57-T)
  //   (2) peaks[i]-peaks[j] is a valid mass jump
  //   If J is empty, then left_jumps[i] is an empty list
  // right_jumps   : right_jumps[i] is a list of all the indices j such that
  //   (1) peaks[i]+57-T < peaks[j] <= peaks[i]+delta2
  //   (2) peaks[j]-peaks[i] is a valid mass jump
  //   If J is empty, then right_jumps[i] is an empty list
  //
  // left_jumps2  : left_jumps2[i] is a list of all the indices j such that
  //   (1) peaks[i]-delta2 <= peaks[j] <= peaks[i]-delta
  //   (2) peaks[i]-peaks[j] is a valid mass jump
  //   If delta >= delta2 then left_jumps2[i] is an empty list (follows from (1))
  // right_jumps2  : right_jumps2[i] is a list of all the indices j such that
  //   (1) peaks[i]+delta < peaks[j] <= peaks[i]+delta2
  //   (2) peaks[j]-peaks[i] is a valid mass jump
  //   If delta >= delta2 then right_jumps2[i] is an empty list (follows from (1))
  //
  // left_neighbors: left_neighbors[i] is a list of all indices j such that
  //   (1) peaks[i]-delta <= peaks[j] < peaks[i]-(57-T)
  //   (2) peaks[i]-peaks[j] is a valid mass jump
  // right_neighbors: right_neighbors[i] is a list of all indices j such that
  //   (1) peaks[i]+57-T < peaks[j] <= peaks[i]+delta
  //   (2) peaks[j]-peaks[i] is a valid mass jump
  //
  // same_vertex_penalty: penalty for using the same vertex twice
  //                      (negative number)
  // ptm_penalty:   penalty for having a PTM (negative number)
  {
    int n0 = common.size();
    int N = n0 - 1;
    int n = (n0 + 1) / 2;

    vector<double> prefix1(n, -infinity);
    vector2 suffix1(n);
    vector2 prefix2(n);
    vector<double> suffix2(n, -infinity);

    vector<double> prefix1_L(n, -infinity);
    vector<double> suffix1max(n, -infinity);
    vector<double> suffix1_R(n, -infinity);
    vector<double> prefix2max(n, -infinity);
    vector<double> prefix2_L(n, -infinity);
    vector<double> suffix2_R(n, -infinity);

    vector2 D1(n, vector<double>(n, -infinity));
    vector3 D2(n, vector2(n));
    vector3 D3(n, vector2(n));
    vector2 D2max(n, vector<double>(n, -infinity));
    vector2 D3max(n, vector<double>(n, -infinity));
    vector2 M1_L(n, vector<double>(n, -infinity));
    vector2 M1_R(n, vector<double>(n, -infinity));
    vector2 M2_L(n, vector<double>(n, -infinity));
    vector3 M2_R(n, vector2(n));
    vector3 M3_L(n, vector2(n));
    vector2 M3_R(n, vector<double>(n, -infinity));

    // compute prefix1
    for (int i = 0; i < n; ++i)
    {
      if (common[i] != None)
      {
        double prev_score = 0.0;
        for (vector<int>::const_iterator it = left_jumps[i].begin();
            it != left_jumps[i].end(); ++it)
        {
          prev_score = max(prev_score, prefix1[*it]);
        }
        int i2 = prev[i];
        if (i2 != None)
        {
          prev_score = max(prev_score, prefix1_L[i2]);
        }
        prefix1[i] = common_scores[i] + prev_score;
      }

      if (i != 0)
        prefix1_L[i] = max(prefix1_L[i - 1], prefix1[i]);
      else
        prefix1_L[i] = prefix1[i];
    }

    // compute suffix1
    for (int j = 0; j < n; ++j)
    {
      int l = right_neighbors[N - j].size();
      suffix1[j].resize(l + 1, -infinity);

      if (common[N - j] != None)
      {
        // s = 0
        double next_score = 0.0;
        for (vector<int>::const_iterator it = right_jumps[N - j].begin();
            it != right_jumps[N - j].end(); ++it)
        {
          next_score = max(next_score, suffix2[N - *it] + ptm_penalty);
        }
        int j2 = next[N - j];
        if (j2 != None)
        {
          next_score = max(next_score, suffix2_R[N - j2] + ptm_penalty);
        }
        for (vector<int>::const_iterator it = right_jumps2[N - j].begin();
            it != right_jumps2[N - j].end(); ++it)
        {
          next_score = max(next_score, suffix1max[N - *it]);
        }
        j2 = next2[N - j];
        if (j2 != None)
        {
          next_score = max(next_score, suffix1_R[N - j2]);
        }

        suffix1[j][0] = common_scores[N - j] + next_score;
        // s > 0
        for (int s = 0; s < l; ++s)
        {
          int j2 = right_neighbors[N - j][s];
          if (common[j2] == None)
            continue;
          double penalty = 0.0;
          double y = peaks2[common[N - j]] + peaks2[common[j2]];
          if (fabs(y - (parent_mass2 + 18)) < mass_tolerance)
            penalty = same_vertex_penalty;
          double next_score = suffix1max[N - j2];
          suffix1[j][s + 1] = common_scores[N - j] + next_score + penalty;
        }
      }

      suffix1max[j] = *max_element(suffix1[j].begin(), suffix1[j].end());
      if (j != 0)
        suffix1_R[j] = max(suffix1_R[j - 1], suffix1max[j]);
      else
        suffix1_R[j] = suffix1max[j];
    }

    // compute prefix2
    for (int i = 0; i < n; ++i)
    {
      int l = left_neighbors[i].size();
      prefix2[i].resize(l + 1, -infinity);

      if (common2[i] != None)
      {
        // s = 0
        double prev_score = 0.0;
        for (vector<int>::const_iterator it = left_jumps[i].begin();
            it != left_jumps[i].end(); ++it)
        {
          prev_score = max(prev_score, prefix1[*it] + ptm_penalty);
        }
        int i2 = prev[i];
        if (i2 != None)
          prev_score = max(prev_score, prefix1_L[i2] + ptm_penalty);
        for (vector<int>::const_iterator it = left_jumps2[i].begin();
            it != left_jumps2[i].end(); ++it)
        {
          prev_score = max(prev_score, prefix2max[*it]);
        }
        i2 = prev2[i];
        if (i2 != None)
        {
          prev_score = max(prev_score, prefix2_L[i2]);
        }
        prefix2[i][0] = common2_scores[i] + prev_score;
        // s > 0
        for (int s = 0; s < l; ++s)
        {
          int i2 = left_neighbors[i][s];
          if (common2[i2] == None)
            continue;
          double penalty = 0.0;
          double y = peaks2[common2[i]] + peaks2[common2[i2]];
          if (fabs(y - (parent_mass2 + 18)) < mass_tolerance)
            penalty = same_vertex_penalty;
          prev_score = prefix2max[i2];
          prefix2[i][s + 1] = common2_scores[i] + prev_score + penalty;
        }
      }

      prefix2max[i] = *max_element(prefix2[i].begin(), prefix2[i].end());
      if (i != 0)
        prefix2_L[i] = max(prefix2_L[i - 1], prefix2max[i]);
      else
        prefix2_L[i] = prefix2max[i];
    }

    // compute suffix2
    for (int j = 0; j < n; ++j)
    {
      if (common2[N - j] != None)
      {
        double next_score = 0.0;
        for (vector<int>::const_iterator it = right_jumps[N - j].begin();
            it != right_jumps[N - j].end(); ++it)
        {
          next_score = max(next_score, suffix2[N - *it]);
        }
        int j2 = next[N - j];
        if (j2 != None)
        {
          next_score = max(next_score, suffix2_R[N - j2]);
        }
        suffix2[j] = common2_scores[N - j] + next_score;
      }
      if (j != 0)
        suffix2_R[j] = max(suffix2_R[j - 1], suffix2[j]);
      else
        suffix2_R[j] = suffix2[j];
    }

    // compute D1/M1_R
    for (int i = 0; i < n; ++i)
    {
      if (common[i] == None)
      {
        if (i != 0)
          M1_L[i] = M1_L[i - 1];
        continue;
      }

      for (int j = 0; j < n; ++j)
      {
        if (common2[N - j] == None)
        {
          if (i != 0)
            M1_L[i][j] = M1_L[i - 1][j];
          if (j != 0)
            M1_R[i][j] = M1_R[i][j - 1];
          continue;
        }

        double penalty = 0.0;
        double x = peaks[i] + peaks[N - j];
        double y = peaks2[common[i]] + peaks2[common2[N - j]];
        if (fabs(x - (parent_mass1 + 18)) < mass_tolerance
            || fabs(y - (parent_mass2 + 18)) < mass_tolerance)
          penalty = same_vertex_penalty;

        if (i >= j)
        {
          double prev_score = suffix2[j] + ptm_penalty;
          // If we choose suffix2[j], then we pay the PTM penalty -
          //  the penalty is paid only once!
          for (vector<int>::const_iterator it = left_jumps[i].begin();
              it != left_jumps[i].end(); ++it)
          {
            prev_score = max(prev_score, D1[*it][j]);
          }
          int i2 = prev[i];
          if (i2 != None)
            prev_score = max(prev_score, M1_L[i2][j]);
          D1[i][j] = common_scores[i] + prev_score + penalty;
        }
        else
        {
          double next_score = prefix1[i] + ptm_penalty;
          for (vector<int>::const_iterator it = right_jumps[N - j].begin();
              it != right_jumps[N - j].end(); ++it)
          {
            next_score = max(next_score, D1[i][N - *it]);
          }
          int j2 = next[N - j];
          if (j2 != None)
          {
            next_score = max(next_score, M1_R[i][N - j2]);
          }
          D1[i][j] = common2_scores[N - j] + next_score + penalty;
        }

        // Compute M1_L
        if (i != 0)
          M1_L[i][j] = max(M1_L[i - 1][j], D1[i][j]);
        else
          M1_L[i][j] = D1[i][j];

        // Compute M1_R
        if (j != 0)
          M1_R[i][j] = max(M1_R[i][j - 1], D1[i][j]);
        else
          M1_R[i][j] = D1[i][j];
      }
    }

    // compute D2/D2max/M2_L
    for (int i = 0; i < n; ++i)
    {
      int l = left_neighbors[i].size();
      for (int j = 0; j < n; ++j)
      {
        D2[i][j].resize(l + 1, -infinity);
        M2_R[i][j].resize(l + 1, -infinity);
      }
      if (common2[i] == None)
      {
        if (i > 0)
          M2_L[i] = M2_L[i - 1];
        continue;
      }

      for (int j = 0; j < n; ++j)
      {
        if (common2[N - j] == None)
        {
          if (i != 0)
            M2_L[i][j] = M2_L[i - 1][j];
          if (j != 0)
            M2_R[i][j] = M2_R[i][j - 1];
          continue;
        }

        double penalty = 0.0;
        double x = peaks[i] + peaks[N - j];
        double y = peaks2[common2[i]] + peaks2[common2[N - j]];
        if (fabs(x - (parent_mass1 + 18)) < mass_tolerance
            || fabs(y - (parent_mass2 + 18)) < mass_tolerance)
          penalty = same_vertex_penalty;

        if (i > j)
        {
          // s = 0
          double prev_score = suffix2[j];
          for (vector<int>::const_iterator it = left_jumps[i].begin();
              it != left_jumps[i].end(); ++it)
          {
            prev_score = max(prev_score, D1[*it][j]);
          }
          int i2 = prev[i];
          if (i2 != None)
            prev_score = max(prev_score, M1_L[i2][j]);
          for (vector<int>::const_iterator it = left_jumps2[i].begin();
              it != left_jumps2[i].end(); ++it)
          {
            prev_score = max(prev_score, D2max[*it][j]);
          }
          i2 = prev2[i];
          if (i2 != None)
            prev_score = max(prev_score, M2_L[i2][j]);
          D2[i][j][0] = common2_scores[i] + prev_score + penalty;

          // s > 0
          for (int s = 0; s < l; ++s)
          {
            int i2 = left_neighbors[i][s];
            if (common2[i2] == None)
              continue;
            double penalty2 = penalty;
            double y = peaks2[common2[i]] + peaks2[common2[i2]];
            if (fabs(y - (parent_mass2 + 18)) < mass_tolerance)
              penalty2 += same_vertex_penalty;
            prev_score = D2max[i2][j];
            D2[i][j][s + 1] = common2_scores[i] + prev_score + penalty2;
          }
        }
        else
        { // i <= j
          // s = 0
          double next_score = prefix2[i][0];
          for (vector<int>::const_iterator it = right_jumps[N - j].begin();
              it != right_jumps[N - j].end(); ++it)
          {
            next_score = max(next_score, D2[i][N - *it][0]);
          }
          int j2 = next[N - j];
          if (j2 != None)
          {
            next_score = max(next_score, M2_R[i][N - j2][0]);
          }
          D2[i][j][0] = common2_scores[N - j] + next_score + penalty;

          // s > 0
          for (int s = 0; s < l; ++s)
          {
            int i2 = left_neighbors[i][s];
            if (common2[i2] == None)
              continue;
            double penalty2 = penalty;
            double y = peaks2[common2[i2]] + peaks2[common2[N - j]];
            if (fabs(y - (parent_mass2 + 18)) < mass_tolerance)
              penalty2 += same_vertex_penalty;
            double next_score = prefix2[i][s + 1];
            for (vector<int>::const_iterator it = right_jumps[N - j].begin();
                it != right_jumps[N - j].end(); ++it)
            {
              next_score = max(next_score, D2[i][N - *it][s + 1]);
            }
            int j2 = next[N - j];
            if (j2 != None)
              next_score = max(next_score, M2_R[i][N - j2][s + 1]);
            D2[i][j][s + 1] = common2_scores[N - j] + next_score + penalty2;
          }
        }

        // Compute D2max
        D2max[i][j] = *max_element(D2[i][j].begin(), D2[i][j].end());

        // Compute M2_L
        if (i != 0)
          M2_L[i][j] = max(M2_L[i - 1][j], D2max[i][j]);
        else
          M2_L[i][j] = D2max[i][j];

        // compute M2_R
        if (j != 0)
        {
          for (int s = 0; s < l + 1; ++s)
          {
            M2_R[i][j][s] = max(M2_R[i][j - 1][s], D2[i][j][s]);
          }
        }
        else
          M2_R[i][j] = D2[i][j];
      }
    }

    // compute D3/M3_R
    for (int i = 0; i < n; ++i)
    {
      for (int j = 0; j < n; ++j)
      {
        int l = right_neighbors[N - j].size();
        D3[i][j].resize(l + 1, -infinity);
        M3_L[i][j].resize(l + 1, -infinity);
      }
      if (common[i] == None)
      {
        if (i > 0)
          M3_L[i] = M3_L[i - 1];
        continue;
      }

      for (int j = 0; j < n; ++j)
      {
        if (common[N - j] == None)
        {
          if (i != 0)
            M3_L[i][j] = M3_L[i - 1][j];
          if (j != 0)
            M3_R[i][j] = M3_R[i][j - 1];
          continue;
        }

        int l = right_neighbors[N - j].size();
        double penalty = 0.0;
        double x = peaks[i] + peaks[N - j];
        double y = peaks2[common[i]] + peaks2[common[N - j]];
        if (fabs(x - (parent_mass1 + 18)) < mass_tolerance
            || fabs(y - (parent_mass2 + 18)) < mass_tolerance)
          penalty = same_vertex_penalty;

        if (i >= j)
        {
          // s = 0
          double prev_score = suffix1[j][0];
          for (vector<int>::const_iterator it = left_jumps[i].begin();
              it != left_jumps[i].end(); ++it)
          {
            prev_score = max(prev_score, D3[*it][j][0]);
          }
          int i2 = prev[i];
          if (i2 != None)
            prev_score = max(prev_score, M3_L[i2][j][0]);
          D3[i][j][0] = common_scores[i] + prev_score + penalty;
          // s > 0
          for (int s = 0; s < l; ++s)
          {
            int j2 = right_neighbors[N - j][s];
            if (common[j2] == None)
              continue;
            double penalty2 = penalty;
            double y = peaks2[common[i]] + peaks2[common[j2]];
            if (fabs(y - (parent_mass2 + 18)) < mass_tolerance)
              penalty2 += same_vertex_penalty;
            double prev_score = suffix1[j][s + 1];
            for (vector<int>::const_iterator it = left_jumps[i].begin();
                it != left_jumps[i].end(); ++it)
            {
              prev_score = max(prev_score, D3[*it][j][s + 1]);
            }
            int i2 = prev[i];
            if (i2 != None)
              prev_score = max(prev_score, M3_L[i2][j][s + 1]);
            D3[i][j][s + 1] = common_scores[i] + prev_score + penalty2;
          }
        }
        else
        {
          // s = 0
          double next_score = prefix1[i];
          for (vector<int>::const_iterator it = right_jumps[N - j].begin();
              it != right_jumps[N - j].end(); ++it)
          {
            next_score = max(next_score, D1[i][N - *it]);
          }
          int j2 = next[N - j];
          if (j2 != None)
          {
            next_score = max(next_score, M1_R[i][N - j2]);
          }
          for (vector<int>::const_iterator it = right_jumps2[N - j].begin();
              it != right_jumps2[N - j].end(); ++it)
          {
            next_score = max(next_score, D3max[i][N - *it]);
          }
          j2 = next2[N - j];
          if (j2 != None)
          {
            next_score = max(next_score, M3_R[i][N - j2]);
          }
          D3[i][j][0] = common_scores[N - j] + next_score + penalty;
          // s > 0
          for (int s = 0; s < l; ++s)
          {
            int j2 = right_neighbors[N - j][s];
            if (common[j2] == None)
              continue;
            double penalty2 = penalty;
            double y = peaks2[common[N - j]] + peaks2[common[j2]];
            if (fabs(y - (parent_mass2 + 18)) < mass_tolerance)
              penalty2 += same_vertex_penalty;
            double next_score = D3max[i][N - j2];
            D3[i][j][s + 1] = common_scores[N - j] + next_score + penalty2;
          }
        }

        // Compute D3max
        D3max[i][j] = *std::max_element(D3[i][j].begin(), D3[i][j].end());

        // Compute M3_L
        if (i != 0)
        {
          for (int s = 0; s < l + 1; ++s)
          {
            M3_L[i][j][s] = max(M3_L[i - 1][j][s], D3[i][j][s]);
          }
        }
        else
          M3_L[i][j] = D3[i][j];

        // Compute M3_R
        if (j != 0)
          M3_R[i][j] = max(M3_R[i][j - 1], D3max[i][j]);
        else
          M3_R[i][j] = D3max[i][j];
      }
    }

    // Find best score
    double best_score = 0.0;
    int best_i = None;
    int best_j = None;
    int best_s = None;
    celltype best_t = INVALID;
    for (int i = 0; i < n; ++i)
    {
      for (vector<int>::const_iterator it = right_jumps[i].begin();
          it != right_jumps[i].end(); ++it)
      {
        int j = N - *it;
        if (j <= n - 1)
        {
          if (best_score < D1[i][j])
          {
            best_score = D1[i][j];
            best_i = i;
            best_j = j;
            best_s = 0;
            best_t = cell_D1;
          }
          if (best_score < D2max[i][j])
          {
            best_score = D2max[i][j];
            best_i = i;
            best_j = j;
            best_s = max_ind(D2[i][j]);
            best_t = cell_D2;
          }
          if (best_score < D3max[i][j])
          {
            best_score = D3max[i][j];
            best_i = i;
            best_j = j;
            best_s = max_ind(D3[i][j]);
            best_t = cell_D3;
          }
        }
      }

      int j0 = next[i];
      if (j0 != None)
      {
        int j = min(n - 1, N - j0);
        if (best_score < M1_R[i][j])
        {
          best_score = M1_R[i][j];
          best_i = i;
          best_j = j;
          best_s = 0;
          best_t = cell_M1_R;
        }

        double tmp = *max_element(M2_R[i][j].begin(), M2_R[i][j].end());
        if (best_score < tmp)
        {
          best_score = tmp;
          best_i = i;
          best_j = j;
          best_s = max_ind(M2_R[i][j]);
          best_t = cell_M2_R;
        }

        if (best_score < M3_R[i][j])
        {
          best_score = M3_R[i][j];
          best_i = i;
          best_j = j;
          best_s = 0;
          best_t = cell_M3_R;
        }
      }

      if (best_score < prefix1[i])
      {
        best_score = prefix1[i];
        best_i = i;
        best_j = 0;
        best_s = 0;
        best_t = cell_prefix1;
      }

      if (best_score < suffix1max[i])
      {
        best_score = suffix1max[i];
        best_j = i;
        best_i = 0;
        best_s = max_ind(suffix1[i]);
        best_t = cell_suffix1;
      }

      if (best_score < prefix2max[i])
      {
        best_score = prefix2max[i];
        best_i = i;
        best_j = 0;
        best_s = max_ind(prefix2[i]);
        best_t = cell_prefix2;
      }

      if (best_score < suffix2[i])
      {
        best_score = suffix2[i];
        best_j = i;
        best_i = 0;
        best_s = 0;
        best_t = cell_suffix2;
      }
    }

    return find_paths(parent_mass1,
                      parent_mass2,
                      common_scores,
                      common2_scores,
                      ptm_penalty,
                      prev,
                      next,
                      prev2,
                      next2,
                      left_jumps,
                      right_jumps,
                      left_jumps2,
                      right_jumps2,
                      left_neighbors,
                      right_neighbors,
                      best_i,
                      best_j,
                      best_s,
                      best_t,
                      best_score,
                      D1,
                      D2,
                      D3,
                      D2max,
                      D3max,
                      M1_L,
                      M1_R,
                      M2_L,
                      M2_R,
                      M3_L,
                      M3_R,
                      prefix1,
                      suffix1,
                      prefix2,
                      suffix2,
                      prefix1_L,
                      suffix1max,
                      suffix1_R,
                      prefix2max,
                      prefix2_L,
                      suffix2_R);
  }

  //////////////////////////////////////////////////////////////////////////////

  pair<double, pair2> find_paths(double parent_mass1,
                                 double parent_mass2,
                                 vector<double> & common_scores,
                                 vector<double> & common2_scores,
                                 double ptm_penalty,
                                 vector<int> & prev,
                                 vector<int> & next,
                                 vector<int> & prev2,
                                 vector<int> & next2,
                                 vector<vector<int> > & left_jumps,
                                 vector<vector<int> > & right_jumps,
                                 vector<vector<int> > & left_jumps2,
                                 vector<vector<int> > & right_jumps2,
                                 vector<vector<int> > & left_neighbors,
                                 vector<vector<int> > & right_neighbors,
                                 int best_i,
                                 int best_j,
                                 int best_s,
                                 celltype best_t,
                                 double best_score,
                                 vector2 & D1,
                                 vector3 & D2,
                                 vector3 & D3,
                                 vector2 & D2max,
                                 vector2 & D3max,
                                 vector2 & M1_L,
                                 vector2 & M1_R,
                                 vector2 & M2_L,
                                 vector3 & M2_R,
                                 vector3 & M3_L,
                                 vector2 & M3_R,
                                 vector<double> & prefix1,
                                 vector2 & suffix1,
                                 vector2 & prefix2,
                                 vector<double> & suffix2,
                                 vector<double> & prefix1_L,
                                 vector<double> & suffix1max,
                                 vector<double> & suffix1_R,
                                 vector<double> & prefix2max,
                                 vector<double> & prefix2_L,
                                 vector<double> & suffix2_R)
  {
    int n0 = prev.size();
    int N = n0 - 1;

    vector<int> path1;
    vector<int> path2;

    if (best_i == None)
      return pair<double, pair2>(best_score / parent_mass1, pair2(path1, path2));

    int i = best_i;
    int j = best_j;
    int s = best_s;
    celltype t = best_t;
    while (1)
    {
      // printf("i/j/s/t %d %d %d %d\n",i,j,s,t);
      if (t == cell_prefix1)
      {
        path1.insert(path1.begin(), i);
        double prev_score = 0.0;
        int index = 0;
        int next_i = 0;
        for (vector<int>::const_iterator it = left_jumps[i].begin();
            it != left_jumps[i].end(); ++it)
        {
          int i2 = *it;
          if (prefix1[i2] > prev_score)
          {
            prev_score = prefix1[i2];
            index = 1;
            next_i = i2;
          }
        }
        int i2 = prev[i];
        if (i2 != None && prefix1_L[i2] > prev_score)
        {
          prev_score = prefix1_L[i2];
          index = 2;
          next_i = i2;
        }

        i = next_i;
        if (index == 0)
        {
          break;
        }
        else if (index == 1)
        {
          // t is unchanged
        }
        else
        {
          t = cell_prefix1_L;
        }

      }
      else if (t == cell_suffix1)
      {
        path1.push_back(N - j);
        if (s == 0)
        {
          double next_score = 0.0;
          int index = 0;
          int next_j = 0;
          for (vector<int>::const_iterator it = right_jumps[N - j].begin();
              it != right_jumps[N - j].end(); ++it)
          {
            int j2 = *it;
            if (suffix2[N - j2] + ptm_penalty > next_score)
            {
              next_score = suffix2[N - j2] + ptm_penalty;
              index = 1;
              next_j = N - j2;
            }
          }
          int j2 = next[N - j];
          if (j2 != None && suffix2_R[N - j2] + ptm_penalty > next_score)
          {
            next_score = suffix2_R[N - j2] + ptm_penalty;
            index = 2;
            next_j = N - j2;
          }
          for (vector<int>::const_iterator it = right_jumps2[N - j].begin();
              it != right_jumps2[N - j].end(); ++it)
          {
            int j2 = *it;
            if (suffix1max[N - j2] > next_score)
            {
              next_score = suffix1max[N - j2];
              index = 3;
              next_j = N - j2;
            }
          }
          j2 = next2[N - j];
          if (j2 != None && suffix1_R[N - j2] > next_score)
          {
            next_score = suffix1_R[N - j2];
            index = 4;
            next_j = N - j2;
          }

          j = next_j;
          if (index == 0)
          {
            break;
          }
          else if (index == 1)
          {
            t = cell_suffix2;
          }
          else if (index == 2)
          {
            t = cell_suffix2_R;
          }
          else if (index == 3)
          {
            // t is unchanged
            s = max_ind(suffix1[j]);
          }
          else
          {
            t = cell_suffix1_R;
          }
        }
        else
        {
          // t is unchanged
          j = right_neighbors[N - j][s - 1];
          j = N - j;
          s = max_ind(suffix1[j]);
        }

      }
      else if (t == cell_prefix2)
      {
        path2.insert(path2.begin(), i);
        if (s == 0)
        {
          double prev_score = 0.0;
          int index = 0;
          int next_i = 0;
          for (vector<int>::const_iterator it = left_jumps[i].begin();
              it != left_jumps[i].end(); ++it)
          {
            int i2 = *it;
            if (prefix1[i2] + ptm_penalty > prev_score)
            {
              prev_score = prefix1[i2] + ptm_penalty;
              index = 1;
              next_i = i2;
            }
          }
          int i2 = prev[i];
          if (i2 != None && prefix1_L[i2] + ptm_penalty > prev_score)
          {
            prev_score = prefix1_L[i2] + ptm_penalty;
            index = 2;
            next_i = i2;
          }
          for (vector<int>::const_iterator it = left_jumps2[i].begin();
              it != left_jumps2[i].end(); ++it)
          {
            int i2 = *it;
            if (prefix2max[i2] > prev_score)
            {
              prev_score = prefix2max[i2];
              index = 3;
              next_i = i2;
            }
          }
          i2 = prev2[i];
          if (i2 != None && prefix2_L[i2] > prev_score)
          {
            prev_score = prefix2_L[i2];
            index = 4;
            next_i = i2;
          }

          i = next_i;
          if (index == 0)
          {
            break;
          }
          else if (index == 1)
          {
            t = cell_prefix1;
          }
          else if (index == 2)
          {
            t = cell_prefix1_L;
          }
          else if (index == 3)
          {
            // t is unchanged
            s = max_ind(prefix2[i]);
          }
          else
          {
            t = cell_prefix2_L;
          }
        }
        else
        {
          // t is unchanged
          i = left_neighbors[i][s - 1];
          s = max_ind(prefix2[i]);
        }

      }
      else if (t == cell_suffix2)
      {
        path2.push_back(N - j);
        double next_score = 0.0;
        int index = 0;
        int next_j = 0;
        for (vector<int>::const_iterator it = right_jumps[N - j].begin();
            it != right_jumps[N - j].end(); ++it)
        {
          int j2 = *it;
          if (suffix2[N - j2] > next_score)
          {
            next_score = suffix2[N - j2];
            index = 1;
            next_j = N - j2;
          }
        }
        int j2 = next[N - j];
        if (j2 != None && suffix2_R[N - j2] > next_score)
        {
          next_score = suffix2_R[N - j2];
          index = 2;
          next_j = N - j2;
        }

        j = next_j;
        if (index == 0)
        {
          break;
        }
        else if (index == 1)
        {
          // t is unchanged
        }
        else
        {
          t = cell_suffix2_R;
        }

      }
      else if (t == cell_prefix1_L)
      {
        if (prefix1_L[i] == prefix1[i])
        {
          t = cell_prefix1;
        }
        else
        {
          i -= 1;
        }

      }
      else if (t == cell_prefix2_L)
      {
        if (prefix2_L[i] == prefix2max[i])
        {
          t = cell_prefix2;
          s = max_ind(prefix2[i]);
        }
        else
        {
          i -= 1;
        }

      }
      else if (t == cell_suffix1_R)
      {
        if (suffix1_R[j] == suffix1max[j])
        {
          t = cell_suffix1;
          s = max_ind(suffix1[j]);
        }
        else
        {
          j -= 1;
        }

      }
      else if (t == cell_suffix2_R)
      {
        if (suffix2_R[j] == suffix2[j])
        {
          t = cell_suffix2;
        }
        else
        {
          j -= 1;
        }

      }
      else if (t == cell_D1)
      {
        if (i >= j)
        {
          path1.insert(path1.begin(), i);
          double prev_score = suffix2[j] + ptm_penalty;
          int index = 0;
          int next_i = 0;
          for (vector<int>::const_iterator it = left_jumps[i].begin();
              it != left_jumps[i].end(); ++it)
          {
            int i2 = *it;
            if (D1[i2][j] > prev_score)
            {
              prev_score = D1[i2][j];
              index = 1;
              next_i = i2;
            }
          }
          int i2 = prev[i];
          if (i2 != None && M1_L[i2][j] > prev_score)
          {
            prev_score = M1_L[i2][j];
            index = 2;
            next_i = i2;
          }

          i = next_i;
          // note that if index=0, then the value of i is irrelevant.
          if (index == 0)
          {
            t = cell_suffix2;
          }
          else if (index == 1)
          {
            // t is unchanged
          }
          else
          {
            t = cell_M1_L;
          }
        }
        else
        {
          path2.push_back(N - j);
          double next_score = prefix1[i] + ptm_penalty;
          int index = 0;
          int next_j = 0;
          for (vector<int>::const_iterator it = right_jumps[N - j].begin();
              it != right_jumps[N - j].end(); ++it)
          {
            int j2 = *it;
            if (D1[i][N - j2] > next_score)
            {
              next_score = D1[i][N - j2];
              index = 1;
              next_j = N - j2;
            }
          }
          int j2 = next[N - j];
          if (j2 != None && M1_R[i][N - j2] > next_score)
          {
            next_score = M1_R[i][N - j2];
            index = 2;
            next_j = N - j2;
          }

          j = next_j;
          if (index == 0)
          {
            t = cell_prefix1;
          }
          else if (index == 1)
          {
            // t is unchanged
          }
          else
          {
            t = cell_M1_R;
          }
        }

      }
      else if (t == cell_D2)
      {
        if (i > j)
        {
          path2.insert(path2.begin(), i);
          if (s == 0)
          {
            double prev_score = suffix2[j];
            int index = 0;
            int next_i = 0;
            for (vector<int>::const_iterator it = left_jumps[i].begin();
                it != left_jumps[i].end(); ++it)
            {
              int i2 = *it;
              if (D1[i2][j] > prev_score)
              {
                prev_score = D1[i2][j];
                index = 1;
                next_i = i2;
              }
            }
            int i2 = prev[i];
            if (i2 != None && M1_L[i2][j] > prev_score)
            {
              prev_score = M1_L[i2][j];
              index = 2;
              next_i = i2;
            }
            for (vector<int>::const_iterator it = left_jumps2[i].begin();
                it != left_jumps2[i].end(); ++it)
            {
              int i2 = *it;
              if (D2max[i2][j] > prev_score)
              {
                prev_score = D2max[i2][j];
                index = 3;
                next_i = i2;
              }
            }
            i2 = prev2[i];
            if (i2 != None && M2_L[i2][j] > prev_score)
            {
              prev_score = M2_L[i2][j];
              index = 4;
              next_i = i2;
            }

            i = next_i;
            if (index == 0)
            {
              t = cell_suffix2;
            }
            else if (index == 1)
            {
              t = cell_D1;
            }
            else if (index == 2)
            {
              t = cell_M1_L;
            }
            else if (index == 3)
            {
              // t is unchanged
              s = max_ind(D2[i][j]);
            }
            else
            {
              t = cell_M2_L;
            }
          }
          else
          {
            // t is unchanged
            i = left_neighbors[i][s - 1];
            s = max_ind(D2[i][j]);
          }
        }
        else
        { // i <= j
          path2.push_back(N - j);
          double next_score = prefix2[i][s];
          int index = 0;
          int next_j = 0;
          for (vector<int>::const_iterator it = right_jumps[N - j].begin();
              it != right_jumps[N - j].end(); ++it)
          {
            int j2 = *it;
            if (D2[i][N - j2][s] > next_score)
            {
              next_score = D2[i][N - j2][s];
              index = 1;
              next_j = N - j2;
            }
          }
          int j2 = next[N - j];
          if (j2 != None && M2_R[i][N - j2][s] > next_score)
          {
            index = 2;
            next_j = N - j2;
          }

          j = next_j;
          if (index == 0)
          {
            t = cell_prefix2;
          }
          else if (index == 1)
          {
            // t is unchanged
          }
          else
          {
            t = cell_M2_R;
          }
        }

      }
      else if (t == cell_D3)
      {
        if (i >= j)
        {
          path1.insert(path1.begin(), i);
          double prev_score = suffix1[j][s];
          int index = 0;
          int next_i = 0;
          for (vector<int>::const_iterator it = left_jumps[i].begin();
              it != left_jumps[i].end(); ++it)
          {
            int i2 = *it;
            if (D3[i2][j][s] > prev_score)
            {
              prev_score = D3[i2][j][s];
              index = 1;
              next_i = i2;
            }
          }
          int i2 = prev[i];
          if (i2 != None && M3_L[i2][j][s] > prev_score)
          {
            index = 2;
            next_i = i2;
          }

          i = next_i;
          if (index == 0)
          {
            t = cell_suffix1;
          }
          else if (index == 1)
          {
            // t is unchanged
          }
          else
          {
            t = cell_M3_L;
          }
        }
        else
        { // i < j
          path1.push_back(N - j);
          if (s == 0)
          {
            double next_score = prefix1[i];
            int index = 0;
            int next_j = 0;
            for (vector<int>::const_iterator it = right_jumps[N - j].begin();
                it != right_jumps[N - j].end(); ++it)
            {
              int j2 = *it;
              if (D1[i][N - j2] > next_score)
              {
                next_score = D1[i][N - j2];
                index = 1;
                next_j = N - j2;
              }
            }
            int j2 = next[N - j];
            if (j2 != None && M1_R[i][N - j2] > next_score)
            {
              next_score = M1_R[i][N - j2];
              index = 2;
              next_j = N - j2;
            }
            for (vector<int>::const_iterator it = right_jumps2[N - j].begin();
                it != right_jumps2[N - j].end(); ++it)
            {
              int j2 = *it;
              if (D3max[i][N - j2] > next_score)
              {
                next_score = D3max[i][N - j2];
                index = 3;
                next_j = N - j2;
              }
            }
            j2 = next2[N - j];
            if (j2 != None && M3_R[i][N - j2] > next_score)
            {
              next_score = M3_R[i][N - j2];
              index = 4;
              next_j = N - j2;
            }

            j = next_j;
            if (index == 0)
            {
              t = cell_prefix1;
            }
            else if (index == 1)
            {
              t = cell_D1;
            }
            else if (index == 2)
            {
              t = cell_M1_R;
            }
            else if (index == 3)
            {
              // t is unchanged
              s = max_ind(D3[i][j]);
            }
            else
            {
              t = cell_M3_R;
            }
          }
          else
          {
            // t is unchanged
            j = right_neighbors[N - j][s - 1];
            j = N - j;
            s = max_ind(D3[i][j]);
          }
        }

      }
      else if (t == cell_M1_L)
      {
        if (M1_L[i][j] == D1[i][j])
        {
          t = cell_D1;
        }
        else
        {
          i -= 1;
        }

      }
      else if (t == cell_M2_L)
      {
        if (M2_L[i][j] == D2max[i][j])
        {
          t = cell_D2;
          s = max_ind(D2[i][j]);
        }
        else
        {
          i -= 1;
        }

      }
      else if (t == cell_M3_L)
      {
        if (M3_L[i][j][s] == D3[i][j][s])
        {
          t = cell_D3;
        }
        else
        {
          i -= 1;
        }

      }
      else if (t == cell_M1_R)
      {
        if (M1_R[i][j] == D1[i][j])
        {
          t = cell_D1;
        }
        else
        {
          j -= 1;
        }

      }
      else if (t == cell_M2_R)
      {
        if (M2_R[i][j][s] == D2[i][j][s])
        {
          t = cell_D2;
        }
        else
        {
          j -= 1;
        }

      }
      else if (t == cell_M3_R)
      {
        if (M3_R[i][j] == D3max[i][j])
        {
          t = cell_D3;
          s = max_ind(D3[i][j]);
        }
        else
        {
          j -= 1;
        }
      }
    }

    return pair<double, pair2>(best_score / parent_mass1, pair2(path1, path2));
  }

}
;
// namespace DEKEL
