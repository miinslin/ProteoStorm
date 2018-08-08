// Header Includes
#include "AlignmentUtils.h"
#include "Logger.h"

// Module Includes

// External Includes
#include "alignment_scoring.h"
//#include "dekel_align.h"

using namespace specnets;

namespace specnets
{
  //---------------------------------------------------------------------------
  // auxUpdateMeanVariance - Updates mean/variance values to also include newValue.
  //
  // IMPORTANT NOTE: sampleSize in NOT UPDATED because multiple different means/variances may
  //                    be updated on the same sample (e.g. num peaks, matched intensities) and
  //                    sampleSize should only be updated once.
  //---------------------------------------------------------------------------
  inline void auxUpdateMeanVariance(float newValue,
                                    float sampleSize,
                                    float &mean,
                                    float &variance)
  {
    float updateRatio = sampleSize / (sampleSize + 1);
    sampleSize++;
    mean = mean * updateRatio + newValue / sampleSize;
    variance = updateRatio * variance + (newValue * newValue) / sampleSize;
  }

  //---------------------------------------------------------------------------
  void mergeVectors(vector<int> &v1, vector<int> &v2, vector<int> &vOut)
  {
    unsigned int i1 = 0, i2 = 0, iOut = 0;
    vOut.resize(v1.size() + v2.size());
    while (i1 < v1.size() or i2 < v2.size())
    {
      if (i1 == v1.size())
      {
        vOut[iOut++] = v2[i2++];
        continue;
      }
      if (i2 == v2.size())
      {
        vOut[iOut++] = v1[i1++];
        continue;
      }
      if (v1[i1] == v2[i2])
      {
        vOut[iOut++] = v1[i1++];
        i2++;
        continue;
      }
      if (v1[i1] < v2[i2])
        vOut[iOut++] = v1[i1++];
      else
        vOut[iOut++] = v2[i2++];
    }
    vOut.resize(iOut);
  }

  //---------------------------------------------------------------------------
  //
  // specSet        - Set of spectra to align
  // baseSpectraIdx - Indices of the spectra to align to others on this run. Note that spectra only align to spectra with higher indices.
  // aaDiff         - Difference between matched spectra's parent masses must be at most the sum of aaDiffs amino acid masses
  //                    or -1 for computation of all shifts between minShift and maxShift
  // minShift, maxShift - Min and max parent mass differences between spectra (to compute ASP shifts)
  // pmTol          - Parent mass error tolerance
  // peakTol        - Peak mass error tolerance
  // minMatchRatio  - minimum ratio of matched score / max score (in both spectra) to retain a match
  // results        - Variable to hold results
  // resolution     - Resolution for parent mass histogram (to determine which spectra to match)
  //
  //---------------------------------------------------------------------------
  void getPairAlignsASP(SpecSet & specSet,
                        vector<int> & baseSpectraIdx,
                        short aaDiff,
                        float minShift,
                        float maxShift,
                        float pmTol,
                        float peakTol,
                        float minMatchRatio,
                        short minNumMatchedPeaks,
                        SpectrumPairSet & results,
                        vector<TwoValues<float> > & ratios,
                        vector<TwoValues<float> > & means,
                        vector<float> & varTerms,
                        ExecTagFilterPairs *tag_filter,
                        float resolution,
                        float symmetryOffset)
  {
    int j, massIdx;
    float pMass1;
    float shiftIncrement = resolution; // Step size used to increment shifts from -pmTol to pmTol
    list<int>::iterator curPair;
    float histResolution = 0.1; // Resolution of the parent mass histogram - only used for binning, independent of overall peak resolution settings
    vector<list<int> > hist;
    specSet.getMassesHistogram(hist, histResolution);
    TwoValues<float> curRatios;
    AAJumps jumps(-1);
    SpectrumPair curResult;

    // Determine parent mass differences for valid matches
	  if (aaDiff > 0)
	  {
		jumps.getjumps(aaDiff, resolution);
		jumps.forceJump(0);
		jumps.forceTolerance(pmTol, resolution);
		jumps.forceDoubleSided();
		jumps.forceUnique(resolution);
		maxShift = jumps.masses[jumps.size()-1] + 3;
	  }
	  DEBUG_VAR(maxShift);

	  {// shiftIncrement should be factor of pmTol.
		int denominator= 2;
		while( pmTol/denominator > peakTol/2 ){
			denominator++;
		}
		while( pmTol/denominator > resolution ){
			denominator++;
		}
		shiftIncrement = pmTol/denominator;
	  }
	  DEBUG_VAR(shiftIncrement);

    // Determine totalScores per spectrum for later computation of match ratios
    vector<float> totalSpecScores;
    totalSpecScores.resize(specSet.size());
    for (unsigned int i = 0; i < specSet.size(); i++)
    {
      totalSpecScores[i] = 0;
      for (unsigned int p = 0; p < specSet[i].size(); p++)
        totalSpecScores[i] += specSet[i][p][1];
    }

    // Initialize means vector
    means.resize(specSet.size());
    varTerms.resize(specSet.size());
    for (unsigned int i = 0; i < means.size(); i++)
    {
      means[i].set(0, 0);
      varTerms[i] = 0;
    }

    // Iterate through pairs using FindMatchPeaksAll and ScoreOverlap6
    //    vector<int> idx1,idx2;                             idx1.reserve(1500);               idx2.reserve(1500);
    vector<int> idxMatched1_zero, idxMatched1_other;
    idxMatched1_zero.reserve(1500);
    idxMatched1_other.reserve(1500);
    vector<int> idxMatched2_zero, idxMatched2_other;
    idxMatched2_zero.reserve(1500);
    idxMatched2_other.reserve(1500);
    vector<int> idxUnion;
    idxUnion.reserve(1500);
    vector<int>::iterator idxUnionEnd, idxUnionStart;
    int spec1, spec2, maxShiftInt = (int)ceil(maxShift / histResolution);
    minShift = max((float)0, minShift);
    short bestNumMatchedPeaks, numMatchedPeaks, numMatchedPeaksU1, numMatchedPeaksU2;
    float pmDiff, score1, score2, shift, bestScore1, bestScore2, bestShift;
    //    clock_t startTime = clock();    double curTime, totTime=0;
    list<int> candidates;
    list<int>::iterator curCandidate;
    for (unsigned int i = 0; i < baseSpectraIdx.size(); i++)
    {
      DEBUG_MSG( "Starting with spectrum " << baseSpectraIdx[i] << "...");

      // Find candidate spectra to compare to
      spec1 = baseSpectraIdx[i];
      pMass1 = specSet[spec1].parentMass;
      int massInt = (int)round(pMass1 / histResolution);
      candidates.clear(); // Make sure candidates is empty from previous run
      if (specSet[spec1].size() == 0 or specSet[spec1].parentMass < 400)
      {
        DEBUG_MSG("(spectrum is empty or parent mass too low)");
        continue;
      }
      // Find candidates in [minShift,maxShift]
      int upperLimit = max(0, (int)ceil((pMass1 - minShift) / histResolution)); // Look at interval below current parent mass
      for (massIdx = (int)floor((pMass1 - min(pMass1, maxShift))
          / histResolution); massIdx < upperLimit; massIdx++)
      {
        for (curPair = hist[massIdx].begin(); curPair != hist[massIdx].end();
            curPair++)
        {
          if (*curPair > spec1
              and fabs(pMass1 - specSet[*curPair].parentMass)
                  >= minShift - 2 * pmTol
              and fabs(pMass1 - specSet[*curPair].parentMass)
                  <= maxShift + 2 * pmTol)
            candidates.push_back(*curPair);
        }
      }
      upperLimit = min((int)hist.size(), massInt + maxShiftInt); // Look at interval above current parent mass
      for (massIdx = (int)floor((pMass1 + minShift) / histResolution);
          massIdx < upperLimit; massIdx++)
      {
        for (curPair = hist[massIdx].begin(); curPair != hist[massIdx].end();
            curPair++)
        {
          if (*curPair > spec1
              and fabs(pMass1 - specSet[*curPair].parentMass)
                  >= minShift - 2 * pmTol
              and fabs(pMass1 - specSet[*curPair].parentMass)
                  <= maxShift + 2 * pmTol)
            candidates.push_back(*curPair);
        }
      }

  /*    if (aaDiff > 0)
      {
        for (j = 0, massIdx = (int)round((pMass1 + jumps[j]) / histResolution);
            j < jumps.size(); j++)
          if ((int)round((pMass1 + jumps[j]) / histResolution) >= 0)
            break; // Skips masses<0
        for (; j < jumps.size(); j++)
        {
          if (fabs(jumps[j]) >= minShift - 2 * pmTol
              and fabs(jumps[j]) <= maxShift + 2 * pmTol)
            continue; // Already added above
          massIdx = (int)round((pMass1 + jumps[j]) / histResolution);
          if (massIdx >= hist.size())
            break;
          for (curPair = hist[massIdx].begin(); curPair != hist[massIdx].end();
              curPair++)
            if (*curPair > spec1
                and fabs(pMass1 + jumps[j] - specSet[*curPair].parentMass)
                    <= 2 * pmTol)
            {
              candidates.push_back(*curPair);
            }
        }
      }//*/

      candidates.sort();
      candidates.unique(); // The same spectrum may occur in multiple bins - make sure the code below only aligns once to each eligible spectrum
      DEBUG_MSG(" num. candidates is " << candidates.size());

	  if (tag_filter->hasTags()) {
		list<int> subsetOfCandidates;
		tag_filter->getPairs(spec1, candidates, subsetOfCandidates);
		candidates.swap(subsetOfCandidates);
		DEBUG_MSG(" (" << candidates.size() << " after tag filtering)");
	  }
	  DEBUG_MSG("..."); //*/

      curCandidate = candidates.begin();
      while (curCandidate != candidates.end())
      {
        spec2 = *curCandidate;
        curCandidate++;

        if (spec1 >= specSet.size() || spec2 >= specSet.size())
        {
          ERROR_MSG("ERROR: Alignment pair (" << spec1 << "," << spec2 << ") is invalid - PRM spectra set has only " << specSet.size() << " spectra!\n");
          exit(-1);
        }

        if (specSet[spec2].size() == 0 or specSet[spec2].parentMass < 400) continue;

        //cerr<<spec2<<":"; cerr.flush();
        pmDiff = specSet[spec1].parentMass - specSet[spec2].parentMass;
        bestScore1 = -10;
        bestScore2 = -10;
        if (abs(pmDiff) > pmTol + 0.00001)
        {
          //	            FindMatchPeaksAll(specSet[spec1], specSet[spec2], 0, peakTol, idx1, idx2);
          //cerr<<"X"; cerr.flush();
          ScoreOverlap6(specSet[spec1],
                        specSet[spec2],
                        0,
                        peakTol,
                        idxMatched1_zero,
                        idxMatched2_zero);
          //cerr<<"X"; cerr.flush();

          for (shift = pmDiff - pmTol; shift <= pmDiff + pmTol + 0.00001;
              shift += shiftIncrement)
          {
            //	                FindMatchPeaksAll(specSet[spec1], specSet[spec2], shift, peakTol, idx1, idx2);
            //cerr<<"<"; cerr.flush();
            ScoreOverlap6(specSet[spec1],
                          specSet[spec2],
                          shift,
                          peakTol,
                          idxMatched1_other,
                          idxMatched2_other);
            //cerr<<">"; cerr.flush();

            idxUnionEnd = set_union(idxMatched1_zero.begin(),
                                    idxMatched1_zero.end(),
                                    idxMatched1_other.begin(),
                                    idxMatched1_other.end(),
                                    idxUnion.begin());
            for (idxUnionStart = idxUnion.begin(), score1 = 0, numMatchedPeaksU1 =
                0; idxUnionStart != idxUnionEnd;
                idxUnionStart++, numMatchedPeaksU1++)
            {
              score1 += specSet[spec1][(*idxUnionStart)][1];
            }

            numMatchedPeaks = numMatchedPeaksU1;

            idxUnionEnd = set_union(idxMatched2_zero.begin(),
                                    idxMatched2_zero.end(),
                                    idxMatched2_other.begin(),
                                    idxMatched2_other.end(),
                                    idxUnion.begin());

            for (idxUnionStart = idxUnion.begin(), score2 = 0, numMatchedPeaksU2 =
                0; idxUnionStart != idxUnionEnd;
                idxUnionStart++, numMatchedPeaksU2++)
            {
              score2 += specSet[spec2][(*idxUnionStart)][1];
            }

            numMatchedPeaks = min(numMatchedPeaks, numMatchedPeaksU2);

            if (((score1 + score2 > bestScore1 + bestScore2)
                || ((score1 + score2 == bestScore1 + bestScore2)
                    && abs(shift - pmDiff) < abs(bestShift - pmDiff))))
            {
              //DEBUG_MSG(*curCandidate << "  " << numMatchedPeaksU1 << "  " << numMatchedPeaksU2 << "  " << numMatchedPeaks);
              bestScore1 = score1;
              bestScore2 = score2;
              bestShift = shift;
              bestNumMatchedPeaks = numMatchedPeaks;
            }
          }
        }
        else
        {
          vector<int> idx1all(0);
          vector<int> idx2all(0);
          vector<int> idxMatched1(0);
          vector<int> idxMatched2(0);

          for (shift = -pmTol; shift <= pmTol + 0.00001; shift +=
              shiftIncrement)
          {
            //cerr<<"{"; cerr.flush();
            FindMatchPeaksAll(specSet[spec1],
                              specSet[spec2],
                              shift,
                              peakTol,
                              idx1all,
                              idx2all);
            //cerr<<"}"; cerr.flush();
            ScoreOverlap7(specSet[spec1],
                          idx1all,
                          specSet[spec2],
                          idx2all,
                          shift,
                          peakTol,
                          idxMatched1,
                          idxMatched2,
                          symmetryOffset);

            numMatchedPeaks = (short)idxMatched1.size();

            score1 = 0;
            score2 = 0;
            for (int i = 0; i < idxMatched1.size(); i++)
              score1 += specSet[spec1][idxMatched1[i]][1];
            for (int i = 0; i < idxMatched2.size(); i++)
              score2 += specSet[spec2][idxMatched2[i]][1];

            if ((score1 + score2 > bestScore1 + bestScore2)
                || (score1 + score2 == bestScore1 + bestScore2
                    && abs(shift) < abs(bestShift)))
            {
              bestScore1 = score1;
              bestScore2 = score2;
              bestShift = shift;
              bestNumMatchedPeaks = numMatchedPeaks;
              /*    cerr<<"Matching peaks with shift "<<bestShift<<", scores ("<<score1<<","<<score2<<"):\nSpectrum "<<spec1<<":\n";
               for(unsigned int i=0; i<idxMatched1.size(); i++) cerr<<specSet[spec1][idxMatched1[i]][0]<<"\t"<<specSet[spec1][idxMatched1[i]][1]<<endl;
               cerr<<"Spectrum "<<spec2<<":\n";
               for(unsigned int i=0; i<idxMatched2.size(); i++) cerr<<specSet[spec2][idxMatched2[i]][0]<<"\t"<<specSet[spec2][idxMatched2[i]][1]<<endl;
               //*/
            }
          }
        }
        //cerr<<"|"; cerr.flush();

        if (bestNumMatchedPeaks >= minNumMatchedPeaks)
        {
          curResult.spec1 = spec1;
          curResult.spec2 = spec2;
          //    curResult.score1 = (bestScore1>0)? bestScore1 : 0;
          //    curResult.score2 = (bestScore2>0)? bestScore2 : 0;
          curResult.score1 = bestScore1;
          curResult.score2 = bestScore2;
          curResult.shift1 = bestShift;

          means[spec1][0] = (means[spec1][0] * means[spec1][1]
              + curResult.score1) / (means[spec1][1] + 1);
          means[spec2][0] = (means[spec2][0] * means[spec2][1]
              + curResult.score2) / (means[spec2][1] + 1);
          means[spec1][1]++;
          means[spec2][1]++;
          float updateRatio = (means[spec1][1] - 1) / means[spec1][1]; // Update variance terms
          varTerms[spec1] = updateRatio * varTerms[spec1]
              + (curResult.score1 * curResult.score1) / means[spec1][1];
          updateRatio = (means[spec2][1] - 1) / means[spec2][1];
          varTerms[spec2] = updateRatio * varTerms[spec2]
              + (curResult.score2 * curResult.score2) / means[spec2][1];

          // Compute match ratios
          curRatios[0] = curResult.score1 / totalSpecScores[spec1];
          curRatios[1] = curResult.score2 / totalSpecScores[spec2];

          if (curRatios[0] >= minMatchRatio && curRatios[1] >= minMatchRatio)
          {
            //! /todo Implement an efficient push_back (or equivalent) for SpectrumPairSet
            results.push_back(curResult);
            ratios.push_back(curRatios);
          }
        }

        //cerr<<"|"; cerr.flush();
      } // while over all candidates for spec1
      DEBUG_MSG(" done. Current number of results in memory: " << results.size());
    } // for i over baseSpectraIdx
    DEBUG_MSG("Done with the alignment!");

  }

  //---------------------------------------------------------------------------
  void getPairAlignsPA(SpecSet & specSet,
                       float pmTol,
                       float peakTol,
                       SpectrumPairSet & results)
  {
    int pair, massOffset, baseShift, shift, shiftSym, middleShift,
        middleTimesTwo, spec1, spec2, i, j, k;
    int intPeakTolMin = -(int)round(peakTol * 10 * 2), intPeakTolMax =
        (int)round(10 * peakTol * 2);
    int intPmTolMin = -(int)round(pmTol * 10), intPmTolMax = (int)round(10
        * pmTol);
    float realShift1, realShift2, maxMass;
    vector<bool> shifts;
    //    vector<int> idx1,idx2;                        idx1.reserve(500);              idx2.reserve(500);
    vector<int> idxMatched1, idxMatched2;
    idxMatched1.reserve(500);
    idxMatched2.reserve(500);
    vector<int> idxMatchedSym1, idxMatchedSym2;
    idxMatchedSym1.reserve(500);
    idxMatchedSym2.reserve(500);
    vector<int> idxUnion;
    idxUnion.reserve(500);
    vector<int>::iterator idxUnionEnd, idxUnionStart;

    ERROR_MSG("getPairAlignsPA() in batch.cpp is not ready for high-accuracy spectra!");

    if (results.size() == 0)
    {
      // Fill results.spec1 and results.spec2 to compute ALL pairwise alignments
      results.resize((specSet.size() * (specSet.size() - 1)) / 2);
      i = 0;
      for (spec1 = 0; spec1 < specSet.size(); spec1++)
        for (spec2 = spec1 + 1; spec2 < specSet.size(); spec2++)
        {
          results[i].spec1 = spec1;
          results[i].spec2 = spec2;
          i++;
        }
    }

    // Pre-allocate vector to keep track of valid shifts - implicit 0.1 Da resolution for shifts
    for (i = 0, maxMass = 0; i <= specSet.size(); i++)
      if (specSet[i].parentMass > maxMass)
        maxMass = specSet[i].parentMass;
    massOffset = (int)ceil(maxMass * 10);
    shifts.resize(2 * massOffset + 1);

    for (pair = 0; pair <= results.size(); pair++)
    {
      spec1 = results[pair].spec1;
      spec2 = results[pair].spec2;

      // Step 1: Get the list of possible shifts
      for (shift = 0; shift < shifts.size(); shift++)
        shifts[shift] = 0;
      for (i = 0; i < specSet[spec1].size(); i++)
        for (j = 0; j < specSet[spec2].size(); j++)
        {
          baseShift = massOffset
              + (int)round(10 * (specSet[spec1][i][0] - specSet[spec2][j][0]));
          for (k = intPeakTolMin; k < intPeakTolMax; k++)
            shifts[baseShift + k] = 1;
        }

      // THERE IS A LOT OF RECOMPUTATION GOING ON HERE!
      // Every symmetric shift is being computed as many times as it is used - once should be enough if the memory requirements stay acceptable

      // Step 2: Compute the summed scores for all pairs of symmetric shifts (use pmTol) + choose the best pair
      float score1, bestScore1, score2, bestScore2, bestShift1, bestShift2;
      middleShift = massOffset
          + (int)floor(5
              * (specSet[spec1].parentMass - specSet[spec2].parentMass));
      middleTimesTwo = (int)round(10
          * (specSet[spec1].parentMass - specSet[spec2].parentMass));

      ERROR_MSG("Spec1: " << spec1 << ", parent mass = " << specSet[spec1].parentMass);
      ERROR_MSG("Spec2: " << spec2 << ", parent mass = " << specSet[spec2].parentMass);
      ERROR_MSG("middleShift = " << middleShift - massOffset << ", middleTimeTwo = " << middleTimesTwo);
      int first = 0;

      for (shift = 0; shift < middleShift; shift++)
      {
        if (shifts[shift] == 0)
          continue;
        realShift1 = (shift - massOffset) / 10;

        //            FindMatchPeaksAll(specSet[spec1], specSet[spec2], realShift1, peakTol, idx1, idx2);
        ScoreOverlap6(specSet[spec1],
                      specSet[spec2],
                      realShift1,
                      peakTol,
                      idxMatched1,
                      idxMatched2);

        shiftSym = massOffset + (middleTimesTwo - (shift - massOffset));
        bestScore1 = 0;
        bestScore2 = 0;
        bestShift1 = 0;
        bestShift2 = 0;

        if (first == 0)
        {
          ERROR_MSG("shift = " << massOffset - shift << ", realShift1 = " << realShift1);
        }

        for (k = intPmTolMin; k < intPmTolMax; k++)
        {
          if (shifts[shiftSym + k] == 0)
            continue;
          realShift2 = (shiftSym - massOffset) / 10;

          //                FindMatchPeaksAll(specSet[spec1], specSet[spec2], realShift2, peakTol, idx1, idx2);
          ScoreOverlap6(specSet[spec1],
                        specSet[spec2],
                        realShift2,
                        peakTol,
                        idxMatchedSym1,
                        idxMatchedSym2);

          // Union of matched peaks
          idxUnionEnd = set_union(idxMatched1.begin(),
                                  idxMatched1.end(),
                                  idxMatchedSym1.begin(),
                                  idxMatchedSym1.end(),
                                  idxUnion.begin());
          for (idxUnionStart = idxUnion.begin(), score1 = 0;
              idxUnionStart != idxUnionEnd; idxUnionStart++)
            score1 += specSet[spec1][(*idxUnionStart)][1];

          idxUnionEnd = set_union(idxMatched2.begin(),
                                  idxMatched2.end(),
                                  idxMatchedSym2.begin(),
                                  idxMatchedSym2.end(),
                                  idxUnion.begin());
          for (idxUnionStart = idxUnion.begin(), score2 = 0;
              idxUnionStart != idxUnionEnd; idxUnionStart++)
            score2 += specSet[spec2][(*idxUnionStart)][1];

          if (first == 0)
          {
            ERROR_MSG("-- shiftSym = " << massOffset - shiftSym << ", realShift2 = " << realShift2 << ", score1 = " << score1 << ", score2 = " << score2);
          }

          if ((score1 + score2 > bestScore1 + bestScore2)
              || (score1 + score2 == bestScore1 + bestScore2
                  && abs(shift + shiftSym - middleTimesTwo - 2 * massOffset)
                      < abs(bestShift1 + bestShift2 - middleTimesTwo
                          - 2 * massOffset)))
          {
            bestScore1 = score1;
            bestScore2 = score2;
            bestShift1 = realShift1;
            bestShift2 = realShift2;
          }
        }
        if (first == 0)
        {
          ERROR_MSG("bestShift1 = " << bestShift1 << ", bestShift2 = " << bestShift2 << ", bestScore1 = " << bestScore1 << ", bestScore2 = " << bestScore2);
          first = 1;
        }
      }

      ERROR_MSG("bestShift1 = " << bestShift1 << ", bestShift2 = " << bestShift2 << ", bestScore1 = " << bestScore1 << ", bestScore2 = " << bestScore2);

      results[pair].score1 = bestScore1;
      results[pair].score2 = bestScore2;
      results[pair].shift1 = bestShift1;
      results[pair].shift2 = bestShift2;
    }
  }

  //---------------------------------------------------------------------------
  //
  //  getPairAlignsPA2 - Like getPairAlignsPA but fixes some problems and adds some extensions
  //    - Each shift is only computed once
  //    - Shift pairs have scores precomputed by computeShifts and are passed to ScoreOverlap6 only if
  //        their potential DP score is higher than the best achieved DP score for the pair
  //    - Only returns pairs that meet the minRatio/minPeakAreaOvlp criteria
  //
  //  alignStats - Alignment statistics per returned pair (each is minimum between the 2 aligned spectra):
  //                  pos 0: num matched peaks (before DP/merging)
  //                  pos 1: matched intensity (before DP/merging)
  //                  pos 2: num matched peaks (after DP/merging)
  //                  pos 3: matched intensity (after DP/merging)
  //  specStats  - Alignment statistics per aligned spectrum: (intensity not considered because intensity distributions should vary per spectrum)
  //                  pos 0: num alignments
  //                  pos 1/2: mean/variance for num matched peaks (before DP/merging)
  //                  pos 3/4: mean/variance num matched peaks (after DP/merging)
  //                  pos 5/6: mean/variance for matched intensity (before DP/merging)
  //                  pos 7/8: mean/variance matched intensity (after DP/merging)
  //---------------------------------------------------------------------------
  void getPairAlignsPA2(SpecSet & specSet,
                        unsigned int startIdx,
                        unsigned int endIdx,
                        float peakTol,
                        float pmTol,
                        float minRatio,
                        float minPeakAreaOvlp,
                        short minNumMatchedPeaks,
                        AAJumps & allowedJumps,
                        float minAbsShift,
                        SpectrumPairSet & results,
                        vector<TwoValues<float> > & ratios,
                        list<TwoValues<int> > & numMatchedPeaks,
                        vector<TwoValues<float> > & means,
                        vector<float> & varTerms,
                        list<vector<float> > & alignStats,
                        vector<vector<float> > & specStats,
                        float resolution,
                        float maxShift)
  {
    unsigned int spec1, spec2;

    // Variables and initializations
    TwoValues<int> res;
    list<float> shiftScores; // Maximum scores for eligible pairs of shifts
    list<TwoValues<unsigned int> > shiftPairs; // List of eligible pairs of shifts (in the same order as shiftScores)
    //  vector<vector<TwoValues<int> > > shiftMatchedPeaks;  // Lists of matched peaks per shift between the two spectra (indexed using
    vector<list<TwoValues<int> > > shiftMatchedPeaks; // Lists of matched peaks per shift between the two spectra (indexed using
    //  a shift offset returned by computeShifts)
    vector<TwoValues<float> > shiftDPscores; // Match score per shift after applying dynamic programming
    vector<float> shiftDPpenalties; // Peak mass mismatch penalties (in DP)
    vector<TwoValues<vector<int> > > shiftDPmatchedPeaks; // Peaks matched by the DP algorithm
    vector<TwoValues<float> > minMaxMatchScores(specSet.size()); // Minimum acceptable / Maximum possible match scores per spectrum in the dataset
    SpectrumPair curResult;
    vector<float> curAlignStats(4);

    // Initialize stats vectors
    means.resize(specSet.size());
    varTerms.resize(specSet.size());
    specStats.resize(specSet.size());
    for (unsigned int i = 0; i < means.size(); i++)
    {
      means[i].set(0, 0);
      varTerms[i] = 0;
      specStats[i].resize(9);
      for (unsigned int j = 0; j < 9; j++)
        specStats[i][j] = 0;
    }

    float maxParentMass = 0;
    unsigned int maxNumPeaks = 0;
    int intPeakTol = (int)round(peakTol / resolution), intPMTol =
        (int)round(pmTol / resolution);
    for (unsigned int i = 0; i < specSet.size(); i++)
    {
      if (specSet[i].parentMass > maxParentMass)
        maxParentMass = specSet[i].parentMass;
      if (specSet[i].size() > maxNumPeaks)
        maxNumPeaks = specSet[i].size();
      minMaxMatchScores[i].set(0, 0);
      for (unsigned int j = 0; j < specSet[i].size(); j++)
        minMaxMatchScores[i][1] += specSet[i][j][1];
      minMaxMatchScores[i][0] = minRatio * minMaxMatchScores[i][1];
    }
    //  shiftMatchedPeaks.resize(1+(int)round(2*(min(maxParentMass,maxShift)+8*peakTol)/resolution)); //for(unsigned int i=0;i<shiftMatchedPeaks.size();i++) shiftMatchedPeaks[i].reserve(10);
    shiftMatchedPeaks.resize(1
        + (int)ceil(2
            * (min(maxParentMass, maxShift) + peakTol + pmTol)
            / resolution)); //for(unsigned int i=0;i<shiftMatchedPeaks.size();i++) shiftMatchedPeaks[i].reserve(10);
    shiftDPscores.reserve(shiftMatchedPeaks.capacity());
    shiftDPpenalties.reserve(shiftMatchedPeaks.capacity());
    shiftDPmatchedPeaks.reserve(shiftMatchedPeaks.capacity());
    vector<int> idx1, idx2;
    idx1.reserve(maxNumPeaks * (2 * intPeakTol + 1) + 1);
    idx2.reserve(idx1.capacity()); // Temporary variables used as input to ScoreOverlap6
    vector<int> idxMerged1, idxMerged2;
    idxMerged1.reserve(idx1.capacity());
    idxMerged2.reserve(idx1.capacity()); // Temporary variables used to merge symmetric results of ScoreOverlap6

    // Compute the pairwise alignments
    time_t curTime = time(0);
    double totTime = 0.0;
    int numResults = 0;
    float bestShiftScore;
    TwoValues<int> bestShiftPair;
    TwoValues<int> bestNumMatchedPeaks;
    TwoValues<float> bestShiftScores, bestCandidateScores;
    list<int> shiftsToCompute;
    for (spec1 = startIdx; spec1 <= endIdx; spec1++)
    {
      DEBUG_MSG("Processing spectrum " << spec1 << "... ");
      if (specSet[spec1].size() == 0 or specSet[spec1].parentMass < 400) continue;

      for (spec2 = spec1 + 1; spec2 < specSet.size(); spec2++)
      {
    	if (specSet[spec2].size() == 0 or specSet[spec2].parentMass < 400) continue;

    	float minPM = specSet[spec1].parentMass, maxPM = specSet[spec2].parentMass;
		if( maxPM < minPM ){
			minPM = specSet[spec2].parentMass;
			 maxPM = specSet[spec1].parentMass;
		}
		if ( minPM < minPeakAreaOvlp*maxPM ) continue;

        for (unsigned int i = 0; i < curAlignStats.size(); i++)
          curAlignStats[i] = 0;

        // Step 1: Get the list of possible shifts
        res = computeShifts(specSet[spec1],
                            specSet[spec2],
                            peakTol,
                            pmTol,
                            minPeakAreaOvlp*maxPM,
                            minMaxMatchScores[spec1][0],
                            minMaxMatchScores[spec2][0],
                            minNumMatchedPeaks,
                            allowedJumps,
                            resolution,
                            shiftScores,
                            shiftPairs,
                            shiftMatchedPeaks,
                            bestCandidateScores,
                            bestNumMatchedPeaks,
                            minAbsShift);
        if (res[1] == 0)
          continue; // No eligible shifts found; could happen when maximum overlap is less than minPeakAreaOvlp

        // Matched intensity before DP/merging
        auxUpdateMeanVariance(bestCandidateScores[0],
                              specStats[spec1][0],
                              specStats[spec1][5],
                              specStats[spec1][6]);
        auxUpdateMeanVariance(bestCandidateScores[1],
                              specStats[spec2][0],
                              specStats[spec2][5],
                              specStats[spec2][6]);
        // Num matched peaks before DP/merging
        auxUpdateMeanVariance(bestNumMatchedPeaks[0],
                              specStats[spec1][0],
                              specStats[spec1][1],
                              specStats[spec1][2]);
        auxUpdateMeanVariance(bestNumMatchedPeaks[1],
                              specStats[spec2][0],
                              specStats[spec2][1],
                              specStats[spec2][2]);
        // Pre-DP alignment stats, in case this alignment gets reported
        curAlignStats[0] = min(bestNumMatchedPeaks[0], bestNumMatchedPeaks[1]);
        curAlignStats[1] = min(bestCandidateScores[0]
                                   / minMaxMatchScores[spec1][1],
                               bestCandidateScores[1]
                                   / minMaxMatchScores[spec2][1]);

        /*          means[spec1][0]=(means[spec1][0]*means[spec1][1]+bestCandidateScores[0])/(means[spec1][1]+1);  // Include accepted pairs in means
         means[spec2][0]=(means[spec2][0]*means[spec2][1]+bestCandidateScores[1])/(means[spec2][1]+1);
         means[spec1][1]++;    means[spec2][1]++;
         float updateRatio = (means[spec1][1]-1)/means[spec1][1];  // Update variance terms
         varTerms[spec1] = updateRatio*varTerms[spec1]+(bestCandidateScores[0]*bestCandidateScores[0])/means[spec1][1];
         updateRatio = (means[spec2][1]-1)/means[spec2][1];
         varTerms[spec2] = updateRatio*varTerms[spec2]+(bestCandidateScores[1]*bestCandidateScores[1])/means[spec2][1];
         */

        //cerr<<"shiftDPscores resized to "<<shiftMatchedPeaks.size()<<endl;
        shiftDPscores.resize(shiftMatchedPeaks.size());
        for (unsigned int i = 0; i < shiftDPscores.size(); i++)
          shiftDPscores[i].set(-1.0, -1.0);
        shiftDPpenalties.resize(shiftMatchedPeaks.size());
        for (unsigned int i = 0; i < shiftDPpenalties.size(); i++)
          shiftDPpenalties[i] = 0.0;
        shiftDPmatchedPeaks.resize(shiftMatchedPeaks.size());
        for (unsigned int i = 0; i < shiftDPmatchedPeaks.size(); i++)
        {
          shiftDPmatchedPeaks[i][0].resize(0);
          shiftDPmatchedPeaks[i][1].resize(0);
        }

        bestShiftScore = 0;
        bestShiftPair.set(-1, -1);
        bestShiftScores.set(0, 0);
        bestCandidateScores.set(0, 0);
        bestNumMatchedPeaks.set(0, 0); // Reuse to store best DP match score / num matched peaks (used to update means/vars)
        list<float>::iterator nextShiftScore = shiftScores.begin();
        list<TwoValues<unsigned int> >::iterator nextShiftPair =
            shiftPairs.begin();
        while (nextShiftScore != shiftScores.end()
            and *nextShiftScore > bestShiftScore)
        {
          int shiftIndex = (*nextShiftPair)[0], shiftSym = (*nextShiftPair)[1];
          //cerr<<"Shift pair "<<shiftIndex<<"/"<<shiftSym<<endl;
          if (shiftDPscores[shiftIndex][0] + shiftDPscores[shiftIndex][1]
              < 0.0001 and shiftMatchedPeaks[shiftIndex].size() > 0)
            shiftsToCompute.push_back(shiftIndex);
          for (int tolIdx = max(0, shiftSym - intPMTol);
              tolIdx <= shiftSym + intPMTol
                  and tolIdx < (int)shiftDPscores.size(); tolIdx++)
            if (shiftDPscores[tolIdx][0] + shiftDPscores[tolIdx][1] < 0.0001
                and shiftMatchedPeaks[tolIdx].size() > 0)
              shiftsToCompute.push_back(tolIdx);

          for (list<int>::iterator curShift = shiftsToCompute.begin();
              curShift != shiftsToCompute.end(); curShift++)
          {
            float shiftMass = ((*curShift) - res[0]) * resolution;

            idx1.resize(shiftMatchedPeaks[*curShift].size());
            idx2.resize(shiftMatchedPeaks[*curShift].size());
            int pivot = 0;
            for (list<TwoValues<int> >::iterator iterPeaks =
                shiftMatchedPeaks[*curShift].begin();
                iterPeaks != shiftMatchedPeaks[*curShift].end(); iterPeaks++)
            {
              idx1[pivot] = (*iterPeaks)[0];
              idx2[pivot] = (*iterPeaks)[1];
              pivot++;
            }

            /*cerr<<"shift = "<<*curShift<<", idx1=(";
             for(unsigned int dbgIdx=0; dbgIdx<idx1.size(); dbgIdx++) cerr<<idx1[dbgIdx]<<",";
             cerr<<"), idx2=(";
             for(unsigned int dbgIdx=0; dbgIdx<idx2.size(); dbgIdx++) cerr<<idx2[dbgIdx]<<",";
             cerr<<")\n";*/
            ScoreOverlap6mp(specSet[spec1],
                            idx1,
                            specSet[spec2],
                            idx2,
                            shiftMass,
                            peakTol,
                            shiftDPmatchedPeaks[*curShift][0],
                            shiftDPmatchedPeaks[*curShift][1],
                            AAJumps::minAAmass,
                            &shiftDPpenalties[*curShift]);
            //          ScoreOverlap6(specSet[spec1], idx1, specSet[spec2], idx2, shiftMass, peakTol, shiftDPmatchedPeaks[*curShift][0], shiftDPmatchedPeaks[*curShift][1], false);

            /*cerr<<"shift = "<<*curShift<<", matched1=(";
             for(unsigned int dbgIdx=0; dbgIdx<shiftDPmatchedPeaks[*curShift][0].size(); dbgIdx++) cerr<<shiftDPmatchedPeaks[*curShift][0][dbgIdx]<<",";
             cerr<<"), matched2=(";
             for(unsigned int dbgIdx=0; dbgIdx<shiftDPmatchedPeaks[*curShift][1].size(); dbgIdx++) cerr<<shiftDPmatchedPeaks[*curShift][1][dbgIdx]<<",";
             cerr<<")\n";*/
            shiftDPscores[*curShift][0] = 0;
            for (unsigned int i = 0;
                i < shiftDPmatchedPeaks[*curShift][0].size(); i++)
              shiftDPscores[*curShift][0] +=
                  specSet[spec1][shiftDPmatchedPeaks[*curShift][0][i]][1];
            shiftDPscores[*curShift][1] = 0;
            for (unsigned int i = 0;
                i < shiftDPmatchedPeaks[*curShift][1].size(); i++)
              shiftDPscores[*curShift][1] +=
                  specSet[spec2][shiftDPmatchedPeaks[*curShift][1][i]][1];
          }
          shiftsToCompute.clear();

          //cerr<<"--- repeat shift pair "<<shiftIndex<<"/"<<shiftSym<<endl;
          for (int tolIdx = max(0, shiftSym - intPMTol);
              tolIdx <= shiftSym + intPMTol
                  and tolIdx < (int)shiftDPmatchedPeaks.size(); tolIdx++)
          {
            //          if(shiftDPscores[shiftIndex][0]+shiftDPscores[tolIdx][0]+shiftDPscores[shiftIndex][1]+shiftDPscores[tolIdx][1]>bestCandidateScores[0]+bestCandidateScores[1])
            //            bestCandidateScores.set(shiftDPscores[shiftIndex][0]+shiftDPscores[tolIdx][0],shiftDPscores[shiftIndex][1]+shiftDPscores[tolIdx][1]);
            //          if(shiftDPscores[shiftIndex][0]+shiftDPscores[tolIdx][0]<minMaxMatchScores[spec1][0] or  // Sole purpose is to avoid unnecessary computations
            //             shiftDPscores[shiftIndex][1]+shiftDPscores[tolIdx][1]<minMaxMatchScores[spec2][0]) continue;
            if (shiftMatchedPeaks[tolIdx].size() == 0)
              continue;

            // merge lists of matched peaks - spec1
            mergeVectors(shiftDPmatchedPeaks[shiftIndex][0],
                         shiftDPmatchedPeaks[tolIdx][0],
                         idxMerged1);
            mergeVectors(shiftDPmatchedPeaks[shiftIndex][1],
                         shiftDPmatchedPeaks[tolIdx][1],
                         idxMerged2);

            /*cerr<<"shifts = "<<shiftIndex<<"/"<<tolIdx<<":\n  idxMerged1: (";
             for(unsigned int dbgIdx=0; dbgIdx<shiftDPmatchedPeaks[shiftIndex][0].size(); dbgIdx++) cerr<<shiftDPmatchedPeaks[shiftIndex][0][dbgIdx]<<",";
             cerr<<") U (";
             for(unsigned int dbgIdx=0; dbgIdx<shiftDPmatchedPeaks[tolIdx][0].size(); dbgIdx++) cerr<<shiftDPmatchedPeaks[tolIdx][0][dbgIdx]<<",";
             cerr<<") -> (";
             for(unsigned int dbgIdx=0; dbgIdx<idxMerged1.size(); dbgIdx++) cerr<<idxMerged1[dbgIdx]<<",";
             cerr<<")\n  idxMerged2: (";
             for(unsigned int dbgIdx=0; dbgIdx<shiftDPmatchedPeaks[shiftIndex][1].size(); dbgIdx++) cerr<<shiftDPmatchedPeaks[shiftIndex][1][dbgIdx]<<",";
             cerr<<") U (";
             for(unsigned int dbgIdx=0; dbgIdx<shiftDPmatchedPeaks[tolIdx][1].size(); dbgIdx++) cerr<<shiftDPmatchedPeaks[tolIdx][1][dbgIdx]<<",";
             cerr<<") -> (";
             for(unsigned int dbgIdx=0; dbgIdx<idxMerged2.size(); dbgIdx++) cerr<<idxMerged2[dbgIdx]<<",";
             cerr<<")\n";*/
            float score1 = 0;
            for (unsigned int i = 0; i < idxMerged1.size(); i++)
              score1 += specSet[spec1][idxMerged1[i]][1];
            float score2 = 0;
            for (unsigned int i = 0; i < idxMerged2.size(); i++)
              score2 += specSet[spec2][idxMerged2[i]][1];
            // Moved (and changed) from above on 2006/09/11
            if (bestShiftScore < 0.0001)
            { // These statements ensure that we get values to update means/variances
              if (score1 + score2
                  > bestCandidateScores[0] + bestCandidateScores[1])
                bestCandidateScores.set(score1, score2); // Used to update means/vars even if this pair is not returned
              if (idxMerged1.size() + idxMerged2.size()
                  > bestNumMatchedPeaks[0] + bestNumMatchedPeaks[1])
                bestNumMatchedPeaks.set(idxMerged1.size(), idxMerged2.size());
            }

            if (score1 < minMaxMatchScores[spec1][0]
                or score2 < minMaxMatchScores[spec2][0])
              continue;
            if (score1 + score2 - shiftDPpenalties[shiftIndex]
                - shiftDPpenalties[tolIdx] > bestShiftScore)
            {
              bestShiftScore = score1 + score2 - shiftDPpenalties[shiftIndex]
                  - shiftDPpenalties[tolIdx];
              bestShiftPair.set(shiftIndex, tolIdx);
              bestShiftScores.set(score1, score2);
              bestNumMatchedPeaks.set(idxMerged1.size(), idxMerged2.size());
            }
          }
          nextShiftScore++;
          nextShiftPair++;
        }

        // Build results for best match
        TwoValues<float> curRatios;
        if (bestShiftScore > 0)
        {
          curResult.spec1 = spec1;
          curResult.spec2 = spec2;
          curResult.score1 = bestShiftScores[0];
          curResult.score2 = bestShiftScores[1];
          curResult.shift1 = (bestShiftPair[0] - res[0])
              * resolution;
          curResult.shift2 = (bestShiftPair[1] - res[0])
              * resolution;

          //! /todo Implement an efficient push_back (or equivalent) for SpectrumPairSet
          results.push_back(curResult);

          curRatios.set(curResult.score1 / minMaxMatchScores[spec1][1],
                        curResult.score2 / minMaxMatchScores[spec2][1]);
          ratios.push_back(curRatios);
          numMatchedPeaks.push_back(bestNumMatchedPeaks);
          bestCandidateScores = bestShiftScores;

          // Post-DP alignment stats since this alignment is reported
          curAlignStats[2] = min(bestNumMatchedPeaks[0],
                                 bestNumMatchedPeaks[1]);
          curAlignStats[3] = min(bestCandidateScores[0]
                                     / minMaxMatchScores[spec1][1],
                                 bestCandidateScores[1]
                                     / minMaxMatchScores[spec2][1]);
          alignStats.push_back(curAlignStats);
        }

        // Regular mean/vars output to text files
        auxUpdateMeanVariance(bestCandidateScores[0],
                              means[spec1][1],
                              means[spec1][0],
                              varTerms[spec1]);
        auxUpdateMeanVariance(bestCandidateScores[1],
                              means[spec2][1],
                              means[spec2][0],
                              varTerms[spec2]);
        // Matched intensity after DP/merging
        auxUpdateMeanVariance(bestCandidateScores[0],
                              specStats[spec1][0],
                              specStats[spec1][7],
                              specStats[spec1][8]);
        auxUpdateMeanVariance(bestCandidateScores[1],
                              specStats[spec2][0],
                              specStats[spec2][7],
                              specStats[spec2][8]);
        // Num matched peaks after DP/merging
        //cerr<<"Num matched peaks: "<<bestNumMatchedPeaks[0]<<"/"<<bestNumMatchedPeaks[1]<<", sampleSize = "<<specStats[spec1][0]<<"/"<<specStats[spec2][0];
        auxUpdateMeanVariance(bestNumMatchedPeaks[0],
                              specStats[spec1][0],
                              specStats[spec1][3],
                              specStats[spec1][4]);
        auxUpdateMeanVariance(bestNumMatchedPeaks[1],
                              specStats[spec2][0],
                              specStats[spec2][3],
                              specStats[spec2][4]);
        //cerr<<", updated means = "<<specStats[spec1][3]<<"/"<<specStats[spec2][3]<<endl;
        specStats[spec1][0]++;
        means[spec1][1]++;
        specStats[spec2][0]++;
        means[spec2][1]++;

        /*      means[spec1][0]=(means[spec1][0]*means[spec1][1]+bestCandidateScores[0])/(means[spec1][1]+1);  // Include accepted pairs in means
         means[spec2][0]=(means[spec2][0]*means[spec2][1]+bestCandidateScores[1])/(means[spec2][1]+1);
         means[spec1][1]++;    means[spec2][1]++;
         float updateRatio = (means[spec1][1]-1)/means[spec1][1];  // Update variance terms
         varTerms[spec1] = updateRatio*varTerms[spec1]+(bestCandidateScores[0]*bestCandidateScores[0])/means[spec1][1];
         updateRatio = (means[spec2][1]-1)/means[spec2][1];
         varTerms[spec2] = updateRatio*varTerms[spec2]+(bestCandidateScores[1]*bestCandidateScores[1])/means[spec2][1];
         */}
      time_t newTime = time(0);
      double ellapsed = difftime(newTime, curTime);
      totTime += ellapsed;
      curTime = newTime;
      int N = endIdx, CUR = spec1 - startIdx + 1;
      DEBUG_MSG(ellapsed << " secs, number of new pair aligns = " << results.size() - numResults << " / " << results.size() << ". Average time = " << totTime / CUR << ", ETA = " << (N - CUR - 1) * (N - CUR - 2) * totTime / (CUR * (N - CUR) + (CUR - 1) * CUR / 2.0));

      numResults = results.size();
    }
  }

//---------------------------------------------------------------------------
//
// specSet        - Set of spectra to align
// baseSpectraIdx - Indices of the spectra to align to others on this run. Note that spectra only align to spectra with higher indices.
// aaDiff         - Difference between matched spectrallowedJumpsa's parent masses must be at most the sum of aaDiffs amino acid masses
//                    or -1 for computation of all shifts between minShift and maxShift
// minShift, maxShift - Min and max parent mass differences between spectra (to compute ASP shifts)
// pmTol          - Parent mass error tolerance
// peakTol        - Peak mass error tolerance
// minMatchRatio  - minimum ratio of matched score / max score (in both spectra) to retain a match
// results        - Variable to hold results
// resolution     - Resolution for parent mass histogram (to determine which spectra to match)
//
//---------------------------------------------------------------------------
void getPairCosines(const SpecSet & specSet,
                    vector<int> & baseSpectraIdx,
                    float maxPMdiff,
                    float pmTol,
                    float peakTol,
                    float minCosine,
                    float minMatchedIntensity,
                    unsigned int minNumMatchedPeaks,
                    SpectrumPairSet & results,
                    vector<TwoValues<float> > & ratios,
                    vector<TwoValues<float> > & means,
                    vector<float> & varTerms,
                    ExecTagFilterPairs *tag_filter,
                    float resolution,
                    float symmetryOffset,
                    bool projectedCosine)
{
    DEBUG_VAR(minCosine);
    DEBUG_VAR(minMatchedIntensity);
    DEBUG_VAR(minNumMatchedPeaks);
    DEBUG_VAR(maxPMdiff);

    int j, massIdx;
    float pMass1;
    float shiftIncrement = 2 * resolution; // Step size used to increment shifts from -pmTol to pmTol
    list<int>::iterator curPair;
    float histResolution = 0.1; // Resolution of the parent mass histogram - only used for binning, independent of overall peak resolution settings
    vector<list<int> > hist;
    specSet.getMassesHistogram(hist, histResolution);
    TwoValues<float> curRatios;
    SpectrumPair curResult;

    unsigned int histCount = 0;
    for (unsigned int i = 0; i < hist.size(); i++)
        histCount += hist[i].size();
    DEBUG_MSG("Parent masses histogram has " << hist.size()
        << " bins with a total of " << histCount << " entries");

    // Initialize means vector
    vector<float> spectrumIntensity; // total intensity per spectrum for later computation of match ratios
    spectrumIntensity.resize(specSet.size());
    means.resize(specSet.size());
    varTerms.resize(specSet.size());
    for (unsigned int i = 0; i < specSet.size(); i++){
        means[i].set(0, 0);
        varTerms[i] = 0;
        spectrumIntensity[i] = 0;
        for (unsigned int p = 0; p < specSet[i].size(); p++)
        spectrumIntensity[i] += specSet[i][p][1];
        //ajustedMinMatchIntensity[i] *= minMatchedIntensity;
    }

    // Iterate through pairs
    //    vector<int> idx1,idx2;                             idx1.reserve(1500);               idx2.reserve(1500);
    vector<int> idxMatched1_zero, idxMatched1_other;
    idxMatched1_zero.reserve(1500);
    idxMatched1_other.reserve(1500);
    vector<int> idxMatched2_zero, idxMatched2_other;
    idxMatched2_zero.reserve(1500);
    idxMatched2_other.reserve(1500);
    vector<int> idxUnion;
    idxUnion.reserve(1500);
    vector<int>::iterator idxUnionEnd, idxUnionStart;

    int spec1, spec2, maxPMdiffInt = (int)ceil(maxPMdiff / histResolution);
    float pmDiff, cosine, score1, score2, shift, bestScore1, bestScore2,
        bestShift;
    vector<GPCAux> peakMatches;
    vector<char> peakUsed1, peakUsed2;
    list<int> candidates;
    list<int>::iterator curCandidate;
    
    DEBUG_VAR(baseSpectraIdx.size());
    DEBUG_VAR(specSet.size());
    
    for (unsigned int i = 0; i < baseSpectraIdx.size(); i++)
    {
        DEBUG_MSG( "Starting with spectrum " << baseSpectraIdx[i] << "...");

        // Find candidate spectra to compare to
        spec1 = baseSpectraIdx[i];
        peakUsed1.resize(specSet[spec1].size());
        pMass1 = specSet[spec1].parentMass;
        int massInt = (int)round(pMass1 / histResolution);
        candidates.clear(); // Make sure candidates is empty from previous run
        if (specSet[spec1].size() < 5)
        {
            DEBUG_MSG("(spectrum has less than 5 peaks)");
            continue;
        }
        // Find candidates within maxPMshift
        int upperLimit = max(0, (int)ceil((pMass1 + maxPMdiff) / histResolution)); // Look at interval around current parent mass
        if (upperLimit >= hist.size())
            upperLimit = hist.size();
        for (massIdx = (int)floor((pMass1 - min(pMass1, maxPMdiff))
            / histResolution); massIdx < upperLimit; massIdx++)
        {
            for (curPair = hist[massIdx].begin(); curPair != hist[massIdx].end();
                curPair++)
            {
                if (*curPair > spec1
                    and fabs(pMass1 - specSet[*curPair].parentMass)
                        <= maxPMdiff + 2 * pmTol)
                candidates.push_back(*curPair);
            }
        }
        
        /*
        for(int i = spec1 + 1; i < specSet.size(); i++){
            candidates.push_back(i);
        }*/
        
        
        
        candidates.sort();
        candidates.unique(); // The same spectrum may occur in multiple bins - make sure the code below only aligns once to each eligible spectrum
        DEBUG_MSG(" num. candidates is " << candidates.size());

        if (tag_filter->hasTags()) {

			list<int> subsetOfCandidates;
			tag_filter->getPairs(spec1, candidates, subsetOfCandidates);
			candidates.swap(subsetOfCandidates);
			DEBUG_MSG(" (" << candidates.size() << " after tag filtering)");
		}
		DEBUG_MSG("..."); //*/

        curCandidate = candidates.begin();
        while (curCandidate != candidates.end()){
            
            spec2 = *curCandidate;
            curCandidate++;
            if (spec1 >= specSet.size() || spec2 >= specSet.size())
            {
                ERROR_MSG("ERROR: Alignment pair (" << spec1 << "," << spec2 << ") is invalid - spectrum set set has only " << specSet.size() << " spectra!\n");
                exit(-1);
            }

            pmDiff = specSet[spec1].parentMass - specSet[spec2].parentMass;
            unsigned int numMatchedPeaks = 0;
            bool projectedCosine2 = projectedCosine;
            float cosine = specSet[spec1].scoreMatch(specSet[spec2],
                            peakTol,
                            numMatchedPeaks,
                            score1,
                            score2,
                            false,
                            projectedCosine2,
                            false);

            // Store results
            curRatios[0] = score1 / spectrumIntensity[spec1];
            curRatios[1] = score2 / spectrumIntensity[spec2];
            curResult.spec1 = spec1;
            curResult.spec2 = spec2;
            curResult.score1 = cosine;
            curResult.score2 = min(curRatios[0], curRatios[1]);
            curResult.shift1 = pmDiff;

            // Update means/variances
            means[spec1][0] = (means[spec1][0] * means[spec1][1] + cosine)
                / (means[spec1][1] + 1);
            means[spec2][0] = (means[spec2][0] * means[spec2][1] + cosine)
                / (means[spec2][1] + 1);
            means[spec1][1]++;
            means[spec2][1]++;
            float updateRatio = (means[spec1][1] - 1) / means[spec1][1]; // Update variance terms
            varTerms[spec1] = updateRatio * varTerms[spec1]
                + (cosine * cosine) / means[spec1][1];
            updateRatio = (means[spec2][1] - 1) / means[spec2][1];
            varTerms[spec2] = updateRatio * varTerms[spec2]
                + (cosine * cosine) / means[spec2][1];

            if (cosine < minCosine or numMatchedPeaks < minNumMatchedPeaks)
                continue;

            if (curRatios[0] >= minMatchedIntensity
                && curRatios[1] >= minMatchedIntensity)
            {
                //! /todo Implement an efficient push_back (or equivalent) for SpectrumPairSet
                results.push_back(curResult);
                ratios.push_back(curRatios);
            }

        } // while over all candidates for spec1
        DEBUG_MSG(" done. Current number of results in memory: " << results.size());
    } // for i over baseSpectraIdx
    DEBUG_MSG("Done with the alignment!");
}

  void getPairAlignGFPValues(SpecSet & specSet,
                             SpectrumPairSet & spectraPairs,
                             vector<TwoValues<double> > & pvalues,
                             float pmTol,
                             float peakTol,
                             bool specTypeMSMS)
  {
    static float shiftIncrement = 0.1; // Step size used to increment shifts from -pmTol to pmTol

    //	  bool specTypeMSMS = true;//ip.getValueBool("SPEC_TYPE_MSMS", false);
    float ionOffset = specTypeMSMS ? AAJumps::massHion : 0;
    for (unsigned int i = 0; i < specSet.size(); i++)
    {
      specSet[i].addZPMpeaks(peakTol, ionOffset, true);
    }

    for (unsigned int i = 0; i < specSet.size(); i++)
    { //convert intensity into integer
      for (unsigned int j = 0; j < specSet[i].size(); j++)
      {
        specSet[i][j][1] = ceil(specSet[i][j][1]);
        // specSet[i][j][1] = floor(specSet[i][j][1]);
        // specSet[i][j][1] = round(specSet[i][j][1]);//
      }
      if (specSet[i].size() < 4)
        continue;
      specSet[i][0][1] = 0;
      specSet[i][specSet[i].size() - 1][1] = 0;
      specSet[i].removePeak(1);
      specSet[i].removePeak(specSet[i].size() - 2);
    }  			 //*/

    vector<float> totalSpecScores;
    totalSpecScores.resize(specSet.size());
    for (unsigned int i = 0; i < specSet.size(); i++)
    {
      totalSpecScores[i] = 0;
      for (unsigned int p = 0; p < specSet[i].size(); p++)
        totalSpecScores[i] += specSet[i][p][1];
    }

    vector<vector<int> > candidates(specSet.size());
    unsigned int numEntries = spectraPairs.size();
    for (unsigned int i = 0; i < numEntries; i++)
    {
      candidates[spectraPairs[i].spec1].push_back(spectraPairs[i].spec2);
    }

    vector<int> idxMatched1_zero, idxMatched1_other;
    idxMatched1_zero.reserve(1500);
    idxMatched1_other.reserve(1500);
    vector<int> idxMatched2_zero, idxMatched2_other;
    idxMatched2_zero.reserve(1500);
    idxMatched2_other.reserve(1500);
    vector<int> idxUnion;
    idxUnion.reserve(1500);
    vector<int>::iterator idxUnionEnd, idxUnionStart;

    for (int spec1 = 0; spec1 < candidates.size(); spec1++)
    {
      if (candidates[spec1].size() == 0)
        continue;
      float pmDiff, score1, score2, shift, bestScore1, bestScore2, bestShift;
      short numMatchedPeaks, numMatchedPeaksU1, numMatchedPeaksU2;

      double ***spec1AlignGFProbabiltyDensity =
          getALGFProbDensity(specSet[spec1], totalSpecScores[spec1], peakTol);
      for (int ipair = 0; ipair < candidates[spec1].size(); ipair++)
      {
        int spec2 = candidates[spec1][ipair];
        pmDiff = specSet[spec1].parentMass - specSet[spec2].parentMass;
        bestScore1 = 0;
        bestScore2 = 0;
        int bestNumMatchPeak = 0;
        int bestFwScore1 = 0, bestFwScore2 = 0, bestRvScore1 = 0, bestRvScore2 =
            0;
        ScoreOverlap6(specSet[spec1],
                      specSet[spec2],
                      0,
                      peakTol,
                      idxMatched1_zero,
                      idxMatched2_zero);

        int fwScore1 = 0, fwScore2 = 0;
        for (int i = 0; i < idxMatched1_zero.size(); i++)
          fwScore1 += specSet[spec1][idxMatched1_zero[i]][1];
        for (int i = 0; i < idxMatched2_zero.size(); i++)
          fwScore2 += specSet[spec2][idxMatched2_zero[i]][1];

        for (shift = pmDiff - pmTol; shift <= pmDiff + pmTol + 0.00001; shift +=
            shiftIncrement)
        {
          ScoreOverlap6(specSet[spec1],
                        specSet[spec2],
                        shift,
                        peakTol,
                        idxMatched1_other,
                        idxMatched2_other);

          int rvScore1 = 0, rvScore2 = 0;
          for (int i = 0; i < idxMatched1_other.size(); i++)
            rvScore1 += specSet[spec1][idxMatched1_other[i]][1];
          for (int i = 0; i < idxMatched2_other.size(); i++)
            rvScore2 += specSet[spec2][idxMatched2_other[i]][1];

          idxUnionEnd = set_union(idxMatched1_zero.begin(),
                                  idxMatched1_zero.end(),
                                  idxMatched1_other.begin(),
                                  idxMatched1_other.end(),
                                  idxUnion.begin());

          for (idxUnionStart = idxUnion.begin(), score1 = 0, numMatchedPeaksU1 =
              0; idxUnionStart != idxUnionEnd;
              idxUnionStart++, numMatchedPeaksU1++)
          {
            score1 += specSet[spec1][(*idxUnionStart)][1];
          }
          numMatchedPeaks = numMatchedPeaksU1;

          idxUnionEnd = set_union(idxMatched2_zero.begin(),
                                  idxMatched2_zero.end(),
                                  idxMatched2_other.begin(),
                                  idxMatched2_other.end(),
                                  idxUnion.begin());

          for (idxUnionStart = idxUnion.begin(), score2 = 0, numMatchedPeaksU2 =
              0; idxUnionStart != idxUnionEnd;
              idxUnionStart++, numMatchedPeaksU2++)
          {
            score2 += specSet[spec2][(*idxUnionStart)][1];
          }
          numMatchedPeaks = min(numMatchedPeaks, numMatchedPeaksU2);

          if (((score1 + score2 > bestScore1 + bestScore2)
              || ((score1 + score2 == bestScore1 + bestScore2)
                  && abs(shift - pmDiff) < abs(bestShift - pmDiff))))
          {
            bestScore1 = score1;
            bestScore2 = score2;
            bestShift = shift;
            bestNumMatchPeak = numMatchedPeaks;
            bestFwScore1 = fwScore1;
            bestFwScore2 = fwScore2;
            bestRvScore1 = rvScore1;
            bestRvScore2 = rvScore2;
          }
        }

        TwoValues<double> pv_pairs;
        double ***spec2AlignGFProbabiltyDensity =
            getALGFProbDensity(specSet[spec2], totalSpecScores[spec2], peakTol);
        pv_pairs[0] = getALGFPValue(spec1AlignGFProbabiltyDensity,
                                    specSet[spec1],
                                    bestFwScore1,
                                    bestRvScore1,
                                    specSet[spec2].parentMass);
        pv_pairs[1] = getALGFPValue(spec2AlignGFProbabiltyDensity,
                                    specSet[spec2],
                                    bestFwScore2,
                                    bestRvScore2,
                                    specSet[spec1].parentMass);
        pvalues.push_back(pv_pairs);
        freeALGFTable(spec2AlignGFProbabiltyDensity, specSet[spec2]);

      }  			 //spec2 loop
      freeALGFTable(spec1AlignGFProbabiltyDensity, specSet[spec1]);
    }  			 //spec1 loop
  }

  static float getAlignScore(Spectrum &spec,
                             Spectrum &prm,
                             float startX,
                             float endX,
                             float peakTol)
  {

    float start = startX - 18, end = endX + 18;
    int sti = 0, edi = prm.size();
    for (int i = 0; i < prm.size(); i++)
    {
      if (prm[i][0] > start)
      {
        sti = i;
        break;
      }
    }
    for (int i = prm.size() - 1; i > -1; i--)
    {
      if (prm[i][0] < end)
      {
        edi = i;
        break;
      }
    }

    float score = 0;
    int cur = 0;
    for (int i = sti; i <= edi; i++)
    {
      float mass = prm[i][0];
      for (int j = cur; j < spec.size(); j++)
      {
        if (spec[j][0] < mass - peakTol)
          continue;
        else if (spec[j][0] > mass + peakTol)
        {
          cur = j;
          break;
        }
        score += spec[j][1];
      }
    }
    return score;
  }

  void getPairAlignGFPValues22(SpecSet & specSet,
                               SpectrumPairSet & spectraPairs,
                               vector<TwoValues<double> > & pvalues,
                               float pmTol,
                               float peakTol,
                               AAJumps & jumps,
                               bool specTypeMSMS)
  {
    float ionOffset = specTypeMSMS ? AAJumps::massHion : 0;
    for (unsigned int i = 0; i < specSet.size(); i++)
    {
      specSet[i].addZPMpeaks(peakTol, ionOffset, true);
    }

    for (unsigned int i = 0; i < specSet.size(); i++)
    { //convert intensity into integer
      for (unsigned int j = 0; j < specSet[i].size(); j++)
      {
        specSet[i][j][1] = ceil(specSet[i][j][1]);
        // specSet[i][j][1] = floor(specSet[i][j][1]);
        // specSet[i][j][1] = round(specSet[i][j][1]);//
      }
      if (specSet[i].size() < 4)
        continue;
      specSet[i][0][1] = 0;
      specSet[i][specSet[i].size() - 1][1] = 0;
      specSet[i].removePeak(1);
      specSet[i].removePeak(specSet[i].size() - 2);
    }    			 //*/

    vector<float> totalSpecScores;
    totalSpecScores.resize(specSet.size());
    for (unsigned int i = 0; i < specSet.size(); i++)
    {
      totalSpecScores[i] = 0;
      for (unsigned int p = 0; p < specSet[i].size(); p++)
        totalSpecScores[i] += specSet[i][p][1];
    }

    vector<vector<TwoValues<int> > > candidates(specSet.size());
    unsigned int numEntries = spectraPairs.size();
    for (unsigned int i = 0; i < numEntries; i++)
    {
      candidates[spectraPairs[i].spec1].push_back(TwoValues<int>(spectraPairs[i].spec2, i));
    }

    pvalues.resize(spectraPairs.size());

    for (int spec1 = 0; spec1 < candidates.size(); spec1++)
    {
      if (candidates[spec1].size() == 0)
        continue;

      float endMass1 = specSet[spec1].parentMass; //specSet[spec1][specSet[spec1].size()-1][0];
      string pept1 = specSet[spec1].psmList.front()->m_annotation;
      Spectrum prm1;
      jumps.getPRMMasses(pept1, prm1, 0, true);
      double ***spec1AlignGFProbabiltyDensity =
          getALGFProbDensity(specSet[spec1], totalSpecScores[spec1], peakTol);

      for (int ipair = 0; ipair < candidates[spec1].size(); ipair++)
      {
        int spec2 = candidates[spec1][ipair][0];

        int pairIndex = candidates[spec1][ipair][1];

        float endMass2 = specSet[spec2].parentMass; //specSet[spec2][specSet[spec2].size()-1][0];
        string pept2 = specSet[spec2].psmList.front()->m_annotation;
        Spectrum prm2;
        jumps.getPRMMasses(pept2, prm2, 0, true);
        double ***spec2AlignGFProbabiltyDensity =
            getALGFProbDensity(specSet[spec2], totalSpecScores[spec2], peakTol);

        TwoValues<double> pv_pairs;
        float shift = spectraPairs[pairIndex].shift1;

        if (shift > 0)
        {
          float overlap = endMass1 - shift;
          if (overlap <= endMass2 + pmTol)
          {
            pv_pairs[0] = getReverseALGFPValue(spec1AlignGFProbabiltyDensity,
                                               specSet[spec1],
                                               getAlignScore(specSet[spec1],
                                                             prm1,
                                                             shift,
                                                             endMass1,
                                                             peakTol),
                                               overlap);
            pv_pairs[1] = getForwardALGFPValue(spec2AlignGFProbabiltyDensity,
                                               specSet[spec2],
                                               getAlignScore(specSet[spec2],
                                                             prm2,
                                                             0,
                                                             overlap,
                                                             peakTol),
                                               overlap);
          }
          else
          {    			 //spec2 is completely included.
            pv_pairs[0] =
                getALGFPValueForInternalRegion(specSet[spec1],
                                               totalSpecScores[spec1],
                                               peakTol,
                                               shift,
                                               shift + endMass2,
                                               getAlignScore(specSet[spec1],
                                                             prm1,
                                                             shift,
                                                             shift + endMass2,
                                                             peakTol));
            pv_pairs[1] = getForwardALGFPValue(spec2AlignGFProbabiltyDensity,
                                               specSet[spec2],
                                               getAlignScore(specSet[spec2],
                                                             prm2,
                                                             0,
                                                             endMass2,
                                                             peakTol),
                                               endMass2);
          }
        }
        else if (shift < 0)
        {
          float overlap = endMass2 + shift;
          if (overlap <= endMass1 + pmTol)
          {
            pv_pairs[0] = getForwardALGFPValue(spec1AlignGFProbabiltyDensity,
                                               specSet[spec1],
                                               getAlignScore(specSet[spec1],
                                                             prm1,
                                                             0,
                                                             overlap,
                                                             peakTol),
                                               overlap);
            pv_pairs[1] = getReverseALGFPValue(spec2AlignGFProbabiltyDensity,
                                               specSet[spec2],
                                               getAlignScore(specSet[spec2],
                                                             prm2,
                                                             -shift,
                                                             endMass2,
                                                             peakTol),
                                               overlap);
          }
          else
          {    			 //spec1 is completely included.
            pv_pairs[0] = getForwardALGFPValue(spec1AlignGFProbabiltyDensity,
                                               specSet[spec1],
                                               getAlignScore(specSet[spec1],
                                                             prm1,
                                                             0,
                                                             endMass1,
                                                             peakTol),
                                               endMass1);
            pv_pairs[1] =
                getALGFPValueForInternalRegion(specSet[spec2],
                                               totalSpecScores[spec2],
                                               peakTol,
                                               -shift,
                                               -shift + endMass1,
                                               getAlignScore(specSet[spec2],
                                                             prm2,
                                                             -shift,
                                                             -shift + endMass1,
                                                             peakTol));
          }
        }
        else
        {
          float overlap = min(endMass1, endMass2);
          pv_pairs[0] = getForwardALGFPValue(spec1AlignGFProbabiltyDensity,
                                             specSet[spec1],
                                             getAlignScore(specSet[spec1],
                                                           prm1,
                                                           0,
                                                           overlap,
                                                           peakTol),
                                             overlap);
          pv_pairs[1] = getForwardALGFPValue(spec2AlignGFProbabiltyDensity,
                                             specSet[spec2],
                                             getAlignScore(specSet[spec2],
                                                           prm2,
                                                           0,
                                                           overlap,
                                                           peakTol),
                                             overlap);
        }
        pvalues[pairIndex] = pv_pairs;
        freeALGFTable(spec2AlignGFProbabiltyDensity, specSet[spec2]);
      }    			 //spec2 loop
      freeALGFTable(spec1AlignGFProbabiltyDensity, specSet[spec1]);
    }    			 //spec1 loop
  }
//---------------------------------------------------------------------------

  void getSingleModPairs(SpecSet & prmSpec,//NO partial overlap, align-gf, tag-alignment
		  	  	  	  	 SpecSet & normMSSpec,//for cosine calculation
						 vector<int> & baseSpectraIdx,
						 short aaDiff,
						 float minShift,
						 float maxShift,
						 float pmTol,
						 float peakTol,
						 float minMatchRatio,
						 short minNumMatchedPeaks,
						 float maxPvalue,
						 bool  highMS2,
						 bool  separateCharge,
						 SpectrumPairSet & results,
						 vector<TwoValues<float> > & ratios,
						 vector<TwoValues<float> > & means,
						 vector<float> & varTerms,
						 ExecTagFilterPairs *tag_filter,
						 float resolution,
						 float symmetryOffset)
  	{
  	  int j, massIdx;
  	  float pMass1;
  	  float shiftIncrement = resolution; // Step size used to increment shifts from -pmTol to pmTol
  	  list<int>::iterator curPair;
  	  float histResolution = 0.1; // Resolution of the parent mass histogram - only used for binning, independent of overall peak resolution settings
  	  vector<list<int> > hist;
  	  TwoValues<float> curRatios;
  	  AAJumps jumps(-1);
  	  SpectrumPair curResult;

  	  SpecSet specSet(prmSpec.size());//creating integer value spectra
  	  for (unsigned int n = 0; n < specSet.size(); n++) {
  		specSet[n] = prmSpec[n];
  		specSet[n].copyNP(prmSpec[n]);
  		specSet[n].rankFilterPeaks(10, 50);//org 10
  	  }

  	  for (unsigned int i = 0; i < specSet.size(); i++) {
		  for (unsigned int j = 0; j < specSet[i].size(); j++) {
			  specSet[i][j][1] = ceil(specSet[i][j][1]);
		//	  specSet[i][j][1] = round(specSet[i][j][1]);
		//	  if( specSet[i][j][1] == 0 ) specSet[i].removePeak(j--);
		  }
		  specSet[i].addZPMpeaks22(peakTol, 0, true, 0.1);
  	  }

  	  specSet.getMassesHistogram(hist, histResolution);

  	  // Determine parent mass differences for valid matches
  	  if (aaDiff > 0)
  	  {
  		jumps.getjumps(aaDiff, resolution);
  		jumps.forceJump(0);
  		jumps.forceTolerance(pmTol, resolution);
  		jumps.forceDoubleSided();
  		jumps.forceUnique(resolution);
  		maxShift = jumps.masses[jumps.size()-1] + 3;
  	  }
  	  DEBUG_VAR(maxShift);

  	  {// shiftIncrement should be factor of pmTol.
  		int denominator= 2;
  		while( pmTol/denominator > peakTol/2 ){
  			denominator++;
  		}
  		while( pmTol/denominator > resolution ){
  			denominator++;
  		}
  		shiftIncrement = pmTol/denominator;
  	  }
  	  DEBUG_VAR(shiftIncrement);

  	  // Determine totalScores per spectrum for later computation of match ratios
  	  vector<float> totalSpecScores(specSet.size());
  	  for (unsigned int i = 0; i < specSet.size(); i++) {

  		totalSpecScores[i] = 0;
  		for (unsigned int p = 0; p < specSet[i].size(); p++){
  		  totalSpecScores[i] += specSet[i][p][1];
  		}
  	  }

  	  vector<vector<float> > alignGF_pvalue_set(specSet.size());
  	  for (unsigned int i = 0; i < specSet.size(); i++)
	  {
		if (specSet[i].size() == 0 or specSet[i].parentMass < 400) continue;

		vector<float> pden;
		getALGFProbDensity(pden, specSet[i], totalSpecScores[i], peakTol, 0);
		alignGF_pvalue_set[i] = pden;
	  }

  	  // Initialize means vector
  	  means.resize(specSet.size());
  	  varTerms.resize(specSet.size());
  	  for (unsigned int i = 0; i < means.size(); i++)
  	  {
  		means[i].set(0, 0);
  		varTerms[i] = 0;
  	  }

  	  // Iterate through pairs using FindMatchPeaksAll and ScoreOverlap6
  	  //    vector<int> idx1,idx2;                             idx1.reserve(1500);               idx2.reserve(1500);
  	  vector<int> idxMatched1_zero, idxMatched1_other;
  	  idxMatched1_zero.reserve(1500);
  	  idxMatched1_other.reserve(1500);
  	  vector<int> idxMatched2_zero, idxMatched2_other;
  	  idxMatched2_zero.reserve(1500);
  	  idxMatched2_other.reserve(1500);
  	  vector<int> idxUnion;
  	  idxUnion.reserve(1500);
  	  vector<int>::iterator idxUnionEnd, idxUnionStart;
  	  int spec1, spec2, maxShiftInt = (int)ceil(maxShift / histResolution);
  	  minShift = max((float)0, minShift);
  	  short numMatchedPeaks, numMatchedPeaksU1, numMatchedPeaksU2;
  	  float pmDiff, score1, score2, shift, bestScore1, bestScore2, bestShift;

  	  //    clock_t startTime = clock();    double curTime, totTime=0;
  	  long befTag=0, aftTag=0;
  	  list<int> candidates;
  	  list<int>::iterator curCandidate;
  	  for (unsigned int i = 0; i < baseSpectraIdx.size(); i++)
  	  {
  		DEBUG_MSG( "Starting with spectrum " << baseSpectraIdx[i] << "...");

  		// Find candidate spectra to compare to
  		spec1 = baseSpectraIdx[i];
  		pMass1 = specSet[spec1].parentMass;

  		int massInt = (int)round(pMass1 / histResolution);
  		candidates.clear(); // Make sure candidates is empty from previous run
  		if (specSet[spec1].size() == 0 or specSet[spec1].parentMass < 400)
  		{
  		  DEBUG_MSG("(spectrum is empty or parent mass too low)");
  		  continue;
  		}
  		// Find candidates in [minShift,maxShift]
  		int upperLimit = max(0, (int)ceil((pMass1 - minShift) / histResolution)); // Look at interval below current parent mass
  		for (massIdx = (int)floor((pMass1 - min(pMass1, maxShift))
  			/ histResolution); massIdx < upperLimit; massIdx++)
  		{
  		  for (curPair = hist[massIdx].begin(); curPair != hist[massIdx].end(); curPair++)
  		  {
  			    if (*curPair > spec1 and fabs(pMass1 - specSet[*curPair].parentMass) >= minShift - 2 * pmTol
  								and fabs(pMass1 - specSet[*curPair].parentMass) <= maxShift + 2 * pmTol)
  				{
  					if( separateCharge && specSet[spec1].parentCharge != specSet[*curPair].parentCharge ) continue;
  					candidates.push_back(*curPair);
  				}
  		  }
  		}
  		upperLimit = min((int)hist.size(), massInt + maxShiftInt); // Look at interval above current parent mass
  		for (massIdx = (int)floor((pMass1 + minShift) / histResolution); massIdx
  			< upperLimit; massIdx++)
  		{
  		  for (curPair = hist[massIdx].begin(); curPair != hist[massIdx].end(); curPair++)
  		  {
  			    if (*curPair > spec1 and fabs(pMass1 - specSet[*curPair].parentMass) >= minShift - 2 * pmTol
  						and fabs(pMass1 - specSet[*curPair].parentMass) <= maxShift + 2 * pmTol)
  				{
  					if( separateCharge && specSet[spec1].parentCharge != specSet[*curPair].parentCharge ) continue;
  					candidates.push_back(*curPair);
  				}
  		  }
  		}

  		candidates.sort();
  		candidates.unique(); // The same spectrum may occur in multiple bins - make sure the code below only aligns once to each eligible spectrum
  		DEBUG_MSG(" num. candidates is " << candidates.size());
  		befTag += candidates.size();

  		if (tag_filter->hasTags())
		{
			list<int> subsetOfCandidates;
			tag_filter->getPairs(spec1, candidates, subsetOfCandidates);
			candidates.swap(subsetOfCandidates);
			aftTag += candidates.size();
			DEBUG_MSG(" (" << candidates.size() << " after tag filtering)");
		}
		DEBUG_MSG("..."); //*/

  		curCandidate = candidates.begin();
  		while (curCandidate != candidates.end())
  		{
  		  spec2 = *curCandidate;
  		  curCandidate++;

  		  if (spec1 >= specSet.size() || spec2 >= specSet.size())
  		  {
  			ERROR_MSG("ERROR: Alignment pair (" << spec1 << "," << spec2
  				<< ") is invalid - PRM spectra set has only " << specSet.size()
  				<< " spectra!\n");
  			exit(-1);
  		  }

  		  if (specSet[spec2].size() == 0 or specSet[spec2].parentMass < 400) continue;

  		  //cerr<<spec2<<":"; cerr.flush();
  		  pmDiff = specSet[spec1].parentMass - specSet[spec2].parentMass;
  		  bestScore1 = 0;
  		  bestScore2 = 0;

  		  int bestNumMatchPeak = 0;
  		  if (abs(pmDiff) > pmTol + 0.00001) {// if modified

  		  //cerr<<"X"; cerr.flush();
  			ScoreOverlap6(specSet[spec1],
  						  specSet[spec2],
  						  0,
  						  peakTol,
  						  idxMatched1_zero,
  						  idxMatched2_zero);//*/

  			//cerr<<"X"; cerr.flush();

  			for (shift = pmDiff - pmTol; shift <= pmDiff + pmTol + 0.00001; shift
  				+= shiftIncrement)
  			{
  			  //FindMatchPeaksAll(specSet[spec1], specSet[spec2], shift, peakTol, idx1, idx2);
  			  //cerr<<"<"; cerr.flush();
  			  ScoreOverlap6(specSet[spec1],
  							specSet[spec2],
  							shift,
  							peakTol,
  							idxMatched1_other,
  							idxMatched2_other);

  			  //cerr<<">"; cerr.flush();

  			 idxUnionEnd = set_union(idxMatched1_zero.begin(),
  									  idxMatched1_zero.end(),
  									  idxMatched1_other.begin(),
  									  idxMatched1_other.end(),
  									  idxUnion.begin());

  			  for (idxUnionStart = idxUnion.begin(), score1 = 0, numMatchedPeaksU1 = 0; idxUnionStart
  				  != idxUnionEnd; idxUnionStart++, numMatchedPeaksU1++) {
  				score1 += specSet[spec1][(*idxUnionStart)][1];
  			  }
  			  numMatchedPeaks = numMatchedPeaksU1;

  			  idxUnionEnd = set_union(idxMatched2_zero.begin(),
  									  idxMatched2_zero.end(),
  									  idxMatched2_other.begin(),
  									  idxMatched2_other.end(),
  									  idxUnion.begin());

  			  for (idxUnionStart = idxUnion.begin(), score2 = 0, numMatchedPeaksU2 = 0; idxUnionStart
  				  != idxUnionEnd; idxUnionStart++, numMatchedPeaksU2++) {
  				score2 += specSet[spec2][(*idxUnionStart)][1];
  			  }
  			  numMatchedPeaks = min(numMatchedPeaks, numMatchedPeaksU2);

  			  if ( (score1 + score2 > bestScore1 + bestScore2) ||
  					 ( (score1 + score2 == bestScore1 + bestScore2) && abs(shift - pmDiff) < abs(bestShift - pmDiff) ) )
  			  {
  				//DEBUG_MSG(*curCandidate << "  " << numMatchedPeaksU1 << "  " << numMatchedPeaksU2 << "  " << numMatchedPeaks);
  				bestScore1 = score1;
  				bestScore2 = score2;
  				bestShift = shift;
  				bestNumMatchPeak = numMatchedPeaks;
  			  }
  			}
  		  }
  		  else {
  			  vector<int> idx1all, idx2all, idxMatched1, idxMatched2;

  			   for (shift = -pmTol; shift <= pmTol + 0.00001; shift
  				   += shiftIncrement)
  			   {
  				 //cerr<<"{"; cerr.flush();
  				 FindMatchPeaksAll(specSet[spec1],
  								   specSet[spec2],
  								   shift,
  								   peakTol,
  								   idx1all,
  								   idx2all);
  				 //cerr<<"}"; cerr.flush();
  				 ScoreOverlap7(specSet[spec1],
  							   idx1all,
  							   specSet[spec2],
  							   idx2all,
  							   shift,
  							   peakTol,
  							   idxMatched1,
  							   idxMatched2,
  							   symmetryOffset);

  				 numMatchedPeaks = (short)idxMatched1.size();

  				 score1 = 0;
  				 score2 = 0;
  				 for (int i = 0; i < idxMatched1.size(); i++){
  				   score1 += specSet[spec1][idxMatched1[i]][1];
  				 }
  				 for (int i = 0; i < idxMatched2.size(); i++){
  				   score2 += specSet[spec2][idxMatched2[i]][1];
  				 }

  				 if( (score1 + score2 > bestScore1 + bestScore2) ||
  					 (score1 + score2 == bestScore1 + bestScore2 && abs(shift) < abs(bestShift) ) )
  				 {
  				   bestScore1 = score1;
  				   bestScore2 = score2;
  				   bestShift = shift;
  				   bestNumMatchPeak = numMatchedPeaks;
  				/*    cerr<<"Matching peaks with shift "<<bestShift<<", scores ("<<score1<<","<<score2<<"):\nSpectrum "<<spec1<<":\n";
  					for(unsigned int i=0; i<idxMatched1.size(); i++) cerr<<specSet[spec1][idxMatched1[i]][0]<<"\t"<<specSet[spec1][idxMatched1[i]][1]<<endl;
  					cerr<<"Spectrum "<<spec2<<":\n";
  					for(unsigned int i=0; i<idxMatched2.size(); i++) cerr<<specSet[spec2][idxMatched2[i]][0]<<"\t"<<specSet[spec2][idxMatched2[i]][1]<<endl;
  					//*/
  				 }
  			   }
  		  }

  		//	cout << bestRvScore1 << " " << totalSpecScores[spec1] << " | " << bestRvScore2 << " " << totalSpecScores[spec2] << endl;
  		//	cout << bestRvScore1 << " ******************** " << endl;
  			curResult.spec1 = spec1;
  			curResult.spec2 = spec2;
  			curResult.score1 = bestScore1;
  			curResult.score2 = bestScore2;
  			curResult.shift1 = bestShift;

  			means[spec1][0]
  				= (means[spec1][0] * means[spec1][1] + curResult.score1)
  					/ (means[spec1][1] + 1);
  			means[spec2][0]
  				= (means[spec2][0] * means[spec2][1] + curResult.score2)
  					/ (means[spec2][1] + 1);
  			means[spec1][1]++;
  			means[spec2][1]++;
  			float updateRatio = (means[spec1][1] - 1) / means[spec1][1]; // Update variance terms
  			varTerms[spec1] = updateRatio * varTerms[spec1] + (curResult.score1
  				* curResult.score1) / means[spec1][1];
  			updateRatio = (means[spec2][1] - 1) / means[spec2][1];
  			varTerms[spec2] = updateRatio * varTerms[spec2] + (curResult.score2
  				* curResult.score2) / means[spec2][1];

  			//  Compute match ratios
  			curRatios[0] = curResult.score1 / totalSpecScores[spec1];
  			curRatios[1] = curResult.score2 / totalSpecScores[spec2];

  			int dynamicNumPeaks = (highMS2)? minNumMatchedPeaks : ( min(specSet[spec1].parentMass, specSet[spec2].parentMass)/300 + minNumMatchedPeaks );
  			if ( curRatios[0] >= minMatchRatio && curRatios[1] >= minMatchRatio && bestNumMatchPeak >= dynamicNumPeaks )
  			{
  			  double bestPValue1 = alignGF_pvalue_set[spec1][bestScore1];
  			  double bestPValue2 = alignGF_pvalue_set[spec2][bestScore2];

  			  if( bestPValue1 < maxPvalue && bestPValue2 < maxPvalue ) {
  				//! /todo Implement an efficient push_back (or equivalent) for SpectrumPairSet
  				curResult.score1 = (bestPValue1==0)? 1000 : -log10(bestPValue1);
  				curResult.score2 = (bestPValue2==0)? 1000 : -log10(bestPValue2);

  				if( normMSSpec.size() != 0 ){
					unsigned int cosMatchedPeaks = 0;
					curResult.shift2 = normMSSpec[spec1].scoreMatch(normMSSpec[spec2],
											peakTol,
											cosMatchedPeaks,
											score1,
											score2,
											false,
											false,
											false);
  				}

  				results.push_back(curResult);
  				ratios.push_back(curRatios);
  			  }
  			}
  		} // while over all candidates for spec1
  		DEBUG_MSG(" done. Current number of results in memory: " << results.size());

  	  } // for i over baseSpectraIdx

  	  if( tag_filter->hasTags() )
  	    DEBUG_MSG("#all possible pairs: " << befTag << " -> #pairs after matching-tag filtering  " << aftTag);

  	  DEBUG_MSG("Done with the alignment!");
  }

  void getPartialOverlapsPairs(SpecSet & prmSpec,//partial overlaps, align-gf, tag-alignment
							   unsigned int startIdx,
							   unsigned int endIdx,
							   float peakTol,
							   float pmTol,
							   float minRatio,
							   float minPeakAreaOvlp,
							   short minNumMatchedPeaks,
							   float maxPvalue,
							   bool  highMS2,
							   bool  separateCharge,
							   AAJumps & allowedJumps,
							   float minAbsShift,
							   SpectrumPairSet & results,
							   vector<TwoValues<float> > & ratios,
							   list<TwoValues<int> > & numMatchedPeaks,
							   vector<TwoValues<float> > & means,
							   vector<float> & varTerms,
							   list<vector<float> > & alignStats,
							   vector<vector<float> > & specStats,
							   ExecTagFilterPairs *tag_filter,
							   float resolution,
							   float maxShift)
    {
      unsigned int spec1, spec2;

      // Variables and initializations
      TwoValues<int> res;
      list<float> shiftScores; // Maximum scores for eligible pairs of shifts
      list<TwoValues<unsigned int> > shiftPairs; // List of eligible pairs of shifts (in the same order as shiftScores)
      //  vector<vector<TwoValues<int> > > shiftMatchedPeaks;  // Lists of matched peaks per shift between the two spectra (indexed using
      vector<list<TwoValues<int> > > shiftMatchedPeaks; // Lists of matched peaks per shift between the two spectra (indexed using
      //  a shift offset returned by computeShifts)
      vector<TwoValues<float> > shiftDPscores; // Match score per shift after applying dynamic programming
      vector<float> shiftDPpenalties; // Peak mass mismatch penalties (in DP)
      vector<TwoValues<vector<int> > > shiftDPmatchedPeaks; // Peaks matched by the DP algorithm
      SpectrumPair curResult;
      vector<float> curAlignStats(4);

      SpecSet specSet(prmSpec.size());//creating integer value spectra
	  for (unsigned int n = 0; n < specSet.size(); n++) {
		specSet[n] = prmSpec[n];
		specSet[n].copyNP(prmSpec[n]);
		specSet[n].rankFilterPeaks(10, 50);//org 10
	  }

	  for (unsigned int i = 0; i < specSet.size(); i++) {
		  for (unsigned int j = 0; j < specSet[i].size(); j++) {
			  specSet[i][j][1] = ceil(specSet[i][j][1]);
		//	  specSet[i][j][1] = round(specSet[i][j][1]);
		//	  if( specSet[i][j][1] == 0 ) specSet[i].removePeak(j--);
		  }
		  specSet[i].addZPMpeaks22(peakTol, 0, true, 0.1);
	  }

	  {
		int denominator = 2;
		if( resolution > peakTol/denominator ) resolution = peakTol/denominator;
		if( resolution > pmTol/denominator )   resolution = pmTol/denominator;
	  }
	  DEBUG_VAR(resolution);//*/

      // Initialize stats vectors
      means.resize(specSet.size());
      varTerms.resize(specSet.size());
      specStats.resize(specSet.size());
      for (unsigned int i = 0; i < means.size(); i++)
      {
        means[i].set(0, 0);
        varTerms[i] = 0;
        specStats[i].resize(9);
        for (unsigned int j = 0; j < 9; j++)
          specStats[i][j] = 0;
      }

      vector<TwoValues<float> > minMaxMatchScores(specSet.size()); // Minimum acceptable / Maximum possible match scores per spectrum in the dataset
      float maxParentMass = 0;
      unsigned int maxNumPeaks = 0;
      int intPeakTol = (int)round(peakTol / resolution), intPMTol = (int)round(pmTol / resolution);
      for (unsigned int i = 0; i < specSet.size(); i++)
      {
        if (specSet[i].parentMass > maxParentMass)
          maxParentMass = specSet[i].parentMass;
        if (specSet[i].size() > maxNumPeaks)
          maxNumPeaks = specSet[i].size();
        minMaxMatchScores[i].set(0, 0);
        for (unsigned int j = 0; j < specSet[i].size(); j++)
          minMaxMatchScores[i][1] += specSet[i][j][1];
        minMaxMatchScores[i][0] = minRatio * minMaxMatchScores[i][1];
      }

      vector<vector<float> > alignGF_pvalue_set(specSet.size());
      for (unsigned int i = 0; i < specSet.size(); i++)
	  {
		if (specSet[i].size() == 0 or specSet[i].parentMass < 400) continue;
		vector<float> pden;
		getALGFProbDensity(pden, specSet[i], minMaxMatchScores[i][1], peakTol, 0);
		alignGF_pvalue_set[i] = pden;
	  }

      //  shiftMatchedPeaks.resize(1+(int)round(2*(min(maxParentMass,maxShift)+8*peakTol)/resolution)); //for(unsigned int i=0;i<shiftMatchedPeaks.size();i++) shiftMatchedPeaks[i].reserve(10);
      shiftMatchedPeaks.resize(1 + (int)ceil( 2*(min(maxParentMass, maxShift)+peakTol+pmTol) / resolution));
      //for(unsigned int i=0;i<shiftMatchedPeaks.size();i++) shiftMatchedPeaks[i].reserve(10);
      shiftDPscores.reserve(shiftMatchedPeaks.capacity());
      shiftDPpenalties.reserve(shiftMatchedPeaks.capacity());
      shiftDPmatchedPeaks.reserve(shiftMatchedPeaks.capacity());
      vector<int> idx1, idx2;
      idx1.reserve(maxNumPeaks * (2 * intPeakTol + 1) + 1);
      idx2.reserve(idx1.capacity()); // Temporary variables used as input to ScoreOverlap6
      vector<int> idxMerged1, idxMerged2;
      idxMerged1.reserve(idx1.capacity());
      idxMerged2.reserve(idx1.capacity()); // Temporary variables used to merge symmetric results of ScoreOverlap6

      // Compute the pairwise alignments
      time_t curTime = time(0);
      double totTime = 0.0;
      int numResults = 0;
      float bestShiftScore;
      TwoValues<int> bestShiftPair;
      TwoValues<int> bestNumMatchedPeaks;
      TwoValues<float> bestShiftScores, bestCandidateScores;
      list<int> shiftsToCompute;

      long numCandidates=0, filteredCands=0;
      for (spec1 = startIdx; spec1 <= endIdx; spec1++)
      {
        DEBUG_MSG("Processing spectrum " << spec1 << "... ");
        if (specSet[spec1].size() == 0 or specSet[spec1].parentMass < 400) continue;

        for (spec2 = spec1 + 1; spec2 < specSet.size(); spec2++)
        {
          if( separateCharge && specSet[spec1].parentCharge != specSet[spec2].parentCharge ) continue;
          if (specSet[spec2].size() == 0 or specSet[spec2].parentMass < 400) continue;

          float minPM = specSet[spec1].parentMass, maxPM = specSet[spec2].parentMass;
          if( maxPM < minPM ){
        	  minPM = specSet[spec2].parentMass;
        	  maxPM = specSet[spec1].parentMass;
          }
          if ( minPM < minPeakAreaOvlp*maxPM ) continue;
          numCandidates++;

          if( tag_filter->hasTags() ){
        	  set<int> possDeltas;
			  tag_filter->getPossibleDeltas(possDeltas, spec1, spec2,
					  specSet[spec1].parentMass, specSet[spec2].parentMass, minPeakAreaOvlp*maxPM);
			  if( possDeltas.size() == 0 ) continue; //NO aligned tags
			  // Step 1: Get the list of possible shifts
			  res = computeShiftsWithFiltered(specSet[spec1],
											  specSet[spec2],
										  	  peakTol,
											  pmTol,
											  minPeakAreaOvlp*maxPM,
											  minMaxMatchScores[spec1][0],
											  minMaxMatchScores[spec2][0],
											  minNumMatchedPeaks,
											  allowedJumps,//not used
											  resolution,
											  possDeltas,
											  shiftScores,
											  shiftPairs,
											  shiftMatchedPeaks,
											  bestCandidateScores,
											  bestNumMatchedPeaks,
											  minAbsShift);//not used

			  if (res[1] == 0) continue; // No eligible deltas found from tag-alignment
			  filteredCands++;
          }
          else{
        	  // Step 1: Get the list of possible shifts
        	  res = computeShifts(specSet[spec1],
								  specSet[spec2],
						   		  peakTol,
								  pmTol,
								  minPeakAreaOvlp*maxPM,
								  minMaxMatchScores[spec1][0],
								  minMaxMatchScores[spec2][0],
								  minNumMatchedPeaks,
								  allowedJumps,
								  resolution,
								  shiftScores,
								  shiftPairs,
								  shiftMatchedPeaks,
								  bestCandidateScores,
								  bestNumMatchedPeaks,
								  minAbsShift);
          }

          for (unsigned int i = 0; i < curAlignStats.size(); i++)
            curAlignStats[i] = 0;

          // Matched intensity before DP/merging
          auxUpdateMeanVariance(bestCandidateScores[0],
                                specStats[spec1][0],
                                specStats[spec1][5],
                                specStats[spec1][6]);
          auxUpdateMeanVariance(bestCandidateScores[1],
                                specStats[spec2][0],
                                specStats[spec2][5],
                                specStats[spec2][6]);
          // Num matched peaks before DP/merging
          auxUpdateMeanVariance(bestNumMatchedPeaks[0],
                                specStats[spec1][0],
                                specStats[spec1][1],
                                specStats[spec1][2]);
          auxUpdateMeanVariance(bestNumMatchedPeaks[1],
                                specStats[spec2][0],
                                specStats[spec2][1],
                                specStats[spec2][2]);
          // Pre-DP alignment stats, in case this alignment gets reported
          curAlignStats[0] = min(bestNumMatchedPeaks[0], bestNumMatchedPeaks[1]);
          curAlignStats[1] = min(bestCandidateScores[0]
                                     / minMaxMatchScores[spec1][1],
                                 bestCandidateScores[1]
                                     / minMaxMatchScores[spec2][1]);

          /*          means[spec1][0]=(means[spec1][0]*means[spec1][1]+bestCandidateScores[0])/(means[spec1][1]+1);  // Include accepted pairs in means
           means[spec2][0]=(means[spec2][0]*means[spec2][1]+bestCandidateScores[1])/(means[spec2][1]+1);
           means[spec1][1]++;    means[spec2][1]++;
           float updateRatio = (means[spec1][1]-1)/means[spec1][1];  // Update variance terms
           varTerms[spec1] = updateRatio*varTerms[spec1]+(bestCandidateScores[0]*bestCandidateScores[0])/means[spec1][1];
           updateRatio = (means[spec2][1]-1)/means[spec2][1];
           varTerms[spec2] = updateRatio*varTerms[spec2]+(bestCandidateScores[1]*bestCandidateScores[1])/means[spec2][1];
           */

          //cerr<<"shiftDPscores resized to "<<shiftMatchedPeaks.size()<<endl;
          shiftDPscores.resize(shiftMatchedPeaks.size());
          for (unsigned int i = 0; i < shiftDPscores.size(); i++)
            shiftDPscores[i].set(-1.0, -1.0);
          shiftDPpenalties.resize(shiftMatchedPeaks.size());
          for (unsigned int i = 0; i < shiftDPpenalties.size(); i++)
            shiftDPpenalties[i] = 0.0;
          shiftDPmatchedPeaks.resize(shiftMatchedPeaks.size());
          for (unsigned int i = 0; i < shiftDPmatchedPeaks.size(); i++)
          {
            shiftDPmatchedPeaks[i][0].resize(0);
            shiftDPmatchedPeaks[i][1].resize(0);
          }

          bestShiftScore = 0;
          bestShiftPair.set(-1, -1);
          bestShiftScores.set(0, 0);
          bestCandidateScores.set(0, 0);
          bestNumMatchedPeaks.set(0, 0); // Reuse to store best DP match score / num matched peaks (used to update means/vars)
          list<float>::iterator nextShiftScore = shiftScores.begin();
          list<TwoValues<unsigned int> >::iterator nextShiftPair =
              shiftPairs.begin();
          while (nextShiftScore != shiftScores.end()
              and *nextShiftScore > bestShiftScore)
          {
            int shiftIndex = (*nextShiftPair)[0], shiftSym = (*nextShiftPair)[1];
            //cerr<<"Shift pair "<<shiftIndex<<"/"<<shiftSym<<endl;
            if (shiftDPscores[shiftIndex][0] + shiftDPscores[shiftIndex][1]
                < 0.0001 and shiftMatchedPeaks[shiftIndex].size() > 0)
              shiftsToCompute.push_back(shiftIndex);
            for (int tolIdx = max(0, shiftSym - intPMTol);
                tolIdx <= shiftSym + intPMTol
                    and tolIdx < (int)shiftDPscores.size(); tolIdx++)
              if (shiftDPscores[tolIdx][0] + shiftDPscores[tolIdx][1] < 0.0001
                  and shiftMatchedPeaks[tolIdx].size() > 0)
                shiftsToCompute.push_back(tolIdx);

            for (list<int>::iterator curShift = shiftsToCompute.begin();
                curShift != shiftsToCompute.end(); curShift++)
            {
              float shiftMass = ((*curShift) - res[0]) * resolution;

              idx1.resize(shiftMatchedPeaks[*curShift].size());
              idx2.resize(shiftMatchedPeaks[*curShift].size());
              int pivot = 0;
              for (list<TwoValues<int> >::iterator iterPeaks =
                  shiftMatchedPeaks[*curShift].begin();
                  iterPeaks != shiftMatchedPeaks[*curShift].end(); iterPeaks++)
              {
                idx1[pivot] = (*iterPeaks)[0];
                idx2[pivot] = (*iterPeaks)[1];
                pivot++;
              }

              /*cerr<<"shift = "<<*curShift<<", idx1=(";
               for(unsigned int dbgIdx=0; dbgIdx<idx1.size(); dbgIdx++) cerr<<idx1[dbgIdx]<<",";
               cerr<<"), idx2=(";
               for(unsigned int dbgIdx=0; dbgIdx<idx2.size(); dbgIdx++) cerr<<idx2[dbgIdx]<<",";
               cerr<<")\n";*/
              ScoreOverlap6mp(specSet[spec1],
                              idx1,
                              specSet[spec2],
                              idx2,
                              shiftMass,
                              peakTol,
                              shiftDPmatchedPeaks[*curShift][0],
                              shiftDPmatchedPeaks[*curShift][1],
                              AAJumps::minAAmass,
                              &shiftDPpenalties[*curShift]);
              //          ScoreOverlap6(specSet[spec1], idx1, specSet[spec2], idx2, shiftMass, peakTol, shiftDPmatchedPeaks[*curShift][0], shiftDPmatchedPeaks[*curShift][1], false);

              /*cerr<<"shift = "<<*curShift<<", matched1=(";
               for(unsigned int dbgIdx=0; dbgIdx<shiftDPmatchedPeaks[*curShift][0].size(); dbgIdx++) cerr<<shiftDPmatchedPeaks[*curShift][0][dbgIdx]<<",";
               cerr<<"), matched2=(";
               for(unsigned int dbgIdx=0; dbgIdx<shiftDPmatchedPeaks[*curShift][1].size(); dbgIdx++) cerr<<shiftDPmatchedPeaks[*curShift][1][dbgIdx]<<",";
               cerr<<")\n";*/
              shiftDPscores[*curShift][0] = 0;
              for (unsigned int i = 0;
                  i < shiftDPmatchedPeaks[*curShift][0].size(); i++)
                shiftDPscores[*curShift][0] +=
                    specSet[spec1][shiftDPmatchedPeaks[*curShift][0][i]][1];
              shiftDPscores[*curShift][1] = 0;
              for (unsigned int i = 0;
                  i < shiftDPmatchedPeaks[*curShift][1].size(); i++)
                shiftDPscores[*curShift][1] +=
                    specSet[spec2][shiftDPmatchedPeaks[*curShift][1][i]][1];
            }
            shiftsToCompute.clear();

            //cerr<<"--- repeat shift pair "<<shiftIndex<<"/"<<shiftSym<<endl;
            for (int tolIdx = max(0, shiftSym - intPMTol);
                tolIdx <= shiftSym + intPMTol
                    and tolIdx < (int)shiftDPmatchedPeaks.size(); tolIdx++)
            {
              //          if(shiftDPscores[shiftIndex][0]+shiftDPscores[tolIdx][0]+shiftDPscores[shiftIndex][1]+shiftDPscores[tolIdx][1]>bestCandidateScores[0]+bestCandidateScores[1])
              //            bestCandidateScores.set(shiftDPscores[shiftIndex][0]+shiftDPscores[tolIdx][0],shiftDPscores[shiftIndex][1]+shiftDPscores[tolIdx][1]);
              //          if(shiftDPscores[shiftIndex][0]+shiftDPscores[tolIdx][0]<minMaxMatchScores[spec1][0] or  // Sole purpose is to avoid unnecessary computations
              //             shiftDPscores[shiftIndex][1]+shiftDPscores[tolIdx][1]<minMaxMatchScores[spec2][0]) continue;
              if (shiftMatchedPeaks[tolIdx].size() == 0)
                continue;

              // merge lists of matched peaks - spec1
              mergeVectors(shiftDPmatchedPeaks[shiftIndex][0],
                           shiftDPmatchedPeaks[tolIdx][0],
                           idxMerged1);
              mergeVectors(shiftDPmatchedPeaks[shiftIndex][1],
                           shiftDPmatchedPeaks[tolIdx][1],
                           idxMerged2);

              /*cerr<<"shifts = "<<shiftIndex<<"/"<<tolIdx<<":\n  idxMerged1: (";
               for(unsigned int dbgIdx=0; dbgIdx<shiftDPmatchedPeaks[shiftIndex][0].size(); dbgIdx++) cerr<<shiftDPmatchedPeaks[shiftIndex][0][dbgIdx]<<",";
               cerr<<") U (";
               for(unsigned int dbgIdx=0; dbgIdx<shiftDPmatchedPeaks[tolIdx][0].size(); dbgIdx++) cerr<<shiftDPmatchedPeaks[tolIdx][0][dbgIdx]<<",";
               cerr<<") -> (";
               for(unsigned int dbgIdx=0; dbgIdx<idxMerged1.size(); dbgIdx++) cerr<<idxMerged1[dbgIdx]<<",";
               cerr<<")\n  idxMerged2: (";
               for(unsigned int dbgIdx=0; dbgIdx<shiftDPmatchedPeaks[shiftIndex][1].size(); dbgIdx++) cerr<<shiftDPmatchedPeaks[shiftIndex][1][dbgIdx]<<",";
               cerr<<") U (";
               for(unsigned int dbgIdx=0; dbgIdx<shiftDPmatchedPeaks[tolIdx][1].size(); dbgIdx++) cerr<<shiftDPmatchedPeaks[tolIdx][1][dbgIdx]<<",";
               cerr<<") -> (";
               for(unsigned int dbgIdx=0; dbgIdx<idxMerged2.size(); dbgIdx++) cerr<<idxMerged2[dbgIdx]<<",";
               cerr<<")\n";*/
              float score1 = 0;
              for (unsigned int i = 0; i < idxMerged1.size(); i++)
                score1 += specSet[spec1][idxMerged1[i]][1];
              float score2 = 0;
              for (unsigned int i = 0; i < idxMerged2.size(); i++)
                score2 += specSet[spec2][idxMerged2[i]][1];
              // Moved (and changed) from above on 2006/09/11
              if (bestShiftScore < 0.0001)
              { // These statements ensure that we get values to update means/variances
                if (score1 + score2
                    > bestCandidateScores[0] + bestCandidateScores[1])
                  bestCandidateScores.set(score1, score2); // Used to update means/vars even if this pair is not returned
                if (idxMerged1.size() + idxMerged2.size()
                    > bestNumMatchedPeaks[0] + bestNumMatchedPeaks[1])
                  bestNumMatchedPeaks.set(idxMerged1.size(), idxMerged2.size());
              }

              if (score1 < minMaxMatchScores[spec1][0]
                  or score2 < minMaxMatchScores[spec2][0])
                continue;
              if (score1 + score2 - shiftDPpenalties[shiftIndex]
                  - shiftDPpenalties[tolIdx] > bestShiftScore)
              {
                bestShiftScore = score1 + score2 - shiftDPpenalties[shiftIndex]
                    - shiftDPpenalties[tolIdx];
                bestShiftPair.set(shiftIndex, tolIdx);
                bestShiftScores.set(score1, score2);
                bestNumMatchedPeaks.set(idxMerged1.size(), idxMerged2.size());
              }
            }
            nextShiftScore++;
            nextShiftPair++;
          }

          // Build results for best match
          TwoValues<float> curRatios;

          float realShift1 = (bestShiftPair[0] - res[0]) * resolution;
          float realShift2 = (bestShiftPair[1] - res[0]) * resolution;
          float overlapedArea = specSet[spec1].parentMass;
          if( 0 < realShift1 ) overlapedArea -= realShift1;
          if( 0 < realShift2 ) overlapedArea -= realShift2;

          int dynamicNumPeaks = (highMS2)? minNumMatchedPeaks : ( overlapedArea/300 + minNumMatchedPeaks );
          if ( bestShiftScore > 0 && min(bestNumMatchedPeaks[0],bestNumMatchedPeaks[1]) >= dynamicNumPeaks )
          {
        	double bestPValue1 = alignGF_pvalue_set[spec1][bestShiftScores[0]];
        	double bestPValue2 = alignGF_pvalue_set[spec2][bestShiftScores[1]];

        	if( bestPValue1 < maxPvalue && bestPValue2 < maxPvalue ) {
				curResult.spec1 = spec1;
				curResult.spec2 = spec2;
			//	curResult.score1 = bestShiftScores[0];
			//	curResult.score2 = bestShiftScores[1];
				curResult.score1 = (bestPValue1==0)? 1000 : -log10(bestPValue1);
				curResult.score2 = (bestPValue2==0)? 1000 : -log10(bestPValue2);
				curResult.shift1 = realShift1;
				curResult.shift2 = realShift2;

				//! /todo Implement an efficient push_back (or equivalent) for SpectrumPairSet
				results.push_back(curResult);

				curRatios.set(curResult.score1 / minMaxMatchScores[spec1][1],
							  curResult.score2 / minMaxMatchScores[spec2][1]);
				ratios.push_back(curRatios);
				numMatchedPeaks.push_back(bestNumMatchedPeaks);
				bestCandidateScores = bestShiftScores;

				// Post-DP alignment stats since this alignment is reported
				curAlignStats[2] = min(bestNumMatchedPeaks[0],
									   bestNumMatchedPeaks[1]);
				curAlignStats[3] = min(bestCandidateScores[0]
										   / minMaxMatchScores[spec1][1],
									   bestCandidateScores[1]
										   / minMaxMatchScores[spec2][1]);
				alignStats.push_back(curAlignStats);
        	}
          }

          // Regular mean/vars output to text files
          auxUpdateMeanVariance(bestCandidateScores[0],
                                means[spec1][1],
                                means[spec1][0],
                                varTerms[spec1]);
          auxUpdateMeanVariance(bestCandidateScores[1],
                                means[spec2][1],
                                means[spec2][0],
                                varTerms[spec2]);
          // Matched intensity after DP/merging
          auxUpdateMeanVariance(bestCandidateScores[0],
                                specStats[spec1][0],
                                specStats[spec1][7],
                                specStats[spec1][8]);
          auxUpdateMeanVariance(bestCandidateScores[1],
                                specStats[spec2][0],
                                specStats[spec2][7],
                                specStats[spec2][8]);
          // Num matched peaks after DP/merging
          //cerr<<"Num matched peaks: "<<bestNumMatchedPeaks[0]<<"/"<<bestNumMatchedPeaks[1]<<", sampleSize = "<<specStats[spec1][0]<<"/"<<specStats[spec2][0];
          auxUpdateMeanVariance(bestNumMatchedPeaks[0],
                                specStats[spec1][0],
                                specStats[spec1][3],
                                specStats[spec1][4]);
          auxUpdateMeanVariance(bestNumMatchedPeaks[1],
                                specStats[spec2][0],
                                specStats[spec2][3],
                                specStats[spec2][4]);
          //cerr<<", updated means = "<<specStats[spec1][3]<<"/"<<specStats[spec2][3]<<endl;
          specStats[spec1][0]++;
          means[spec1][1]++;
          specStats[spec2][0]++;
          means[spec2][1]++;

          /*      means[spec1][0]=(means[spec1][0]*means[spec1][1]+bestCandidateScores[0])/(means[spec1][1]+1);  // Include accepted pairs in means
           means[spec2][0]=(means[spec2][0]*means[spec2][1]+bestCandidateScores[1])/(means[spec2][1]+1);
           means[spec1][1]++;    means[spec2][1]++;
           float updateRatio = (means[spec1][1]-1)/means[spec1][1];  // Update variance terms
           varTerms[spec1] = updateRatio*varTerms[spec1]+(bestCandidateScores[0]*bestCandidateScores[0])/means[spec1][1];
           updateRatio = (means[spec2][1]-1)/means[spec2][1];
           varTerms[spec2] = updateRatio*varTerms[spec2]+(bestCandidateScores[1]*bestCandidateScores[1])/means[spec2][1];
           */}
        time_t newTime = time(0);
        double ellapsed = difftime(newTime, curTime);
        totTime += ellapsed;
        curTime = newTime;
        int N = endIdx, CUR = spec1 - startIdx + 1;
        DEBUG_MSG(ellapsed << " secs, number of new pair aligns = " << results.size() - numResults << " / " << results.size() << ". Average time = " << totTime / CUR << ", ETA = " << (N - CUR - 1) * (N - CUR - 2) * totTime / (CUR * (N - CUR) + (CUR - 1) * CUR / 2.0));

        numResults = results.size();

      }//spec1 loop

      if( tag_filter->hasTags() )
    	  DEBUG_MSG("#all possible pairs: " << numCandidates << " -> #pairs after matching-tag filtering  " << filteredCands);
    }

} //namespace specnets
