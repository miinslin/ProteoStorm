/*
 * IsoEnvelope.cpp
 *
 *  Created on: Jan 18, 2012
 *      Author: aguthals
 */

#include "IsoEnvelope.h"

using namespace std;

namespace specnets
{

  //
  //  Isotopic envelope classes & functions
  //

  // findEnvelope - finds the index of the envelope in envelopes with the
  //   monoisotopic mass closest to mass
  inline unsigned int IsoEnvelope::findEnvelope(float mass)
  {
    unsigned int lowerLim, upperLim, middle;

    if (envelopes.size() == 0)
    {
      cerr << "ERROR in findEnvelope (ms1.cpp): No envelope masses loaded!\n";
      exit(-1);
    }

    lowerLim = 0;
    upperLim = envelopes.size() - 1;
    while (upperLim > lowerLim + 1)
    { // Set lowerLim to the index of the largest monoisotopic mass < mass
      middle = (lowerLim + upperLim) / 2;
      if (envelopes[middle][0] > mass)
        upperLim = middle;
      else
        lowerLim = middle;
    }
    if (fabs(envelopes[lowerLim][0] - mass) < fabs(envelopes[upperLim][0]
        - mass))
      return lowerLim;
    else
      return upperLim;
  }

  //
  //  makeStrict - converts a vector [mass p1 p2 p3] to [mass 0 p1 0 p2 0 p3]
  //
  inline void IsoEnvelope::makeStrict(unsigned int idxEnv, vector<float> &out)
  {
    out.resize(1 + envelopeSize + envelopeSize);
    //out[0] = envelopes[idxEnv][0];
    for (unsigned int pivot = 1, counter = 1; counter <= envelopeSize; pivot
        += 2, counter++)
    {
      out[pivot] = 0;
      out[pivot + 1] = envelopes[idxEnv][counter];
    }
    out[0] = 0;
    normalize(out, true); // Need non-zero values in every bin to avoid division errors in the hypothesis test
    out[0] = envelopes[idxEnv][0];
  }

  float IsoEnvelope::GetMonoisotopicMass(float in_mass, unsigned short charge)
  {
    return (in_mass * ((float)((short)charge))) - ((float)((short)charge))
        + AAJumps::massHion;
  }

  //
  //  LoadModel - loads a set of monoisotopic masses and associated average isotopic envelopes.
  //    File is generated in Matlab (getIsotopeTable.m/saveIsotopeTable.m) and the format is:
  //     number of lines, cols in data (2x int32)
  //     isotopic mass, propensity of isotopic peak at isotopic mass + [0 1 2 3 4] (6x float)
  //
  bool IsoEnvelope::LoadModel(const char *filename)
  {
    if (Load_binArray(filename, envelopes) > 0)
      if (envelopes.size() > 0)
      {
        envelopeSize = envelopes[0].size() - 1;
        for (unsigned int i = 0; i < envelopes.size(); i++)
          normalize(envelopes[i], true, 1);
        return true;
      }
    return false;
  }

  float IsoEnvelope::ScoreEnvelope(float monoisotopicMass,
                                   vector<float> &massEnvelope,
                                   bool strictMode)
  {
    vector<float> *curEnvelope;
    unsigned int idxEnv = findEnvelope(monoisotopicMass), curEnvelopeSize;
    float matchScore = 0;

    if (strictMode)
    {
      curEnvelope = new vector<float> ;
      makeStrict(idxEnv, *curEnvelope);
      curEnvelopeSize = curEnvelope->size() - 1;
    }
    else
    {
      curEnvelope = &envelopes[idxEnv];
      curEnvelopeSize = envelopeSize;
    }

    //cerr<<"curEnvelope: "; for(unsigned int i=0; i<curEnvelope->size(); i++) cerr<<(*curEnvelope)[i]<<" "; cerr<<"\n"; cerr.flush();

    if (massEnvelope.size() != curEnvelopeSize)
    {
      cerr
          << "ERROR in ScoreEnvelope (ms1.cpp): Envelope dimensions don't match!\n";
      exit(-1);
    }
    for (unsigned int pivot = 0; pivot < curEnvelopeSize; pivot++)
      matchScore += massEnvelope[pivot] * log(massEnvelope[pivot]
          / (*curEnvelope)[pivot + 1]);

    if (strictMode)
      delete curEnvelope;

    return matchScore;
  }

  void IsoEnvelope::ExtractEnvelope(float monoisotopicMass,
                                    unsigned short charge,
                                    Spectrum &spec,
                                    float peakWidth,
                                    vector<float> &intensities,
                                    bool strictMode)
  {
    if (strictMode)
      intensities.resize(envelopeSize + envelopeSize);
    else
      intensities.resize(envelopeSize);
    for (unsigned int pivot = 0; pivot < intensities.size(); pivot++)
      intensities[pivot] = 0;

    // Find monoisotopic mass and get its summed intensity
    float increment = AAJumps::massIsotopeSpace / ((float)charge);
    vector<int> matches;
    int curSpecIdx, massEnvIdx = 0;
    if (strictMode)
    {
      // Collect intensity between monoisotopicMass-increment-peakWidth and monoisotopicMass-peakWidth (should total roughly zero for the correct Z/PM)
      curSpecIdx = spec.findClosest(monoisotopicMass - increment - peakWidth);
      if (curSpecIdx < 0)
        return;
      while (curSpecIdx > 0 and spec[curSpecIdx][0] >= monoisotopicMass
          - increment - peakWidth)
        curSpecIdx--;
      curSpecIdx++;
      for (; curSpecIdx < (int)spec.size() and spec[curSpecIdx][0]
          < monoisotopicMass - peakWidth; curSpecIdx++)
        intensities[massEnvIdx] += spec[curSpecIdx][1];
      massEnvIdx++;
    }
    else
    {
      curSpecIdx = spec.findClosest(monoisotopicMass);
      if (curSpecIdx < 0)
        return;
    }
    spec.findMatches(monoisotopicMass, peakWidth, matches, -1.0);
    for (unsigned int pivot = 0; pivot < matches.size(); pivot++)
      intensities[massEnvIdx] += spec[matches[pivot]][1];
    massEnvIdx++;
    if (matches.size() > 0)
      curSpecIdx = matches[matches.size() - 1] + 1; // Set curSpecIdx to first peak with
    else if (spec[curSpecIdx][0] < monoisotopicMass)
      curSpecIdx++; //  mass > monoisotopicMass+peakWidth

    // Get intensities for the remaining peaks in the envelope
    float nextMass = monoisotopicMass + increment;
    unsigned int envIdx;
    for (envIdx = 1; envIdx < envelopeSize; envIdx++)
    {
      // Advance to the next peak in the envelope
      for (; curSpecIdx < (int)spec.size() and spec[curSpecIdx][0] < nextMass
          - peakWidth; curSpecIdx++)
        if (strictMode)
          intensities[massEnvIdx] += spec[curSpecIdx][1]; // Add the intensities of the peaks between the current and next isotopic peaks
      if (strictMode)
        massEnvIdx++;

      // Get the intensity for the next peak in the envelope
      for (; curSpecIdx < (int)spec.size() and spec[curSpecIdx][0] <= nextMass
          + peakWidth; curSpecIdx++)
        intensities[massEnvIdx] += spec[curSpecIdx][1];
      massEnvIdx++;
      nextMass += increment;
    }
  }

  void IsoEnvelope::ExtractEnvelope(float monoisotopicMass,
                                    unsigned short charge,
                                    Spectrum &spec,
                                    set<unsigned int>& ignorePeaks,
                                    vector<float> &intensities,
                                    vector<list<unsigned int> >& indicesUsed)
  {
    set<unsigned int> used(ignorePeaks);
    intensities.resize(envelopeSize);
    indicesUsed.resize(envelopeSize);

    for (unsigned int pivot = 0; pivot < intensities.size(); pivot++)
    {
      intensities[pivot] = 0;
      indicesUsed[pivot].clear();
    }

    float increment = AAJumps::massIsotopeSpace / ((float)charge);
    list<int> matches;
    unsigned int massEnvIdx = 0;

    int peakIdx = spec.findPeaks(monoisotopicMass, 0, &matches);

    if (peakIdx < 0)
    {
      return;
    }

    float tolerance = spec.getTolerance(peakIdx);

    for (list<int>::iterator pivot = matches.begin(); pivot != matches.end(); pivot++)
    {
      if (used.count(*pivot) > 0 || abs(spec[*pivot][0] - (monoisotopicMass
          + increment)) < abs(spec[*pivot][0] - monoisotopicMass))
      {
        continue;
      }
      intensities[massEnvIdx] += spec[*pivot][1];
      indicesUsed[massEnvIdx].push_back(*pivot);
      used.insert(*pivot);
    }

    if (indicesUsed[massEnvIdx].size() == 0)
    {
      return;
    }

    float nextMass = monoisotopicMass + increment;
    for (massEnvIdx = 1; massEnvIdx < envelopeSize; massEnvIdx++)
    {
      peakIdx = spec.findPeaks(nextMass, tolerance, &matches);
      if (peakIdx < 0)
      {
        nextMass += increment;
        continue;
      }
      for (list<int>::iterator pivot = matches.begin(); pivot != matches.end(); pivot++)
      {
        if (used.count(*pivot) > 0 || (massEnvIdx < envelopeSize - 1
            && abs(spec[*pivot][0] - (nextMass + increment))
                < abs(spec[*pivot][0] - nextMass)))
        {
          continue;
        }
        intensities[massEnvIdx] += spec[*pivot][1];
        indicesUsed[massEnvIdx].push_back(*pivot);
        used.insert(*pivot);
      }
      nextMass += increment;
    }
  }

  void IsoEnvelope::normalize(vector<float> &values,
                              bool avoidZeros,
                              unsigned int startIdx)
  {
    float totalIntensity = 0;
    for (unsigned int pivot = startIdx; pivot < values.size(); pivot++)
    {
      if (avoidZeros and values[pivot] <= 0.0001)
        values[pivot] = 0.0001;
      totalIntensity += values[pivot];
    }
    for (unsigned int pivot = startIdx; pivot < values.size(); pivot++)
      values[pivot] /= totalIntensity;
  }

}

