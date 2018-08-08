/*
 * DeconvSpectrum.cpp
 *
 *  Created on: Jun 13, 2011
 *      Author: aguthals
 */

#include "DeconvSpectrum.h"

namespace specnets
{

  DeconvSpectrum::DeconvSpectrum(void)
  {
    scoredEnvelopes.resize(0);
  }

  DeconvSpectrum::DeconvSpectrum(DeconvSpectrum* other) :
    Spectrum((Spectrum)(*other))
  {
    scoredEnvelopes = other->scoredEnvelopes;
  }

  DeconvSpectrum& DeconvSpectrum::operator=(const DeconvSpectrum &other)
  {
    Spectrum::operator =((Spectrum&)other);
    scoredEnvelopes = other.scoredEnvelopes;
    return *this;
  }

  DeconvSpectrum& DeconvSpectrum::operator=(const Spectrum &other)
  {
    Spectrum::operator =(other);
    return *this;
  }

  unsigned int DeconvSpectrum::resize(unsigned int newSize)
  {
    Spectrum::resize(newSize);
    scoredEnvelopes.resize(newSize);
    return size();
  }

  /**
   * Assigns peak charges based on the KL divergence between observed and expected isotopic envelopes.
   *   The score of each charge assignment is 1/KL divergence
   * @param env pointer to IsoEnvelope instance holding an expected isotopic envelope for each charge
   * @param theshold maximum allowable divergence between observed and expected isotopic envelopes in order to assign a charge
   * @return
   */
  void DeconvSpectrum::AssignChargesKLDiv(IsoEnvelope* env, float threshold)
  {
    scoredEnvelopes.resize(size());

    vector<list<unsigned int> > indicesUsed;
    vector<float> massEnvelope;
    float mass, monoisotopicmass, intensity, divergence;
    set<unsigned int> assignedPeaks;
    for (int index = 0; index < size(); index++)
    {
      mass = peakList[index][0];
      intensity = peakList[index][1];

      if (assignedPeaks.count(index) > 0)
      {
        continue;
      }

      // clear out old assigned charges
      scoredEnvelopes[index].charge = 0;
      scoredEnvelopes[index].score = 1.1;
      scoredEnvelopes[index].binnedIntensity = 0.0;
      scoredEnvelopes[index].peakEnvelope.clear();

      for (int charge = parentCharge; charge > 1; charge--)
      {
        monoisotopicmass = GetMonoisotopicMass(mass, charge);

        // fragment peak can't be greater than parent mass
        if (monoisotopicmass >= parentMass)
        {
          continue;
        }
        massEnvelope.resize(0);
        indicesUsed.resize(0);

        // get the binned peak intensities (1/charge apart) that might match this envelope
        env->ExtractEnvelope(mass,
                             charge,
                             *this,
                             assignedPeaks,
                             massEnvelope,
                             indicesUsed);

        // ignore this envelope if there are no peaks in the first or second bin
        if (massEnvelope[0] == 0 || massEnvelope[1] == 0)
        {
          continue;
        }

        // loop over subsets of the envelope
        for (int i = massEnvelope.size(); i > 1; i--)
        {

          if (i < massEnvelope.size() && indicesUsed[i].size() == 0)
          {
            continue;
          }
          if (i < massEnvelope.size())
          {
            massEnvelope[i] = 0.0;
            indicesUsed[i].clear();
          }

          // try not to include a peak into the envelope with a relatively high intensity
          bool badcharge = false;
          for (int j = 0; j < indicesUsed.size(); j++)
          {
            for (list<unsigned int>::iterator indiciesIt =
                indicesUsed[j].begin(); indiciesIt != indicesUsed[j].end(); indiciesIt++)
            {
              if (peakList[*indiciesIt][1] > 12.0 * intensity)
              {
                badcharge = true;
                break;
              }
            }
            if (badcharge)
              break;
          }
          if (badcharge)
          {
            continue;
          }

          // normalize and score envelope
          vector<float> normalizedEnvelope(massEnvelope);
          env->normalize(normalizedEnvelope, true);
          divergence = env->ScoreEnvelope(monoisotopicmass,
                                          normalizedEnvelope,
                                          false);

          // overwrite previous envelopes for this peak with higher divergence
          if (divergence < scoredEnvelopes[index].score && divergence
              <= threshold)
          {
            scoredEnvelopes[index].charge = charge;
            scoredEnvelopes[index].score = divergence;
            scoredEnvelopes[index].peakEnvelope.clear();
            for (int j = 0; j < indicesUsed.size(); j++)
            {
              scoredEnvelopes[index].peakEnvelope.insert(indicesUsed[j].begin(),
                                                         indicesUsed[j].end());
            }

          }
        }
      }

      // if this peak was assigned an envelope, add up cummulative intensity and cascade assignment to all peaks in the isotopic envelope
      if (scoredEnvelopes[index].charge > 0)
      {
        //fprintf(stdout, "Charge at %.6f is %d with score %.6f : ", mass, scoredEnvelopes[index].charge, scoredEnvelopes[index].score);

        DeconvPeak* thisPeak = &scoredEnvelopes[index];
        for (set<int>::iterator peakIt = thisPeak->peakEnvelope.begin(); peakIt
            != thisPeak->peakEnvelope.end(); peakIt++)
        {
          scoredEnvelopes[*peakIt].charge = thisPeak->charge;
          scoredEnvelopes[*peakIt].score = thisPeak->score;
          scoredEnvelopes[*peakIt].peakEnvelope = thisPeak->peakEnvelope;

          //fprintf(stdout, "%.6f, ", peakList[*peakIt][0]);

          // don't consider isotopic peaks for membership in another isotopic envelope
          assignedPeaks.insert(*peakIt);
          scoredEnvelopes[index].binnedIntensity += peakList[*peakIt][1];
        }
        //fprintf(stdout, "\n");
        for (set<int>::iterator peakIt = thisPeak->peakEnvelope.begin(); peakIt
            != thisPeak->peakEnvelope.end(); peakIt++)
        {
          scoredEnvelopes[*peakIt].binnedIntensity
              = scoredEnvelopes[index].binnedIntensity;
        }
      }
    }
  }

  /**
   * Converts all peaks to charge one. Any peak that was assigned a charge is converted to a
   *   charge one peak in the new spectrum. If a peak was not assigned a charge, it is also
   *   added to the new spectrum.
   * @param putSpec put the new spectrum here. If 0, modify this spectrum
   * @return
   */
  void DeconvSpectrum::ConvertPeaksChargeOne(Spectrum& putSpec)
  {
    putSpec = *this;
    DeconvPeak* thisPeak;

    putSpec.resize(size());
    int peakIdxUse = 0;

    set<int> assignedPeaks;
    for (int i = 0; i < size(); i++)
    {
      if (assignedPeaks.count(i) > 0)
      {
        continue;
      }
      thisPeak = &(scoredEnvelopes[i]);

      // add un-assigned peaks to updated spectrum
      if (thisPeak->charge == 0)
      {
        putSpec[peakIdxUse] = peakList[i];
        putSpec.setTolerance(peakIdxUse, peakTols[i]);
        ++peakIdxUse;
      }
      else
      { // or add peaks that have assigned charge
        putSpec[peakIdxUse][0] = GetMonoisotopicMass(peakList[i][0],
                                                     thisPeak->charge);
        putSpec[peakIdxUse][1] = thisPeak->binnedIntensity;
        putSpec.setTolerance(peakIdxUse, peakTols[i]
            * ((float)thisPeak->charge));
        assignedPeaks.insert(thisPeak->peakEnvelope.begin(),
                             thisPeak->peakEnvelope.end());

        // if the monoisotopic peak is already present in the spectrum, just merge the peaks;
        float pkTolUse = putSpec.getTolerance(peakIdxUse);
        int closestIdx = findClosest(putSpec[peakIdxUse][0]);
        if (MZRange::EqualWithinRange(peakList[closestIdx][0],
                                      putSpec[peakIdxUse][0],
                                      pkTolUse))
        {
          putSpec[peakIdxUse][0] = (peakList[closestIdx][0]
              + putSpec[peakIdxUse][0]) / 2.0;
          putSpec[peakIdxUse][1] += peakList[closestIdx][1];
          assignedPeaks.insert(closestIdx);
        }

        ++peakIdxUse;
      }
    }
    putSpec.resize(peakIdxUse);
    putSpec.sortPeaks();
  }
}
