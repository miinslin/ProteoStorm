/*
 * DeconvSpectrum.h
 *
 *  Created on: Jun 13, 2011
 *      Author: aguthals
 */

#ifndef DECONVSPECTRUM_H_

#include "Logger.h"
#include "IsoEnvelope.h"

#include <vector>
#include <list>
//#include <set>

using namespace std;

namespace specnets
{
  struct DeconvPeak
  {
    int charge; //peak's assigned charge, 0 if no assignment
    float score; //confidence score of charge assignment
    float binnedIntensity; //summed intensity of all peaks in this peak's isotopic envelope
    set<int> peakEnvelope; //all peak indices in this peak's isotopic envelope
  };

  class DeconvSpectrum : public Spectrum
  {
  public:

    /**
     * Computes the singly charged monoisotopic mass of a given peak mass and charge
     * @param mass peak mass in Da
     * @param charge peak charge
     * @return singly charged monoisotopic mass
     */
    static float GetMonoisotopicMass(float mass, int charge)
    {
      if (charge < 1)
      { // currently only works for charges >= 1
        ERROR_MSG("Computing the monoisopic mass for a peak with an invalid charge " << charge);
      }
      return (mass * ((float)charge)) - (((float)(charge - 1))
          * AAJumps::massHion);
    }

    // charge assignments parallel to Spectrum peaks
    vector<DeconvPeak> scoredEnvelopes;

    /**
     * Default constructor, implicitly calls Spectrum constructor
     */
    DeconvSpectrum(void);

    /**
     * Copy constructor, also calls Spectrum copy constructor
     * @param other pointer to template DeconvSpectrum
     */
    DeconvSpectrum(DeconvSpectrum* other);

    /**
     * Assigns this DeconvSpectrum to a copy of other
     * @param other template DeconvSpectrum
     */
    DeconvSpectrum &operator=(const DeconvSpectrum &other);

    /**
     * Assigns this DeconvSpectrum to a copy of other
     * @param other template Spectrum
     */
    DeconvSpectrum &operator=(const Spectrum &other);

    /**
     * Resizes scoredEnvelopes and calls Spectrum's resize method
     * @param newSize new requested size of DeconvSpectrum
     * @return size of modified DeconvSpectrum
     */
    unsigned int resize(unsigned int newSize);

    /**
     * Assigns peak charges based on the KL divergence between observed and expected isotopic envelopes.
     *   The score of each charge assignment is 1/KL divergence
     * @param env pointer to IsoEnvelope instance holding an expected isotopic envelope for each charge
     * @param theshold maximum allowable divergence between observed and expected isotopic envelopes in order to assign a charge
     * @return
     */
    void AssignChargesKLDiv(IsoEnvelope* env, float threshold);

    /**
     * Converts all peaks to charge one. Any peak that was assigned a charge is converted to a
     *   charge one peak in the new spectrum. If a peak was not assigned a charge, it is also
     *   added to the new spectrum.
     * @param putSpec put the new spectrum here
     * @return
     */
    void ConvertPeaksChargeOne(Spectrum& putSpec);
  };
}

#define DECONVSPECTRUM_H_

#endif /* DECONVSPECTRUM_H_ */
