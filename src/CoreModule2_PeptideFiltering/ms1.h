#ifndef MS1_H
#define MS1_H

#include "spectrum.h"
#include "SpecSet.h"
#include "IsoEnvelope.h"
#include <set>

namespace specnets
{
  using namespace std;

  /**
   * TODO: add description
   *
   *@param spec1
   *@param spec2
   *@return
   */
  float BinnedMatchProduct(Spectrum &spec1, Spectrum &spec2);

  /**
   * TODO: add description
   *
   *@param run1
   *@param run2
   *@param peakTol
   *@param matchedScans
   *@param binSpectra
   */
  void AlignChromatography(SpecSet &run1, SpecSet &run2, float peakTol, vector<
      TwoValues<unsigned int> > &matchedScans, bool binSpectra = true);

}
#endif
