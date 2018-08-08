#include "mzrange.h"

using namespace std;

namespace specnets
{

  const float PPM_FACTOR = 1.0e-6;

  /**
   * Checks if two 1D values are equal within a range
   * @param center1
   * @param center2
   * @param range
   * @return true if both centers are within radius of each other, false otherwise
   */
  bool MZRange::EqualWithinRange(float center1, float center2, float radius)
  {
    return (center1 >= center2 - radius) && (center1 <= center2 + radius);
  }

  /**
   * Computes the radius of a mass with PPM tolerance
   * @param in_mass mass
   * @param PPM_tol PPM tolerance
   * @return radius of mass error
   */
  float MZRange::GetPPMRadius(float in_mass, float PPM_tol)
  {
    return in_mass * PPM_tol * PPM_FACTOR;
  }

  /**
   * Computes the weighted average of a sequence of MZRanges. Weight of each
   *   MZRange is the intensity.
   * @param ranges sequence of MZRanges.
   * @return a MZRange with weighted average center, weighted average tolerance,
   *   and summed intensity over all in ranges
   */
  MZRange MZRange::WeightedAverage(list<MZRange>* ranges)
  {

    if (ranges->size() == 0) {
      return MZRange(0, 0, 0);
    }

    list<MZRange>::iterator rangeIt = ranges->begin();

    float new_average = 0;
    float new_tolerance = 0;
    float offset = 0;

    float min_score = rangeIt->getIntensity();
    float total_score = min_score;

    rangeIt++;

    for (; rangeIt != ranges->end(); rangeIt++) {

      total_score += rangeIt->getIntensity();
      min_score = min(min_score, rangeIt->getIntensity());
    }

    //check for negative scores
    if (min_score < 0) {
      //make all scores positive
      offset = 0 - min_score;
      total_score += offset * ((float) ranges->size());
    }
    else if (EqualWithinRange(total_score, 0, 0.00001)) {

      total_score = (float) ranges->size();
      for (rangeIt = ranges->begin(); rangeIt != ranges->end(); rangeIt++) {
        new_average += rangeIt->getMass() / total_score;
        new_tolerance += rangeIt->getTolerance() / total_score;
      }

      return MZRange(new_average, 0, new_tolerance);
    }

    //add up weighted average
    for (rangeIt = ranges->begin(); rangeIt != ranges->end(); rangeIt++) {
      new_average += rangeIt->getMass() * ((offset + rangeIt->getIntensity())
          / total_score);
      new_tolerance += rangeIt->getTolerance() * ((offset
          + rangeIt->getIntensity()) / total_score);
    }

    return MZRange(new_average, total_score, new_tolerance);
  }

  /**
   * Inserts a mzrange into existing mzrange_bins
   * @param mzrange_bins each coalesced mzrange mapped to a sequence of mzranges it is derived from
   * @param new_mzrange mzrange to insert
   * @param mergeType same as in MergeMZRanges
   * @param ranges_to_insert specifies a bin of MZRanges to associate with the merged. Used during
   *   recursive calls.
   * @return coalesced MZRange of overlapping ranges
   */
  /*
   MZRange MZRange::InsertMZRange(
   map<MZRange, list<MZRange>*>* mzrange_bins, MZRange& new_mzrange,
   short mergeType, list<MZRange>* ranges_to_insert) {

   if (mzrange_bins->count(new_mzrange) > 0) {

   list<MZRange>* new_bin;

   if (ranges_to_insert == (list<MZRange>*) 0) {
   new_bin = new list<MZRange> ;
   new_bin->push_back(new_mzrange);
   } else {
   new_bin = ranges_to_insert;
   }

   (*mzrange_bins)[new_mzrange] = new_bin;

   return MZRange(new_mzrange);
   }

   list<MZRange>* existing_bin = (*mzrange_bins)[new_mzrange];
   mzrange_bins->erase(new_mzrange);

   if (ranges_to_insert == 0) {

   existing_bin->push_back(new_mzrange);

   //compute average mzrange
   MZRange average_mzrange;
   average_mzrange.MergeMZRanges(existing_bin, mergeType);

   return InsertMZRange(mzrange_bins, average_mzrange, mergeType,
   existing_bin);
   }

   ranges_to_insert->splice(ranges_to_insert->begin(), *existing_bin);
   delete existing_bin;

   //compute average mzrange
   MZRange average_mzrange;
   average_mzrange.MergeMZRanges(ranges_to_insert, mergeType);

   return InsertMZRange(mzrange_bins, average_mzrange, mergeType,
   ranges_to_insert);

   }
   */

  /**
   * Sets the center and radius of this MZRange to the center MZRange of
   *   a sequence of MZRanges.
   * @param peaks sequence of MZRanges
   * @param mergeType specifies how to compute center MZRange
   *   0: MZRange with lowest tolerance is taken as center
   *   1: Take weighted average of centers and tolerances (intensity is the "weight")
   *   2: MZRange with the highest intensity is taken as the center
   *   3: the first MZRange in the list is taken as the center
   * @return
   */
  void MZRange::MergeMZRanges(list<MZRange>* peaks, short mergeType)
  {

    if (peaks->size() == 0) {
      return;
    }

    float new_tolerance, new_mass, new_intensity, highest_intensity;
    list<MZRange>::iterator peakIt;
    MZRange weightedCenter;
    switch (mergeType) {
      case 0:
        peakIt = peaks->begin();
        new_tolerance = peakIt->getTolerance();
        new_mass = peakIt->getMass();
        new_intensity = peakIt->getIntensity();

        peakIt++;
        for (; peakIt != peaks->end(); peakIt++) {
          new_intensity += peakIt->getIntensity();
          if (peakIt->getTolerance() < new_tolerance) {
            new_tolerance = peakIt->getTolerance();
            new_mass = peakIt->getMass();
          }
        }

        set(new_mass, new_intensity, new_tolerance);

      case 1:
        weightedCenter = WeightedAverage(peaks);

        set(weightedCenter.getMass(),
            weightedCenter.getIntensity(),
            weightedCenter.getTolerance());

      case 2:
        peakIt = peaks->begin();
        new_tolerance = peakIt->getTolerance();
        new_mass = peakIt->getMass();
        new_intensity = peakIt->getIntensity();
        highest_intensity = peakIt->getIntensity();

        peakIt++;
        for (; peakIt != peaks->end(); peakIt++) {
          new_intensity += peakIt->getIntensity();
          if (peakIt->getIntensity() > highest_intensity) {
            new_tolerance = peakIt->getTolerance();
            new_mass = peakIt->getMass();
            highest_intensity = peakIt->getTolerance();
          }
        }

        set(new_mass, new_intensity, new_tolerance);

      case 3:
        peakIt = peaks->begin();
        new_intensity = peakIt->getIntensity();
        new_tolerance = peakIt->getTolerance();
        new_mass = peakIt->getMass();

        peakIt++;
        for (; peakIt != peaks->end(); peakIt++) {
          new_intensity += peakIt->getIntensity();
        }

        set(new_mass, new_intensity, new_tolerance);
    }

  }

}
