#ifndef _AlignmentUtils_H_
#define _AlignmentUtils_H_

// Module Includes
#include "SpectrumPairSet.h"
#include "ExecTagFilterPairs.h"
// External Module Includes
#include "PepnovoTags.h"
#include "spectrum.h"
#include "SpecSet.h"

// System Includes
#include <list>
#include <vector>

namespace specnets
{

  /**
   * Computes Almost Same Peptide matches between spectra.
   *
   *@param specSet
   *@param baseSpectraIdx
   *@param aaDiff
   *@param minShift
   *@param maxShift
   *@param pmTol
   *@param peakTol
   *@param minMatchRatio
   *@param results
   *@param ratios
   *@param means
   *@param tagsMatchFlank
   *@param tagsMatchCount
   *@param resolution
   *@param symmetryOffset
   */
  void getPairAlignsASP(SpecSet & specSet,
                        vector<int> & baseSpectraIdx,
                        short aaDiff,
                        float minShift,
                        float maxShift,
                        float pmTol,
                        float peakTol,
                        float minMatchRatio,
                        short minNumMatchedPeaks,
                        SpectrumPairSet &results,
                        std::vector<TwoValues<float> > & ratios,
                        vector<TwoValues<float> > & means,
                        vector<float> & varTerms,
                        ExecTagFilterPairs *tag_filter,
                        float resolution = 0.1,
                        float symmetryOffset = 0);

  /**
   *
   * Computes Pairwise Alignments between spectra (TODO: NEEDS DEBUGGING).
   *
   *@param specSet
   *@param pmTol
   *@param peakTol
   *@param results
   */
  void getPairAlignsPA(SpecSet & specSet,
                       float pmTol,
                       float peakTol,
                       SpectrumPairSet & results);

  /**
   * Computes Pairwise Alignments between spectra (TODO: NEEDS DEBUGGING).
   *
   *@param specSet
   *@param startIdx
   *@param endIdx
   *@param pmTol
   *@param minRatio
   *@param minPeakAreaOvlp
   *@param minNumMatchedPeaks
   *@param allowedJumps
   *@param minAbsShift
   *@param results
   *@param ratios
   *@param numMatchedPeaks
   *@param means
   *@param varTerms
   *@param alignStats
   *@param specStats
   */
  void getPairAlignsPA2(SpecSet &specSet,
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
                        std::vector<TwoValues<float> > & ratios,
                        list<TwoValues<int> > & numMatchedPeaks,
                        vector<TwoValues<float> > & means,
                        vector<float> &varTerms,
                        list<vector<float> > & alignStats,
                        vector<vector<float> > & specStats,
                        float resolution,
                        float maxShift);

  /**
    * Finds pairs of spectra with high cosines after correcting peak masses
    * for the parent mass differences.
    *
    *@param specSet
    *@param baseSpectraIdx
    *@param aaDiff
    *@param minShift
    *@param maxShift
    *@param pmTol
    *@param peakTol
    *@param minMatchRatio
    *@param results
    *@param ratios
    *@param means
    *@param tagsMatchFlank
    *@param tagsMatchCount
    *@param resolution
    *@param symmetryOffset
    */
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
                       float resolution = 0.1,
                       float symmetryOffset = 0,
                       bool projectedCosine = false);

   void getPairAlignGFPValues(SpecSet & specSet,
      		  	  	  	  	  SpectrumPairSet & spectraPairs,
      		  	  	  	  	  vector<TwoValues<double> > & pvalues,
      		  	  	  	  	  float pmTol,
                              float peakTol,
                              bool specTypeMSMS = false);

   void getPairAlignGFPValues22(SpecSet & specSet,
         		  	  	  	  	SpectrumPairSet & spectraPairs,
         		  	  	  	  	vector<TwoValues<double> > & pvalues,
         		  	  	  	  	float pmTol,
                                float peakTol,
                                AAJumps & jumps,
                                bool specTypeMSMS = false);

   void getSingleModPairs(SpecSet & prmSpec,//No partial overlap, align-gf, tag-alignment
		                  SpecSet & normMSSpec,//for cosine
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
 						  float symmetryOffset);

   void getPartialOverlapsPairs(SpecSet & prmSpec,//partial overlap, align-gf, tag-alignment
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
   								float maxShift);

} //namespace specnets

#endif // _AlignmentUtils_H_
