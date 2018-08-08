#ifndef _SpectrumPairSet_H_
#define _SpectrumPairSet_H_

// Module Includes
#include "SpectrumPair.h"
#include "SpecSet.h"

// System Includes
#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace specnets
{
  /**
   * TODO: add description
   */
  class SpectrumPairSet
  {
  public:
    SpectrumPairSet(void);
    SpectrumPairSet(unsigned int nsize);

    void copy(SpectrumPairSet & that);

    int loadFromBinaryFile(const std::string & filename);
    bool saveToBinaryFile(const std::string & filename);

    unsigned int size(void) const;
    void resize(unsigned int newSize);
    SpectrumPair const & operator[](unsigned int index) const;
    SpectrumPair & operator[](unsigned int index);
    void push_back(const SpectrumPair & newPair);

    void sort_pairs();
    void sort_pairs_by_index();
    void sort_descending_by_score();
    
    bool filter_by_component_size(unsigned int max_component_size);
    bool filter_by_unique_mass_number_in_component(SpecSet & specSet, unsigned int max_unique_mass_number);
    bool filter_by_max_spec_per_variant(SpecSet & specSet, unsigned int max_per_variant);

    bool filter_by_edge_fdr22(SpecSet & specSet, PeptideSpectrumMatchSet & psms, float fdrCut, float pmTol, float fragTol, bool canonicalPSMForm);
    bool filter_by_edge_fdr33(SpecSet & specSet, PeptideSpectrumMatchSet & psms, float fdrCut, float pmTol, float fragTol, const std::string & filename, bool canonicalPSMForm);

    bool getSpecComponentID(SpecSet & specSet, vector<int> & specCompID, vector<int> & compSize);

    bool getModificationFrequencies(float resolution, 
                                    float maxDiffMass, 
                                    std::map<float, float> & modFreqs);

    float getFDR() { return m_fdr; }

    bool cosine_precision_recall(SpecSet & specSet, PeptideSpectrumMatchSet & psms, const std::string & filename, float maxDelta, float fragTol, bool seperateCharge, bool canonicalPSMForm);
    bool cosine_precision_recall(PeptideSpectrumMatchSet & psmSet, const std::string & filename, float minOverlap, bool seperateCharge, bool canonicalPSMForm);
    bool cosine_precision_recall_curve(PeptideSpectrumMatchSet & psmSet, float minOverlap, bool seperateCharge, bool canonicalPSMForm);

  protected:
    std::vector<SpectrumPair> thePairs;
    float m_fdr;
  };

} //namespace specnets

#endif // _SpectrumPairSet_H_
