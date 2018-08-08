#ifndef _SpectrumPair_H_
#define _SpectrumPair_H_

// System Includes
#include <fstream>
#include <string>
#include <vector>

namespace specnets
{
  /**
   * TODO: add description
   */
  class SpectrumPair
  {
  public:

    /**
     * TODO: add description
     *
     */
    SpectrumPair();

    /**
     * TODO: add description
     *
     */
    ~SpectrumPair();

    /**
     * TODO: add description
     *
     *@param other
     *@return
     */
    virtual SpectrumPair &operator=(const SpectrumPair & that);

    /**
     * TODO: add description
     *
     *@paramout
     *@param separator
     */
    void output(std::ostream &out, char separator);

    /**
     * TODO: add description
     *
     *@param v
     */
    void serialize(std::vector<float> &v);

    /**
     * TODO: add description
     *
     *@return
     */
    unsigned int loadSz();

    /**
     * TODO: add description
     *
     *@param v
     */
    void load(std::vector<float> &v);

    /**
     * Index of the first spectrum to align.
     */
    int spec1;
    /**
     * Index of the second spectrum to align.
     */
    int spec2;

    /**
     * Match scores in spectrum1.
     */
    float score1;

    /**
     * Match scores in spectrum2.
     */
    float score2;

    /**
     * Value of the symmetric shift in the Almost Same Peptide
     * problem (first shift is always zero).
     */
    float shift1;

    /**
     * Value of the second shift in the Pairwise Alignment problem
     * (first shift is stored in shift1).
     */
    float shift2;

    /**
     *
     * Value of the second shift in the Pairwise Alignment problem
     * (first shift is stored in shift1).
     */
    int specC;

    /**
     * Boolean value indicating whether spec2 was matched as-is (==false)or
     * reversed(==true).
     */
    bool spec2rev;

  };

} //namespace specnets

#endif // _SpectrumPair_H_
