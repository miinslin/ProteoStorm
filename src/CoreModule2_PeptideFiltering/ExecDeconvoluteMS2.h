/*
 * ExecDeconvoluteMS2.h
 *
 *  Created on: Feb 28, 2011
 *      Author: aguthals
 */

#ifndef EXECDECONVOLUTEMS2_H_
#define EXECDECONVOLUTEMS2_H_

// Module Includes
#include "ExecBase.h"
#include "Logger.h"
#include "FileUtils.h"
#include "ParameterList.h"
#include "ExecMergeConvert.h"

// External Includes
#include "DeconvSpectrum.h"
#include "ms1.h"
#include "utils.h"

// System Includes
#include <string>
#include <vector>

namespace specnets
{
  class ExecMergeConvert;

  class ExecDeconvoluteMS2 : public ExecBase
  {
  public:
    ExecDeconvoluteMS2(void);

    ExecDeconvoluteMS2(const ParameterList & inputParams);

    ExecDeconvoluteMS2(const ParameterList & inputParams,
                       SpecSet * inputSpectra,
                       IsoEnvelope * m_inputEnv);

    ExecDeconvoluteMS2(const ParameterList & inputParams,
                       SpecSet * inputSpectra,
                       IsoEnvelope * m_inputEnv,
                       SpecSet * outputSpectra);

    virtual ~ExecDeconvoluteMS2(void);

    virtual void setOwnInput(bool _ownInput);

    virtual ExecBase * clone(const ParameterList & input_params) const;

    virtual bool invoke(void);

    virtual bool loadInputData(void);

    virtual bool saveOutputData(void);

    virtual bool saveInputData(std::vector<std::string> & filenames);

    virtual bool loadOutputData(void);

    virtual std::vector<ExecBase *> const & split(int numSplit);

    virtual bool merge(void);

    virtual bool validateParams(std::string & error);

  private:
    SpecSet * m_inputSpectra; //! The set of input spectra
    IsoEnvelope * m_inputEnv; //! Input reference isotopic envelope
    SpecSet * m_outputSpectra; //! The set of output spectra

    ExecMergeConvert* m_loader; // Universal loader/converter
    vector<pair<int, int> >* m_loadedIndices;

    bool ownInput; //! Does this object "own" the input data structures (and hence have to free them)
    bool ownOutput; //! Does this object "own" the output data pointers (and hence have to free them)
  };
}

#endif /* EXECDECONVOLUTEMS2_H_ */
