/*
 * ExecDeconvoluteMS2.cpp
 *
 *  Created on: Feb 28, 2011
 *      Author: aguthals
 */

// Header includes
#include "ExecDeconvoluteMS2.h"

using namespace specnets;
using namespace std;

namespace specnets
{

  ExecDeconvoluteMS2::ExecDeconvoluteMS2(void) :
    m_inputSpectra(0x0), m_inputEnv(0x0), m_outputSpectra(0x0), ownInput(true),
        ownOutput(true), m_loader(0x0), m_loadedIndices(0x0)
  {
    m_name = "ExecDeconvoluteMS2";
    m_type = "ExecDeconvoluteMS2";
  }

  ExecDeconvoluteMS2::ExecDeconvoluteMS2(const ParameterList & inputParams) :

    ExecBase(inputParams), m_inputSpectra(0x0), m_inputEnv(0x0),
        m_outputSpectra(0x0), ownInput(true), ownOutput(true), m_loader(0x0),
        m_loadedIndices(0x0)
  {

    m_name = "ExecDeconvoluteMS2";
    m_type = "ExecDeconvoluteMS2";
  }

  ExecDeconvoluteMS2::ExecDeconvoluteMS2(const ParameterList & inputParams,
                                         SpecSet * inputSpectra,
                                         IsoEnvelope * m_inputEnv) :
    ExecBase(inputParams), m_inputSpectra(inputSpectra),
        m_inputEnv(m_inputEnv), m_outputSpectra(0x0), ownInput(false),
        ownOutput(true), m_loader(0x0), m_loadedIndices(0x0)
  {

    m_name = "ExecDeconvoluteMS2";
    m_type = "ExecDeconvoluteMS2";
  }

  ExecDeconvoluteMS2::ExecDeconvoluteMS2(const ParameterList & inputParams,
                                         SpecSet * inputSpectra,
                                         IsoEnvelope * m_inputEnv,
                                         SpecSet * outputSpectra) :

    ExecBase(inputParams), m_inputSpectra(inputSpectra),
        m_inputEnv(m_inputEnv), m_outputSpectra(outputSpectra),
        ownInput(false), ownOutput(false), m_loader(0x0), m_loadedIndices(0x0)
  {

    m_name = "ExecDeconvoluteMS2";
    m_type = "ExecDeconvoluteMS2";
  }

  ExecDeconvoluteMS2::~ExecDeconvoluteMS2(void)
  {
    if (ownInput)
    {
      if (m_inputSpectra != 0)
      {
        delete m_inputSpectra;
      }
      if (m_inputEnv != 0)
      {
        delete m_inputEnv;
      }
    }
    if (ownOutput)
    {
      if (m_outputSpectra != 0)
      {
        delete m_outputSpectra;
      }
    }

    if (m_loader != 0)
    {
      delete m_loader;
    }

    if (m_loadedIndices != 0)
    {
      delete m_loadedIndices;
    }
  }

  void ExecDeconvoluteMS2::setOwnInput(bool _ownInput)
  {
    ownInput = _ownInput;
  }

  ExecBase * ExecDeconvoluteMS2::clone(const ParameterList & inputParams) const
  {
    return new ExecDeconvoluteMS2(inputParams);
  }

  bool ExecDeconvoluteMS2::invoke(void)
  {

    if (m_inputSpectra == 0)
    {
      ERROR_MSG("No input spectra!");
      return false;
    }

    if (ownOutput && m_outputSpectra == 0)
    {
      m_outputSpectra = new SpecSet;
    }

    if (m_loader == 0)
    {
      m_loadedIndices = new vector<pair<int, int> > ;
      m_loader = new ExecMergeConvert(m_params,
                                      m_loadedIndices,
                                      m_inputSpectra,
                                      m_outputSpectra);
    }

    if (!m_loader->invoke())
    {
      return false;
    }

    float threshold = m_params.getValueFloat("MAX_KLDiv", 0.50);

    int start = (m_params.exists("IDX_START"))
        ? m_params.getValueInt("IDX_START") : 0;
    int end = (m_params.exists("IDX_END")) ? m_params.getValueInt("IDX_END")
        : m_inputSpectra->size() - 1;

    DEBUG_VAR(m_outputSpectra->size());
    DEBUG_MSG("Deconvoluting spectra with maximum KL-divergence = " << threshold);

    DeconvSpectrum workingSpec;
    for (int i = start; i <= end; i++)
    {
      workingSpec = (*m_outputSpectra)[i];
      workingSpec.AssignChargesKLDiv(m_inputEnv, threshold);
      workingSpec.ConvertPeaksChargeOne((*m_outputSpectra)[i]);
      //(*m_outputSpectra)[i].output(cout);
      //break;
    }
    return true;
  }

  bool ExecDeconvoluteMS2::loadInputData(void)
  {
    if (ownInput && m_inputSpectra == 0)
    {
      m_inputSpectra = new SpecSet();
    }
    if (ownOutput && m_outputSpectra == 0)
    {
      m_outputSpectra = new SpecSet;
    }

    if (m_loader == 0)
    {
      m_loadedIndices = new vector<pair<int, int> > ;
      m_loader = new ExecMergeConvert(m_params,
                                      m_loadedIndices,
                                      m_inputSpectra,
                                      m_outputSpectra);
    }
    m_inputEnv = new IsoEnvelope();

    if (!m_loader->loadInputData())
    {
      return false;
    }

    if (m_inputSpectra->size() == 0)
    {
      ERROR_MSG("Input spectra size is 0!, did you specify INPUT_SPECTRA?");
      return false;
    }

    if (m_params.exists("INPUT_ISO_ENV"))
    {
      DEBUG_MSG("Loading input iso envelope from <" << m_params.getValue("INPUT_ISO_ENV") << "> ...");
      if (!m_inputEnv->LoadModel(m_params.getValue("INPUT_ISO_ENV").c_str()))
      {
        ERROR_MSG("Could not load " << m_params.getValue("INPUT_ISO_ENV"));
        return false;
      }
    }
    else if (m_params.exists("EXE_DIR"))
    {
      string exeDir = m_params.getValue("EXE_DIR");
      string rsrcDir = getPath(exeDir, "resources", false);
      string modelPath = getPath(rsrcDir, "model_isoenv.bin", false);

      DEBUG_MSG("Loading input iso envelope from <" << modelPath << "> ...");
      if (!m_inputEnv->LoadModel(modelPath.c_str()))
      {
        ERROR_MSG("Could not load " << modelPath);
        return false;
      }
    }
    else
    {
      ERROR_MSG("Must specify INPUT_ISO_ENV");
      return false;
    }

    return true;
  }

  bool ExecDeconvoluteMS2::saveInputData(std::vector<std::string> & filenames)
  {
    //SpecSet m_inputContigs; // the input spectra
    std::string spectraFilename = getName() + "_spectra.pklbin";
    m_params.setValue("INPUT_SPECTRA_PKLBIN", spectraFilename);
    if (!fileExists(spectraFilename))
    {
      m_inputSpectra->savePklBin(spectraFilename.c_str());
    }

    // No method for saving IsoEnvelope

    // Have to set up the output files also so the params will be correct on reload
    m_params.setValue("OUTPUT_SPECTRA_PKLBIN", getName() + "_spectra_z.pklbin");

    std::string paramFilename = getName() + ".params";
    m_params.writeToFile(paramFilename);

    filenames.push_back(paramFilename); // Parameter file MUST be first in vector
    filenames.push_back(spectraFilename);

    return true;
  }

  bool ExecDeconvoluteMS2::saveOutputData(void)
  {
    if (ownOutput && m_outputSpectra == 0)
    {
      m_outputSpectra = new SpecSet;
    }

    if (m_loader == 0)
    {
      m_loadedIndices = new vector<pair<int, int> > ;
      m_loader = new ExecMergeConvert(m_params,
                                      m_loadedIndices,
                                      m_inputSpectra,
                                      m_outputSpectra);
    }

    if (!m_loader->saveOutputData())
    {
      return false;
    }
    return true;
  }

  bool ExecDeconvoluteMS2::loadOutputData(void)
  {
    if (m_params.exists("OUTPUT_SPECTRA"))
    {
      if (!ExecMergeConvert::loadSpecset(m_params.getValue("EXE_DIR"),
                                         m_params.getValue("OUTPUT_SPECTRA"),
                                         m_outputSpectra))
      {
        return false;
      }
    }
    return true;
  }

  std::vector<ExecBase *> const & ExecDeconvoluteMS2::split(int numSplit)
  {

    m_subModules.resize(0);

    if (numSplit < 2)
    {
      DEBUG_MSG("Number split [" << numSplit << "] must be at least 2");
      return m_subModules;
    }

    if (m_inputSpectra->size() == 0)
    {
      DEBUG_MSG("Must have at least one spectrum");
      return m_subModules;
    }

    m_subModules.resize(numSplit);

    long num_ops = 0;
    for (int i = 0; i < m_inputSpectra->size(); i++)
    {
      num_ops += (*m_inputSpectra)[i].size();
    }

    long numOpsPerChild = num_ops / ((long)numSplit);
    int globalIdx = 0, childIdx = 0;
    int startIdx, endIdx;

    DEBUG_MSG("Splitting into " << numSplit << " children");
    for (int i = 0; i < numSplit; i++)
    {

      num_ops = 0;
      startIdx = globalIdx;
      endIdx = globalIdx;

      while (globalIdx < m_inputSpectra->size() && (num_ops <= numOpsPerChild
          || i == numSplit - 1))
      {

        num_ops += (*m_inputSpectra)[globalIdx].size();

        ++globalIdx;
      }
      endIdx = globalIdx - 1;

      ParameterList childParams(m_params);
      childParams.setValue("IDX_START", parseInt(startIdx));
      childParams.setValue("IDX_END", parseInt(endIdx));

      ExecBase * theClone = new ExecDeconvoluteMS2(childParams,
                                                   m_inputSpectra,
                                                   m_inputEnv);

      theClone ->setName(makeName(m_name, i));
      m_subModules[i] = theClone;
    }

    DEBUG_MSG("Splitting success");
    DEBUG_TRACE;
    return m_subModules;
  }

  bool ExecDeconvoluteMS2::merge(void)
  {
    if (m_subModules.size() == 0)
    {
      DEBUG_MSG("No children found when merging");
      return false;
    }

    int num_specs = m_inputSpectra->size();

    m_outputSpectra->resize(num_specs);
    int specIdx = 0;

    DEBUG_MSG("Merging " << m_subModules.size() << " children");
    for (int child = 0; child < m_subModules.size(); child++)
    {
      ExecDeconvoluteMS2* theChild = (ExecDeconvoluteMS2*)m_subModules[child];
      SpecSet* childSpecs = theChild->m_outputSpectra;

      for (int i = 0; i < childSpecs->size(); i++)
      {
        (*m_outputSpectra)[specIdx] = (*childSpecs)[i];
        ++specIdx;
      }
    }
    m_outputSpectra->resize(specIdx);

    DEBUG_MSG("Merging success");
    return true;
  }

  bool ExecDeconvoluteMS2::validateParams(std::string & error)
  {
    VALIDATE_PARAM_EXIST("INPUT_SPECTRA");
    VALIDATE_PARAM_EXIST("OUTPUT_SPECTRA");
    return true;
  }
}

