/*
 * ExecMergeConvert.cpp
 *
 *  Created on: Oct 5, 2011
 *      Author: aguthals
 */

#include <algorithm>
#include <cmath>

#include "ExecMergeConvert.h"

#define DIR_SEP       '/'
#define DIR_SEP_OTHER '\\'

using namespace std;

namespace specnets
{

  bool ExecMergeConvert::loadSpecsetMultiple(const string &exeDir,
                                             const string& specFilePaths,
                                             SpecSet* spectra,
                                             vector<string>* fileNames,
                                             vector<pair<int, int> >* fileSpecIdxs,
                                             const char* separator)
  {
    list<string> fileList;
    splitText(specFilePaths.c_str(), fileList, separator);
    if (fileNames)
    {
      fileNames->assign(fileList.begin(), fileList.end());
    }
    if (fileSpecIdxs)
    {
      fileSpecIdxs->resize(fileList.size());
    }

    if (fileList.size() == 1)
    {
      if (!loadSpecset(exeDir, fileList.front(), spectra))
      {
        return false;
      }
      if (fileNames)
      {
        (*fileNames)[0] = fileList.front();
      }
      if (fileSpecIdxs)
      {
        (*fileSpecIdxs)[0].first = 0;
        (*fileSpecIdxs)[0].second = spectra->size() - 1;
      }
      return true;
    }

    SpecSet tempSpecs;
    int idxUse = 0;
    unsigned int lastIdx = 0, thisIdx = 0;

    for (list<string>::iterator fileIt = fileList.begin();
        fileIt != fileList.end(); fileIt++)
    {
      if (!loadSpecset(exeDir, *fileIt, &tempSpecs))
      {
        return false;
      }
      if (fileSpecIdxs)
      {
        thisIdx = tempSpecs.size() - 1 + lastIdx;
        (*fileSpecIdxs)[idxUse].first = lastIdx;
        (*fileSpecIdxs)[idxUse].second = thisIdx;

        lastIdx = thisIdx + 1;
        ++idxUse;
      }
      spectra->swapAppendSpecSet(tempSpecs, false);
    }

    return true;
  }

  bool ExecMergeConvert::convertSpecset(const string& exeDir,
                                        const string& inFile,
                                        const string& outFile)
  {
    // convert command
    std::string convertCommand;
    // Path
    convertCommand = exeDir;
    // replace directory separators
    replace(convertCommand.begin(),
            convertCommand.end(),
            DIR_SEP_OTHER,
            DIR_SEP);
    // Executable name
    convertCommand += DIR_SEP;
    convertCommand += "convert";
    // input file
    convertCommand += " \"";
    convertCommand += inFile;
    // outpur file
    convertCommand += "\" \"";
    convertCommand += outFile;
    // Force convertCommand output into a file
    convertCommand += "\" >> .";
    convertCommand += DIR_SEP;
    convertCommand += "log_convert.txt";
    //Send to sdtout DEBUG
    DEBUG_MSG("call convert: " << convertCommand);
    // call the command
    bool res = spsSystem(convertCommand.c_str());

    return res;
  }

  bool ExecMergeConvert::loadSingleSpecset(const string& exeDir,
                                           const string &specFilePath,
                                           SpecSet* spectra,
                                           FilenameManager &fm)
  {
    if (fm.checkExtension(false))
    {
      if (!spectra->Load(specFilePath.c_str(), NULL))
      {
        ERROR_MSG("Failed to load spectra from \'" << specFilePath << "\'");
        //return false;
      }
      return true;
    }

    string outFile = "TEMP__";
    //outFile += parseInt(spsGetPID());
    outFile += parseInt(rand() % 100000);
    outFile += ".pklbin";

    if (convertSpecset(exeDir, specFilePath, outFile))
    {
      DEBUG_TRACE;
      spectra->Load(outFile.c_str(), NULL);
      DEBUG_TRACE;
      int remRes = remove(outFile.c_str());
      //unlink(outFile.c_str());
    }
    else
    {
      ERROR_MSG("convert exited in error for file: " << specFilePath);
      //return false;
    }

    return true;
  }

  bool ExecMergeConvert::loadSpecset(const string &exeDir,
                                     const string &specFilePath,
                                     SpecSet* spectra,
                                     vector<string>* fileNames,
                                     vector<pair<int, int> >* fileSpecIdxs)
  {
    // split the filename into components
    FilenameManager fm(specFilePath.c_str());

    DEBUG_MSG("Loading spectra from \'" << fm.filename << "\' ...");

    if (fm.checkExtension())
    {
      //      if(!spectra->Load(specFilePath.c_str()))
      if (!loadSingleSpecset(exeDir, specFilePath, spectra, fm))
        goto load_fail;

    }
    else
    {
      // try loading a directory of pkl or mic files if all else fails
      vector<string> possibleDirFiles = directoryContents(specFilePath,
                                                          "",
                                                          "",
                                                          true);

      if (possibleDirFiles.size() == 0)
      {
        ERROR_MSG("/'" << specFilePath << "\' is an empty directory, " << "is an un-reachable directory, or is not named by a " << "supported spectrum-set file format (.mgf, .pklbin, .prms, or .pkl");
        goto load_fail;
      }
      sort(possibleDirFiles.begin(), possibleDirFiles.end());

      int idxUse = 0;
      for (int i = 0; i < possibleDirFiles.size(); i++)
      {
        string fullPath = specFilePath + "/" + possibleDirFiles[i];
        if (possibleDirFiles[i].find(".pkl") == string::npos
            and possibleDirFiles[i].find(".mic") == string::npos)
        {
          WARN_MSG("Skipping file \'" << fullPath << "\'");
          continue;
        }

        SpecSet tempSpecs;
        if (!tempSpecs.LoadSpecSet_pkl_mic(fullPath.c_str()))
        {
          ERROR_MSG("Failed to load file \'" << fullPath << "\' in pkl or mic format!");
          goto load_fail;
        }

        spectra->resize(idxUse + tempSpecs.size());
        if (fileNames)
        {
          fileNames->push_back(fullPath);
        }
        if (fileSpecIdxs)
        {
          fileSpecIdxs->resize(fileSpecIdxs->size() + 1);
          (*fileSpecIdxs)[fileSpecIdxs->size() - 1].first = spectra->size();
          (*fileSpecIdxs)[fileSpecIdxs->size() - 1].second = spectra->size()
              + tempSpecs.size() - 1;
        }
        for (int j = 0; j < tempSpecs.size(); j++)
        {
          (*spectra)[idxUse] = tempSpecs[j];
          ++idxUse;
        }
      }
    }

    DEBUG_MSG("Loaded " << spectra->size() << " spectra from \'" << specFilePath << "\'");

    return true;

    load_fail:

    ERROR_MSG("Failed to load from \'" << specFilePath << "\'!");

    return false;

  }

  bool ExecMergeConvert::saveSpecset(const string& specFilePath,
                                     SpecSet* spectra,
                                     bool outputScans)
  {

    if (!spectra)
    {
      DEBUG_MSG("SpecSet object is NULL!");
      return false;
    }

    int ret = spectra->SaveSpecSet(specFilePath.c_str(), outputScans);

    DEBUG_VAR(ret);

    if (ret == 1)
    {
      DEBUG_MSG("Saved " << spectra->size() << " spectra to \'" << specFilePath << "\'");
      return true;
    }

    FilenameManager fm(specFilePath.c_str());
    if (ret == -2)
    {
      ERROR_MSG("Unrecognized file extension \'" << fm.extension << "\'");
    }

    ERROR_MSG("Failed to save to \'" << specFilePath << "\'!");
    return false;
  }

  bool ExecMergeConvert::saveSpecsetMultiple(const string& dirName,
                                             const string& specFilePaths,
                                             SpecSet* spectra,
                                             vector<string>* fileNames,
                                             const char* separator)
  {
    list<string> fileList;
    splitText(specFilePaths.c_str(), fileList, separator);
    if (fileNames) {
      fileNames->assign(fileList.begin(), fileList.end());
    }
    for (list<string>::iterator fileIt = fileList.begin();
        fileIt != fileList.end(); fileIt++)
    {
      string fileName = *fileIt;
      if (!dirName.empty()) {
        fileName = dirName + "/" + fileName;
      }
      if (!saveSpecset(fileName, spectra))
      {
        return false;
      }
    }

    return true;
  }

  bool ExecMergeConvert::saveSpecsetMultiple(const string& specFilePaths,
                                             SpecSet* spectra,
                                             vector<string>* fileNames,
                                             const char* separator)
  {
    return saveSpecsetMultiple("",
                               specFilePaths,
                               spectra,
                               fileNames,
                               separator);
  }

  pair<bool, list<string> > ExecMergeConvert::getFileList(const string& listFilePath)
  {
    BufferedLineReader blr;
    list<string> fileList;

    if (blr.Load(listFilePath.c_str()) <= 0)
    {
      ERROR_MSG("Could not load list of files from \'" << listFilePath << "\'");
      return pair<bool, list<string> >(false, fileList);
    }

    for (int lineIdx = 0; lineIdx < blr.size(); lineIdx++)
    {
      string file = blr.getline(lineIdx);
      file = replaceAll(file, "\n", "");
      if (file.length() > 0)
      {
        fileList.push_back(file);
      }
    }
    blr.reset();

    return pair<bool, list<string> >(true, fileList);
  }

  ExecMergeConvert::ExecMergeConvert(void) :
      m_inputSpecsetIdx(0x0), ownInput(true), m_outputSpecset(0x0),
          m_inputSpecset(0x0), ownOutput(true), m_recordedInput(0x0),
          m_recordedOutput(0x0)
  {
    m_name = "ExecMergeConvert";
    m_type = "ExecMergeConvert";
    m_recordedInput = new vector<pair<string, int> >;
    m_recordedOutput = new vector<pair<string, int> >;
    m_separator = ";";

    if (!initialize())
    {
      abort();
    }
  }

  ExecMergeConvert::ExecMergeConvert(const ParameterList & inputParams) :

      ExecBase(inputParams), m_inputSpecsetIdx(0x0), ownInput(true),
          m_outputSpecset(0x0), ownOutput(true), m_inputSpecset(0x0),
          m_recordedInput(0x0), m_recordedOutput(0x0)
  {

    m_name = "ExecMergeConvert";
    m_type = "ExecMergeConvert";
    m_recordedInput = new vector<pair<string, int> >;
    m_recordedOutput = new vector<pair<string, int> >;
    m_separator = ";";

    if (!initialize())
    {
      abort();
    }
  }

  ExecMergeConvert::ExecMergeConvert(const ParameterList & inputParams,
                                     vector<pair<int, int> >* inputSpecsetIdx,
                                     SpecSet * inputSpecset) :

      ExecBase(inputParams), m_inputSpecsetIdx(inputSpecsetIdx),
          m_inputSpecset(inputSpecset), ownInput(false), m_outputSpecset(0x0),
          ownOutput(true), m_recordedInput(0x0), m_recordedOutput(0x0)
  {

    m_name = "ExecMergeConvert";
    m_type = "ExecMergeConvert";
    m_recordedInput = new vector<pair<string, int> >;
    m_recordedOutput = new vector<pair<string, int> >;
    m_separator = ";";

    if (!initialize())
    {
      abort();
    }
  }

  ExecMergeConvert::ExecMergeConvert(const ParameterList & inputParams,
                                     SpecSet * outputSpecset) :
      ExecBase(inputParams), m_inputSpecsetIdx(0x0), m_inputSpecset(0x0),
          ownInput(true), m_outputSpecset(outputSpecset), ownOutput(false),
          m_recordedInput(0x0), m_recordedOutput(0x0)
  {

    m_name = "ExecMergeConvert";
    m_type = "ExecMergeConvert";
    m_recordedInput = new vector<pair<string, int> >;
    m_recordedOutput = new vector<pair<string, int> >;
    m_separator = ";";

    if (!initialize())
    {
      abort();
    }
  }

  ExecMergeConvert::ExecMergeConvert(const ParameterList & inputParams,
                                     vector<pair<int, int> >* inputSpecsetIdx,
                                     SpecSet * inputSpecset,
                                     SpecSet * outputSpecset) :

      ExecBase(inputParams), m_inputSpecsetIdx(inputSpecsetIdx),
          m_inputSpecset(inputSpecset), ownInput(false),
          m_outputSpecset(outputSpecset), ownOutput(false),
          m_recordedInput(0x0), m_recordedOutput(0x0)
  {

    m_name = "ExecMergeConvert";
    m_type = "ExecMergeConvert";
    m_recordedInput = new vector<pair<string, int> >;
    m_recordedOutput = new vector<pair<string, int> >;
    m_separator = ";";

    if (!initialize())
    {
      abort();
    }
  }

  ExecMergeConvert::~ExecMergeConvert(void)
  {
    if (ownInput && ownOutput && m_inputSpecset == m_outputSpecset
        && m_outputSpecset != 0)
    {
      delete m_outputSpecset;
      m_inputSpecset = 0;
      m_outputSpecset = 0;
    }
    else if (m_inputSpecset == m_outputSpecset && m_outputSpecset != 0)
    {
      m_inputSpecset = 0;
      m_outputSpecset = 0;
    }

    if (ownInput)
    {
      if (m_inputSpecsetIdx != 0)
      {
        delete m_inputSpecsetIdx;
      }
      if (m_inputSpecset != 0)
      {
        delete m_inputSpecset;
      }
    }
    if (ownOutput)
    {
      if (m_outputSpecset != 0)
      {
        delete m_outputSpecset;
      }
    }

    if (m_recordedInput != 0x0)
    {
      delete m_recordedInput;
    }

    if (m_recordedOutput != 0x0)
    {
      delete m_recordedOutput;
    }
  }

  ExecBase * ExecMergeConvert::clone(const ParameterList & inputParams) const
  {
    return new ExecMergeConvert(inputParams);
  }

  bool ExecMergeConvert::invoke(void)
  {

    if (m_inputSpecset == 0)
    {
      ERROR_MSG("No input spectra !!");
      return false;
    }

    if (!initialize())
    {
      return false;
    }

    const char* separator = m_separator.c_str();

    if (m_params.exists("TOLERANCE_PEAK")
        && !m_params.exists("TOLERANCE_PEAK_PPM"))
    {
      float pktol = m_params.getValueFloat("TOLERANCE_PEAK");
      m_outputSpecset->setPeakTolerance(pktol, false);
      m_recordedProcessing << "Setting peak tolerance to " << pktol << " Da\n";
      DEBUG_MSG("Setting peak tolerance to " << pktol << " Da");
    }
    else if (m_params.exists("TOLERANCE_PEAK_PPM"))
    {
      float pktol = m_params.getValueFloat("TOLERANCE_PEAK_PPM");
      m_outputSpecset->setPeakTolerance(pktol, true);
      m_recordedProcessing << "Setting peak tolerance to " << pktol << " ppm\n";
      DEBUG_MSG("Setting peak tolerance to " << pktol << " ppm");
    }

    if (m_params.exists("TOLERANCE_PM") && !m_params.exists("TOLERANCE_PM_PPM"))
    {
      float pktol = m_params.getValueFloat("TOLERANCE_PM");
      m_outputSpecset->setParentMassTolerance(pktol, false);
      m_recordedProcessing << "Setting parent mass tolerance to " << pktol
          << " Da\n";
      DEBUG_MSG("Setting parent mass tolerance to " << pktol << " Da");
    }
    else if (m_params.exists("TOLERANCE_PM_PPM"))
    {
      float pktol = m_params.getValueFloat("TOLERANCE_PM_PPM");
      m_outputSpecset->setParentMassTolerance(pktol, true);
      m_recordedProcessing << "Setting parent mass tolerance to " << pktol
          << " ppm\n";
      DEBUG_MSG("Setting parent mass tolerance to " << pktol << " ppm");
    }

    if (m_outputSpecset->size() == 0)
    {
      WARN_MSG("Detected 0 spectra, skipping ExecMergeConvert::invoke()");
      return true;
    }

    if (m_params.exists("FILTER_CHARGE"))
    {
      if (m_params.exists("FILTER_OUTPUT_CHARGES"))
      {
        WARN_MSG("Ignoring FILTER_CHARGE, using FILTER_OUTPUT_CHARGES instead");
      }
      m_params.addIfDoesntExist("FILTER_OUTPUT_CHARGES",
                                m_params.getValue("FILTER_CHARGE"));
    }

    if (m_params.exists("FILTER_INPUT_CHARGES"))
    {
      list<string> chargeRanges;
      splitText(m_params.getValue("FILTER_INPUT_CHARGES").c_str(),
                chargeRanges,
                separator);

      if (chargeRanges.size() != m_inputSpecsetIdx->size())
      {
        ERROR_MSG("Found " << chargeRanges.size() << " input charge filters for " << m_inputSpecsetIdx->size() << " input files, must have same for both");
        return false;
      }

      int rangeIdx = 0;
      set<short> locRangeVals;
      for (list<string>::iterator rangeIt = chargeRanges.begin();
          rangeIt != chargeRanges.end(); rangeIt++)
      {

        locRangeVals.clear();
        if (!getRanges(rangeIt->c_str(), locRangeVals))
        {
          ERROR_MSG("Failed to parse charge range \'" << *rangeIt << "\'");
          return false;
        }

        m_recordedProcessing << "Keeping spectra of precursor charge "
            << *rangeIt << " from indices "
            << (*m_inputSpecsetIdx)[rangeIdx].first << " to "
            << (*m_inputSpecsetIdx)[rangeIdx].second << "\n";
        DEBUG_MSG("Keeping spectra of precursor charge " << *rangeIt << " from indices " << (*m_inputSpecsetIdx)[rangeIdx].first << " to " << (*m_inputSpecsetIdx)[rangeIdx].second);

        for (int i = (*m_inputSpecsetIdx)[rangeIdx].first;
            i <= (*m_inputSpecsetIdx)[rangeIdx].second; i++)
        {
          if (locRangeVals.count((*m_outputSpecset)[i].parentCharge) == 0)
          {
            (*m_outputSpecset)[i].resize(0);
          }
        }
        ++rangeIdx;
      }
    }

    if (m_params.exists("SET_FILENAMES"))
    {
      vector<string> filenames;
      splitText(m_params.getValue("SET_FILENAMES").c_str(),
                filenames,
                separator);

      if (filenames.size() != m_inputSpecsetIdx->size())
      {
        ERROR_MSG("SET_FILENAMES length (" << filenames.size() << ") must match INPUT_SPECTRA length (" << m_inputSpecsetIdx->size() << ")");
        return false;
      }

      DEBUG_MSG("Setting filenames to \'" << m_params.getValue("SET_FILENAMES") << "\'");
      m_recordedProcessing << "Setting filenames to \'"
          << m_params.getValue("SET_FILENAMES") << "\'\n";

      for (int fIdx = 0; fIdx < m_inputSpecsetIdx->size(); fIdx++)
      {
        for (int i = (*m_inputSpecsetIdx)[fIdx].first;
            i <= (*m_inputSpecsetIdx)[fIdx].second; i++)
        {
          (*m_outputSpecset)[i].setFilename(filenames[fIdx]);
        }
      }
    }

    if (m_params.exists("INSTRUMENT_TYPE"))
    {
      string instType = m_params.getValue("INSTRUMENT_TYPE");
      std::transform(instType.begin(),
                     instType.end(),
                     instType.begin(),
                     ::tolower);
      if (instType == "it")
      {
        instType = "ion trap";
      }
      else if (instType == "ft")
      {
        instType = "orbitrap";
      }

      DEBUG_MSG("Setting mass analyzer to \'" << instType << "\' in all input spectra");
      m_recordedProcessing << "Setting mass analyzer to \'" << instType
          << "\' in all input spectra" << endl;
      for (int i = 0; i < m_outputSpecset->size(); i++)
      {
        (*m_outputSpecset)[i].msMassAnalyzerType =
            Spectrum::parseMassAnalyzer(instType);
      }
    }

    if (m_params.getValueBool("CONVERT_FT_CID_TO_HCD", false))
    {
      float minCollisonEnergy = 39.9;
      int numConvertedHCD = 0;
      for (int i = 0; i < m_outputSpecset->size(); i++)
      {
        if ((*m_outputSpecset)[i].msFragType == Spectrum::FragType_CID
            && (*m_outputSpecset)[i].msMassAnalyzerType
                == Spectrum::MassAnalyzer_ORBI_TRAP
            && (*m_outputSpecset)[i].collision_energy > minCollisonEnergy)
        {
          (*m_outputSpecset)[i].msFragType = Spectrum::FragType_HCD;
          numConvertedHCD++;
        }
      }
      DEBUG_MSG("Set " << numConvertedHCD << " CID spectra acquired in an orbitrap mass analyzer with collision enrgy > " << minCollisonEnergy << " to HCD");
      m_recordedProcessing << "Set " << numConvertedHCD
          << " CID spectra acquired in an orbitrap mass analyzer with collision enrgy > "
          << minCollisonEnergy << " to HCD" << endl;
    }

    if (m_params.exists("FILTER_ACTIVATION"))
    {
      string activsFilter = m_params.getValue("FILTER_ACTIVATION");
      int numKept = 0;

      list<string> activList;
      splitText(activsFilter.c_str(), activList, ";");

      set<Spectrum::FragType> activSet;

      for (list<string>::const_iterator aIt = activList.begin();
          aIt != activList.end(); aIt++)
      {
        Spectrum::FragType activType = Spectrum::parseActivation(*aIt);
        activSet.insert(activType);
      }

      SpecSet newSpecs(m_outputSpecset->size());
      unsigned int idxUse = 0;
      for (unsigned int i = 0; i < m_outputSpecset->size(); i++)
      {
        if (activSet.count((*m_outputSpecset)[i].msFragType) > 0)
        {
          newSpecs[idxUse++] = (*m_outputSpecset)[i];
          numKept++;
        }
      }
      newSpecs.resize(idxUse);
      m_outputSpecset->operator =(newSpecs);

      DEBUG_MSG("Only keeping " << numKept << " " << activsFilter << " spectra");
      m_recordedProcessing << "Only keeping " << numKept << " " << activsFilter
          << " spectra\n";
    }

    if (m_params.exists("FILTER_MASS_ANALYZER"))
    {
      string activsFilter = m_params.getValue("FILTER_MASS_ANALYZER");
      int numKept = 0;

      list<string> activList;
      splitText(activsFilter.c_str(), activList, ";");

      set<Spectrum::MassAnalyzerType> activSet;

      for (list<string>::const_iterator aIt = activList.begin();
          aIt != activList.end(); aIt++)
      {
        Spectrum::MassAnalyzerType activType =
            Spectrum::parseMassAnalyzer(*aIt);
        activSet.insert(activType);
      }

      SpecSet newSpecs(m_outputSpecset->size());
      unsigned int idxUse = 0;
      for (unsigned int i = 0; i < m_outputSpecset->size(); i++)
      {
        if (activSet.count((*m_outputSpecset)[i].msMassAnalyzerType) > 0)
        {
          newSpecs[idxUse++] = (*m_outputSpecset)[i];
          numKept++;
        }
      }
      newSpecs.resize(idxUse);
      m_outputSpecset->operator =(newSpecs);

      DEBUG_MSG("Only keeping " << numKept << " " << activsFilter << " spectra");
      m_recordedProcessing << "Only keeping " << numKept << " " << activsFilter
          << " spectra\n";
    }

    if (m_params.getValueInt("FIX_CHARGE_ZEROS", 0) == 1)
    {
      DEBUG_MSG("Guessing the charge for any charge 0 spectrum where the parent mass is less than the last mass in the spectrum");
      for (int i = 0; i < m_outputSpecset->size(); i++)
      {
        if ((*m_outputSpecset)[i].parentCharge == 0
            && (*m_outputSpecset)[i].parentMass + AAJumps::massMH
                < (*m_outputSpecset)[i][(*m_outputSpecset)[i].size() - 1][0])
        {

          while ((*m_outputSpecset)[i].parentMass + AAJumps::massMH
              < (*m_outputSpecset)[i][(*m_outputSpecset)[i].size() - 1][0])
          {
            (*m_outputSpecset)[i].parentCharge += 1;
            (*m_outputSpecset)[i].parentMass = ((*m_outputSpecset)[i].parentMZ
                * (float)(*m_outputSpecset)[i].parentCharge)
                - (AAJumps::massHion
                    * ((float)(*m_outputSpecset)[i].parentCharge - 1.0));
          }

          m_recordedProcessing << "Assigned charge "
              << (*m_outputSpecset)[i].parentCharge << " to spectrum " << i
              << "\n";

          DEBUG_MSG("Assigned charge " << (*m_outputSpecset)[i].parentCharge << " to spectrum " << i);

        }
      }
    }

    if (m_params.exists("CONVERT_CHARGE"))
    {

      list<string> chargeRanges;
      splitText(m_params.getValue("CONVERT_CHARGE").c_str(),
                chargeRanges,
                separator);

      map<short, short> chargeRangeVals;
      const char* rangeSep = ">";

      for (list<string>::iterator rangeIt = chargeRanges.begin();
          rangeIt != chargeRanges.end(); rangeIt++)
      {
        vector<string> rangeConvert;
        splitText(rangeIt->c_str(), rangeConvert, rangeSep);

        if (rangeConvert.size() != 2)
        {
          ERROR_MSG("Failed to parse charge range conversion \'" << rangeIt->c_str() << "\', must be something like \'4-7" << rangeSep << "3\'");
          return false;
        }

        set<short> rangeVals;
        if (!getRanges(rangeConvert[0].c_str(), rangeVals))
        {
          ERROR_MSG("Failed to parse charge range \'" << rangeConvert[0] << "\'");
          return false;
        }
        short newCharge = (short)getInt(rangeConvert[1].c_str());
        m_recordedProcessing << "Converting all charge " << rangeConvert[0]
            << " spectra to charge " << newCharge << "\n";
        DEBUG_MSG("Converting all charge " << rangeConvert[0] << " spectra to charge " << newCharge);

        for (set<short>::iterator cIt = rangeVals.begin();
            cIt != rangeVals.end(); cIt++)
        {
          chargeRangeVals[*cIt] = newCharge;
        }
      }

      for (int i = 0; i < m_outputSpecset->size(); i++)
      {
        if (chargeRangeVals.count((*m_outputSpecset)[i].parentCharge) > 0)
        {
          (*m_outputSpecset)[i].setCharge(chargeRangeVals[(*m_outputSpecset)[i].parentCharge]);
        }
      }
    }

    if (m_params.getValueInt("WEAVE", 0) > 1)
    {
      int numWeave = m_params.getValueInt("WEAVE");

      m_recordedProcessing << "Weaving spectra from the first " << numWeave
          << " files\n";
      DEBUG_MSG("Weaving spectra from the first " << numWeave << " files");

      if (m_inputSpecsetIdx->size() < numWeave)
      {
        ERROR_MSG("Only loaded " << m_inputSpecsetIdx->size() << " files, need at least " << numWeave << " to weave!");
        return false;
      }

      SpecSet mergedSpecs(m_outputSpecset->size());
      int numSpecWeave = ((*m_inputSpecsetIdx)[0].second
          - (*m_inputSpecsetIdx)[0].first) + 1;

      DEBUG_VAR(numSpecWeave);
      int weaveIdx = -1;

      for (int i = 0; i < numSpecWeave; i++)
      {
        weaveIdx += 1;
        mergedSpecs[weaveIdx] = (*m_outputSpecset)[i];
        int parallelIdx = i + numSpecWeave;
        int numWeaved = 1;
        while (numWeaved < numWeave)
        {

          if (!MZRange::EqualWithinRange((*m_outputSpecset)[i].parentMass,
                                         (*m_outputSpecset)[parallelIdx].parentMass,
                                         0.001))
          {
            WARN_MSG("ERROR: Spectra at parallel index " << i << " have different parent masses (" << (*m_outputSpecset)[i].parentMass << ") and (" << (*m_outputSpecset)[parallelIdx].parentMass << ")");
            DEBUG_VAR(parallelIdx);
          }

          weaveIdx += 1;
          mergedSpecs[weaveIdx] = (*m_outputSpecset)[parallelIdx];
          parallelIdx += numSpecWeave;
          ++numWeaved;
        }
      }
      *m_outputSpecset = mergedSpecs;

    }

    //float samePrecTol = 0.001;

    if (m_params.getValueInt("MERGE_CONSECUTIVE", 0) > 1)
    {
      int allowedConsec = m_params.getValueInt("MERGE_CONSECUTIVE");
      SpecSet mergedSpecs(m_outputSpecset->size());

      m_recordedProcessing << "Merging every " << allowedConsec
          << " consecutive spectra from the same precursor\n";
      DEBUG_MSG("Merging every " << allowedConsec << " consecutive spectra from the same precursor");

      int rootSpecIdx = 0;
      mergedSpecs[rootSpecIdx] = (*m_outputSpecset)[0];
      float lastPM = (*m_outputSpecset)[0].parentMass;
      int curConsec = 0;
      for (int i = 1; i < m_outputSpecset->size(); i++)
      {
        if (curConsec < allowedConsec)
        {

          mergedSpecs[rootSpecIdx].mergeClosestPeaks((*m_outputSpecset)[i], 2);
          ++curConsec;
        }
        else
        {
          /*
           if (curConsec < allowedConsec) {
           WARN_MSG("Only merged " << curConsec << " consecutive spectra from idx "
           << rootSpecIdx << " to " << i - 1 << ". Parent masses not equal ("
           << lastPM << " and " << (*m_outputSpecset)[i].parentMass << ")");
           }
           */
          curConsec = 0;
          ++rootSpecIdx;
          mergedSpecs[rootSpecIdx] = (*m_outputSpecset)[i];
          lastPM = (*m_outputSpecset)[i].parentMass;
        }
      }
      mergedSpecs.resize(rootSpecIdx + 1);

      *m_outputSpecset = mergedSpecs;
    }

    if (m_params.getValueInt("MERGE_PARALLEL", 0) > 1)
    {
      int allowedMerged = m_params.getValueInt("MERGE_PARALLEL");

      m_recordedProcessing
          << "Merging spectra at parallel indices in the first "
          << allowedMerged << " files\n";
      DEBUG_MSG("Merging spectra at parallel indices in the first " << allowedMerged << " files");

      int expectedParallel = (*m_inputSpecsetIdx)[0].second + 1;

      for (int i = 1; i < allowedMerged; i++)
      {
        if (i >= m_inputSpecsetIdx->size())
        {
          ERROR_MSG("Only loaded " << m_inputSpecsetIdx->size() << " files, need at least " << allowedMerged << " to merge parallel indices!");
          return false;
        }

        int detectedParallel = ((*m_inputSpecsetIdx)[i].second
            - (*m_inputSpecsetIdx)[i].first) + 1;

        if (detectedParallel != expectedParallel)
        {
          ERROR_MSG("Detected " << expectedParallel << " spectra in the first file, need exactly " << expectedParallel << " in the first " << allowedMerged << " files to merge parallel indices (found " << detectedParallel << " spectra in file " << i + 1 << ")");
          return false;
        }
      }

      SpecSet mergedSpecs(expectedParallel);

      for (int i = 0; i < expectedParallel; i++)
      {
        mergedSpecs[i] = (*m_outputSpecset)[i];
        int parallelIdx = i + expectedParallel;
        int numMerged = 1;
        while (numMerged < allowedMerged)
        {
          /*
           if (!MZRange::EqualWithinRange(mergedSpecs[i].parentMass,
           allSpecs[parallelIdx].parentMass,
           0.001)) {
           cerr << "ERROR: Spectra at parallel index " << i
           << " have different parent masses ("
           << mergedSpecs[i].parentMass << ") and ("
           << allSpecs[parallelIdx].parentMass << ")\n";
           return -1;
           }
           */

          mergedSpecs[i].mergeClosestPeaks((*m_outputSpecset)[parallelIdx], 2);
          parallelIdx += expectedParallel;
          ++numMerged;
        }
      }
      *m_outputSpecset = mergedSpecs;
    }

    if (m_params.exists("FILTER_OUTPUT_CHARGES"))
    {
      list<string> chargeRanges;
      splitText(m_params.getValue("FILTER_OUTPUT_CHARGES").c_str(),
                chargeRanges,
                separator);
      set<short> rangeVals;

      for (list<string>::iterator rangeIt = chargeRanges.begin();
          rangeIt != chargeRanges.end(); rangeIt++)
      {

        set<short> locRangeVals;
        if (!getRanges(rangeIt->c_str(), rangeVals))
        {
          ERROR_MSG("Failed to parse charge range \'" << *rangeIt << "\'");
          return false;
        }
        rangeVals.insert(locRangeVals.begin(), locRangeVals.end());
      }

      m_recordedProcessing << "Keeping spectra of precursor charge "
          << m_params.getValue("FILTER_OUTPUT_CHARGES") << "\n";
      DEBUG_MSG("Keeping spectra of precursor charge " << m_params.getValue( "FILTER_OUTPUT_CHARGES"));

      for (int i = 0; i < m_outputSpecset->size(); i++)
      {
        if (rangeVals.count((*m_outputSpecset)[i].parentCharge) == 0)
        {
          (*m_outputSpecset)[i].resize(0);
        }
      }
    }

    if (m_params.exists("MASS_OFFSET"))
    {
      float massOffset = m_params.getValueFloat("MASS_OFFSET");

      DEBUG_MSG("Adding " << massOffset << " to every mass");
      m_recordedProcessing << "Adding " << massOffset << " to every mass\n";

      for (int i = 0; i < m_outputSpecset->size(); i++)
      {
        for (int j = 0; j < (*m_outputSpecset)[i].size(); j++)
        {
          (*m_outputSpecset)[i][j][0] += massOffset;
        }
      }
    }

    if (m_params.exists("ACTIVATION"))
    {
      string activStr = m_params.getValue("ACTIVATION");
      Spectrum::FragType msFragType = Spectrum::parseActivation(activStr);

      m_recordedProcessing << "Setting activation method to \'" << activStr
          << "\' for all spectra\n";
      DEBUG_MSG("Setting activation method to \'" << activStr << "\' for all spectra");

      for (int i = 0; i < m_outputSpecset->size(); i++)
      {
        (*m_outputSpecset)[i].msFragType = msFragType;
      }
    }

    //Removing Peaks Below STD DEV of Noise
    if (m_params.getValueFloat("FILTER_STDDEV_PEAK_INT", 0.0) > 0.0)
    {
      DEBUG_MSG("Removing peaks below " << m_params.getValueFloat("FILTER_STDDEV_PEAK_INT", 0.0) <<" STD DEV ABOVE MEAN");

      for (int i = 0; i < m_outputSpecset->size(); i++)
      {

        list<int> peaksToRemove;

        //Removing peaks lower than min_peak_intensity
        peaksToRemove.clear();

        //Calculating Mean
        float total_int = 0.0;
        int total_count = 0;
        vector<float> peak_intensity_list;

        for (int peak_idx = 0; peak_idx < (*m_outputSpecset)[i].size();
            peak_idx++)
        {
          peak_intensity_list.push_back((*m_outputSpecset)[i][peak_idx][1]);
        }

        //sorting and finding bottom 25%
        sort(peak_intensity_list.begin(), peak_intensity_list.end());
        if (peak_intensity_list.size() < 10)
        {
          continue;
        }

        //float bottom_threshold = peak_intensity_list[peak_intensity_list.size()/4];
        for (int peak_idx = 0; peak_idx < peak_intensity_list.size() / 4;
            peak_idx++)
        {
          total_int += peak_intensity_list[peak_idx];
          total_count++;
        }

        float mean = total_int / total_count;
        float sum_variance = 0.0;
        for (int peak_idx = 0; peak_idx < peak_intensity_list.size() / 4;
            peak_idx++)
        {
          sum_variance += (mean - peak_intensity_list[peak_idx])
              * (mean - peak_intensity_list[peak_idx]);
        }
        float variance = sum_variance / total_count;
        float std_dev = sqrt(variance);

        float upper_threshold = m_params.getValueFloat("FILTER_STDDEV_PEAK_INT",
                                                       0.0) * std_dev + mean;

        //DEBUG_MSG("FILTERING OUT INTENSITIES BELOW : "<<  upper_threshold <<" BECAUSE OF MEAN " << mean << " AND STD DEV " << std_dev << " SCAN: "<<(*m_outputSpecset)[i].scan << " " << (*m_outputSpecset)[i].fileName);

        for (int peak_idx = 0; peak_idx < (*m_outputSpecset)[i].size();
            peak_idx++)
        {
          if ((*m_outputSpecset)[i][peak_idx][1] < upper_threshold)
          {
            peaksToRemove.push_back(peak_idx);
          }
        }

        (*m_outputSpecset)[i].removePeaks(peaksToRemove);

      }
    }

    if (m_params.getValueInt("RANK_FILTER", 0) > 0)
    {
      int rankFiltK = m_params.getValueInt("RANK_FILTER");
      float windowRadius = m_params.getValueFloat("RANK_FILTER_RADIUS",
                                                  AAJumps::minAAmass - 1.0);

      m_recordedProcessing << "Applying rank filtering, only keeping the top "
          << rankFiltK << " peaks(s) in every +/- " << windowRadius
          << " Da window\n";

      DEBUG_MSG("Applying rank filtering, only keeping the top " << rankFiltK << " peaks(s) in every +/- " << windowRadius << " Da window");

      for (int i = 0; i < m_outputSpecset->size(); i++)
      {
        (*m_outputSpecset)[i].rankFilterPeaks(rankFiltK, windowRadius);
      }

    }

    if (m_params.getValueInt("TOP_PEAKS", 0) > 0)
    {
      int topK = m_params.getValueInt("TOP_PEAKS");

      m_recordedProcessing << "Keeping the top " << topK
          << " peaks in every spectrum\n";

      DEBUG_MSG("Keeping the top " << topK << " peaks in every spectrum");

      for (int i = 0; i < m_outputSpecset->size(); i++)
      {
        vector<float> scores((*m_outputSpecset)[i].size());
        for (int j = 0; j < (*m_outputSpecset)[i].size(); j++)
        {
          scores[j] = (*m_outputSpecset)[i][j][1];
        }
        sort(scores.begin(), scores.end());

        if (scores.size() > topK)
        {
          (*m_outputSpecset)[i].filterLowIntensity(scores[scores.size() - topK]);
        }
      }
    }

    if (m_params.exists("REVERSE_OFFSET"))
    {
      float revOffset = m_params.getValueFloat("REVERSE_OFFSET");

      for (int i = 0; i < m_outputSpecset->size(); i++)
      {
        (*m_outputSpecset)[i].reverse(revOffset);
      }
    }

    if (m_params.getValueInt("SET_SCAN_NUMS", 0) > 0)
    {
      DEBUG_MSG("Setting the scan number of each spectrum to its one-based index");
      m_recordedProcessing
          << "Setting the scan number of each spectrum to its index\n";
      for (int i = 0; i < m_outputSpecset->size(); i++)
      {
        (*m_outputSpecset)[i].scan = (unsigned int)(i + 1);
      }
    }

    if (m_params.getValueInt("COMPRESS", 0) > 0)
    {
      DEBUG_MSG("Compressing spectrum indices w/o peaks");
      m_recordedProcessing << "Compressing spectrum indices w/o peaks\n";
      unsigned int idxUse = 0;
      DEBUG_VAR(m_outputSpecset->size());
      for (int i = 0; i < m_outputSpecset->size(); i++)
      {
        if ((*m_outputSpecset)[i].size() > 0 && i > idxUse)
        {
          (*m_outputSpecset)[idxUse++] = (*m_outputSpecset)[i];
        }
      }
      m_outputSpecset->resize(idxUse);
      DEBUG_VAR(m_outputSpecset->size());
    }

    if (m_params.getValueInt("ABS_INTENSITIES", 0) > 0)
    {
      DEBUG_MSG("Setting the absolute value of every peak intensity");
      m_recordedProcessing
          << "Setting the absolute value of every peak intensity";
      for (int i = 0; i < m_outputSpecset->size(); i++)
      {
        for (int j = 0; j < (*m_outputSpecset)[i].size(); j++)
        {
          (*m_outputSpecset)[i][j][1] = abs((*m_outputSpecset)[i][j][1]);
        }
        //(*m_outputSpecset)[i].scan = (unsigned int) (i + 1);
      }
    }

    if (m_params.getValueInt("FILTER_PRECURSOR_WINDOW", 0) > 0)
    {
      DEBUG_MSG("Removing peaks near precursor");
      for (int i = 0; i < m_outputSpecset->size(); i++)
      {
        //Removing peaks around precursor
        float precursor_mz = (*m_outputSpecset)[i].parentMZ;
        list<int> peaksToRemove;
        for (int peak_idx = 0; peak_idx < (*m_outputSpecset)[i].size();
            peak_idx++)
        {
          if ((*m_outputSpecset)[i][peak_idx][0] > (precursor_mz - 17)
              && (*m_outputSpecset)[i][peak_idx][0] < (precursor_mz + 15))
          {
            //(*m_outputSpecset)[i][peak_idx][1] = 0.f;
            peaksToRemove.push_back(peak_idx);
          }
        }

        (*m_outputSpecset)[i].removePeaks(peaksToRemove);

      }
    }

    if (m_params.getValueFloat("MIN_PEAK_INT", 0.0) > 0.0)
    {
      DEBUG_MSG("Removing peaks below " << m_params.getValueFloat("MIN_PEAK_INT", 0.0));
      for (int i = 0; i < m_outputSpecset->size(); i++)
      {

        list<int> peaksToRemove;

        //Removing peaks lower than min_peak_intensity
        peaksToRemove.clear();
        float min_peak_intensity = m_params.getValueFloat("MIN_PEAK_INT", 50.0);
        for (int peak_idx = 0; peak_idx < (*m_outputSpecset)[i].size();
            peak_idx++)
        {
          if ((*m_outputSpecset)[i][peak_idx][1] < min_peak_intensity)
          {
            //(*m_outputSpecset)[i][peak_idx][1] = 0.f;
            peaksToRemove.push_back(peak_idx);
          }
        }

        (*m_outputSpecset)[i].removePeaks(peaksToRemove);

      }
    }

    if (m_params.exists("FILENAME"))
    {
      string newFilename = m_params.getValue("FILENAME");

      m_recordedProcessing << "Setting each spectrum's filename to \'"
          << newFilename << "\'\n";

      DEBUG_MSG("Setting each spectrum's filename to \'" << newFilename << "\'");

      for (int i = 0; i < m_outputSpecset->size(); i++)
      {
        (*m_outputSpecset)[i].fileName = newFilename;
      }
    }

    if (m_params.getValueBool("REMOVE_ETD_REDUCED_Z", false))
    {
      int totalPeaksRemoved = 0;
      for (int i = 0; i < m_outputSpecset->size(); i++)
      {
        if ((*m_outputSpecset)[i].msFragType != Spectrum::FragType_ETD)
        {
          continue;
        }
        totalPeaksRemoved +=
            (*m_outputSpecset)[i].removeChargeReducedPrecursors();
      }

      m_recordedProcessing << "Removed " << totalPeaksRemoved
          << " charge-reduced and neutral-loss precursors from MS/MS spectra\n";

      DEBUG_MSG("Removed " << totalPeaksRemoved << " charge-reduced and neutral-loss precursors from MS/MS spectra\n");
    }

    return true;
  }

  bool ExecMergeConvert::loadInputData(void)
  {

    m_inputSpecsetIdx->resize(0);
    m_inputSpecset->resize(0);
    m_recordedInput->resize(0);

    const char* separator = m_separator.c_str();

    bool loadSpecResult;
    int idxUse = 0;

    pair<string, int> nextPair;

    // Load spectra from text files of paths
    if (m_params.exists("INPUT_SPECTRA_LISTS"))
    {
      SpecSet* locSpecs = new SpecSet;
      list<string> listFilesList;
      pair<bool, list<string> > fileListResult;
      splitText(m_params.getValue("INPUT_SPECTRA_LISTS").c_str(),
                listFilesList,
                separator);
      for (list<string>::iterator fileListIt = listFilesList.begin();
          fileListIt != listFilesList.end(); fileListIt++)
      {
        fileListResult = getFileList(*fileListIt);
        if (!fileListResult.first)
        {
          return false;
        }
        for (list<string>::iterator fileIt = fileListResult.second.begin();
            fileIt != fileListResult.second.end(); fileIt++)
        {
          loadSpecResult = loadSpecset(m_params.getValue("EXE_DIR"),
                                       *fileIt,
                                       locSpecs);

          if (!loadSpecResult)
          {
            return false;
          }
          m_recordedInput->push_back(pair<string, int>(*fileIt,
                                                       locSpecs->size()));
          m_inputSpecsetIdx->push_back(pair<int, int>(idxUse,
                                                      idxUse + locSpecs->size()
                                                          - 1));

          idxUse += locSpecs->size();

          m_inputSpecset->swapAppendSpecSet(*locSpecs, false);
        }
      }

      delete locSpecs;
    }

    // Load spectra from direct paths
    if (m_params.exists("INPUT_SPECTRA"))
    {
      vector<string> fileNames;

      vector<TwoValues<unsigned int> > fileSpecIdxs;
      if (!loadSpecsetMultiple(m_params.getValue("EXE_DIR"),
                               m_params.getValue("INPUT_SPECTRA"),
                               m_inputSpecset,
                               &fileNames,
                               m_inputSpecsetIdx))
      {
        return false;
      }

      for (unsigned int i = 0; i < fileNames.size(); i++)
      {
        m_recordedInput->push_back(pair<string, int>(fileNames[i],
                                                     (*m_inputSpecsetIdx)[i].second
                                                         - (*m_inputSpecsetIdx)[i].first
                                                         + 1));
      }
    }

    int numSpecs = 0;
    for (int i = 0; i < m_inputSpecset->size(); i++)
    {
      if ((*m_inputSpecset)[i].parentCharge >= 2)
      {
        numSpecs++;
      }
    }

    DEBUG_MSG("Have " << numSpecs << " spectra of charge >= 2");

    return true;
  }

  bool ExecMergeConvert::saveInputData(std::vector<std::string> & filenames)
  {
    return true;
  }

  bool ExecMergeConvert::saveOutputData(void)
  {
    const char* separator = m_separator.c_str();

    bool splitOut = m_params.getValueInt("SPLIT_OUTPUT", 0) > 0;

    int curIdx = 0, specsPerFile = 0;

    // Save merged spectra as single SpecSet
    string outputSpectraNames = m_params.getValue("OUTPUT_SPECTRA");
    // This is for ProteoSAFe
    if (outputSpectraNames.empty()) {
      outputSpectraNames = m_params.getValue("MERGE_CONVERT_OUTPUT_SPECTRA");
    }
    DEBUG_VAR(outputSpectraNames);

    string outDir = m_params.getValue("OUTPUT_DIR");
    if (outDir.empty()) {
      outDir = ".";
    }
    DEBUG_VAR(outDir);

    if (!outputSpectraNames.empty())
    {
      list<string> listFiles;
      splitText(outputSpectraNames.c_str(),
                listFiles,
                separator);

      if (splitOut)
      {
        int numOut = (int)listFiles.size();
        specsPerFile = ((int)m_outputSpecset->size()) / (numOut);
        int remainder = ((int)m_outputSpecset->size()) % (numOut);
        if (remainder > 0)
        {
          ++specsPerFile;
        }
        DEBUG_MSG("Splitting output spectra into " << numOut << " sets of ~" << specsPerFile << " spectra");
        m_recordedProcessing << "Splitting output spectra into " << numOut
            << " sets of ~" << specsPerFile << " spectra\n";
      }
      for (list<string>::iterator fileIt = listFiles.begin();
          fileIt != listFiles.end(); fileIt++)
      {

        if (splitOut)
        {
          int startIdx = curIdx;
          int endIdx = min(curIdx + specsPerFile - 1,
                           ((int)m_outputSpecset->size()) - 1);
          SpecSet tempSpecs(endIdx - startIdx + 1);
          int locIdx = 0;
          for (int i = startIdx; i <= endIdx; i++)
          {
            tempSpecs[locIdx] = (*m_outputSpecset)[i];
            ++locIdx;
          }
          string saveName = outDir + "/" + *fileIt;
          DEBUG_VAR(saveName);
          if (!saveSpecset(saveName, &tempSpecs))
          {
            return false;
          }
          m_recordedOutput->push_back(pair<string, int>(saveName,
                                                        tempSpecs.size()));
          curIdx = endIdx + 1;
        }
        else
        {
          if (!saveSpecset(*fileIt, m_outputSpecset))
          {
            return false;
          }
          m_recordedOutput->push_back(pair<string, int>(*fileIt,
                                                        m_outputSpecset->size()));
        }
      }
    }

    if (m_params.getValueBool("ORIGIN_SPLIT", false))
    {
      for (int i = 0; i < m_inputSpecsetIdx->size(); i++)
      {
        const int idxStart = (*m_inputSpecsetIdx)[i].first;
        const int idxEnd = (*m_inputSpecsetIdx)[i].second;
        SpecSet tempSpecs(idxEnd - idxStart + 1);

        int tempIdx = 0;
        for (int j = idxStart; j <= idxEnd; j++)
        {
          tempSpecs[tempIdx++] = (*m_outputSpecset)[j];
        }

        string filename = (*m_recordedInput)[i].first;

        FilenameManager fm(filename);

        string outputFilename = fm.filename;
        outputFilename += ".mgf";

        if (!saveSpecset(outputFilename, &tempSpecs))
        {
          return false;
        }
        m_recordedOutput->push_back(pair<string, int>(outputFilename,
                                                      tempSpecs.size()));
      }
    }

    // Save a record of what was loaded, merged, filtered, and/or saved
    if (ownInput and ownOutput)
    {
      stringstream filePath;
      filePath << m_params.m_inputFilePath << ".info";

      DEBUG_MSG("Saving activity record to \'" << filePath.str() << "\' ...");
      bool success = saveActivity(filePath.str());
    }

    return true;
  }

  bool ExecMergeConvert::loadOutputData(void)
  {

    return true;
  }

  std::vector<ExecBase *> const & ExecMergeConvert::split(int numSplit)
  {

    m_subModules.resize(0);

    return m_subModules;
  }

  bool ExecMergeConvert::merge(void)
  {

    return true;
  }

  bool ExecMergeConvert::validateParams(std::string & error)
  {
    m_isValid = false;

    if (!m_params.exists("INPUT_SPECTRA_LISTS"))
    {
      VALIDATE_PARAM_EXIST("INPUT_SPECTRA");
    }
    else
    {
      VALIDATE_PARAM_EXIST("INPUT_SPECTRA_LISTS");
    }

    m_isValid = true;
    return true;
  }

  bool ExecMergeConvert::saveActivity(const string& infoFilePath)
  {

    FILE* output = fopen(infoFilePath.c_str(), "wb");
    if (output == NULL)
    {
      ERROR_MSG("Failed to save activity record to file \'" << infoFilePath << "\'");
      return false;
    }

    fprintf(output, "<Input File>\t<# Spectra Loaded>\n");
    for (vector<pair<string, int> >::const_iterator inIt =
        m_recordedInput->begin(); inIt != m_recordedInput->end(); inIt++)
    {
      fprintf(output, "%s\t%d\n", inIt->first.c_str(), inIt->second);
    }

    fprintf(output, "\n%s\n", m_recordedProcessing.str().c_str());

    fprintf(output, "<Output File>\t<# Spectra Saved>\n");
    for (vector<pair<string, int> >::const_iterator inIt =
        m_recordedOutput->begin(); inIt != m_recordedOutput->end(); inIt++)
    {
      fprintf(output, "%s\t%d\n", inIt->first.c_str(), inIt->second);
    }

    fclose(output);

    return true;
  }

  bool ExecMergeConvert::initialize()
  {
    // If we own the input, there is no point keeping the input
    //   spectra separate from the output spectra (saves memory w/ big files)

    if (ownInput && m_inputSpecsetIdx == 0)
    {
      m_inputSpecsetIdx = new vector<pair<int, int> >;
    }

    if (ownInput && ownOutput)
    {
      if (m_outputSpecset == 0)
      {
        m_outputSpecset = new SpecSet;
      }
      m_inputSpecset = m_outputSpecset;
    }
    else if (ownInput && !ownOutput && m_outputSpecset != 0)
    {
      m_inputSpecset = m_outputSpecset;
    }
    else if (!ownInput && ownOutput && m_inputSpecset != 0)
    {
      m_outputSpecset = m_inputSpecset;
    }
    else if (!ownInput && !ownOutput && m_inputSpecset != 0
        && m_outputSpecset != 0)
    {
      m_outputSpecset->operator =(*m_inputSpecset);
    }
    else
    {
      ERROR_MSG("Input/output specsets not properly set!!");
      DEBUG_VAR(ownInput);
      DEBUG_VAR(ownOutput);
      DEBUG_VAR(m_inputSpecset);
      DEBUG_VAR(m_outputSpecset);
      return false;
    }
    return true;
  }

}
