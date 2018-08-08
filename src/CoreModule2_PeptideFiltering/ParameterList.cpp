// Header Include
#include "Logger.h"
#include "ParameterList.h"
#include "utils.h"

// System Include
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>

//#include <sys/types.h>
//#include <pwd.h>

using namespace std;

namespace specnets
{
  const unsigned int NUM_PROTEOSAFE_PARAMS = 3;
  const string PROTEOSAFE_PARAMS[NUM_PROTEOSAFE_PARAMS][2] = {
      { "tolerance.PM_tolerance", "TOLERANCE_PM" }, { "tolerance.Ion_tolerance",
                                                      "TOLERANCE_PEAK" },
      { "instrument.instrument", "INSTRUMENT_TYPE" } };

  const string PROTEOSAFE_INFILE_PARAM = "upload_file_mapping";
  const string SPECNETS_INFILE_PARAM = "INPUT_SPECS_MS";

  string ParameterList::getInputSpectraC13Param(int c13Offset)
  {
    string baseParam = "INPUT_SPECTRA";
    if (c13Offset == 0)
    {
      return baseParam;
    }
    else if (c13Offset < 0)
    {
      baseParam += "_C13";
    }
    else
    {
      baseParam += "_C13+";
    }
    baseParam += parseInt(c13Offset);
    return baseParam;
  }

  string ParameterList::getOutputSpectraC13Param(int c13Offset)
  {
    string baseParam = "OUTPUT_SPECTRA";
    if (c13Offset == 0)
    {
      return baseParam;
    }
    else if (c13Offset < 0)
    {
      baseParam += "_C13";
    }
    else
    {
      baseParam += "_C13+";
    }
    baseParam += parseInt(c13Offset);
    return baseParam;
  }

  //---------------------------------------------------------------------------
  ParameterList::ParameterList(void)
  {
    // EMPTY
  }

  //---------------------------------------------------------------------------
  ParameterList::ParameterList(const ParameterList & that)
  {
    this->parameters = that.parameters;
    m_inputFilePath = that.m_inputFilePath;
  }

  //---------------------------------------------------------------------------
  ParameterList::~ParameterList(void)
  {
    clear();
  }

  //---------------------------------------------------------------------------
  bool ParameterList::exists(const string & paramName) const
  {
    return parameters.find(paramName) != parameters.end();
  }

  //---------------------------------------------------------------------------
  bool ParameterList::exists(vector<string> paramList) const
  {
    for (int r = 0; r < paramList.size(); r++)
    {
      if (!exists(paramList[r]))
      {
        return false;
      }
    }
    return true;
  }

  //---------------------------------------------------------------------------
  std::string ParameterList::getValue(const std::string & paramName,
                                      const std::string & defaultValue) const
  {
    if (parameters.find(paramName) != parameters.end())
      return parameters.find(paramName)->second.c_str();
    else
      return defaultValue;
  }

  //---------------------------------------------------------------------------
  bool ParameterList::getValueBool(const std::string & paramName,
                                   bool defaultValue) const
  {
    if (exists(paramName)) {
      string stringValue = getValue(paramName);
      std::transform(stringValue.begin(), stringValue.end(), 
                     stringValue.begin(), ::tolower);
      if (stringValue == "0" || 
          stringValue == "no" || 
          stringValue == "off" || 
          stringValue == "false") {
        return false;
      } else if (stringValue == "1" || 
                 stringValue == "yes" || 
                 stringValue == "on" || 
                 stringValue == "true") {
        return true;
      }
    }
    return defaultValue;
  }

  //---------------------------------------------------------------------------
  double ParameterList::getValueDouble(const std::string & paramName,
                                       double defaultValue) const
  {
    if (exists(paramName))
    {
      return atof(getValue(paramName).c_str());
    }
    return defaultValue;
  }

  //---------------------------------------------------------------------------
  float ParameterList::getValueFloat(const std::string & paramName,
                                     float defaultValue) const
  {
    if (exists(paramName))
    {
      return atof(getValue(paramName).c_str());
    }
    return defaultValue;
  }

  //---------------------------------------------------------------------------
  int ParameterList::getValueInt(const std::string & paramName,
                                 int defaultValue) const
  {
    if (exists(paramName))
    {
      return atoi(getValue(paramName).c_str());
    }
    return defaultValue;
  }

  //---------------------------------------------------------------------------
  void ParameterList::print(ostream & ostr) const
  {
    for (map<string, string>::const_iterator i = parameters.begin();
        i != parameters.end(); i++)
      ostr << i->first << " = " << i->second << endl;
  }

  //---------------------------------------------------------------------------
  unsigned int ParameterList::size() const
  {
    return parameters.size();
  }

  //---------------------------------------------------------------------------
  void ParameterList::getGroups(map<std::string, std::string> &groups,
                                string &groupName)
  {
    for (map<std::string, std::string>::iterator it = parameters.begin();
        it != parameters.end(); it++)
    {
      if (it->first.compare(0, groupName.length(), groupName) == 0)
      {
        groups[it->first] = it->second;
      }
    }
  }

  //---------------------------------------------------------------------------
  bool ParameterList::readFromFile(const string & filename)
  {
    ifstream input(filename.c_str(), ios::binary);

    if (!input)
    {
      ERROR_MSG("Error opening " << filename)
      return false;
    }

    // default
    parameters.clear();
    //parameters["USER"] = (char*)getpwuid(geteuid()); // cuserid(0);

    for (size_t npos[3]; input.good();)
    {
      string buffer;

      getline(input, buffer);

      //if this is a DOS file, strip off carriage return
      if ('\r' == buffer[buffer.size() - 1])
      { //strip off carriage return for dos
        buffer.resize(buffer.size() - 1);
      }

      npos[0] = buffer.find('=');
      if (npos[0] != string::npos)
      {
        string svalue = buffer.substr(npos[0] + 1);

        // substitute variables
        for (npos[1] = svalue.find('$'); npos[1] != string::npos; npos[1] =
            svalue.find('$', npos[1]))
        {
          for (npos[2] = npos[1] + 1; npos[2] < svalue.size(); npos[2]++)
            if (!isalpha(svalue[npos[2]]))
              break;

          svalue.replace(npos[1],
                         npos[2] - npos[1],
                         parameters[svalue.substr(npos[1] + 1,
                                                  npos[2] - npos[1] - 1)]);
        }

        parameters[buffer.substr(0, npos[0])] = svalue;
      }
    }

    input.close();

    if (parameters.empty())
      return false;

    m_inputFilePath = filename;

    return true;
  }

  //---------------------------------------------------------------------------
  bool ParameterList::readFromProteosafeXMLFile(const string & filename)
  {
    ifstream input(filename.c_str(), ios::binary);

    if (!input)
    {
      ERROR_MSG("Error opening " << filename)
      return false;
    }

    // default
    parameters.clear();
    //parameters["USER"] = (char*)getpwuid(geteuid()); // cuserid(0);

    int upload_mapping_count = 0;

    for (size_t npos[3]; input.good();)
    {
      string buffer;

      getline(input, buffer);

      //This is from proteosafe so we are stripping out stuff to make it compatible with the current parsing
      if (buffer.find("<parameter name=\"") == string::npos)
        continue;
      std::string param_first = "<parameter name=\"";
      buffer.replace(buffer.find(param_first), param_first.length(), "");
      std::string param_second = "\">";
      buffer.replace(buffer.find(param_second), param_second.length(), "=");
      std::string param_third = "</parameter>";
      buffer.replace(buffer.find(param_third), param_third.length(), "");

      //if this is a DOS file, strip off carriage return
      if ('\r' == buffer[buffer.size() - 1])
      { //strip off carriage return for dos
        buffer.resize(buffer.size() - 1);
      }

      npos[0] = buffer.find('=');
      if (npos[0] != string::npos)
      {
        string svalue = buffer.substr(npos[0] + 1);
        // substitute variables
        for (npos[1] = svalue.find('$'); npos[1] != string::npos; npos[1] =
            svalue.find('$', npos[1]))
        {
          for (npos[2] = npos[1] + 1; npos[2] < svalue.size(); npos[2]++)
            if (!isalpha(svalue[npos[2]]))
              break;

          svalue.replace(npos[1],
                         npos[2] - npos[1],
                         parameters[svalue.substr(npos[1] + 1,
                                                  npos[2] - npos[1] - 1)]);
        }
        string key = buffer.substr(0, npos[0]);
        if (key == PROTEOSAFE_INFILE_PARAM)
        {
          char upload_mapping_buf[100];
          sprintf(upload_mapping_buf,
                  "%s%d",
                  key.c_str(),
                  upload_mapping_count);
          key = upload_mapping_buf;
          upload_mapping_count++;
        }

        // It is possible to have multiple same keys in Proteosafe
        //   if this happens we need to add a unique modifier
        string originalKey = key;
        int i = 1;
        while (parameters.find(key) != parameters.end()) {
          char buf[128];
          sprintf(buf, "%d", i);
          key = originalKey + "." + string(buf);
          i++;
        }
        parameters[key] = svalue;
      }
    }

    input.close();

    if (parameters.empty())
      return false;

    ParameterList copy(*this);
    //addProteosafeParams(copy);
    m_inputFilePath = filename;

    return true;
  }

  void ParameterList::addProteosafeParams(const ParameterList & that)
  {
    for (unsigned i = 0; i < NUM_PROTEOSAFE_PARAMS; i++)
    {
      string proteoSafeParam = PROTEOSAFE_PARAMS[i][0];
      string specnetsParam = PROTEOSAFE_PARAMS[i][1];

      if (that.exists(proteoSafeParam))
      {
        //replace ProteoSAFe params with specnets params of the same value
        setValue(specnetsParam, that.getValue(proteoSafeParam));
      }
    }

    unsigned int upload_mapping_count = 0;
    char upload_mapping_buf[100];

    // compute initial ProteoSAFe param for input spectra
    sprintf(upload_mapping_buf,
            "%s%d",
            PROTEOSAFE_INFILE_PARAM.c_str(),
            upload_mapping_count);
    upload_mapping_count++;
    string fileParam = upload_mapping_buf;

    string fileNames = "";

    while (that.exists(fileParam))
    {
      // combine input spectra into one parameter value delimited by ";"
      fileNames += that.getValue(fileParam);
      fileNames += ";";

      // compute next ProteoSAFe param for input spectra
      sprintf(upload_mapping_buf,
              "%s%d",
              PROTEOSAFE_INFILE_PARAM.c_str(),
              upload_mapping_count);
      upload_mapping_count++;
      string fileParam = upload_mapping_buf;
    }

    if (fileNames.length() > 0)
    {
      // add specnets parameter for input spectra
      fileNames.erase(fileNames.length() - 1, 1);
      setValue(SPECNETS_INFILE_PARAM, fileNames);
    }

  }

  //---------------------------------------------------------------------------
  bool ParameterList::writeToFile(const std::string & filename) const
  {
    ofstream of(filename.c_str(), ios::out | ios::binary);
    if (!of)
    {
      return false;
    }

    map<string, string>::const_iterator itr = parameters.begin();
    map<string, string>::const_iterator itr_end = parameters.end();

    for (; itr != itr_end; itr++)
    {
      of << itr->first << "=" << itr->second << endl;
    }
    of.close();

    return true;
  }

  bool addList(const ParameterList & that);

  //---------------------------------------------------------------------------
  void ParameterList::setValue(const std::string & paramName,
                               const std::string & paramValue)
  {
    parameters[paramName] = paramValue;
    return;
  }

  //---------------------------------------------------------------------------
  bool ParameterList::addList(const ParameterList & that, bool overwrite)
  {
    map<string, string>::const_iterator itr = that.parameters.begin();
    map<string, string>::const_iterator itr_end = that.parameters.end();

    for (; itr != itr_end; itr++)
    {
      if (overwrite)
      {
        parameters[itr->first] = itr->second;
      }
      else
      {
        // If we are not overwriting then check for existence of parameter
        if (!exists(itr->first))
        {
          parameters[itr->first] = itr->second;
        }
      }
    }

    return true;
  }

  //---------------------------------------------------------------------------
  bool ParameterList::addIfExists(const ParameterList & that,
                                  const std::string & paramName)
  {
    if (that.exists(paramName))
    {
      parameters[paramName] = that.getValue(paramName);
      return true;
    }
    return false;
  }

  //---------------------------------------------------------------------------
  bool ParameterList::addIfDoesntExist(const std::string & paramName,
                                       const std::string & paramValue)
  {
    if (!exists(paramName))

    {
      parameters[paramName] = paramValue;
      return true;
    }
    return false;
  }

  //---------------------------------------------------------------------------
  void ParameterList::removeParam(const std::string & paramName)
  {
    std::map<std::string, std::string>::iterator itrFind =
        parameters.find(paramName);
    if (itrFind != parameters.end())
    {
      parameters.erase(itrFind);
    }
    return;
  }

  //---------------------------------------------------------------------------
  ParameterList & ParameterList::operator=(const ParameterList & that)
  {
    if (this != &that)
    {
      this->parameters = that.parameters;
      this->m_inputFilePath = that.m_inputFilePath;
    }
    return *this;
  }

} // namespace specnets
