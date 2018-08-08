#include "inputParams.h"

#if !defined(__MINGW32__)
#include "stdio.h"
#endif

#include <cctype>
#include <fstream>
#include <iostream>


#if defined(__MINGW32__)
namespace
{
inline char * cuserid(char *)
{
  static char sbuffer[32767];
  long ncount = sizeof(sbuffer);

  //GetUserName(sbuffer, & ncount);

  return sbuffer;
}
}
#endif

const char *InputParams::emptyString = "";
float InputParams::Resolution = 0.1;
float InputParams::PeakTol    = 0.5;
float InputParams::PMTol      = 1.0;
float InputParams::MaxShift   = 3000.0;
float InputParams::PPM        = 0.001;

bool InputParams::readParams(const char *filename)
{
  ifstream input(filename, ios::binary);

  if (!input)
  {
    cerr << "Error opening " << filename << endl;
    return false;
  }

  // default
  parameters.clear();
  //parameters["USER"] = cuserid(0);

  for (size_t npos[3]; input.good(); )
  {
    string buffer;

    getline(input, buffer);

    npos[0] = buffer.find('=');
    if (npos[0] != string::npos)
    {
      string svalue = buffer.substr(npos[0] + 1);

      // substitute variables
      for (npos[1] = svalue.find('$'); npos[1] != string::npos; npos[1] = svalue.find('$', npos[1]))
      {
        for (npos[2] = npos[1] + 1; npos[2] < svalue.size(); npos[2] ++)
          if (! isalpha(svalue[npos[2]]))
            break;

        svalue.replace(npos[1], npos[2] - npos[1], parameters[svalue.substr(npos[1] + 1, npos[2] - npos[1] - 1)]);
      }

      parameters[buffer.substr(0, npos[0])] = svalue;
    }
  }

  input.close();

  if (parameters.empty())
    return false;

  if(paramPresent("TOLERANCE_PEAK"))
    PeakTol = (float)getValueDouble("TOLERANCE_PEAK");
  if(paramPresent("TOLERANCE_PM"))
    PMTol = (float)getValueDouble("TOLERANCE_PM");
  if(paramPresent("RESOLUTION"))
    Resolution = (float)getValueDouble("RESOLUTION");
  if(paramPresent("MAX_SHIFT"))
    MaxShift = (float)getValueDouble("MAX_SHIFT");
  if(paramPresent("TOLERANCE_PPM"))
	PPM = (float)getValueDouble("TOLERANCE_PPM");

  return true;
}

bool InputParams::addParam(string param_str, string value_str) {
	if (parameters.count(param_str) > 0) return false;
	parameters[param_str] = value_str;
	return true;
}

bool InputParams::addParam(const char* param, const char* value) {
	string param_str = param;
	string value_str = value;
	return addParam(param_str, value_str);
}

bool InputParams::addParams(vector<TwoValues<string> >& params, vector<bool>& paramCnfrm) {
	paramCnfrm.resize(params.size());
	bool one_added = false;
	for (int i = 0; i < params.size(); i++) {
		paramCnfrm[i] = addParam(params[i][0], params[i][1]);
		if (! one_added && paramCnfrm[i]) one_added = true;
	}
	return one_added;
}

const char * InputParams::getValue(const char *paramName) const
{
  static const char * snull = "";

  if (parameters.find(paramName) != parameters.end())
    return const_cast<char *>(parameters.find(paramName)->second.c_str());
  else
    return snull;
}

bool InputParams::paramPresent(const char *paramName)
{
  return parameters.find(paramName) != parameters.end();
}

//
//  Checks if all required parameters are present
//
bool InputParams::confirmParams(vector<const char *> required)
{
  for (int r = 0; r < required.size(); r ++)
    if (! paramPresent(required[r]))
      return false;

  return true;
}

bool InputParams::dumpParams(FILE* output, vector<string>* paramOrder) {
	if (paramOrder == 0) {
		for (map<string, string>::iterator paramIt = parameters.begin(); paramIt != parameters.end(); paramIt ++) {
			fprintf(output, "%s=%s\n", (paramIt->first).c_str(), (paramIt->second).c_str());
		}
	} else {
		for (vector<string>::iterator paramIt = paramOrder->begin(); paramIt != paramOrder->end(); paramIt ++) {
			if ( ! paramPresent(paramIt->c_str())) {
				cerr << "ERROR: Missing parameter " << paramIt->c_str() << "\n";
				return false;
			}
			fprintf(output, "%s=%s\n", paramIt->c_str(), getValue(paramIt->c_str()));
		}
	}
	return true;
}
