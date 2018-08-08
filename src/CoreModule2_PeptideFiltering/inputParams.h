#ifndef INPUT_PARAMS
#define INPUT_PARAMS

#include <stdlib.h>
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <stdio.h>

#include "twovalues.h"

using namespace std;

/**
 * TODO: add description
 *
 *
 * Parameters file is of the format
 * <num_parameters_in_file>
 * <param_name>=<param_value>  (note no spaces around =)
 *
 */
class InputParams {
    map<string, string> parameters;
    static const char *emptyString;

public:

    /**
    * TODO: add description
    */
    static float Resolution;

    /**
    * TODO: add description
    */
    static float PeakTol;

    /**
    * TODO: add description
    */
    static float PMTol;

    /**
    * TODO: add description
    */
    static float MaxShift;
	
	static float PPM;

    InputParams() {
    }

    /**
    * Reads parameters from a given file.
    *
    *@param filename the filename to read the parameters from.
    *@return
    */
    bool readParams(const char *filename);
	
	/**
	* Adds a parameter with corresponding value.
	*
	*@param param_str the parameter's name
	*@param value_str the parameter's value
	*@return true if parameter was added, false if not added because a parameter with the same name already exists
	*/
	bool addParam(string param_str, string value_str = "");
	
	/**
	* Adds a parameter with corresponding value.
	*
	*@param param the parameter's name
	*@param value the parameter's value
	*@return true if parameter was added, false if not added because a parameter with the same name already exists
	*/
	bool addParam(const char* param, const char* value = "");
	
	/**
	* Adds a sequence of parameters with corresponding values.
	*
	*@param params sequence of parameters' name paired with corresponding values
	*@param paramCnfrm for each parameter in params, true if it was added by "addParam"
	*@return true if one of params was added by addParam, false if all were not added
	*/
	bool addParams(vector<TwoValues<string> >& params, vector<bool>& paramCnfrm);

    /**
    * Gets the value of the requested parameter.
    *
    *@param paramName the name of the parameter to read.
    *@return a char * to the value of the parameter.
    */
    const char *getValue(const char *paramName) const;

    /**
    * Gets the value of the requested parameter as an int.
    *
    *@param paramName the name of the parameter to read.
    *@return the int value of the given parameter. Wraps atoi, so if the given
    *parameter does not exist, this function will return accordingly.
    */
    int getValueInt(const char *paramName)
    {
        return atoi(getValue(paramName));
    }

    /**
    * Gets the value of of the parameter at the specified index.
    *
    *@param index the index of the parameter to get.
    *@return the int value of the requested parameter or 0 if the parameter
    *doesn't exist.
    */
    int getValueInt(unsigned int index)
    {
      unsigned i = 0;
      for (map<string, string>::iterator j = parameters.begin(); j != parameters.end(); i ++, j ++)
        if (i == index)
          return atoi(j->second.c_str());

      return 0;
    }

    /**
    * Gets the value of the requested parameter as a double.
    *
    *@param paramName the name of the parameter to get.
    *@return the requested parameter as a double, or the value of atof
    *if the requested parameter does not exist.
    */
    double getValueDouble(const char *paramName)
    {
        return atof(getValue(paramName));
    }

    /**
    * Gets the value of the requested parameter as a double at a specified
    * index.
    *
    *@param index the index of the parameter to get.
    *@return the value of the requested parameter as a double.
    */
    double getValueDouble(unsigned int index)
    {
      unsigned i = 0;
      for (map<string, string>::iterator j = parameters.begin(); j != parameters.end(); i ++, j ++)
        if (i == index)
          return atof(j->second.c_str());

      return 0.0;
    }

    /**
    * Gets the size of the parameters collection found.
    *
    *@return the size of the parameters collection - zero if empty.
    */
    unsigned int size()
    {
      return parameters.size();
    }

    /**
    * Gets the name of the parameter at the given index.
    *
    *@param index the index of the parameter to get.
    *@return a pointer to the name of the parameter at the given index or
    *emptyString if one doesn't exist.
    */
    const char *getParamName(unsigned int index)
    {
      unsigned i = 0;
      for (map<string, string>::iterator j = parameters.begin(); j != parameters.end(); i ++, j ++)
        if (i == index)
          return const_cast<char *>(j->first.c_str());

      return emptyString;
    }

    /**
    * Determines whether the specified parameter is present.
    *
    *@param paramName the name of the parameter to test for.
    *@return true if the parameter is present, false otherwise.
    */
    bool paramPresent(const char *paramName);

    /**
    * Tests to see whether the given required set of parameters is present.
    *
    *@param required the required set of parameters to verify.
    *@return true if the parameters are present and false otherwise.
    */
    bool confirmParams(vector<const char *> required);
	
	/**
	* Ouptuts all parameters with values to file
	*
	*@param output opened file output stream to write to
	*@param paramOrder optional order of parameter names
	*@return true if parameters were written to file successfully, false if not
	*/
	bool dumpParams(FILE* output, vector<string>* paramOrder = 0);
	
    /**
    * Prints all currently loaded parameters to cerr.
    */
    void showAll()
    {
      for (map<string, string>::iterator i = parameters.begin(); i != parameters.end(); i++)
        cerr << i->first << " = " << i->second << endl;
    }
};

#endif
