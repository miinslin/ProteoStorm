#ifndef _Paramters_H_
#define _Paramters_H_

#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;
namespace specnets
{

  class ParameterList
  {
  public:

    // This ParameterList's input parameters file.
    // Always set to the last file path that was successfully loaded by readFromFile.
    std::string m_inputFilePath;

    static string getInputSpectraC13Param(int c13Offset);

    static string getOutputSpectraC13Param(int c13Offset);

    //! \name CONSTRUCTORS
    //@{
    ParameterList();
    ParameterList(const ParameterList & that);
    //@}

    //! \name DESTRUCTOR
    //@{
    virtual ~ParameterList();
    //@}

    //! \name ACCESSORS
    //@{
    /*! \brief Determines whether the specified parameter is present.

     @param paramName the name of the parameter to test for.
     @return True if the parameter is present, False otherwise.
     */
    bool exists(const std::string & paramName) const;

    /*! \brief Determines whether the given required set of parameters is present.

     @param required the required set of parameters to verify.
     @return true if the parameters are present and false otherwise.
     */
    bool exists(std::vector<std::string> paramList) const;

    /*! \brief Gets the value of the requested parameter.

     @param paramName the name of the parameter to read.
     @param defaultValue default value if parameter is not found
     @return a char * to the value of the parameter.
     */
    std::string
    getValue(const std::string & paramName, const std::string & defaultValue =
                 std::string("")) const;

    /*! \brief Gets the value of the requested parameter as a boolean.

     @param paramName the name of the parameter to read.
     @param defaultValue default value if parameter is not found
     @return Boolean value of the requested parameter, or 0 if is not present.
     */
    bool
    getValueBool(const std::string & paramName, bool defaultValue = 0) const;

    /*! \brief Gets the value of the requested parameter as a double.

     @param paramName the name of the parameter to get.
     @param defaultValue default value if parameter is not found
     @return Double value of the requested parameter, or 0.0 if is not present.
     */
    double getValueDouble(const std::string & paramName, double defaultValue =
                              0.0) const;

    /*! \brief Gets the value of the requested parameter as a float.

     @param paramName the name of the parameter to get.
     @param defaultValue default value if parameter is not found
     @return Double value of the requested parameter, or 0.0 if is not present.
     */
    float
    getValueFloat(const std::string & paramName,
                  float defaultValue = 0.0) const;

    /*! \brief Gets the value of the requested parameter as an integer.

     @param paramName the name of the parameter to read.
     @param defaultValue default value if parameter is not found
     @return Integer value of the requested parameter, or 0 if is not present.
     */
    int getValueInt(const std::string & paramName, int defaultValue = 0) const;

    /*! \brief Prints all currently loaded parameters to the given stream
     */
    void print(std::ostream & strm) const;

    /*! \brief Gets the size of the parameters collection found.

     @return the size of the parameters collection - zero if empty.
     */
    unsigned int size(void) const;

    /*! \brief Writes the parameters to the given file.

     @param filename the filename to write the parameters to.
     @return
     */
    bool writeToFile(const std::string & filename) const;

    //@}

    //! \name MODIFIERS
    //@{
    /*! /brief Adds a parameter it is found in another list.
     *
     * If the parameter specified by paramName exists in the list being passed
     * to the method then add it to our list. Otherwise do nothing.
     @param that ParamterList to check for existence
     @param paramName the parameter we wish to add (if it exists)
     @return true if the parameter was added. False if not
     */
    bool addIfExists(const ParameterList & that, const std::string & paramName);

    /*! /brief Adds a parameter if it does not already exist in list.
     *
     * If the parameter specified by paramName does not exist in our list already
     * then add it to our list. Otherwise do nothing.
     *@param paramName the name of the parameter to set.
     *@param paramValue the value of the named parameter.
     @return true if the parameter was added. False if not
     */
    bool addIfDoesntExist(const std::string & paramName,
                          const std::string & paramValue);

    /*! /brief Adds another list of parameters to this one.
     *
     * If the paramter list being added contains parameters
     * with the same name as the existing list, the existing
     * parameters may be overwritten or kept as indicated by
     * the overwrite flag.
     @param that another ParamterList to combine with this one.
     @param overwirte should existing parameters be overwritten.
     @return
     */
    bool addList(const ParameterList & that, bool overwrite);

    /*! /brief Reads parameters from a given file.

     @param filename the filename to read the parameters from.
     @return
     */
    bool readFromFile(const std::string & filename);

    /*! /brief Reads parameters from a given file that is in an XML provided by proteosafe

     @param filename the filename to read the parameters from.
     @return
     */
    bool readFromProteosafeXMLFile(const std::string & filename);

    /**
     * copies any pre-defined ProteoSafe parameters to pre-defined specnets parameters
     * @param that proteoSafe parameter list
     * @return true
     */
    void addProteosafeParams(const ParameterList & that);

    /*! /brief Sets the value of a parameter.
     *
     *@param paramName the name of the parameter to set.
     *@param paramValue the value of the named parameter.
     */
    void
    setValue(const std::string & paramName, const std::string & paramValue);

    /*! /brief Removes a parameter from the list.
     *
     *@param paramName the name of the parameter to remove.
     */
    void removeParam(const std::string & paramName);
    //@}

    //! \name MODIFIERS
    //@{
    ParameterList & operator=(const ParameterList & that);
    //@}

    void getGroups(std::map<std::string, std::string> &goups,
                   std::string &groupName);

    inline void clear()
    {
      parameters.clear();
    }

  protected:

    std::map<std::string, std::string> parameters;

  };

} // namespace specnets

#endif // _Paramters_H_
