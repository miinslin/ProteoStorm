#ifndef __ExecBase_H__
#define __ExecBase_H__

// External Includes
#include "ParameterList.h"

// System Includes
#include <string>
#include <vector>

#define VALIDATE_PARAM_EXIST(value)  \
    if (!m_params.exists(value))                \
    {                                         \
      error = value;                          \
      error += " not present";                \
      return false;                           \
    }

namespace specnets
{
  /*! \brief Base class for all execution framework modules.

   Declares a set of virtual functions required to execute an Specnets module
   within the specnets execution framework. Implements some very basic
   functions such as the isValid() method and holds the parameters structure
   that all modules will require as well as the member data necessary for
   splitting into sub-modules and re-merging the results.

   An execution module of this type may be used in two basic ways:<br>
   1) The module is created using the constructor that takes an InputParams
   parameter and then run using the invoke() method.<br>
   2) The module is created using the constructor that takes an InputParams
   parameter and then one of the classes derived from ParallelExecution is
   created to split the module into sub-modules which will be run in parallel.
   */
  class ExecBase
  {
  public:
    //! \name CONSTRUCTORS
    //@{
    /*! \brief The exemplar constructor.

     Generally this constructor should not be used. It is only used by the
     module execution factory in order to create an exemplar object (without
     valid parameters) which is then used to create a real (valid) object
     using the clone() method.
     @sa clone()
     */
    ExecBase(void);

    /*! \brief The default constructor

     This is the default constructor. A valid set of parameters must be
     supplied in the input parameter. The parameters can then be verified
     using the validateParams() method.
     @sa validateParams()
     @param inputParams structure containing all input parameters necessary for execution
     */
    ExecBase(const ParameterList & inputParams);
    //@}

    //! \name DESTRUCTOR
    //@{
    virtual ~ExecBase(void);
    //@}

    //! \name ACCESSORS
    //@{
    /*! \brief Creates a new module of the virtual class with the given params.

     The parameters should be sufficient for the derived class to be invoked properly.

     @return A pointer to the newly created object
     */
    virtual ExecBase * clone(const ParameterList & input_params) const = 0;

    /*! \brief Returns the class name of the module.

     By default the name of the object is same as the class itself, however it may
     be set to a different value using the setName() method. Or in the case of
     modules that are children or other modules, the name may be a variation of
     the parent's name.

     @sa setName()
     @return The name of the class
     */
    std::string getName(void) const;

    /*! \brief Returns a filename unique to this process ID and ExecModule

     Concatenates this module's name with the process ID between a supplied prefix/suffix

     @sa getTemporaryFilename()
     @return The name of the class
     */
    std::string getTemporaryFilename(const string &prefix, const string &suffix, const bool addPID = true) const;

    /*! \brief Returns the type name of the module.

     The type name is automatically set by the constructor to be the same name
     as the class. This method may be used to determine the true class from
     a pointer to the base (this) class.

     @return The type name of the class
     */
    std::string getType(void) const;

    /*! \brief Returns validity of the module.

     This is a non-virtual method that simply returns the value of the internal
     validity flag. This flag should have been set before calling this method
     by calling the validate() method.

     @return The value of the internal object validity flag.
     @sa validateParams(std::string & error)
     */
    bool isValid(void) const;
    //@}

    //! \name MODIFIERS
    //@{
    /*! \brief Executes the module.

     In order to call this method succesfully all the necessary data for
     execution must already be loaded into memory (data members). This can
     be accomplished using the loadInputData() method.

     @return True if execution finished successfully, false otherwise.
     @sa loadInputData()
     */
    virtual bool invoke(void) = 0;

    /*! \brief Loads the input data from the files specified in the params.

     Loads all the data from files into the data members. This method is
     primarily used by the execution module to load necessary data when
     executing in a separate process..

     @return True if data was loaded successfully, false otherwise.
     @sa ExecBase(const ParameterList & input_params), saveOutputData()
     */
    virtual bool loadInputData(void) = 0;

    /*! \brief Saves all the result data to files specified in the params.

     Saves all the result data generated by the module to files specified
     by values in the params. This is used to either save the data permanantly,
     or to be loaded back in after remote execution is finished and results
     need to be merged by the merge() method.

     @param filenames A list of file names that contain the data necessary to run the module
     @return True if data was saved successfully, false otherwise.
     @sa ExecBase(const ParameterList & input_params), loadInputData(), merge()
     */
    virtual bool saveOutputData(void) = 0;

    /*! \brief Saves all the internal data into files specified in the params.

     Saves all the data required for an external process to execute the
     module. The external process would call loadInputData() to reload the
     data into the members before calling invoke(). The user passes a vector
     that will contain the names of all files necessary to run the module as
     a separate process. The first file in this list will always be the main
     parameter file.

     @param filenames A list of file names that contain the data necessary to run the module
     @return True if data was saved successfully, false otherwise.
     @sa ExecBase(const ParameterList & input_params), loadInputData()
     */
    virtual bool saveInputData(std::vector<std::string> & filenames) = 0;

    /*! \brief Loads the output data from the files specified in the params.

     Loads all the data from output files into the data members. The purpose
     of this is to ready "child" modules to be merged back together after
     being executed separately.

     @return True if data was loaded successfully, false otherwise.
     @sa ExecBase(const ParameterList & input_params), saveOutputData()
     */
    virtual bool loadOutputData(void) = 0;

    /*! \brief Splits the module into multiple "children" for parallel execution.

     Divides the work required by the module into a vector of sub-modules
     that can be executed in parallel (by the the ParallelExecution() class.

     @param numSplit Number of separate modules to split into
     @return The set of sub-modules (children) that the original module has been split into.
     An empty vector implies an error.
     @sa merge()
     */
    virtual std::vector<ExecBase *> const & split(int numSplit) = 0;

    /*! \brief Merges the child modules back into a complete result

     This method is only called when split is used for parallel execution
     and after the execution of all children has been completed. This method
     will merge the results generated by each of the child modules into
     one cohesive result as if the module had been run as a single entity.

     @return True if merge could be performed successfully, false otherwise.
     @sa split(int numSplit)
     */
    virtual bool merge(void) = 0;

    /*! \brief Sets the name of the module.

     By default the name of the object is same as the class itself, however this
     method may be used to set the name to a different value.

     @sa getName()
     @return The name of the class
     */
    void setName(std::string name);

    /*! \brief Performs validation of the input parameters.

     Checks the parameters structure provided in the constructor to see if they
     are sufficient and correct to invoke the module. Also sets the internal
     validity flag so that isValid() will return the correct result.

     @param error A description of the error (if any occurs)
     @return True if the parameters for the module are valid, false otherwise.
     @sa isValid()
     */
    virtual bool validateParams(std::string & error) = 0;
    //@}

    ParameterList m_params;
  protected:
    /*! \brief Creates a unique file name for a sub_module

     Convenience method that creates a filename based on the base name given,
     along with the number of the split.

     @param baseFilename the base name of the file on which will be appended the
     numSplit
     @return The full filename
     */
    std::string makeName(std::string baseFilename, int numSplit) const;

    bool m_isValid;
    std::string m_name;
    std::string m_type;

    int m_numNodes;
    int m_numCpus;
    std::vector<ExecBase *> m_subModules;

  private:
    //! \name NOT IMPLEMENTED
    //@{
    ExecBase(const ExecBase & that);
    //@}
  };

} // namespace specnets

#endif // __ExecBase_H__
