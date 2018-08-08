#ifndef AMINOACID_H
#define AMINOACID_H

#include "utils.h"
#include "inputParams.h"
#include "mzrange.h"

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>

namespace specnets
{

  extern const double AAmasses[];
  extern const unsigned int AAcount;
  extern const char AAletters[];
  extern const char* AAnames[];

  extern const double AA_ROUNDING;

  class Spectrum;
// defined in spectrum.h

//deprecated, use AAJumps::getPRMMasses
  void getMassesCummulative(const char *sequence,
                            vector<float> &masses,
                            float offset);

//deprecated, use AAJumps::getPRMMasses
  void getMassesCummulativeNoShift(char *sequence, vector<float> &masses);

  /**
   * Transform zero terminated string into vector of amino acid masses.
   *
   *@param sequence
   *@param masses
   */
  void getMasses(char *sequence, vector<float> &masses);

  /**
   * Transform zero terminated string into vector of amino acid masses.
   *
   *@param sequence
   *@param masses
   *
   */
  void getMasses(vector<char> &sequence, vector<float> &masses);

  /**
   * From a c string get all characters that correspond to the mass of an amino acid
   *@param sequence char* to input
   *@param destination char* to output
   *@return
   */
  void getPepSeq(const char* sequence, char* destination);

  /**
   * Transform zero terminated string into vector of amino acid masses.
   *
   *@param sequence
   *@param masses
   *
   */
  void getMasses(vector<char> &sequence, vector<float> &masses);

  /**
   * Get mass of a single char
   *
   *@param aa
   */
  float getMass(char aa);

  /**
   * Define sets of valid mass jumps based on amino acid and modifications masses.
   */
  class AAJumps
  {
  protected:

    /**
     * Global reference amino acid masses, initialized from const AAmasses
     * (see .cpp).
     *
     */
    static vector<double> glbMasses;

    /**
     * Global reference amino acid letters, initialized from const AAletters
     * (see .cpp).
     *
     */
    static vector<char> glbLetters;

    /**
     * Global amino acid modifications.
     */
    static vector<vector<double> > glbMods;

    /**
     * Reference amino acid masses, initialized from const AAmasses (see .cpp)
     */
    vector<double> refMasses;

    /**
     * TODO: add description
     *
     */
    vector<char> refLetters;

    vector<vector<double> > refMods;

    bool m_generateLabels;

    double m_resolution;
  public:

    /*
     static const short NO_MODS=0, USE_MODS=1;
     static const double massHion =  1.0072763,
     massH2O  = 18.010564686,
     minAAmass= 57.0214637230,
     massMH = 18.010564686+1.0072763;
     */

    /**
     * TODO: add description
     *
     */
    static const short NO_MODS;

    /**
     * TODO: add description
     *
     */
    static const short USE_MODS;

    /**
     * TODO: add description
     *
     */
    static const double massHion;

    /**
     * TODO: add description
     *
     */
    static const double minAAmass;

    /**
     * TODO: add description
     *
     */
    static const double massMH;

    /**
     * TODO: add description
     *
     */
    static const double massH2O;

    /**
     * TODO: add description
     *
     */
    static const double massNH3;

    /**
     * TODO: add description
     *
     */
    static const double massCO;

    static const double massNH;

    static const double massC;

    static const double massC13;

    static const double massC_Iso;

    static const double massIsotopeSpace;

    /**
     * Generates a modification string of an amino acid with a modification
     * FORMAT:
     *   (<AA character>,<+/- mod mass>)
     * @param aa
     * @param modMass
     * @return mod string
     */
    static string getModString(char aa, double modMass);

    /**
     * Replaces all modified residues in a peptide string with their unmodified residues
     * @param peptideStr
     * @return un-modified string
     */
    static string stripMods(const string& peptideStr);

    /**
     * Returns true if the peptide string contains a modification or mass gap, false otherwise
     */
    inline static bool isModified(const string& peptideStr)
    {
      return (peptideStr.find('(') != string::npos)
          || (peptideStr.find('[') != string::npos);
    }

    /**
     * Computes the set of all suffix jumps w/ mods given a modified peptide label
     * @param peptideStr
     * @param outputSuffixes
     * @return
     */
    static void getSuffixJumps(const string& peptideStr,
                               vector<string>& outputSuffixes);

    /**
     * Computes the set of all prefix jumps w/ mods given a modified peptide label
     * @param peptideStr
     * @param outputPrefixes
     * @return
     */
    static void getPrefixJumps(const string& peptideStr,
                               vector<string>& outputPrefixes);

    /**
     * Converts a modified peptide label to a vector of single jumps
     * @param peptideStr
     * @param outputJumps
     * @return
     */
    static void getSingleJumps(const string& peptideStr,
                               vector<string>& outputJumps);

    /**
     * Computes the reversed version of a peptide string w/ mods
     * @param peptideStr
     * @return reversed peptide
     */
    static string reversePeptide(const string& peptideStr);

    /**
     * Computes the number of jumps in a peptide (a modified residue equals 1 jump)
     * @param peptideStr
     * @return number of jumps
     */
    static string::size_type getNumJumps(const string& peptideStr);

    /**
     * Initializes the global AAJumps object. This is usually only necessary if generating
     *   AAJumps with labels, in which case the computation can get intensive if done many
     *   times.
     * @param maxJumpSize
     * @param resolution
     * @param useMods
     */
    static void initializeGlobalJumps(short maxJumpSize, double resolution =
                                          0.01,
                                      short useMods = USE_MODS);

    /**
     * Returns a reference to the global AAJumps object, which must have been initialized by initializeGlobalJumps
     */
    static AAJumps& getGlobalJumps();

    /**
     * Returns true if the global jumps have been initialized, false if not
     */
    static bool globalJumpsInitialized();

    /**
     * De-allocates the global AAJumps object
     */
    static void cleanupGlobalJumps();

    /**
     * TODO: add description
     *
     */
    short modsUsed;

    /**
     * List of valid amino acid jumps' masses.
     */
    vector<float> masses;

    /**
     * Replacement for masses when labels are generated
     *   first - label
     *   second - mass
     * This is only filled if generateLabels = true.
     */
    vector<pair<string, double> > jumpLabels;

    //    vector< vector<double> > massesMods; // List of modification masses that can be added to entries in masses
    //    vector<list<string> > letters;

    /**
     * Amino acid letters for the single-residue jumps.
     */
    vector<char> aaLetters;

    /**
     * Position i contains the start/end indices of jumps within tolerance of mass
     * i*resolution (see computeIndex() below)
     *
     * Note: index is not valid for negative-mass jumps
     */
    vector<TwoValues<unsigned int> > index;

    /**
     * Maintains a mapping of every single AA jump to its index in 'masses' or 'jumpLabels'
     */
    map<string, int> massLookup;

    /**
     * Default constructor, equivalent to calling other constructor as AAJumps(0)
     */
    AAJumps();

    /**
     * TODO: add description
     *
     *@param maxJumpSize
     *@param resolution
     *@param peakTol
     *@param useMods
     */
    AAJumps(short maxJumpSize,
            double resolution = 0.01,
            double peakTol = -1,
            short useMods = NO_MODS,
            bool uniqueMasses = true,
            bool generateLabels = false);

    /**
     * TODO: add description
     *
     *@param maxJumpSize
     *@param resolution
     *@param peakTol
     *@param useMods
     */
    void getjumps(short maxJumpSize,
                  double resolution = 0.01,
                  double peakTol = -1,
                  short useMods = NO_MODS,
                  bool uniqueMasses = true,
                  bool generateLabels = false);

    /**
     * Quick test function to check AA refmass
     */
    /**
     * Fills the vector with all reference AA's
     */
    bool getAAref(const char aa, double &mass);

    bool getAllAArefs(vector<pair<char, float> > & returnAAs);

    void getEquivalentIndices(map<char, vector<vector<char> > > & equivalences);

    /**
     * TODO: add description
     *
     *@param maxJumpMass
     *@param resolution
     *@param peakTol
     *@param useMods
     */
    void alljumps(double maxJumpMass, double resolution = 0.01, double peakTol =
                      -1,
                  short useMods = NO_MODS);

    /**
     * TODO: add description
     *
     *@param newJumps
     *@param newNames
     */
    void addJumps(vector<double> &newJumps, vector<char> *newNames = 0);

    /**
     * TODO: add description
     *
     *@param filename
     *@param setGlobal
     *@return
     */
    bool loadJumps(const char *filename, bool setGlobal = false);

    /**
     * Loads jumps possibly containing modifications
     * FORMAT:
     *  <AA character>*=<+/- modification mass>
     * @param filename
     * @param setGlobal if true, set the loaded jumps to the static global jumps
     * @param clearOldAA if true, ignore all old jumps that are not specified in the file
     * @return true if file was loaded succcessfully, false if not
     */
    bool loadJumpsWMods(const char *filename,
                        bool setGlobal = false,
                        bool clearOldAA = false);

    /**
     * Saves the amino acid masses to a file
     *
     *@param filename
     *@return
     */
    bool saveJumps(const char *filename);

    /**
     * TODO: add description
     *
     *@param i
     *@return
     */
    double operator[](const int i) const
    {
      if (m_generateLabels)
      {
        return jumpLabels[i].second;
      }
      else
      {
        return masses[i];
      }
    }

    /**
     * Returns the jump label of the mass at index i
     */
    string getLabel(const unsigned int i) const
    {
      if (m_generateLabels)
      {
        return jumpLabels[i].first;
      }
      else
      {
        cerr << "did not initialize with label generation\n";
        abort();
      }
    }

    /**
     * TODO: add description
     *
     *@return
     */
    const unsigned int size() const
    {
      if (m_generateLabels)
      {
        return jumpLabels.size();
      }
      else
      {
        return masses.size();
      }
    }

    void multiplyMasses(double coefficient);

    /**
     * TODO: add description
     *
     *@param resolution
     */
    void forceUnique(double resolution = 0.01)
    {
      masses = Utils::unique(masses, resolution);
      index.resize(0);
    }

    /**
     * Makes masses = [-masses; masses];
     */
    void forceDoubleSided();

    /**
     * Adds a given mass to the set of valid jumps.
     *
     *@param mass
     */
    void forceJump(double mass)
    {
      masses.push_back(mass);
      sort(masses.begin(), masses.end());
      index.resize(0);
    }

    /**
     * Adds a set of given masses to the set of valid jumps.
     *
     *@param newMasses
     */
    void forceJumps(vector<float> newMasses);

    /**
     * Every mass m in masses is replaced by a set of masses
     * m+[-tolerance:resolution:tolerance].
     *
     *@param tolerance
     *@param resolution
     */
    void forceTolerance(float tolerance,
                        float resolution = 0.01,
                        bool resetRefMasses = false);

    /**
     * TODO: add description
     *
     *@param largestJump
     */
    void removeHigherJumps(float largestJump);

    /**
     * Test whether a given mass is a valid jump in the current set.
     *
     *@param mass
     *@param tolerance
     *@return
     */
    bool isValid(float mass, float tolerance);

    /**
     * TODO: add description
     *
     *@param mass
     *@param tolerance
     *@param idxBounds
     *@return
     */
    unsigned int find(float mass,
                      float tolerance,
                      TwoValues<unsigned int> &idxBounds);

    /**
     * TODO: add description
     *
     *@param peakTol
     *@param resolution
     *@param maxMass
     */
    void computeIndex(double peakTol, double resolution, double maxMass = -1);

    /**
     * Locates AAJumps matching a query mass within tolerance
     * @param mass query
     * @param tolerance tolerance of query
     * @param outputJumps output data structure of matched jumps (first = label, second = mass)
     * @param maxNumMods maximum allowable number of mods to allow per query match
     * @param maxNumJumps maximum allowable number of jumps to allow per query match
     * @param strictNumJumps only allow a certain number of jumps per query match
     * @return
     */
    void findJumpsWLabels(double mass,
                          double tolerance,
                          list<pair<string, double> >& outputJumps,
                          int strictNumJumps = -1,
                          int maxNumJumps = -1,
                          int maxNumMods = -1);

    /** AAJumps::getPRMMasses
     * Calculates PRM masses for input annotation in SpecNets format
     * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
     * @param masses: 	Vector of masses to be modified.
     * @param offset: offset for prefix masses
     * @param addZeroMass: if true add offset mass to beginning of masses
     */
    bool getPRMMasses(const string sequence,
                      vector<float> &masses,
                      const float offset = 0.0,
                      vector<string> * tokens = (vector<string> *)0,
                      const bool addZeroMass = false) const;

    /** AAJumps::getPRMMasses
     * Calculates PRM masses for input annotation in SpecNets format. Same as
     *   getPRMMasses(string,vector<float> &,float,bool) but returning the theoretical
     *   PRM masses in a Spectrum object instead of in a vector<float>
     * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
     * @param prmMassesAsSpectrum:   Output Spectrum object with vector of masses to be modified to contain theoretical PRM masses.
     * @param offset: offset for prefix masses
     * @param addZeroMass: if true add offset mass to beginning of masses
     */
    bool getPRMMasses(const string sequence,
                      Spectrum &spec,
                      const float offset = 0.0,
                      bool addZeroMass = false) const;

    /**
     * Calculates PRM masses as sums of rounded masses
     * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
     * @param masses:   Vector of masses to be modified.
     * @param offset: offset for prefix masses
     * @param addZeroMass: if true add offset mass to beginning of masses
     */
    bool getRoundedPRMMasses(const string &sequence,
                             vector<float> &masses,
                             const float offset = 0.0,
                             vector<string> * tokens = 0,
                             const bool addZeroMass = false) const;

    bool getRoundedPRMMasses(const string sequence,
                             Spectrum &spec,
                             const float offset = 0.0,
                             bool addZeroMass = false) const;

    /**
     * Checks if two peptides overlap using the specified shift
     * @param pep1
     * @param pep2
     * @param shift Da shift of pep2 in relation to pep1
     * @param strictMods If true, AA jumps must match exactly, if false then modified forms of
     *                   the same residue can overlap
     * @return >=0 if both peptides overlap at the specified shift, otherwise -1. 0 means they are the same peptide, 1 means prefix/suffix
     */
    int comparePeptideOverlap(const string& pep1,
                              const string& pep2,
                              float shift,
                              bool strictMods = true) const;

    /** AAJumps::getSRMMasses
     * Calculates SRM masses for input annotation in SpecNets format *.SEQ[-17]UENCE.*
     * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
     * @param masses: 	Vector of masses to be modified.
     * @param offset: offset for suffix masses
     * @param addZeroMass: if true add offset mass to beginning of masses
     */
    bool getSRMMasses(const string sequence,
                      vector<float> &masses,
                      const float offset = 0.0,
                      bool const addZeroMass = false) const;

    /** AAJumps::getPRMandSRMMasses
     * Calculates SRM masses for input annotation in SpecNets format *.SEQ[-17]UENCE.*
     * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
     * @param prm_masses: 	Vector of prefix masses to be modified.
     * @param srm_masses: 	Vector of suffix masses to be modified.
     * @param addZeroMass: if true add offset mass to beginning of masses
     *
     */
    bool getPRMandSRMMasses(const string &sequence,
                            vector<float> &prm_masses,
                            vector<float> &srm_masses,
                            float &peptide_mass,
                            const bool addZeroMass = false) const;

    /** AAJumps::getPeptideFromSpectrum
     * Converts a spectrum to a peptide sequence by instantiating consecutive mass differences
     *   to amino acids; mass differences are reported as "[mass]" whenever indistinguishable between amino acids (e.g., I/L)
     * @param spec: spectrum with the sequence represented as a series of cumulative prefix masses
     * @param sequence: output sequence
     * @param peakTol: mass error tolerance when deciding when a mass difference matches an amino acid (defaults to 0.45 Da)
     * @param offset: offset for prefix masses from theoretical masses (set to 0.0 for PRMs)
     */
    void getPeptideFromSpectrum(const Spectrum &spec,
                                string &sequence,
                                float peakTol = 0.45,
                                float offset = 0.0);

    /** AAJumps::getPeptideMass
     * Calculates sum of all aa masses for input annotation in SpecNets format *.SEQ[-17]UENCE.* NOT PARENT MASS
     * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
     *
     */
    double getPeptideMass(const string &sequence) const;

    /**
     * Calculates the sum of all AA and modified AA residues in sequence
     * @param sequence
     * @return cummulative mass
     */
    double getModPeptideMass(const string &sequence) const;

    /** AAJumps::getPeptideLength
     * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
     */
    int getPeptideLength(const string &sequence) const;

    bool getPeptideShift(const string &pep1,
                         const string &pep2,
                         int minAAOverlap,
                         int *putDaShift,
                         float *putPeakShift,
                         int *putResShift) const;

    string getPeptideSuffix(const string &pep, const int prefix) const;

    string getPeptidePrefix(const string &pep, const int suffix) const;

    /** AAJumps::checkSequence
     * Checks if a sequence is consistent with the SpecNets format *.SEQ[-17]UENCE.*
     * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
     * @return true if sequence is consistent with the SpecNets format
     */
    bool checkSequence(const string &sequence) const;

  private:

    void init(short maxJumpSize,
              double resolution = 0.01,
              double peakTol = -1,
              short useMods = NO_MODS,
              bool uniqueMasses = true,
              bool generateLabels = false);
  };

  /**
   * Finds the longest amino acid tag contained in a sequence of masses
   * @param masses sequence of cummulative masses
   * @param jumps AAJumps class
   * return longest amino acid found by checking consecutive mass differences
   */
  string getLongestTag(vector<MZRange>& masses, AAJumps& jumps);
}
#endif
