#include "aminoacid.h"
#include "Logger.h"
#include "spectrum.h"

namespace specnets
{

  struct SortLabeledJumps : public std::binary_function<pair<string, double>,
      pair<string, double>, bool>
  {
    bool operator()(pair<string, double> left, pair<string, double> right) const
    {
      return left.second < right.second;
    }
    ;
  };

/*
// Original monoisotopic masses
  const double AAmasses[] = {71.037113805, 156.10111105, 115.026943065,
							114.04292747, 103.009184505, 129.042593135,
							128.05857754, 57.021463735, 137.058911875,
							113.084064015, 113.084064015, 128.09496305,
							131.040484645, 147.068413945, 97.052763875,
							87.032028435, 101.047678505, 186.07931298,
							163.063328575, 99.068413945};
*/
 
// block cysteine 57.021463735, 160.03064824
  const double AAmasses[] = {71.037113805, 156.10111105, 115.026943065,
							114.04292747, 160.03064824, 129.042593135,
							128.05857754, 57.021463735, 137.058911875,
							113.084064015, 113.084064015, 128.09496305,
							131.040484645, 147.068413945, 97.052763875,
							87.032028435, 101.047678505, 186.07931298,
							163.063328575, 99.068413945};

/*    // block cysteine + TMT
    const double AAmasses[] = {71.037113805, 156.10111105, 115.026943065,
                               114.04292747, 160.03064824, 129.042593135,
                               128.05857754, 57.021463735, 137.058911875,
                               113.084064015, 113.084064015, 357.257895228,
                               131.040484645, 147.068413945, 97.052763875,
                               87.032028435, 101.047678505, 186.07931298,
                               163.063328575, 99.068413945};*/


  const unsigned int AAcount = 20;
  const char AAletters[] = { 'A', 'R', 'D', 'N', 'C', 'E', 'Q', 'G', 'H', 'I',
                             'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' };
  const char* AAnames[] = { "Alanine", "Arginine", "Aspartic acid",
                            "Asparagine", "Cysteine", "Glutamic acid",
                            "Glutamine", "Glycine", "Histidine", "Isoleucine",
                            "Leucine", "Lysine", "Methionine", "Phenylalanine",
                            "Proline", "Serine", "Threonine", "Tryptophan",
                            "Tyrosine", "Valine" };

  const double AA_ROUNDING = 0.9995;

  static AAJumps* globalJumps = 0;

  void AAJumps::initializeGlobalJumps(short maxJumpSize,
                                      double resolution,
                                      short useMods)
  {
    if (globalJumps != 0)
    {
      delete globalJumps;
    }
    globalJumps = new AAJumps(maxJumpSize,
                              resolution,
                              -1,
                              useMods,
                              false,
                              true);
  }

  AAJumps& AAJumps::getGlobalJumps()
  {
    if (globalJumps != 0)
    {
      return *globalJumps;
    }
    else
    {
      ERROR_MSG("Global jumps have not been initialized");
      abort();
    }
  }

  bool AAJumps::globalJumpsInitialized()
  {
    return globalJumps != 0;
  }

  void AAJumps::cleanupGlobalJumps()
  {
    if (globalJumps != 0)
    {
      delete globalJumps;
      globalJumps = 0;
    }
  }

  string AAJumps::getModString(char aa, double modMass)
  {
    string modStr("(");
    modStr.push_back(aa);
    modStr += ",";
    modStr += parseDouble(modMass, 2);
    modStr += ")";
    return modStr;
  }

  string AAJumps::stripMods(const string& peptideStr)
  {
    string unModStr("");
    string::size_type i = 0;
    while (i < peptideStr.length())
    {
      char c = peptideStr.at(i);

      if (c != '(')
      {
        unModStr += c;
        i++;
        continue;
      }
      else if (i == peptideStr.length() - 1)
      {
        break;
      }

      i++;
      string::size_type j = 0;
      while (j + i < peptideStr.length() && peptideStr.at(j + i) != ',')
      {
        j++;
      }

      unModStr += peptideStr.substr(i, j);

      i = i + j;

      while (i < peptideStr.length() && peptideStr.at(i) != ')')
      {
        i++;
      }
      i++;
    }
    return unModStr;
  }

  void AAJumps::getSuffixJumps(const string& peptideStr,
                               vector<string>& outputSuffixes)
  {
    outputSuffixes.resize(getNumJumps(peptideStr));
    string::size_type i = 0;
    unsigned int idxUse = 0;
    while (i < peptideStr.length())
    {
      char c = peptideStr.at(i);

      if (c != '(')
      {
        outputSuffixes[idxUse++] = peptideStr.substr(i,
                                                     peptideStr.length() - i);
        i++;
        continue;
      }
      else if (i == peptideStr.length() - 1)
      {
        break;
      }

      outputSuffixes[idxUse++] = peptideStr.substr(i, peptideStr.length() - i);

      while (i < peptideStr.length() && peptideStr.at(i) != ')')
      {
        i++;
      }
      i++;
    }
  }

  void AAJumps::getPrefixJumps(const string& peptideStr,
                               vector<string>& outputPrefixes)
  {
    outputPrefixes.resize(getNumJumps(peptideStr));
    string::size_type i = 0;
    unsigned int idxUse = 0;
    while (i < peptideStr.length())
    {
      char c = peptideStr.at(i);

      if (c != '(')
      {
        outputPrefixes[idxUse++] = peptideStr.substr(0, i + 1);
        i++;
        continue;
      }
      else if (i == peptideStr.length() - 1)
      {
        break;
      }

      while (i < peptideStr.length() && peptideStr.at(i) != ')')
      {
        i++;
      }
      i++;

      outputPrefixes[idxUse++] = peptideStr.substr(0, i);
    }
  }

  void AAJumps::getSingleJumps(const string& peptideStr,
                               vector<string>& outputJumps)
  {
    outputJumps.resize(getNumJumps(peptideStr));
    string::size_type i = 0, j = 0;
    unsigned int idxUse = 0;
    while (i < peptideStr.length())
    {
      char c = peptideStr.at(i);

      if (c != '(')
      {
        outputJumps[idxUse++] = peptideStr.substr(i, 1);
        i++;
        continue;
      }
      else if (i == peptideStr.length() - 1)
      {
        break;
      }

      j = 1;
      while (i + j < peptideStr.length() && peptideStr.at(i + j) != ')')
      {
        j++;
      }

      outputJumps[idxUse++] = peptideStr.substr(i, j + 1);
      i += j + 1;
    }
  }

  string AAJumps::reversePeptide(const string& peptideStr)
  {
    string revStr("");
    int i = peptideStr.length() - 1;
    while (i >= 0)
    {
      char c = peptideStr.at(i);

      if (c != ')')
      {
        revStr += c;
        i--;
        continue;
      }
      else if (i == 0)
      {
        break;
      }

      int j = i;
      while (i >= 0 && peptideStr.at(j) != '(')
      {
        j--;
      }

      revStr += peptideStr.substr(j, i - j + 1);

      i = j - 1;
    }
    return revStr;
  }

  string::size_type AAJumps::getNumJumps(const string& peptideStr)
  {
    string::size_type numJumps = 0;
    string::size_type i = 0;
    while (i < peptideStr.length())
    {
      char c = peptideStr.at(i);

      if (c != '(')
      {
        numJumps++;
        i++;
        continue;
      }
      else if (i == peptideStr.length() - 1)
      {
        break;
      }

      i++;
      while (i < peptideStr.length() && peptideStr.at(i) != ')')
      {
        i++;
      }

      numJumps++;
      i++;
    }
    return numJumps;
  }

  vector<double> AAJumps::glbMasses;
  vector<char> AAJumps::glbLetters;
  vector<vector<double> > AAJumps::glbMods;
  const short AAJumps::NO_MODS = 0, AAJumps::USE_MODS = 1;
  //  const double AAJumps::massHion = 1.007825035; // this is H, not ion
  const double AAJumps::massHion = 1.007276035;
  const double AAJumps::minAAmass = 57.0214637230;
  const double AAJumps::massH2O = 18.010564686;
  const double AAJumps::massMH = AAJumps::massH2O + AAJumps::massHion;
  const double AAJumps::massNH3 = 17.026549105;
  const double AAJumps::massCO = 27.99491463;
  const double AAJumps::massNH = 15.010899035;
  const double AAJumps::massC = 12;
  const double AAJumps::massC13 = 13.00335483;
  const double AAJumps::massC_Iso = AAJumps::massC13 - AAJumps::massC;
  const double AAJumps::massIsotopeSpace = 1.00235;

  void getMasses(char *sequence, vector<float> &masses)
  {
    masses.resize(strlen(sequence));
    int m;
    for (unsigned int i = 0; i < masses.size(); i++)
    {
      m = 0;
      while (m < AAcount and AAletters[m] != sequence[i])
        m++;
      if (m >= AAcount)
        masses[i] = 0;
      else
        masses[i] = AAmasses[m];
    }
  }

  void getPepSeq(const char* sequence, char* destination)
  {
    int m;
    int dest_index = 0;
    for (unsigned int i = 0; sequence[i] != 0; i++)
    {
      m = 0;
      while (m < AAcount and AAletters[m] != sequence[i])
        m++;
      if (m < AAcount)
      {
        destination[dest_index] = sequence[i];
        dest_index++;
      }
    }
    destination[dest_index] = 0;
  }

  void getMassesCummulative(const char *sequence,
                            vector<float> &masses,
                            float offset)
  {
    int size = strlen(sequence);
    masses.resize(size + 1);
    float f;
    int m;
    int float_i;
    char* modStringBuffer = (char*)malloc(size + 1);
    //char float_seq[25];
    int spec_index = 1;
    masses[0] = 0.0;
    float total = offset;
    for (unsigned int i = 0; sequence[i] != 0; i++)
    {
      m = 0;
      while (m < AAcount and AAletters[m] != sequence[i])
        m++;
      if (m >= AAcount)
      {
        if (sequence[i] == '[')
        {
          float_i = i + 1;
          for (; sequence[i] != ']'; i++)
            ;
          strncpy(modStringBuffer, &sequence[float_i], i - float_i);
          modStringBuffer[i - float_i] = '\0';
          f = getFloat(modStringBuffer);
          masses[spec_index - 1] = masses[spec_index - 1] + f;
          total += f;
        }
      }
      else
      {
        masses[spec_index] = AAmasses[m] + total;
        total += AAmasses[m];
        spec_index++;
      }
    }
    free(modStringBuffer);
    masses.resize(spec_index);
  }

  void getMassesCummulativeNoShift(char *sequence, vector<float> &masses)
  {
    masses.resize(strlen(sequence) + 1);
    int m;
    int count = 1;
    float total = 0.0;
    masses[0] = 0.0;
    for (unsigned int i = 0; sequence[i] != 0; i++)
    {
      m = 0;
      while (m < AAcount and AAletters[m] != sequence[i])
        m++;
      if (m < AAcount)
      {
        masses[count] += AAmasses[m] + total;
        total += AAmasses[m];
        count++;
      }
    }
    masses.resize(count);
  }

  void getMasses(vector<char> &sequence, vector<float> &masses)
  {
    masses.resize(sequence.size());
    int m;
    for (unsigned int i = 0; i < masses.size(); i++)
    {
      m = 0;
      while (m < AAcount and AAletters[m] != sequence[i])
        m++;
      if (m >= AAcount)
        masses[i] = 0;
      else
        masses[i] = AAmasses[m];
    }
  }

  float getMass(char aa)
  {
    int m = 0;
    while (m < AAcount and AAletters[m] != aa) {
      m++;
    }
    if (m >= AAcount) {
      return 0.0;
    } else {
      return AAmasses[m];
    }
  }

//
// ************************************************************************************************
//  AAJumps::AAJumps(short maxJumpSize, float resolution, short useMods)
// ************************************************************************************************
//
  AAJumps::AAJumps()
  {
    init(0);
  }
  AAJumps::AAJumps(short maxJumpSize,
                   double resolution,
                   double peakTol,
                   short useMods,
                   bool uniqueMasses,
                   bool generateLabels)
  {
    init(maxJumpSize,
         resolution,
         peakTol,
         useMods,
         uniqueMasses,
         generateLabels);
  }

  void AAJumps::init(short maxJumpSize,
                     double resolution,
                     double peakTol,
                     short useMods,
                     bool uniqueMasses,
                     bool generateLabels)
  {
    m_generateLabels = generateLabels;
    m_resolution = resolution;
    unsigned int aaIndex;
    if (glbMasses.empty())
    {
      glbMasses.resize(AAcount);
      glbLetters.resize(AAcount);
      glbMods.resize(AAcount);
      for (aaIndex = 0; aaIndex < AAcount; aaIndex++)
      {
        glbMasses[aaIndex] = AAmasses[aaIndex];
        glbMods[aaIndex].resize(0);
      }
      for (aaIndex = 0; aaIndex < AAcount; aaIndex++)
        glbLetters[aaIndex] = AAletters[aaIndex];
    }
    getjumps(maxJumpSize,
             resolution,
             peakTol,
             useMods,
             uniqueMasses,
             generateLabels);
  }

// ************************************************************************************************
  void AAJumps::getEquivalentIndices(map<char, vector<vector<char> > > & equivalences)
  {
    for (unsigned int m = 0; m < masses.size(); m++)
    {
      int iMass = (int)masses[m];
      vector<vector<char> > vecTemp;
      equivalences[m] = vecTemp;
      for (unsigned int i = 0; i < masses.size(); i++)
      {
        for (unsigned int j = 0; j < masses.size(); j++)
        {
          int iTotal = (int)masses[i] + (int)masses[j];
          if (iTotal == iMass)
          {
            vector<char> vecIndexes(2);
            vecIndexes[0] = i;
            vecIndexes[1] = j;
            equivalences[m].push_back(vecIndexes);
          }
        }
      }
    }

#if 0
    map<char, vector<vector<char> > >::iterator itr = equivalences.begin();
    map<char, vector<vector<char> > >::iterator itr_end = equivalences.end();
    for (; itr != itr_end; itr++)
    {
      for (unsigned int i = 0; i < itr->second.size(); i++)
      {
        DEBUG_MSG((int)itr->first << "  " << aaLetters[itr->first] << "  " <<
            aaLetters[itr->second[i][0]] << "  " << aaLetters[itr->second[i][1]] << "  " <<
            (int)itr->first << "  " << masses[itr->first] << "  " <<
            masses[itr->second[i][0]] << "  " << masses[itr->second[i][1]])
      }
    }
#endif

    return;
  }

//
// ************************************************************************************************
//  void AAJumps::getjumps(short maxJumpSize, double resolution, short useMods)
// ************************************************************************************************
//
  void AAJumps::getjumps(short maxJumpSize,
                         double resolution,
                         double peakTol,
                         short useMods,
                         bool uniqueMasses,
                         bool generateLabels)
  {
    unsigned int step = 0, aaIndex = 0, modIndex = 0, numMasses;
    m_generateLabels = generateLabels;

    if (generateLabels && uniqueMasses)
    {
      ERROR_MSG("Getting labels with unique masses is not yet supported")
      masses.resize(0);
      abort();
    }

    if (maxJumpSize == 0)
    {
      aaLetters.resize(1);
      aaLetters[0] = '0';
      if (m_generateLabels)
      {
        jumpLabels.resize(1);
        jumpLabels[0].first = "";
        jumpLabels[0].second = 0;
      }
      else
      {
        masses.resize(1);
        masses[0] = 0;
      }
      return;
    }

    if (maxJumpSize < 0)
    {
      aaLetters.resize(1);
      if (m_generateLabels)
      {
        jumpLabels.resize(0);
      }
      else
      {
        masses.resize(0);
      }
      return;
    }
    /*
     if (useMods == USE_MODS)
     {
     ERROR_MSG("AAjumps::getjumps: USE_MODS not implemented yet");
     masses.resize(0);
     exit(-1);
     }
     else
     {
     modsUsed = useMods;
     } */

    modsUsed = useMods;

    if (refMasses.empty())
    {
      refMasses.resize(glbMasses.size());
      refLetters.resize(glbMasses.size());
      for (aaIndex = 0; aaIndex < glbMasses.size(); aaIndex++)
      {
        refMasses[aaIndex] = glbMasses[aaIndex];
        refLetters[aaIndex] = glbLetters[aaIndex];
        string letter("");
        letter += refLetters[aaIndex];
      }

      refMods.resize(glbMods.size());
      for (aaIndex = 0; aaIndex < glbMods.size(); aaIndex++)
      {
        refMods[aaIndex] = glbMods[aaIndex];
        for (modIndex = 0; modIndex < refMods[aaIndex].size(); modIndex++)
        {
          string modStr = getModString(refLetters[aaIndex],
                                       refMods[aaIndex][modIndex]);
        }
      }
    }

    int totalMasses = refMasses.size();
    if (useMods == USE_MODS)
    {
      for (aaIndex = 0; aaIndex < refMods.size(); aaIndex++)
      {
        totalMasses += refMods[aaIndex].size();
      }
    }

    int numMassesPerStep = totalMasses, totalNumMasses = numMassesPerStep;
    for (step = 2; step <= maxJumpSize; step++)
    {
      numMassesPerStep *= totalMasses;
      totalNumMasses += numMassesPerStep;
    }

    if (m_generateLabels)
    {
      jumpLabels.resize(totalNumMasses);
    }
    else
    {
      masses.resize(totalNumMasses);
    }
    //cout << "getjumps: masses.size() = " << masses.size() << ", " << endl;

    massLookup.clear();

    // This takes care of step 1
    for (aaIndex = 0; aaIndex < refMasses.size(); aaIndex++)
    {
      if (m_generateLabels)
      {
        jumpLabels[aaIndex].first = refLetters[aaIndex];
        jumpLabels[aaIndex].second = refMasses[aaIndex];
        massLookup[jumpLabels[aaIndex].first] = aaIndex;
      }
      else
      {
        masses[aaIndex] = refMasses[aaIndex];
        string label = "";
        label += refLetters[aaIndex];
        massLookup[label] = aaIndex;
      }
    }
    int idxUse = aaIndex;
    vector<pair<string, double> > allMods(totalMasses - refMasses.size());

    if (useMods == USE_MODS)
    {
      for (aaIndex = 0; aaIndex < refMods.size(); aaIndex++)
      {
        for (modIndex = 0; modIndex < refMods[aaIndex].size(); modIndex++)
        {
          if (m_generateLabels)
          {
            jumpLabels[idxUse].first = getModString(refLetters[aaIndex],
                                                    refMods[aaIndex][modIndex]);
            allMods[idxUse - refMasses.size()].first = jumpLabels[idxUse].first;
            jumpLabels[idxUse].second = refMasses[aaIndex]
                + refMods[aaIndex][modIndex];

            massLookup[jumpLabels[idxUse].first] = idxUse;

            allMods[idxUse - refMasses.size()].second =
                jumpLabels[idxUse].second;
          }
          else
          {
            masses[idxUse] = refMasses[aaIndex] + refMods[aaIndex][modIndex];
            allMods[idxUse - refMasses.size()].second = masses[idxUse];
          }
          idxUse++;
        }
      }
    }

    int prevStepStart = 0, // Start of previous step's jumps
        prevStepEnd = totalMasses - 1, // End of previous step's jumps
        prevIter, // Iterator for previous step's jumps
        curStepStart = totalMasses, // Start of current step's jumps
        curStepIter = curStepStart; // Iterates through current step's jumps

    for (step = 2; step <= maxJumpSize; step++)
    {
      //cerr << "Step = " << step << ": ["<<prevStepStart<<","<<prevStepEnd<<"] -> " << curStepStart << endl;
      for (prevIter = prevStepStart; prevIter <= prevStepEnd; prevIter++)
      {
        for (aaIndex = 0; aaIndex < totalMasses; aaIndex++)
        {
          double addMass =
              (aaIndex < refMasses.size()) ? refMasses[aaIndex] :
                  allMods[aaIndex - refMasses.size()].second;

          if (m_generateLabels)
          {
            jumpLabels[curStepIter].first = jumpLabels[prevIter].first;
            string addLabel;
            if (aaIndex < refMasses.size())
            {
              addLabel = refLetters[aaIndex];
            }
            else
            {
              addLabel = allMods[aaIndex - refMasses.size()].first;
            }
            jumpLabels[curStepIter].first += addLabel;
            jumpLabels[curStepIter].second = jumpLabels[prevIter].second
                + addMass;
          }
          else
          {
            masses[curStepIter] = masses[prevIter] + addMass;
          }
          curStepIter++;
        }
      }
      prevStepStart = curStepStart;
      prevStepEnd = curStepIter - 1;
      curStepStart = curStepIter;
    }

    if (maxJumpSize == 1 && (!m_generateLabels) && useMods == NO_MODS)
    {
      if (uniqueMasses)
      {
        vector<unsigned int> idx;
        Utils::unique(masses, resolution, &idx);
        aaLetters.resize(idx.size());
        for (unsigned int aaIdx = 0; aaIdx < idx.size(); aaIdx++)
          aaLetters[aaIdx] = refLetters[idx[aaIdx]];
      }
      else
      {
        aaLetters.resize(refLetters.size());
        for (unsigned int aaIdx = 0; aaIdx < refLetters.size(); aaIdx++)
          aaLetters[aaIdx] = refLetters[aaIdx];
      }
    }
    else if (!m_generateLabels)
    {
      aaLetters.resize(0);
      Utils::unique(masses, resolution);
    }
    else
    {
      sort(jumpLabels.begin(), jumpLabels.end(), SortLabeledJumps());
      for (unsigned int i = 0; i < jumpLabels.size(); i++)
      {
        if (getNumJumps(jumpLabels[i].first) == 1)
        {
          massLookup[jumpLabels[i].first] = i;
        }
      }
    }

    //cout << "getjumps: masses.size() = " << masses.size() << endl;
    if (m_generateLabels)
    {
      computeIndex(peakTol, resolution);
    }
    else
    {
      if (peakTol < 0)
        index.resize(0);
      else
        computeIndex(peakTol, resolution);
    }
  }

//
// Test function
  bool AAJumps::getAAref(const char aa, double &mass)
  {
    for (int i = 0; i < refLetters.size(); i++)
    {
      if (aa == refLetters[i])
      {
        mass = refMasses[i];
        return true;
      }
    }
    return false;
  }

  bool AAJumps::getAllAArefs(vector<pair<char, float> > & returnAAs)
  {
    for (int i = 0; i < refLetters.size(); i++)
    {
      pair<char, float> newPair = pair<char, float>(refLetters[i],
                                                         refMasses[i]);
      returnAAs.push_back(newPair);
    }
    return (returnAAs.size() != 0);
  }

//
// ************************************************************************************************
//  void AAJumps::alljumps(short maxJumpMass, double resolution, short useMods)
// ************************************************************************************************
//
  void AAJumps::alljumps(double maxJumpMass,
                         double resolution,
                         double peakTol,
                         short useMods)
  {
    unsigned int vecSize = (unsigned int)ceil(maxJumpMass / resolution) + 1,
        jumpIdx, destIdx;
    unsigned int aaIdx;
    vector<bool> *jumpOk = new vector<bool>(vecSize);

    for (jumpIdx = 0; jumpIdx < vecSize; jumpIdx++)
      (*jumpOk)[jumpIdx] = false;
    getjumps(4, resolution, useMods); // Compute exact jump masses up to 4 amino acid jumps
    for (jumpIdx = 0; jumpIdx < masses.size(); jumpIdx++) // and initialize valid jumps
      (*jumpOk)[(int)round(masses[jumpIdx] / resolution)] = true;

    // Add new valid jumps from every valid jump
    unsigned int countOk = 0;
    for (jumpIdx = 0; jumpIdx < vecSize; jumpIdx++)
      if ((*jumpOk)[jumpIdx])
      {
        countOk++;
        for (aaIdx = 0; aaIdx < refMasses.size(); aaIdx++)
        {
          destIdx = jumpIdx
              + (unsigned int)round(refMasses[aaIdx] / resolution);
          if (destIdx < vecSize)
            (*jumpOk)[destIdx] = true;
        }
      }

    masses.resize(countOk);
    destIdx = 0;
    for (jumpIdx = 0; jumpIdx < vecSize; jumpIdx++)
      if ((*jumpOk)[jumpIdx])
        masses[destIdx++] = jumpIdx * resolution;
    if (peakTol < 0)
      index.resize(0);
    else
      computeIndex(peakTol, resolution, maxJumpMass);
    delete jumpOk;
  }

//
// ************************************************************************************************
//  void AAJumps::addJumps(vector<double> &newJumps, vector<char> *newNames = 0)
//
//       Adds new jumps to masses (set of current jumps)
// ************************************************************************************************
//
  void AAJumps::addJumps(vector<double> &newJumps, vector<char> *newNames)
  {
    unsigned int numJumps = masses.size();
    masses.resize(numJumps + newJumps.size());
    aaLetters.resize(masses.size());
    index.resize(0);
    for (unsigned int i = 0; i < newJumps.size(); i++)
    {
      masses[numJumps + i] = newJumps[i];
      if (newNames != (vector<char> *)0 and i < newNames->size())
        aaLetters[numJumps + i] = (*newNames)[i];
      else
        aaLetters[numJumps + i] = 'X';
    }

    // Sort by increasing mass
    vector<pair<double, unsigned int> > massesIdx(masses.size());
    vector<double> oldMasses(masses.size());
    vector<char> oldLetters(masses.size());
    for (unsigned int i = 0; i < masses.size(); i++)
    {
      oldMasses[i] = masses[i];
      oldLetters[i] = aaLetters[i];
      massesIdx[i].first = masses[i];
      massesIdx[i].second = i;
    }
    sort(massesIdx.begin(), massesIdx.end());

    for (unsigned int i = 0; i < masses.size(); i++)
    {
      masses[i] = oldMasses[massesIdx[i].second];
      aaLetters[i] = oldLetters[massesIdx[i].second];
    }
  }

//
// ************************************************************************************************
//  void AAJumps::loadJumps(char *filename)
//
//       Adds new jumps to masses (set of current jumps)
// ************************************************************************************************
//
  bool AAJumps::loadJumps(const char *filename, bool setGlobal)
  {
    InputParams aaMassParams;
    if (not aaMassParams.readParams(filename))
      return false;
    //cerr<<" -- Loaded "<<filename<<"\n"; cerr.flush();

    refMasses.resize(aaMassParams.size());
    refLetters.resize(aaMassParams.size());
    refMods.resize(aaMassParams.size());
    const char *paramName;
    for (unsigned int i = 0; i < aaMassParams.size(); i++)
    {
      refMods[i].resize(0);
      refMasses[i] = aaMassParams.getValueDouble(i);
      paramName = aaMassParams.getParamName(i);
      if (strlen(paramName) > 0)
        refLetters[i] = paramName[0];
      else
        refLetters[i] = '-';
      //cerr<<" == Got -"<<refLetters[i]<<"- with mass "<<refMasses[i]<<"\n"; cerr.flush();
      //printf(" == Got - %c - with mass %.9f\n",refLetters[i],refMasses[i]);
    }
    //cerr<<" -- Done reading file\n"; cerr.flush();
    if (setGlobal)
    {
      glbMasses = refMasses;
      glbLetters = refLetters;
      glbMods = refMods;
    }
    addJumps(refMasses, &refLetters);
    //cerr<<" -- Returning\n"; cerr.flush();
    return true;
  }

  bool AAJumps::loadJumpsWMods(const char *filename,
                               bool setGlobal,
                               bool clearOldAA)
  {
    InputParams aaMassParams;
    if (not aaMassParams.readParams(filename))
      return false;

    map<char, double> aaMasses;

    if (!clearOldAA)
    {
      for (unsigned int i = 0; i < refMasses.size(); i++)
      {
        aaMasses[refLetters[i]] = refMasses[i];
      }
    }

    map<char, list<double> > aaMods;
    list<double> modList;

    const char *paramName;
    for (unsigned int i = 0; i < aaMassParams.size(); i++)
    {
      double massJump = aaMassParams.getValueDouble(i);
      paramName = aaMassParams.getParamName(i);
      if (strlen(paramName) == 1)
      {
        aaMasses[paramName[0]] = massJump;
      }
      else if (strlen(paramName) == 2 && paramName[1] == '*')
      {
        if (aaMods.count(paramName[0]) == 0)
        {
          modList.clear();
          modList.push_back(massJump);
          aaMods[paramName[0]] = modList;
        }
        else
        {
          aaMods[paramName[0]].push_back(massJump);
        }
      }
      //cerr<<" == Got -"<<refLetters[i]<<"- with mass "<<refMasses[i]<<"\n"; cerr.flush();
      //printf(" == Got - %c - with mass %.9f\n",refLetters[i],refMasses[i]);
    }
    refMasses.resize(aaMasses.size());
    refLetters.resize(aaMasses.size());
    refMods.resize(aaMasses.size());
    unsigned int idxUse = 0;
    map<char, unsigned int> charToIdx;
    for (map<char, double>::iterator massIt = aaMasses.begin();
        massIt != aaMasses.end(); massIt++)
    {
      refMasses[idxUse] = massIt->second;
      refLetters[idxUse] = massIt->first;
      refMods[idxUse].resize(0);
      charToIdx[massIt->first] = idxUse;
      idxUse++;
    }
    for (map<char, list<double> >::iterator massIt = aaMods.begin();
        massIt != aaMods.end(); massIt++)
    {
      unsigned int refIdx = charToIdx[massIt->first];
      refMods[refIdx].resize(massIt->second.size());
      idxUse = 0;
      for (list<double>::iterator modIt = massIt->second.begin();
          modIt != massIt->second.end(); modIt++)
      {
        refMods[refIdx][idxUse++] = *modIt;
      }
    }

    //cerr<<" -- Done reading file\n"; cerr.flush();
    if (setGlobal)
    {
      glbMasses = refMasses;
      glbLetters = refLetters;
      glbMods = refMods;
    }

    //cerr<<" -- Returning\n"; cerr.flush();
    return true;

  }

// ************************************************************************************************
//
//  AAJumps::saveJumps(const char *filename)
//
// ************************************************************************************************
//
  bool AAJumps::saveJumps(const char *filename)
  {
    ofstream ofs(filename, ios_base::out | ios_base::binary);

    if (!ofs)
    {
      ERROR_MSG("Error opening " << filename);
      return false;
    }
    if (refMasses.size() == 0)
    {
      ERROR_MSG("No amino acid masses to save in saveJumps()");
      return false;
    }

    ofs << refMasses.size() << endl;
    for (int i = 0; i < refMasses.size(); i++)
    {
      ofs << refLetters[i] << "=" << refMasses[i] << endl;
    }

    return true;
  }

  void AAJumps::multiplyMasses(double coefficient)
  {
    for (unsigned int i = 0; i < refMasses.size(); i++)
    {
      refMasses[i] *= coefficient;
    }
    for (unsigned int i = 0; i < refMods.size(); i++)
    {
      for (unsigned int j = 0; j < refMods[i].size(); j++)
      {
        refMods[i][j] *= coefficient;
      }
    }
    if (m_generateLabels)
    {
      for (unsigned int i = 0; i < jumpLabels.size(); i++)
      {
        jumpLabels[i].second *= coefficient;
      }
    }
    else
    {
      for (unsigned int i = 0; i < masses.size(); i++)
      {
        masses[i] *= coefficient;
      }
    }
  }

//
// ************************************************************************************************
//  void AAJumps::forceDoubleSided()
//
//       Makes masses = [-masses; masses];
// ************************************************************************************************
//
  void AAJumps::forceDoubleSided()
  {
    vector<float> tmp(masses);

    masses.resize(2 * masses.size());
    index.resize(0);
    for (unsigned int i = 0; i < tmp.size(); i++)
    {
      masses[i] = -tmp[i];
      masses[i + tmp.size()] = tmp[i];
    }
  }

//
// ************************************************************************************************
//  void AAJumps::forceJump(double mass)
//
//      Adds a given mass to the set of valid jumps
// ************************************************************************************************
//
  void AAJumps::forceJumps(vector<float> newMasses)
  {
    masses.reserve(masses.size() + newMasses.size());
    index.resize(0);
    for (unsigned int i = 0; i < newMasses.size(); i++)
      masses.push_back(newMasses[i]);
    sort(masses.begin(), masses.end());
  }

//
// ************************************************************************************************
//  void AAJumps::forceTolerance(double tolerance, double resolution)
//
//      Every mass m in masses is replaced by a set of masses m+[-tolerance:resolution:tolerance]
// ************************************************************************************************
//
  void AAJumps::forceTolerance(float tolerance,
                               float resolution,
                               bool resetRefMasses)
  {
    double iter;

    int szRange = 0;
    for (iter = -tolerance; iter <= tolerance + 0.00001; iter += resolution)
      szRange++;
    if (szRange == 0)
      return;

    vector<float> tmp(masses);
    masses.resize(masses.size() * szRange);
    index.resize(0);
    for (unsigned int i = 0, curMass = 0; i < tmp.size(); i++)
      for (iter = -tolerance; iter <= tolerance + 0.00001; iter += resolution)
        masses[curMass++] = tmp[i] + iter;

    if (resetRefMasses)
    {
      for (unsigned int i = 0; i < refMasses.size(); i++)
        refMasses[i] = round(refMasses[i] / resolution) * resolution;
    }
  }

//
// ************************************************************************************************
//  void AAJumps::removeHigherJumps(double highestMass)
//
//      Removes all jumps of mass higher than largestJump
// ************************************************************************************************
//
  void AAJumps::removeHigherJumps(float largestJump)
  {
    unsigned int idx;
    for (idx = 0; idx < masses.size() and masses[idx] <= largestJump; idx++)
      ;
    masses.resize(idx);
  }

//
// ************************************************************************************************
//  bool AAJumps::isValid(double mass, double tolerance)
//
//      Test whether a given mass is a valid jump in the current set
// ************************************************************************************************
//
  bool AAJumps::isValid(float mass, float tolerance)
  {
    if (masses.size() == 0)
      return false;
    unsigned int curIdx = (masses.size()) / 2, high = masses.size(), low = 0;

    while (low < high)
    {
      if (fabs(masses[curIdx] - mass) <= tolerance + 0.000001)
        return true;
      if (masses[curIdx] < mass)
        low = curIdx + 1;
      else
        high = curIdx; // High is first non-checked position, low is lowest eligible position
      curIdx = (int)(high + low) / 2;
      //		szIdxJump = (int)abs(prev-curIdx)/2;   prev = curIdx;
      //		if(masses[curIdx]<mass) curIdx+=szIdxJump; else curIdx-=szIdxJump;
    }
    if (curIdx < masses.size()
        and fabs(masses[curIdx] - mass) <= tolerance + 0.0001)
      return true;
    //cerr<<"Closest was "<<masses[curIdx]<<endl;
    return false;
  }

//
// ************************************************************************************************
//  unsigned int AAJumps::find(double mass, double tolerance, TwoValues<unsigned int> &idxBounds)
//
//      Find the indices of all jumps within tolerance of mass. On exit, idxBounds[0] is the
//    index of the first matching mass and idxBounds[1] is the index of the last match plus one.
// ************************************************************************************************
//
  unsigned int AAJumps::find(float mass,
                             float tolerance,
                             TwoValues<unsigned int> &idxBounds)
  {
    idxBounds.set(0, 0);
    if (masses.size() == 0)
      return 0;

    unsigned int curIdx = (masses.size()) / 2, high = masses.size(), low = 0;
    bool found = false;
    while (low < high and not found)
    {
      if (fabs(masses[curIdx] - mass) <= tolerance + 0.000001)
      {
        found = true;
        break;
      }
      if (masses[curIdx] < mass)
        low = curIdx + 1;
      else
        high = curIdx; // High is first non-checked position, low is lowest eligible position
      curIdx = (int)(high + low) / 2;
    }
    if (not found)
      return 0;

    for (; curIdx > 0 and fabs(masses[curIdx] - mass) <= tolerance + 0.000001;
        curIdx--)
      ; // Find first match
    if (fabs(masses[curIdx] - mass) > tolerance + 0.000001)
      curIdx++;
    idxBounds[0] = curIdx;
    for (;
        curIdx < masses.size() and masses[curIdx] - mass <= tolerance + 0.000001;
        curIdx++)
      ; // Find first index after last match
    idxBounds[1] = curIdx;

    return idxBounds[1] - idxBounds[0];
  }

  void AAJumps::computeIndex(double peakTol, double resolution, double maxMass)
  {
    m_resolution = resolution;
    if (!m_generateLabels)
    {
      if (maxMass < 0)
        maxMass = masses[masses.size() - 1];
      maxMass += peakTol + resolution;
      unsigned int jStart, jEnd; // start/end indices of the jumps within tolerance of the current mass
      index.resize((unsigned int)ceil(maxMass / resolution));

      // Initialize jStart/jEnd
      for (jStart = 0; jStart < masses.size() and masses[jStart] <= 0; jStart++)
        ;
      if (jStart == masses.size())
      {
        index.resize(0);
        return;
      }
      jEnd = jStart;

      // Fill-out index
      float curMass;
      for (unsigned int idxMass = 0; idxMass < index.size(); idxMass++)
      {
        curMass = ((float)idxMass) * resolution;
        while (jStart < masses.size() and masses[jStart] < curMass - peakTol)
          jStart++;
        jEnd = max(jEnd, jStart);
        while (jEnd < masses.size() and masses[jEnd] <= curMass + peakTol)
          jEnd++;
        index[idxMass].set(jStart, jEnd);
      }
    }
    else
    {
      if (maxMass < 0)
      {
        maxMass = jumpLabels[jumpLabels.size() - 1].second;
      }
      maxMass += m_resolution;
      unsigned int jStart, jEnd; // start/end indices of the jumps within tolerance of the current mass
      index.resize((unsigned int)ceil(maxMass / m_resolution));

      // Initialize jStart/jEnd
      for (jStart = 0;
          jStart < jumpLabels.size() and jumpLabels[jStart].second <= 0;
          jStart++)
        ;
      if (jStart == jumpLabels.size())
      {
        index.resize(0);
        return;
      }
      jEnd = jStart;

      for (unsigned int idxMass = 0; idxMass < index.size(); idxMass++)
      {
        while (jStart < jumpLabels.size()
            and ((unsigned int)(doubleToInt(jumpLabels[jStart].second
                / m_resolution))) < idxMass)
        {
          jStart++;
        }

        jEnd = max(jEnd, jStart);
        while (jEnd < jumpLabels.size()
            and ((unsigned int)(doubleToInt(jumpLabels[jEnd].second
                / m_resolution))) <= idxMass)
          jEnd++;
        jEnd = (jEnd > 0) ? jEnd - 1 : jEnd;
        index[idxMass].set(jStart, jEnd);
      }

    }
  }

  void AAJumps::findJumpsWLabels(double mass,
                                 double tolerance,
                                 list<pair<string, double> >& outputJumps,
                                 int strictNumJumps,
                                 int maxNumJumps,
                                 int maxNumMods)
  {
    if (mass < 0 || tolerance < 0)
    {
      ERROR_MSG("both mass and tolerance must be positive");
      abort();
    }
    if (!m_generateLabels)
    {
      ERROR_MSG("AAJumps was not initialized to generate labels");
      abort();
    }

    outputJumps.clear();

    int startMass = floatToInt((mass - tolerance) / m_resolution);
    if (startMass < 0)
    {
      startMass = 0;
    }

    if (startMass >= index.size()
        || jumpLabels[index[startMass][0]].second > mass + tolerance)
    {
      return;
    }

    int endMass = floatToInt((mass + tolerance) / m_resolution);
    if (endMass >= index.size())
    {
      endMass = index.size() - 1;
    }

    pair<string, double> nextJump;
    for (unsigned int i = index[startMass][0]; i <= index[endMass][1]; i++)
    {
      if (maxNumMods >= 0)
      {
        int numMods = 0;
        size_t found = jumpLabels[i].first.find(",");
        while (found != string::npos)
        {
          numMods++;
          found = jumpLabels[i].first.find(",", found + 1);
        }

        if (numMods > maxNumMods)
        {
          continue;
        }
      }

      if (maxNumJumps >= 0 && getNumJumps(jumpLabels[i].first) > maxNumJumps)
      {
        continue;
      }

      if (strictNumJumps >= 0
          && getNumJumps(jumpLabels[i].first) != strictNumJumps)
      {
        continue;
      }

      outputJumps.push_back(jumpLabels[i]);
    }

  }

  /** AAJumps::checkSequence
   * Checks if a sequence is consistent with the SpecNets format *.SEQ[-17]UENCE.*
   * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
   * @return true if sequence is consistent with the SpecNets format
   */
  bool AAJumps::checkSequence(const string &sequence) const
  {
    int nPar = 0, nRect = 0;
    int nCharsInMass = 0, nCharsInSeq = 0, nChars = 0;
    bool hasComma = true;

    string numbers = "0123456789";
    string aas = "ACDEFGHIKLMNPQRSTVWY";

    for (int i = 0; i < sequence.length(); i++)
    {

      switch (sequence[i])
      {

      case '[':
        if (nPar > 0 || nRect > 0)
          return false;
        nRect++;
        nCharsInMass = 0;
        nCharsInSeq = 0;
        break;

      case ']':
        if (nPar > 0 || nRect != 1 || nCharsInMass == 0 || nCharsInSeq != 0)
          return false;
        nRect--;
        break;

      case '(':
        if (nPar > 0 || nRect > 0)
          return false;
        nPar++;
        nCharsInMass = 0;
        nCharsInSeq = 0;
        hasComma = false;
        break;

      case ')':
        if (nPar != 1 || nRect > 0)
          return false;
        if (hasComma && (nCharsInMass == 0 || nCharsInSeq != 0))
          return false;
        if (!hasComma && (nCharsInMass != 0 || nCharsInSeq == 0))
          return false;
        nPar--;
        break;

      case ',':
        if (hasComma)
          return false;
        if (nCharsInSeq == 0)
          return false;
        if (nCharsInMass != 0)
          return false;
        nCharsInSeq = 0;
        hasComma = true;
        break;

      case '+':
      case '-':
        if (!((hasComma && nPar > 0) || nRect > 0))
          return false;
        if (nCharsInMass != 0)
          return false;
        break;

      case '.':
        break;

      default:

        int found = numbers.find_first_of(sequence[i]);
        if (found != std::string::npos)
        {
          nCharsInMass++;
          break;
        }

        char aux = toupper((unsigned char)sequence[i]);
        found = aas.find_first_of(aux);
        if (found != std::string::npos)
        {
          nCharsInSeq++;
          break;
        }

        return false;
      }

    }

    if (nPar > 0 || nRect > 0)
      return false;

    return true;
  }

  /** AAJumps::getPRMMasses
   * Calculates PRM masses for input annotation in SpecNets format *.SEQ[-17]UENCE.*
   * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
   * @param masses: 	Vector of masses to be modified.
   * @param offset: offset for prefix masses
   */
  bool AAJumps::getPRMMasses(const string sequence,
                             vector<float> &masses,
                             const float offset,
                             vector<string> * tokens,
                             const bool addZeroMass) const
  {
    bool sequenceOK = checkSequence(sequence);

    const char * sequence_ptr = sequence.c_str();

    int size = strlen(sequence_ptr);

    masses.resize(size + 2);

    float f;
    int float_i;
    char * float_seq = new char[sequence.size() + 1];
    int spec_index = 0;
    float total = offset;
    int start = 2;

    float nTermMod = 0;

    if (addZeroMass)
    {
      masses[spec_index] = total;
      spec_index++;
    }

    if (sequence_ptr[1] != '.')
    {
      start = 0;
    }

    //note: we include the total peptide mass in this vector, which is a bogus
    //"prefix mass", but is useful for generating suffix ions. Never use the
    //last ion for any prefix ion calculations!
    for (unsigned int i = start; i < size - start; i++)
    {
      // find global matching aa for curr sequence aa
      int max_aa = refLetters.size();

      // Square brackets: an unknown mass
      if ('[' == sequence_ptr[i])
      {
        int tokenStart = i;

        float_i = i + 1;
        for (; sequence_ptr[i] != ']'; i++)
          ; //iterate to closing bracket
        strncpy(float_seq, &sequence_ptr[float_i], i - float_i);
        float_seq[i - float_i] = '\0';

        f = getFloat(float_seq);
        //if (spec_index > 0) {
        // if we have a modification at the very beginning of the peptide, we don't want it to appear as separate mass
        //masses[spec_index - 1] = masses[spec_index - 1] + f;
        //}
        // Square brackets define a fixed mass for which we don't know the sequence
        masses[spec_index] = f + total; // added
        spec_index++; // added
        total += f;

        if (tokens)
        {
          string token(sequence.substr(tokenStart, i - tokenStart + 1));
          tokens->push_back(token);
        }
        // round brackets: a modification
      }
      else if (sequence_ptr[i] == '(')
      {
        int tokenStart = i;

        int j;
        if (i + 1 < size && sequence_ptr[i + 1] == '(')
        { // check if we have a N-term modification inside another modification
          float_i = i + 1;
          for (; sequence_ptr[i] != ')'; i++)
            ; //iterate to closing bracket
          for (j = i; j > 0 && sequence_ptr[j] != ','; j--)
            ;
          if (j >= float_i && i - j >= 2)
          {
            // define space for value string
            char * float_seq2 = new char[sequence.size() + 1];
            // get the string
            strncpy(float_seq2, &sequence_ptr[j + 1], i - j - 1);
            // add string terminator
            float_seq2[i - j - 1] = '\0';
            // get the value
            nTermMod += getFloat(float_seq2);
            delete[] float_seq2;
          }
          continue;
        }

        float_i = i + 1;

        for (; sequence_ptr[i] != ')'; i++)
          ; //iterate to closing bracket

        strncpy(float_seq, &sequence_ptr[float_i], i - float_i);
        float_seq[i - float_i] = '\0';

        // find comma
        for (j = 0; (j < i - float_i) && (float_seq[j] != ','); j++)
          ;

        // Parse aa's in sequence
        float massAux = 0.0;
        // Cycle through all amino acids before the comma
        for (int k = 0; k < j; k++)
        {
          // skip spaces in the middle of the sequence
          if (float_seq[k] == ' ' || float_seq[k] == '\t')
            continue;
          // initial aa index in aa table
          int aa_index = 0;
          // search for the aa in aa table
          while (aa_index < max_aa and refLetters[aa_index] != float_seq[k])
            aa_index++;
          // if not found, issue message and skip
          if (aa_index >= max_aa)
          {
            WARN_MSG("Warning1: unknown aa! " << float_seq[k]);
            // otherwise, add it's mass
          }
          else
          {
            massAux += refMasses[aa_index];
          }
        }
        // Now process the modification value - number after the comma
        // initialize modification mass
        f = 0.0;
        // if there is a comma, we have a modification
        if (j < i - float_i)
        {
          // define space for value string
          char * float_seq2 = new char[sequence.size() + 1];
          // get the string
          strncpy(float_seq2, &float_seq[j + 1], i - float_i - j);
          // add string terminator
          float_seq2[i - float_i - j] = '\0';
          // get the value
          f = getFloat(float_seq2);
          delete[] float_seq2;
        }
        // add aa sequence mass value plus modification
        masses[spec_index] = massAux + f + total + nTermMod; // added
        total += massAux + f;
        spec_index++;

        //cout << "i2 " << i << endl;

        if (tokens)
        {
          string token(sequence.substr(tokenStart, i - tokenStart + 1));
          tokens->push_back(token);
        }

        // An aminoacid.
      }
      else
      {
        //cout << "i3 " << i << endl;

        int aa_index = 0;
        while (aa_index < max_aa and refLetters[aa_index] != sequence[i])
          aa_index++;
        if (aa_index >= max_aa)
        {
          ERROR_MSG( "Warning2: unknown aa! " << sequence[i] << " in sequence " << sequence);
        }
        else
        {
          masses[spec_index] = refMasses[aa_index] + total;
          total += refMasses[aa_index];
          spec_index++;

          if (tokens)
          {
            string token("X");
            token[0] = sequence[i];
            tokens->push_back(token);
          }

        }
      }
    }
    masses.resize(spec_index);
    delete[] float_seq;
    return sequenceOK;
  }

// -------------------------------------------------------------------------

  bool AAJumps::getPRMMasses(const string sequence,
                             Spectrum &prmMassesAsSpectrum,
                             float offset,
                             bool addZeroMass) const
  {
    vector<float> masses;
    bool sequenceOK = getPRMMasses(sequence,
                                   masses,
                                   offset,
                                   (vector<string> *)0,
                                   addZeroMass);

    prmMassesAsSpectrum.resize(masses.size());
    if (masses.empty())
      return sequenceOK;
    for (unsigned int idxPeak = 0; idxPeak < prmMassesAsSpectrum.size();
        idxPeak++)
    {
      prmMassesAsSpectrum[idxPeak].set(masses[idxPeak], 1.0);
    }
    float pm = masses[masses.size() - 1] + massMH;
    if (prmMassesAsSpectrum.parentMZ > massMH)
      prmMassesAsSpectrum.parentCharge = (short)round(pm
          / prmMassesAsSpectrum.parentMZ);
    prmMassesAsSpectrum.parentMass = pm;

    return sequenceOK;
  }

  bool AAJumps::getRoundedPRMMasses(const string &sequence,
                                    vector<float> &masses,
                                    const float offset,
                                    vector<string> * tokens,
                                    const bool addZeroMass) const
  {
    vector<float> prmMasses;
    vector<string> singleAAJumps;

    if (!getPRMMasses(sequence, prmMasses, offset, &singleAAJumps, addZeroMass))
    {
      return false;
    }

    masses.assign(prmMasses.size(), 0);
    int prmIdx = 0;

    if (addZeroMass)
    {
      masses[prmIdx] = prmMasses[prmIdx];
      prmIdx++;
    }
    float total = 0;
    for (int jumpIdx = 0; jumpIdx < sequence.length(); jumpIdx++)
    {
      const string &curJumpStr = singleAAJumps[jumpIdx];
      const unsigned int &aa = massLookup.at(curJumpStr);

      masses[prmIdx] = total + round(this->operator [](aa));
      total = masses[prmIdx];
      prmIdx++;
    }

    if (tokens != 0)
    {
      tokens->operator =(singleAAJumps);
    }
    return true;
  }

  bool AAJumps::getRoundedPRMMasses(const string sequence,
                                    Spectrum &spec,
                                    const float offset,
                                    bool addZeroMass) const
  {
    vector<float> masses;
    bool sequenceOK = getRoundedPRMMasses(sequence,
                                          masses,
                                          offset,
                                          (vector<string> *)0,
                                          addZeroMass);

    spec.resize(masses.size());
    if (masses.empty())
      return sequenceOK;
    for (unsigned int idxPeak = 0; idxPeak < spec.size(); idxPeak++)
    {
      spec[idxPeak].set(masses[idxPeak], 1.0);
    }
    float pm = masses[masses.size() - 1] + massMH;
    if (spec.parentMZ > massMH)
      spec.parentCharge = (short)round(pm / spec.parentMZ);
    spec.setParentMass(pm);

    return sequenceOK;
  }

  int AAJumps::comparePeptideOverlap(const string& pep1,
                                     const string& pep2,
                                     float shift,
                                     bool strictMods) const
  {
    if (pep1.size() == 0 || pep2.size() == 0)
    {
      return -1;
    }

    const string &pepUse1 = (shift > -0.01) ? pep1 : pep2;
    const string &pepUse2 = (shift > -0.01) ? pep2 : pep1;
    float shiftUse = (shift > -0.01) ? shift : 0.0 - shift;

    vector<float> prmMasses1;
    vector<string> aaJumps1;

    vector<float> prmMasses2;
    vector<string> aaJumps2;

    getPRMMasses(pepUse1, prmMasses1, 0, &aaJumps1, true);
    getPRMMasses(pepUse2, prmMasses2, 0, &aaJumps2, true);

    unsigned int startAAIdx = 0;
    while (startAAIdx < prmMasses1.size()
        && prmMasses1[startAAIdx] + 2.0 < shiftUse)
    {
      startAAIdx++;
    }

    if (startAAIdx >= aaJumps1.size())
    {
      return -1;
    }

    unsigned int overlapLen = min((unsigned int)aaJumps2.size(),
                                  (unsigned int)aaJumps1.size() - startAAIdx);

    string overlapPep1("");
    string overlapPep2("");
    for (unsigned int i = 0; i < overlapLen; i++)
    {
      overlapPep1 += aaJumps1[i + startAAIdx];
      overlapPep2 += aaJumps2[i];
    }

    if (!strictMods)
    {
      overlapPep1 = AAJumps::stripMods(overlapPep1);
      overlapPep2 = AAJumps::stripMods(overlapPep2);
    }

    if (overlapPep1.compare(overlapPep2) != 0)
    {
      return -1;
    }

    if (startAAIdx == 0 && aaJumps1.size() == aaJumps2.size())
    {
      return 0;
    }
    else
    {
      return 1;
    }
  }

// -------------------------------------------------------------------------

  /** AAJumps::getSRMMasses
   * Calculates SRM masses for input annotation in SpecNets format *.SEQ[-17]UENCE.*
   * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
   * @param masses: 	Vector of masses to be modified.
   * @param offset: offset for suffix masses
   */
  bool AAJumps::getSRMMasses(const string sequence,
                             vector<float> &masses,
                             float offset,
                             bool addZeroMass) const
  {
    bool sequenceOK = checkSequence(sequence);

    vector<float> prm_masses; //this is a bit of a cheat, but I'm lazy. Calculate SRM from PRM masses

    getPRMMasses(sequence, prm_masses, addZeroMass);

    int prm_length = prm_masses.size();

    masses.resize(prm_length - 1);
    int i;
    int mass_i = 0;

    for (i = prm_length - 2; i >= 0; i--)
    {
      masses[mass_i] = prm_masses[prm_length - 1] - prm_masses[i] + offset; //total mass of peptide minus masses starting from end.
      mass_i++;
    }

    return sequenceOK;
  }

  /** AAJumps::getPRMandSRMMasses
   * Calculates SRM masses for input annotation in SpecNets format *.SEQ[-17]UENCE.*
   * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
   * @param prm_masses: 	Vector of prefix masses to be modified.
   * @param srm_masses: 	Vector of suffix masses to be modified.
   * @param peptide_mass: Summed mass of all amino acids. NOT THE SAME AS PRECURSOR MASS.
   *
   */
  bool AAJumps::getPRMandSRMMasses(const string &sequence,
                                   vector<float> &prm_masses,
                                   vector<float> &srm_masses,
                                   float &peptide_mass,
                                   bool addZeroMass) const
  {
    bool sequenceOK = getPRMMasses(sequence, prm_masses, addZeroMass);

    int prm_length = prm_masses.size();

    peptide_mass = prm_masses[prm_length - 1];

    srm_masses.resize(prm_length - 1);
    int i;
    int mass_i = 0;

    for (i = prm_length - 2; i >= 0; i--)
    {
      srm_masses[mass_i] = prm_masses[prm_length - 1] - prm_masses[i]; //total mass of peptide minus masses starting from end.
      mass_i++;
    }

    return sequenceOK;
  }

  /** AAJumps::getPeptideFromSpectrum
   * Converts a spectrum to a peptide sequence by instantiating consecutive mass differences
   *   to amino acids; I/L are always reported as L, Q/K are reported as "[mass]" whenever indistinguishable
   * @param spec: spectrum with the sequence represented as a series of cumulative prefix masses
   * @param sequence: output sequence
   * @param peakTol: mass error tolerance when deciding when a mass difference matches an amino acid
   * @param offset: offset for prefix masses from theoretical masses (set to 0.0 for PRMs)
   */
  void AAJumps::getPeptideFromSpectrum(const Spectrum &spec,
                                       string &sequence,
                                       float peakTol,
                                       float offset)
  {
    float massDiff;
    TwoValues<unsigned int> idxAAMasses; // Index of lowest/highest AA mass matched to massDiff
    ostringstream newSeq;
    Spectrum epSpec; // Version of spec guaranteed to have EndPoints peaks at mass zero and parent mass
    epSpec = spec;
    epSpec.addZPMpeaks(peakTol, offset, false);

    sequence = "";
    for (unsigned int idxPeak = 1; idxPeak < epSpec.size(); idxPeak++)
    {
      massDiff = epSpec[idxPeak][0] - epSpec[idxPeak - 1][0];
      find(massDiff, peakTol, idxAAMasses);
      if (idxAAMasses[1] - idxAAMasses[0] != 1)
        newSeq << "[" << massDiff << "]";
      else
        newSeq << aaLetters[idxAAMasses[0]];
    }
    sequence = newSeq.str();
  }

  /** AAJumps::getPeptideMass
   * Calculates sum of all aa masses for input annotation in SpecNets format *.SEQ[-17]UENCE.* NOT PARENT MASS
   * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
   */
  double AAJumps::getPeptideMass(const string &sequence) const
  {
    vector<float> masses;
    float offset = 0;
    getPRMMasses(sequence, masses, 0);
    return masses[masses.size() - 1];
  }

  double AAJumps::getModPeptideMass(const string &sequence) const
  {
    vector<string> singleJumps;
    getSingleJumps(sequence, singleJumps);
    double totalMass = 0;
    for (unsigned int i = 0; i < singleJumps.size(); i++)
    {
      if (masses.size() > 0)
      {
        totalMass += masses[massLookup.at(singleJumps[i])];
      }
      else
      {
        totalMass += jumpLabels[massLookup.at(singleJumps[i])].second;
      }
    }
    return totalMass;
  }

  int AAJumps::getPeptideLength(const string &sequence) const
  {
    const char * sequence_ptr = sequence.c_str();

    int size = strlen(sequence_ptr);
    int float_i;
    char* modStringBuffer = (char*)malloc(size + 1);
    //char float_seq[25];

    int spec_index = 0;
    int start = 2;

    if (sequence_ptr[1] != '.')
    {
      start = 0;
    }

    for (unsigned int i = start; i < size - start; i++)
    {
      if ('[' == sequence_ptr[i])
      {
        for (; sequence_ptr[i] != ']'; i++)
          ;
        spec_index++;
      }
      else if (sequence_ptr[i] == '(')
      {
        float_i = i + 1;
        for (; sequence_ptr[i] != ')'; i++)
          ; //iterate to closing bracket
        strncpy(modStringBuffer, &sequence_ptr[float_i], i - float_i);
        modStringBuffer[i - float_i] = '\0';

        // find comma
        int j;
        for (j = 0; (j < i - float_i) && (modStringBuffer[j] != ','); j++)
          //j = # of aas before comma
          ;
        spec_index += j;
      }
      else
      {
        int aa_index = 0;
        while (aa_index < AAcount and refLetters[aa_index] != sequence[i])
          aa_index++;
        if (aa_index >= AAcount)
        {
          WARN_MSG("Warning3: unknown aa! " << sequence[i]);
        }
        else
        {
          spec_index++;
        }
      }
    }
    free(modStringBuffer);
    return spec_index;
  }

  bool AAJumps::getPeptideShift(const string &pep1,
                                const string &pep2,
                                int minAAOverlap,
                                int *putDaShift,
                                float *putPeakShift,
                                int *putResShift) const
  {
    if (pep1 == pep2 && pep1.length() >= minAAOverlap)
    {
      *putDaShift = 0;
      *putResShift = 0;
      return true;
    }

    int bestStart1 = -1, bestStart2 = -1, maxOverlap = -1;

    bool foundOverlap = false;
    int maxStart2 = pep2.length() - minAAOverlap;
    for (int start2 = 0; start2 <= maxStart2; start2++)
    {
      int remain2 = pep2.length() - start2;
      remain2 = min(remain2, (int)pep1.length());

      bool res = true;

      for (int strI = 0; strI < remain2; strI++)
      {
        if (pep1[strI] != pep2[strI + start2])
        {
          res = false;
          break;
        }
      }

      if (res && remain2 > maxOverlap)
      {
        maxOverlap = remain2;
        bestStart1 = 0;
        bestStart2 = start2;
      }
    }

    int maxStart1 = pep1.length() - minAAOverlap;
    for (int start1 = 0; start1 <= maxStart1; start1++)
    {
      int remain1 = pep1.length() - start1;
      remain1 = min(remain1, (int)pep2.length());

      bool res = true;

      for (int strI = 0; strI < remain1; strI++)
      {
        if (pep1[strI + start1] != pep2[strI])
        {
          res = false;
          break;
        }
      }

      if (res && remain1 > maxOverlap)
      {
        maxOverlap = remain1;
        bestStart1 = start1;
        bestStart2 = 0;
      }
    }

    if (maxOverlap >= minAAOverlap)
    {
      if (bestStart1 == 0 && bestStart2 == 0)
      {
        *putDaShift = 0;
        *putResShift = 0;
        *putPeakShift = 0;
      }
      else
      {
        string overlap("");
        int shift = 0;
        float unRoundedShift = 0;
        if (bestStart2 > 0)
        {
          overlap = pep2.substr(0, bestStart2);
          *putResShift = 0 - bestStart2;
        }
        else if (bestStart1 > 0)
        {
          overlap = pep1.substr(0, bestStart1);
          *putResShift = bestStart1;
        }

        if (overlap.length() > 0)
        {
          vector<float> masses;
          getRoundedPRMMasses(overlap, masses);

          shift = floatToInt(masses[masses.size() - 1]);

          getPRMMasses(overlap, masses);

          unRoundedShift = masses[masses.size() - 1];
        }

        if (bestStart2 > 0)
        {
          shift = 0 - shift;
          unRoundedShift = 0.0 - unRoundedShift;
        }

        *putDaShift = shift;
        *putPeakShift = unRoundedShift;
      }

      return true;
    }

    return false;
  }

  string AAJumps::getPeptideSuffix(const string &pep, const int prefix) const
  {
    if (pep.length() == 0)
    {
      return string("");
    }
    vector<float> masses;
    getRoundedPRMMasses(pep, masses);

    if (masses.size() == 0)
    {
      return string("");
    }

    vector<int> massesInt(masses.size());
    for (int p = 0; p < masses.size(); p++)
    {
      massesInt[p] = floatToInt(masses[p]);
    }

    for (int p = 0; p < pep.length(); p++)
    {
      if (massesInt[p] > prefix)
      {
        return pep.substr(p, string::npos);
      }
    }
    return string("");
  }

  string AAJumps::getPeptidePrefix(const string &pep, const int suffix) const
  {
    if (pep.length() == 0)
    {
      return string("");
    }
    vector<float> masses;
    getRoundedPRMMasses(pep, masses);

    if (masses.size() == 0)
    {
      return string("");
    }

    vector<int> massesInt(masses.size());
    for (int p = 0; p < masses.size(); p++)
    {
      massesInt[p] = floatToInt(masses[p]);
    }

    int pepMass = massesInt[masses.size() - 1];

    for (int p = pep.length() - 1; p >= 0; p--)
    {
      if ((pepMass - massesInt[p]) >= suffix)
      {
        return pep.substr(0, p + 1);
      }
    }
    return string("");
  }

  /**
   * Finds the longest amino acid tag contained in a sequence of masses
   * @param masses sequence of cummulative masses
   * @param jumps AAJumps class
   * return longest amino acid found by checking consecutive mass differences
   */
  string getLongestTag(vector<MZRange>& masses, AAJumps& jumps)
  {
    string next_tag;
    string max_tag;

    if (masses.size() == 0)
    {
      return string(max_tag);
    }

    MZRange prev_mass = masses[0];
    MZRange next_mass;
    MZRange diff_mass;
    int AAIdx;

    for (int i = 1; i < masses.size(); i++)
    {
      next_mass = masses[i];
      diff_mass = next_mass - prev_mass;

      AAIdx = -1;
      for (int j = 0; j < jumps.masses.size(); j++)
      {

        if (diff_mass == jumps.masses[j])
        {
          AAIdx = j;
          break;
        }
      }

      if (AAIdx >= 0)
      {

        next_tag.append(1, jumps.aaLetters[AAIdx]);
      }
      else
      {
        if (next_tag.length() > max_tag.length())
        {
          max_tag = next_tag;
        }
        next_tag.clear();
      }

      prev_mass = next_mass;
    }
    if (next_tag.length() > max_tag.length())
    {
      max_tag = next_tag;
    }
    return string(max_tag);
  }
}
