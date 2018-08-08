#include "spectrum_scoring.h"
#include "spectrum.h"
#include "aminoacid.h"
#include <cstring>
#include <cmath>

namespace specnets
{

  MS2ScoringModel::MS2ScoringModel()
  {
    peakTolPPM = 0;
    peakTolDA = 0;
    noiseProb = 0.02;
    probs.clear();
    probs_ci.clear();
    probs_co.clear();
  }

  void MS2ScoringModel::getParentMassIonNames(vector<string> &parentIons) const
  {
    parentIons.clear();
    for (int i = 0; i < probs.size(); i++)
    {
      if (!probs[i].isIF)
      {
        parentIons.push_back(probs[i].name);
      }
    }
  }

  void MS2ScoringModel::getBreakIonsNames(vector<string> &bIons,
                                          vector<string> &yIons) const
  {
    bIons.clear();
    yIons.clear();
    for (int i = 0; i < probs.size(); i++)
      if (probs[i].isBreakIon)
        if (probs[i].isNTerm)
          bIons.push_back(probs[i].name);
        else
          yIons.push_back(probs[i].name);
  }

  bool MS2ScoringModel::ParseParameter(char *parameterLine)
  {
    struct ftIonFragment tmpFrag;
    struct ftCorrelatedIntensity tmpCI;
    struct ftCorrelatedOccurrence tmpCO;
    char *token;

    token = strtok(parameterLine, ";");
    if (strcmp(token, "MTD") == 0)
    { // peak Mass Tolerance in Daltons
      token = strtok(NULL, ";");
      if (token == NULL)
        return false; //we're missing the MTD token
      peakTolDA = atof(token);
      return true;
    }

    if (strcmp(token, "MT") == 0)
    { // peak Mass Tolerance in ppm
      token = strtok(NULL, ";");
      if (token == NULL)
        return false;
      peakTolPPM = atof(token);
      return true;
    }

    if (strcmp(token, "N") == 0)
    { // probability of Noise
      token = strtok(NULL, ";");
      if (token == NULL)
        return false;
      noiseProb = atof(token);
      return true;
    }

    // Ion Fragment: isNTerm (bool) ;
    // ion mass offset (float) ;
    // ion probability (float) ;
    // ion name (string) ;
    // charge (int)

    if (strcmp(token, "IF") == 0)
    {
      //required parameters:
      token = strtok(NULL, ";");
      if (token == NULL)
        return false;
      tmpFrag.isNTerm = (bool)atoi(token);

      token = strtok(NULL, ";");
      if (token == NULL)
        return false;
      tmpFrag.massOffset = (float)atof(token);

      token = strtok(NULL, ";");
      if (token == NULL)
        return false;
      tmpFrag.prob = (float)atof(token);

      token = strtok(NULL, ";");
      if (token == NULL)
        return false;
      tmpFrag.name.assign(token);

      token = strtok(NULL, ";");
      if (token == NULL)
        return false;
      tmpFrag.charge = (int)atoi(token);

      //optional parameters

      // check for isMainIon flag
      token = strtok(NULL, ";");
      if (token != NULL)
        tmpFrag.isMainIon = (int)atoi(token);

      // check for isBreakIon flag
      token = strtok(NULL, ";");
      if (token != NULL)
        tmpFrag.isBreakIon = (int)atoi(token);

      // check for hasLabel flag
      token = strtok(NULL, ";");
      if (token != NULL)
        tmpFrag.hasLabel = (int)atoi(token);

      // check for label alternative representation
      token = strtok(NULL, ";");
      if (token != NULL)
        tmpFrag.label.assign(token);

      tmpFrag.isIF = true;

      tmp_probs.push_back(tmpFrag);
      return true;
    }

    //Parent mass type. I.e., mass offset from total aa mass of peptide. ;
    // isNterm (bool) - Always not applicable, but set to null. ;
    // ion mass offset (float) ;
    // ion probability (float) ;
    // ion name (string) ;
    // charge (int)

    if (strcmp(token, "PM") == 0)
    {
      //required parameters:
      token = strtok(NULL, ";");
      tmpFrag.isNTerm = false; //ignore set parameter, means nothing

      token = strtok(NULL, ";");
      tmpFrag.massOffset = (float)atof(token);

      token = strtok(NULL, ";");
      if (token == NULL)
        return false;
      tmpFrag.prob = (float)atof(token);

      token = strtok(NULL, ";");
      if (token == NULL)
        return false;
      tmpFrag.name.assign(token);

      token = strtok(NULL, ";");
      if (token == NULL)
        return false;
      tmpFrag.charge = (int)atoi(token);

      //optional parameters

      // check for hasLabel flag
      token = strtok(NULL, ";");
      if (token != NULL)
        tmpFrag.hasLabel = (int)atoi(token);

      // check for label alternative representation
      token = strtok(NULL, ";");
      if (token != NULL)
        tmpFrag.label.assign(token);

      tmpFrag.isIF = false;

      tmp_probs.push_back(tmpFrag);
      return true;
    }

    if (strcmp(token, "CI") == 0)
    { // Correlated Ions: 1-based index of the first ion ; 1-based index of the second ion ; probability of observing the second ion given that the first ion is observed
      token = strtok(NULL, ";");
      tmpCI.fragIndex1 = (unsigned int)atoi(token) - 1;
      token = strtok(NULL, ";");
      tmpCI.fragIndex2 = (unsigned int)atoi(token) - 1;
      token = strtok(NULL, ";");
      tmpCI.prob = (float)atof(token);
      probs_ci.push_back(tmpCI);
      return true;
    }

    if (strcmp(token, "CO") == 0)
    { // Correlated iOns: 1-based index of the first ion ; 1-based index of the second ion ; probability of observing the second ion given that the first ion is observed ; probability of observing the second ion given that the first ion is NOT observed
      token = strtok(NULL, ";");
      tmpCO.fragIndex1 = (unsigned int)atoi(token) - 1;
      token = strtok(NULL, ";");
      tmpCO.fragIndex2 = (unsigned int)atoi(token) - 1;
      token = strtok(NULL, ";");
      tmpCO.probIfPresent = (float)atof(token);
      token = strtok(NULL, ";");
      tmpCO.probIfAbsent = (float)atof(token);
      if (tmpCO.fragIndex1 >= tmp_probs.size() or tmpCO.fragIndex2
          >= tmp_probs.size())
      {
        cerr << "ERROR reading parameter CO for Ion Fragment indices "
            << tmpCO.fragIndex1 + 1 << ", " << tmpCO.fragIndex2 + 1
            << " - only have " << tmp_probs.size()
            << " ion fragment (IF) types defined so far!\n";
        return false;
      }
      else
      {
        list<ftIonFragment>::iterator iter = tmp_probs.begin();
        for (unsigned int listPos = 0; listPos < tmpCO.fragIndex1; listPos++)
          iter++;
        (*iter).prob = 0; // Makes sure that only CO will be used to score this mass offset
      }
      probs_co.push_back(tmpCO);
      return true;
    }

    return false;
  }

  bool MS2ScoringModel::LoadModel(const char *filename)
  {
    unsigned int lineIdx;
    char *curLine;
    bool parseOk;
    BufferedLineReader input;

    if (input.Load(filename) < 0)
    {
      cerr << "ERROR loading model from " << filename
          << " (in MS2ScoringModel::LoadModel())\n";
      return false;
    }
    for (lineIdx = 0; lineIdx < input.size(); lineIdx++)
    {
      curLine = input.getline(lineIdx);
      parseOk = true;

      if (strlen(curLine) > 0 and curLine[0] != '#')
      {
        parseOk = ParseParameter(curLine);
      }

      if (!parseOk)
      {
        cerr << "ERROR parsing " << filename << " on line " << lineIdx + 1
            << ": " << curLine << endl;
        return false;
      }
      else
      {
#ifdef DEBUG
        cerr << " --> got |" << curLine << "|\n";
#endif
      }
    }

    probs.resize(tmp_probs.size());

    unsigned int pIdx = 0;
    for (list<ftIonFragment>::iterator iter = tmp_probs.begin(); iter
        != tmp_probs.end(); iter++)
    {
      probs[pIdx++] = (*iter);
    }

    tmp_probs.clear();

    return true;
  }

  bool MS3ScoringModel::ParseParameter(char *parameterLine)
  {
    struct ftCorrelatedOccurrence tmpCO;
    char *token, *lineBuffer;

    lineBuffer = (char *)malloc(strlen(parameterLine) + 1);
    strcpy(lineBuffer, parameterLine);
    if (MS2ScoringModel::ParseParameter(lineBuffer))
    {
      free(lineBuffer);
      return true;
    }
    free(lineBuffer);

    token = strtok(parameterLine, ";");
    if (strcmp(token, "CO23") == 0)
    { // Correlated iOns between ms2/ms3-dependent spectra: 1-based index of the first ion ; 1-based index of the second ion ; probability of observing the second ion given that the first ion is observed ; probability of observing the second ion given that the first ion is NOT observed
      token = strtok(NULL, ";");
      tmpCO.fragIndex1 = (unsigned int)atoi(token) - 1;
      token = strtok(NULL, ";");
      tmpCO.fragIndex2 = (unsigned int)atoi(token) - 1;
      token = strtok(NULL, ";");
      tmpCO.probIfPresent = (float)atof(token);
      token = strtok(NULL, ";");
      tmpCO.probIfAbsent = (float)atof(token);
      if (tmpCO.fragIndex1 >= tmp_probs.size() or tmpCO.fragIndex2
          >= tmp_probs.size())
      {
        cerr << "ERROR reading parameter CO23 for Ion Fragment indices "
            << tmpCO.fragIndex1 + 1 << ", " << tmpCO.fragIndex2 + 1
            << " - only have " << tmp_probs.size()
            << " ion fragment (IF) types defined so far!\n";
        return false;
      }
      else
      {
        list<ftIonFragment>::iterator iter = tmp_probs.begin();
        for (unsigned int listPos = 0; listPos < tmpCO.fragIndex1; listPos++)
          iter++;
        (*iter).prob = 0; // Makes sure that only CO23 will be used to score this mass offset
      }
      probs_co23.push_back(tmpCO);
      return true;
    }

    return false;
  }

  void ScoreSpectrum(Spectrum &in,
                     MS2ScoringModel &model,
                     bool removeNegativeScores,
                     Spectrum *out)
  {
    vector<float> ifIntensity; // Keeps the ion current for each mass offset
    ifIntensity.resize(model.probs.size());

    // If input/output spectra are the same then use temporary to hold the scored spectrum
    bool sameInOut = (out == NULL or &in == out);
    if (sameInOut)
      out = new Spectrum;
    out->copyNP(in);

    // Find all masses that need to be scored (every peak interpreted as b or y-ion)
    Spectrum masses = in;
    masses.makeSymmetric(2, 0); // Create closer complementary masses. Corresponds to guaranteeing
    //   that every peak is interpreted as b _and_ as y.
    //	if(model.peakTolDA>0) masses.makeSymmetric(2,model.peakTolDA);
    //	else masses.makeSymmetric(2,model.peakTolPPM*1500/1000000);  // Da tolerance at m/z=1500
    out->resize(masses.size());

    float massToScore, curScore, peakTol, ionMass;
    unsigned int peakIdx, ionIdx, matchesIdx, outIdx = 0;
    vector<int> matches;
    for (peakIdx = 0; peakIdx < masses.size(); peakIdx++)
    {
      if (masses[peakIdx][0] <= 1.0072763)
        continue;
      massToScore = masses[peakIdx][0] - 1.0072763; // Assume each peak as a possible b-ion (or y-ion, by symmetry in masses)
      //cerr<<"Current massToScore = "<<massToScore<<"\n";
      if (model.peakTolDA > 0)
        peakTol = model.peakTolDA;
      else
        peakTol = model.peakTolPPM * massToScore / 1000000;

      curScore = 0;
      ionIdx = 0;
      vector<ftIonFragment>::iterator curIonFrag;
      for (curIonFrag = model.probs.begin(); curIonFrag != model.probs.end(); curIonFrag++, ionIdx++)
      {
        // Assuming MH+ parent mass - should be the case for MS2 spectra
        if (curIonFrag->isNTerm)
          ionMass = massToScore + curIonFrag->massOffset;
        else
          ionMass = in.parentMass - AAJumps::massMH - massToScore
              + curIonFrag->massOffset;
        in.findMatches(ionMass, peakTol, matches);

        ifIntensity[ionIdx] = 0;
        for (matchesIdx = 0; matchesIdx < matches.size(); matchesIdx++)
          ifIntensity[ionIdx] += in[matches[matchesIdx]][1];

        if (curIonFrag->prob > 0)
        {
          if (ifIntensity[ionIdx] > 0)
          {
            //cerr<<" -- IF("<<ionMass<<"): intensity = "<<ifIntensity[ionIdx]<<", score change = "<<log(curIonFrag->prob/model.noiseProb)<<" (found peak)\n";
            curScore += log(curIonFrag->prob / model.noiseProb);
          }
          else
          {
            //cerr<<" -- IF("<<ionMass<<"): intensity = "<<ifIntensity[ionIdx]<<",score change = "<<log((1-curIonFrag->prob)/(1-model.noiseProb))<<" (peak not found)\n";
            curScore += log((1 - curIonFrag->prob) / (1 - model.noiseProb));
          }
        }
      }

      // Score correlated intensities
      list<ftCorrelatedIntensity>::iterator curCI;
      for (curCI = model.probs_ci.begin(); curCI != model.probs_ci.end(); curCI++)
        if (ifIntensity[curCI->fragIndex1] > 0
            and ifIntensity[curCI->fragIndex2] > 0)
        {
          if (ifIntensity[curCI->fragIndex1] >= ifIntensity[curCI->fragIndex2])
            curScore += log(curCI->prob / 0.5);
          else
            curScore += log((1 - curCI->prob) / 0.5);
        }

      // Score correlated occurrences
      list<ftCorrelatedOccurrence>::iterator curCO;
      for (curCO = model.probs_co.begin(); curCO != model.probs_co.end(); curCO++)
        if (ifIntensity[curCO->fragIndex2] > 0)
        {
          if (ifIntensity[curCI->fragIndex1] > 0)
            curScore += log(curCO->probIfPresent / model.noiseProb);
          else
            curScore += log((1 - curCO->probIfPresent) / (1 - model.noiseProb));
        }
        else
        {
          if (ifIntensity[curCI->fragIndex1] > 0)
            curScore += log(curCO->probIfAbsent / model.noiseProb);
          else
            curScore += log((1 - curCO->probIfAbsent) / (1 - model.noiseProb));
        }

      (*out)[outIdx++].set(massToScore, curScore);
    }
    out->resize(outIdx);

    if (removeNegativeScores)
    {
      for (outIdx = 0, peakIdx = 0; peakIdx < out->size(); peakIdx++)
        if ((*out)[peakIdx][1] > 0)
          (*out)[outIdx++] = (*out)[peakIdx];
      out->resize(outIdx);
    }

    if (sameInOut)
    {
      in = (*out);
      delete out;
    }
  }

  void ScoreSpectrumMS3(Spectrum &ms3in,
                        Spectrum &ms2in,
                        MS3ScoringModel &model,
                        float ms3shift,
                        bool removeNegativeScores,
                        Spectrum *ms3out)
  {
    vector<float> ifIntensity; // Keeps the ion current for each mass offset
    ifIntensity.resize(model.probs.size());

    // If input/output spectra are the same then use temporary to hold the scored spectrum
    bool sameInOut = (ms3out == NULL or &ms3in == ms3out);
    if (sameInOut)
      ms3out = new Spectrum;
    ms3out->copyNP(ms3in);

    // Find all masses that need to be scored (every peak interpreted as b or y-ion)
    Spectrum masses = ms3in;
    masses.makeSymmetric(2, 0); // Create closer complementary masses. Corresponds to guaranteeing
    //   that every peak is interpreted as b _and_ as y.
    //	if(model.peakTolDA>0) masses.makeSymmetric(2,model.peakTolDA);
    //	else masses.makeSymmetric(2,model.peakTolPPM*1500/1000000);  // Da tolerance at m/z=1500
    ms3out->resize(masses.size());

    float massToScore, curScore, peakTol, ionMass;
    unsigned int peakIdx, ionIdx, matchesIdx, outIdx = 0;
    vector<int> matches;
    for (peakIdx = 0; peakIdx < masses.size(); peakIdx++)
    {
      massToScore = masses[peakIdx][0] - 1.0072763; // Assume each peak as a possible b-ion (or y-ion, by symmetry in masses)
      if (massToScore <= 0 or massToScore >= ms3in.parentMass)
        continue;
      if (model.peakTolDA > 0)
        peakTol = model.peakTolDA;
      else
        peakTol = model.peakTolPPM * massToScore / 1000000;

      curScore = 0;
      ionIdx = 0;
      for (ionIdx = 0; ionIdx < model.probs.size(); ionIdx++)
      {
        // Assuming MH+ parent mass - should be the case for MS2 spectra
        if (model.probs[ionIdx].isNTerm)
          ionMass = massToScore + model.probs[ionIdx].massOffset;
        else
          ionMass = ms3in.parentMass - AAJumps::massMH - massToScore
              + model.probs[ionIdx].massOffset;
        ms3in.findMatches(ionMass, peakTol, matches);

        ifIntensity[ionIdx] = 0;
        for (matchesIdx = 0; matchesIdx < matches.size(); matchesIdx++)
          ifIntensity[ionIdx] += ms3in[matches[matchesIdx]][1];

        if (model.probs[ionIdx].prob > 0)
        {
          if (ifIntensity[ionIdx] > 0)
            curScore += log(model.probs[ionIdx].prob / model.noiseProb);
          else
            curScore += log((1 - model.probs[ionIdx].prob) / (1
                - model.noiseProb));
        }
      }

      // Score correlated intensities
      list<ftCorrelatedIntensity>::iterator curCI;
      for (curCI = model.probs_ci.begin(); curCI != model.probs_ci.end(); curCI++)
        if (ifIntensity[curCI->fragIndex1] > 0
            and ifIntensity[curCI->fragIndex2] > 0)
        {
          if (ifIntensity[curCI->fragIndex1] >= ifIntensity[curCI->fragIndex2])
            curScore += log(curCI->prob / 0.5);
          else
            curScore += log((1 - curCI->prob) / 0.5);
        }

      // Score correlated occurrences
      list<ftCorrelatedOccurrence>::iterator curCO;
      for (curCO = model.probs_co.begin(); curCO != model.probs_co.end(); curCO++)
        if (ifIntensity[curCO->fragIndex2] > 0)
        {
          if (ifIntensity[curCI->fragIndex1] > 0)
            curScore += log(curCO->probIfPresent / model.noiseProb);
          else
            curScore += log((1 - curCO->probIfPresent) / (1 - model.noiseProb));
        }
        else
        {
          if (ifIntensity[curCI->fragIndex1] > 0)
            curScore += log(curCO->probIfAbsent / model.noiseProb);
          else
            curScore += log((1 - curCO->probIfAbsent) / (1 - model.noiseProb));
        }

      // Score correlated occurrences in MS2/MS3
      for (curCO = model.probs_co23.begin(); curCO != model.probs_co23.end(); curCO++)
      {
        if (model.probs[curCO->fragIndex2].isNTerm)
          ionMass = (massToScore + ms3shift)
              + model.probs[curCO->fragIndex2].massOffset;
        else
          ionMass = ms2in.parentMass - AAJumps::massMH - (massToScore
              + ms3shift) + model.probs[curCO->fragIndex2].massOffset;
        ms2in.findMatches(ionMass, peakTol, matches);

        if (matches.size() > 0)
        {
          if (ifIntensity[curCO->fragIndex1] > 0)
            curScore += log(curCO->probIfPresent / model.noiseProb);
          else
            curScore += log((1 - curCO->probIfPresent) / (1 - model.noiseProb));
        }
        else
        {
          if (ifIntensity[curCO->fragIndex1] > 0)
            curScore += log(curCO->probIfAbsent / model.noiseProb);
          else
            curScore += log((1 - curCO->probIfAbsent) / (1 - model.noiseProb));
        }
      }

      (*ms3out)[outIdx++].set(massToScore, curScore);
    }
    ms3out->resize(outIdx);

    if (removeNegativeScores)
    {
      for (outIdx = 0, peakIdx = 0; peakIdx < ms3out->size(); peakIdx++)
        if ((*ms3out)[peakIdx][1] > 0)
          (*ms3out)[outIdx++] = (*ms3out)[peakIdx];
      ms3out->resize(outIdx);
    }

    if (sameInOut)
    {
      ms3in = (*ms3out);
      delete ms3out;
    }
  }

}
