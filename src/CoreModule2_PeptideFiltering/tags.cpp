#include "tags.h"
#include "aminoacid.h"
#include<cmath>

#define DEBUG_TAGS 0
#define DEBUG_TAGS2 0
#define DEBUG_EQUIV_TAGS 0

namespace specnets
{
  string getTagString(AAJumps & jumps, Tag & tagIn)
  { 
    string seq;
    for (int c = 0; c < tagIn.sequence.size(); c++) {
      seq += jumps.aaLetters[tagIn.sequence[c]];
    }
    return seq;
  }

  bool cmp_score(Tag & t1, Tag & t2)
  {
    return t1.score > t2.score;
  } // Used to sort tags by descending tag score

  bool cmp_seq(Tag & t1, Tag & t2)
  {
    return t1.strSequence > t2.strSequence;
  }

  bool equal_seq(Tag & t1, Tag & t2)
  {
    return t1.strSequence == t2.strSequence;
  }

  // FindMassNeighbors - indices[i] contains the indices of the peaks in spec
  //    such that abs( mass - jumps[i] - spec[any peak in indices[i]][0] ) <= peakTol
  //    idxSupremum - index of a peak in spec with a mass higher than mass-min(jumps) (optional)
  void FindMassNeighbors(float mass, Spectrum &spec, AAJumps &jumps, vector<
      list<short> > &indices, float peakTol, int idxSupremum = -1)
  {
    if (idxSupremum < 0 or idxSupremum >= (int) spec.size())
      idxSupremum = (int) spec.size() - 1;
    if (idxSupremum < 0)
      return;

    unsigned int jumpIdx;
    int peakIdx = idxSupremum;
    indices.resize(jumps.size());
    for (jumpIdx = 0; jumpIdx < jumps.size(); jumpIdx++)
      indices[jumpIdx].clear();
    if (spec.size() == 0)
      return;
    float curMass;

    for (jumpIdx = 0; jumpIdx < jumps.size(); jumpIdx++) {
      curMass = mass - jumps[jumpIdx];
      if (curMass < -peakTol)
        break;
      if (peakIdx < 0)
        peakIdx = 0;
      while (peakIdx < (int) spec.size() and spec[peakIdx][0] <= curMass
          + peakTol + 0.0001)
        peakIdx++;
      if (peakIdx >= (int) spec.size())
        peakIdx = (int) spec.size() - 1;
      while (peakIdx >= 0 and spec[peakIdx][0] > curMass + peakTol + 0.0001)
        peakIdx--;
      while (peakIdx >= 0 and spec[peakIdx][0] > curMass - peakTol - 0.0001)
        indices[jumpIdx].push_back(peakIdx--);
    }
  }

  float findPeakFromMass(Spectrum &spec, float mass)
  {
    //DEBUG_VAR(mass);
    int peakIndex = -1;
    for (int i = 0; i < spec.size(); i++) {
      //DEBUG_VAR(spec[i][0]);
      if (abs(spec[i][0] - mass) < 0.3) {
        peakIndex = i;
        break;
      }
    }
    //DEBUG_VAR(peakIndex);
    return peakIndex;
  }

  void createEquivalentTags(AAJumps & jumps, 
			    Spectrum &spec,
                            Tag & tagIn, 
                            list<Tag> & tagsOut, 
                            map<char, vector<vector<char> > > & equivalences)
  {
    //DEBUG_VAR(tagIn.strSequence);
    // Loop over all characters (actually indexes) in the tag sequence
    for (int s = 0; s < tagIn.sequence.size(); s++) {
      map<char, vector<vector<char> > >::iterator itrFound = equivalences.find(tagIn.sequence[s]);
      //DEBUG_VAR(jumps.aaLetters[tagIn.sequence[s]]);
      if (itrFound != equivalences.end()) {
        // Loop over all the equivalent combinations
        for (int i = 0; i < itrFound->second.size(); i++) {
          Tag newTag(tagIn);
          newTag.sequence.clear();
          newTag.score = 0.0;  // We set the score to 0.0 so we can check against it later
          float additonalMass = 0.0;
          for (int s2 = 0; s2 < s; s2++) {
            if (DEBUG_EQUIV_TAGS) DEBUG_MSG((int)tagIn.sequence[s2] << "  " << jumps.aaLetters[tagIn.sequence[s2]]);
            newTag.sequence.push_back(tagIn.sequence[s2]);
            additonalMass += jumps.aaLetters[tagIn.sequence[s2]];
          }
          for (int s3 = 0; s3 < itrFound->second[i].size(); s3++) {
            if (DEBUG_EQUIV_TAGS) DEBUG_MSG((int)itrFound->second[i][s3] << "  " << jumps.aaLetters[itrFound->second[i][s3]]);
            newTag.sequence.push_back(itrFound->second[i][s3]);
            if (s3 < itrFound->second[i].size() - 1) {
              additonalMass += jumps.aaLetters[tagIn.sequence[s3]];
              if (DEBUG_EQUIV_TAGS) DEBUG_VAR(tagIn.flankingPrefix);
              if (DEBUG_EQUIV_TAGS) DEBUG_VAR(additonalMass);
              float mass = tagIn.flankingPrefix + additonalMass;
              if (DEBUG_EQUIV_TAGS) DEBUG_VAR(mass);
              int peakIndex = findPeakFromMass(spec, mass);
              if (DEBUG_EQUIV_TAGS) DEBUG_VAR(peakIndex);
              if (peakIndex > 0) {
                newTag.score += spec[peakIndex][1];
              }
            }
          }
          for (int s2 = s + 1; s2 < tagIn.sequence.size(); s2++) {
            if (DEBUG_EQUIV_TAGS) DEBUG_MSG((int)tagIn.sequence[s2] << "  " << jumps.aaLetters[tagIn.sequence[s2]]);
            newTag.sequence.push_back(tagIn.sequence[s2]);
          }
          // We want the new score to be low unless it would be better (so it will get tossed if it isn't better)
          // If we don't do this the score would be equal to the original and we don't want that
          if (newTag.score > 0.0) {
            newTag.score += tagIn.score;
            newTag.strSequence = getTagString(jumps, newTag);
            if (DEBUG_EQUIV_TAGS) DEBUG_VAR(newTag.strSequence);
            if (DEBUG_EQUIV_TAGS) DEBUG_MSG("PUSH");
            tagsOut.push_back(newTag);
          }
        }
      }
    }
    return;
  }

  void getTopNTags(list<Tag> &tags, unsigned int maxNumTags) {
    // Chop list down to the top N tags as requested
    tags.sort(cmp_score);
    list<Tag>::iterator iter = tags.begin();
    list<Tag>::iterator iter_end = tags.end();
    float lastScore = 0.0;
    if (maxNumTags > 0 and tags.size() > maxNumTags) {
      // Keep top N tags
      for (unsigned int count = 0; count < maxNumTags && iter->score != 0; count++) {
        lastScore = iter->score;
        iter++;
      }
      // And keep all with the same score
      while (iter->score == lastScore) {
        iter++;
      }
      tags.erase(iter, tags.end());
    }
    return;
  }


  void splitTagsToLength(AAJumps & jumps, list<Tag> &tagsIn, list<Tag> &tagsOut, unsigned int tagLen)
  {
    // Chop tags down to tagLen
    list<Tag>::iterator iter = tagsIn.begin();
    list<Tag>::iterator iter_end = tagsIn.end();
    for ( ; iter != iter_end; iter++) {
      int addLength = iter->sequence.size() - tagLen;
      if (addLength != 0) {
        for (int i = 0; i <= addLength; i++) {
          Tag newTag(*iter);
          newTag.sequence.clear();
          for (int j = 0; j < tagLen; j++) {
            newTag.sequence.push_back(iter->sequence[j + i]);
          }
          newTag.strSequence = getTagString(jumps, newTag);
          tagsOut.push_back(newTag);
        }
      } else {
        tagsOut.push_back(*iter);
      }
    }
    return;
  }

  unsigned int ExtractTagsAllJumps(Spectrum &spec,
                           list<Tag> &tags,
                           float peakTol,
                           unsigned int tagLen,
                           unsigned int maxNumJumps,
                           unsigned int maxNumTags,
                           float        peakSkipPenalty)
  {
    AAJumps jumps(1);
    tags.clear();
    list<Tag> tagsOneJumpSize;
    for (int j = 0; j <= maxNumJumps; j++) {
     ExtractTags(spec,
                 tagsOneJumpSize,
                 peakTol,
                 tagLen,
                 j,
                 maxNumTags,
                 peakSkipPenalty);
      list<Tag>::iterator iter = tagsOneJumpSize.begin();
      list<Tag>::iterator iter_end = tagsOneJumpSize.end();
      for ( ; iter != iter_end; iter++) {
        iter->strSequence = getTagString(jumps, *iter);
        tags.push_back(*iter);
      }
    }

    map<char, vector<vector<char> > > equivalences;
    jumps.getEquivalentIndices(equivalences);

    list<Tag> tagsDone;
    list<Tag> tagsEquiv;
    list<Tag> & tagsIn = tags;
    list<Tag> & tagsOut = tagsEquiv;
    int newTags = 0;

    list<Tag>::iterator iter = tagsIn.begin();
    list<Tag>::iterator iter_end = tagsIn.end();

    do {
      //DEBUG_VAR(tagsIn.size());
      iter = tagsIn.begin();
      // Create all the equivalent tags for the input set of tags
      for ( ; iter != iter_end; iter++) {
        createEquivalentTags(jumps, spec, *iter, tagsOut, equivalences);
      }
      newTags = tagsOut.size();
      //DEBUG_VAR(tagsOut.size());

      // Take the input tags and add them to the "done" list
      tagsDone.insert(tagsDone.begin(), tagsIn.begin(), tagsIn.end());
      //DEBUG_VAR(tagsDone.size());

      // Move the output list to the input
      tagsIn.clear();
      list<Tag> & tagsTemp = tagsIn;
      tagsIn = tagsOut;
      tagsOut = tagsTemp;
      tagsOut.clear();
    } while (newTags != 0);

    tagsDone.insert(tagsDone.begin(), tagsOut.begin(), tagsOut.end());
    tags.clear();
    splitTagsToLength(jumps, tagsDone, tags, tagLen);

    // Make sure all tags have their string sequences (for comparisons)
    list<Tag>::iterator itrTag = tags.begin();
    list<Tag>::iterator itrTagEnd = tags.end();
    for ( ; itrTagEnd != itrTag; itrTag++) {
      if (itrTag->strSequence.empty()) itrTag->strSequence = getTagString(jumps, *itrTag);
    }
    
    tags.sort(cmp_seq);
    tags.unique(equal_seq);

    if (DEBUG_TAGS) {
      DEBUG_MSG("All Tags");
      for (iter = tags.begin(), iter_end = tags.end(); iter != iter_end; iter++) {
        DEBUG_MSG(iter->strSequence << "  " << iter->score);
      }
    }


    // Chop list down to the top N tags as requested
    getTopNTags(tags, maxNumTags);

    if (DEBUG_TAGS) {
      DEBUG_MSG("TopN Tags");
      for (iter = tags.begin(), iter_end = tags.end(); iter != iter_end; iter++) {
        DEBUG_MSG(iter->strSequence << "  " << iter->score);
      }
    }

    return (tags.size());
  }

  unsigned int ExtractTags(Spectrum &spec,
                           list<Tag> &tags,
                           float peakTol,
                           unsigned int tagLen,
                           unsigned int maxNumJumps,
                           unsigned int maxNumTags,
                           float        peakSkipPenalty)
  {
    if (DEBUG_TAGS) {
      DEBUG_VAR(peakTol);
      DEBUG_VAR(tagLen);
      DEBUG_VAR(maxNumJumps);
      DEBUG_VAR(maxNumTags);
      DEBUG_VAR(peakSkipPenalty);
    }
    // Position [i,j,k] contains all the tags of length k
    //   with j di-peptide jumps and ending on the i-th spectrum peak
    // dimension 1 - Index of the spectrum peak where the tags end
    // dimension 2 - Number of used di-peptide jumps (usually 0-2)
    // dimension 3 - Tag length
    AAJumps jumps(1);
    unsigned int diIdx, peakIdx, lenIdx, aaIdx, numPeaks = spec.size(),
        numJumps = jumps.size();
    vector < vector<vector<list<Tag> > > > partialTags; // Table containing all partial tags
    partialTags.resize(numPeaks);
    Tag curTag(tagLen);
    tags.clear();

    vector<list<short> > aaNeighs(numJumps), // Peaks with mass one aa to the left of current peak
        diNeighs(numJumps); // Peaks with mass two aa to the left of current peak

    if (DEBUG_TAGS) {
        for (int i = 0; i < spec.size(); i++) {
          DEBUG_MSG(i << "  " << spec[i][0] << "  " << spec[i][1]);
        }
    }

    // Initializations
    for (peakIdx = 0; peakIdx < numPeaks; peakIdx++) {
      partialTags[peakIdx].resize(maxNumJumps + 1);
      for (diIdx = 0; diIdx <= maxNumJumps; diIdx++) {
        partialTags[peakIdx][diIdx].resize(tagLen);
        for (lenIdx = 0; lenIdx < tagLen; lenIdx++)
          partialTags[peakIdx][diIdx][lenIdx].clear();
      }
    }

    // Generate tags
    if (DEBUG_TAGS) DEBUG_VAR(numPeaks);
    for (peakIdx = 0; peakIdx < numPeaks; peakIdx++) {
      FindMassNeighbors(spec[peakIdx][0],
                        spec,
                        jumps,
                        aaNeighs,
                        peakTol,
                        (int) peakIdx);
      unsigned int lastAAneigh = peakIdx; // Index of the lowest-mass detected aa neighbor of peakIdx (used for diNeighs searches)

      if (DEBUG_TAGS2) DEBUG_VAR(numJumps);
      for (aaIdx = 0; aaIdx < numJumps; aaIdx++) {

        // No neighbors 1-aa away, look for neighbors 2 aa masses away
        if (DEBUG_TAGS2) DEBUG_VAR(aaIdx);
        if (DEBUG_TAGS2) DEBUG_VAR(aaNeighs[aaIdx].size());
        if (DEBUG_TAGS2) DEBUG_VAR(maxNumJumps);
        if (aaNeighs[aaIdx].size() == 0 and maxNumJumps > 0) {
          FindMassNeighbors(spec[peakIdx][0] - jumps[aaIdx],
                            spec,
                            jumps,
                            diNeighs,
                            peakTol,
                            (int) lastAAneigh);

          for (unsigned int aaDiIdx = 0; aaDiIdx < numJumps; aaDiIdx++) {
            if (DEBUG_TAGS) DEBUG_VAR(diNeighs[aaDiIdx].size());
            for (list<short>::iterator neigh = diNeighs[aaDiIdx].begin(); neigh
                != diNeighs[aaDiIdx].end(); neigh++) {

              // Initialize tags of length 2 with no middle peak
              curTag.sequence[0] = (char) aaDiIdx;
              curTag.sequence[1] = (char) aaIdx;
              float skipPenalty = 0.0;
              if (DEBUG_TAGS) DEBUG_MSG(*neigh << "  " << peakIdx);
              for (int iSkipped = *neigh + 1; iSkipped < peakIdx; iSkipped++) {
                skipPenalty += spec[iSkipped][1] * peakSkipPenalty;
              }
              if (DEBUG_TAGS) DEBUG_VAR(skipPenalty);
              curTag.score = spec[*neigh][1] + spec[peakIdx][1] - skipPenalty;
              if (DEBUG_TAGS) DEBUG_VAR(curTag.score);
              curTag.flankingPrefix = spec[*neigh][0];
              curTag.flankingSuffix = spec.parentMass - AAJumps::massMH
                  - spec[peakIdx][0];
              if (DEBUG_TAGS) DEBUG_MSG(peakIdx << "  " << jumps.aaLetters[curTag.sequence[0]] << "  " << jumps.aaLetters[curTag.sequence[1]]);
              partialTags[peakIdx][1][1].push_back(curTag);
              if (tagLen == 2) {
                tags.push_back(curTag);
              }

              // Extend to tags of length 3+
              for (diIdx = 1; diIdx <= maxNumJumps; diIdx++) {
                for (lenIdx = 2; lenIdx < tagLen; lenIdx++) {
                  if (partialTags[*neigh][diIdx - 1][lenIdx - 2].size() == 0)
                    continue;
                  for (list<Tag>::iterator tagIter = partialTags[*neigh][diIdx
                      - 1][lenIdx - 2].begin(); tagIter
                      != partialTags[*neigh][diIdx - 1][lenIdx - 2].end(); tagIter++) {
                    curTag = *tagIter;
                    curTag.sequence[lenIdx - 1] = (char) aaDiIdx;
                    curTag.sequence[lenIdx] = (char) aaIdx;
                    if (DEBUG_TAGS) DEBUG_MSG(*neigh << "  " << peakIdx);
                    for (int iSkipped = *neigh + 1; iSkipped < peakIdx; iSkipped++) {
                      skipPenalty += spec[iSkipped][1] * peakSkipPenalty;
                    }
                    if (DEBUG_TAGS) DEBUG_VAR(skipPenalty);
                    curTag.score += spec[peakIdx][1] - skipPenalty;
                    if (DEBUG_TAGS) DEBUG_VAR(curTag.score);
                    curTag.flankingSuffix = spec.parentMass - AAJumps::massMH
                        - spec[peakIdx][0];
                    partialTags[peakIdx][diIdx][lenIdx].push_back(curTag);
                    if (DEBUG_TAGS) DEBUG_MSG(peakIdx << "  " << jumps.aaLetters[curTag.sequence[lenIdx - 1]] << "  " << jumps.aaLetters[curTag.sequence[lenIdx]]);
                    if (lenIdx == tagLen - 1)
                      tags.push_back(curTag);
                  }
                }
              }
            }
          }
          continue;
        }

        // Handle neighbors 1-aa away
        for (list<short>::iterator neigh = aaNeighs[aaIdx].begin(); neigh
            != aaNeighs[aaIdx].end(); neigh++) {
          lastAAneigh = *neigh;

          // Initialize tags of length 1
          curTag.sequence[0] = (char) aaIdx;
          curTag.score = spec[*neigh][1] + spec[peakIdx][1];
          if (DEBUG_TAGS) DEBUG_VAR(curTag.score);
          curTag.flankingPrefix = spec[*neigh][0];
          curTag.flankingSuffix = spec.parentMass - AAJumps::massMH
              - spec[peakIdx][0];
          if (DEBUG_TAGS) DEBUG_MSG(peakIdx << "  " << jumps.aaLetters[curTag.sequence[0]]);
          partialTags[peakIdx][0][0].push_back(curTag);
          if (tagLen == 1)
            tags.push_back(curTag);

          // Extend to tags of length 2+
          for (diIdx = 0; diIdx <= maxNumJumps; diIdx++) {
            for (lenIdx = 1; lenIdx < tagLen; lenIdx++) {
              if (partialTags[*neigh][diIdx][lenIdx - 1].size() == 0)
                continue;
              for (list<Tag>::iterator tagIter =
                  partialTags[*neigh][diIdx][lenIdx - 1].begin(); tagIter
                  != partialTags[*neigh][diIdx][lenIdx - 1].end(); tagIter++) {
                curTag = *tagIter;
                curTag.sequence[lenIdx] = (char) aaIdx;
                curTag.score += spec[peakIdx][1];
                if (DEBUG_TAGS) DEBUG_VAR(curTag.score);
                curTag.flankingSuffix = spec.parentMass - AAJumps::massMH
                    - spec[peakIdx][0];
                if (DEBUG_TAGS) DEBUG_MSG(peakIdx << "  " << jumps.aaLetters[curTag.sequence[lenIdx]]);
                partialTags[peakIdx][diIdx][lenIdx].push_back(curTag);
                if (lenIdx == tagLen - 1)
                  tags.push_back(curTag);
              }
            }
          }
        }
      }
    }

    // Chop list down to the top N tags as requested
    getTopNTags(tags, maxNumTags);

    return (tags.size());
  }

  unsigned int ExtractTags(char *sequence,
                           vector<Tag> &tags,
                           unsigned int tagLen)
  {
    tags.resize(0);
    if (!sequence or strlen(sequence) < tagLen)
      return 0;
    unsigned int maxStart = strlen(sequence) - tagLen;
    tags.resize(maxStart + 1);
    for (unsigned int tagStart = 0; tagStart <= maxStart; tagStart++) {
      tags[tagStart].sequence.resize(tagLen);
      for (unsigned int pivot = 0; pivot < tagLen; pivot++)
        tags[tagStart].sequence[pivot] = sequence[tagStart + pivot];
    }
    return maxStart + 1;
  }

  /*
   * MatchTagsToSpectra - Matches a set of amino acid sequence tags to each of a set of spectra
   *   and returns the indices and match-score of the matched tags.
   *
   * specs     - set of spectra
   * tags      - set of tags
   * peakTol   - peak mass tolerance when matching tags to the spectrum
   * maxCharge - maximum fragment charge to consider when matching tags to the spectrum
   * missingPeakPenalty - score penalty for missing a peak in the tag
   * noisePeakPenalty   - score penalty factor (times score of noise peaks) applied to unmatched
   *                       spectrum peaks between first and last tag masses
   * maxNumMissedPeaks  - maximum number of missed tag peaks
   * matchedTagsIdx     - indices of matched tags per spectrum
   * matchedTagsScore   - scores of matched tags per spectrum
   * matchedTagsPos     - index of the leftmost spectrum peak where the tag was matched
   */
  void MatchTagsToSpectra(SpecSet &specs,
                          list<Tag> &tags,
                          float peakTol,
                          float maxCharge,
                          float missingPeakPenalty,
                          float noisePeakPenalty,
                          unsigned int maxNumMissedPeaks,
                          vector<list<unsigned int> > &matchedTagsIdx,
                          vector<list<float> > &matchedTagsScores,
                          vector<list<unsigned int> > &matchedTagsPos)
  {
    matchedTagsIdx.resize(specs.size());
    matchedTagsScores.resize(specs.size());
    matchedTagsPos.resize(specs.size());
    for (unsigned int i = 0; i < specs.size(); i++) {
      matchedTagsIdx[i].clear();
      matchedTagsScores[i].clear();
      matchedTagsPos[i].clear();
    }
    if (specs.size() == 0 or tags.size() == 0)
      return;
    vector < vector<float> > tagMasses(tags.size()); // Mass-vector representation of all tags
    vector<float> curTagMasses;
    unsigned int specIdx, peakIdx, tagsIdx, tagMassIdx;
    float tagMatchScore, cumTagMass;

    //	for(tagsIdx=0; tagsIdx<tags.size(); tagsIdx++) getMasses(tags[tagsIdx].sequence, tagMasses[tagsIdx]);
    list<Tag>::iterator tagsIter = tags.begin();
    for (tagsIdx = 0; tagsIter != tags.end(); tagsIter++, tagsIdx++)
      getMasses(tagsIter->sequence, tagMasses[tagsIdx]);

    for (specIdx = 0; specIdx < specs.size(); specIdx++) {
      for (float charge = 1.0; charge < maxCharge + 0.0001; charge
          = (float) round(charge + 1.0)) {
        for (tagsIdx = 0; tagsIdx < tagMasses.size(); tagsIdx++) {
          //			for(tagsIdx=6; tagsIdx<7; tagsIdx++) {
          curTagMasses.resize(tagMasses[tagsIdx].size() + 1);
          curTagMasses[0] = 0;
          for (tagMassIdx = 0; tagMassIdx < tagMasses[tagsIdx].size(); tagMassIdx++)
            curTagMasses[tagMassIdx + 1] = tagMasses[tagsIdx][tagMassIdx]
                / charge;

          for (peakIdx = 0; peakIdx < specs[specIdx].size(); peakIdx++) {
            vector<int> matchesIdx;
            unsigned int pivot, lastMatchedPeak = peakIdx, missedPeaks = 0;
            cumTagMass = specs[specIdx][peakIdx][0];
            tagMatchScore = 0;
            for (tagMassIdx = 0; tagMassIdx < curTagMasses.size()
                and missedPeaks <= maxNumMissedPeaks; tagMassIdx++) {
              cumTagMass += curTagMasses[tagMassIdx];

              // Look for a matching tag peak
              if (specs[specIdx].findMatches(cumTagMass,
                                             peakTol,
                                             matchesIdx,
                                             lastMatchedPeak)) {
                for (pivot = 0; pivot < matchesIdx.size(); pivot++)
                  tagMatchScore += specs[specIdx][matchesIdx[pivot]][1];
                if (tagMassIdx == 0) {
                  lastMatchedPeak = matchesIdx[matchesIdx.size() - 1];
                  continue;
                }
              }
              else {
                tagMatchScore += missingPeakPenalty;
                missedPeaks++;
                //cerr<<"["<<specIdx<<","<<peakIdx<<","<<charge<<"]: missed cumTagMass = "<<cumTagMass<<", tagMassIdx = "<<tagMassIdx<<endl; cerr.flush();
              }

              // Penalize for noise peaks
              for (pivot = lastMatchedPeak + 1; pivot < specs[specIdx].size()
                  and specs[specIdx][pivot][0] < cumTagMass - peakTol - 0.0001; pivot++)
                tagMatchScore += (noisePeakPenalty * specs[specIdx][pivot][1]);
              if (matchesIdx.size())
                lastMatchedPeak = matchesIdx[matchesIdx.size() - 1];
              else {
                for (lastMatchedPeak = pivot - 1; lastMatchedPeak
                    < specs[specIdx].size()
                    and specs[specIdx][lastMatchedPeak][0] < cumTagMass
                        + peakTol + 0.0001; lastMatchedPeak++)
                  ;
                lastMatchedPeak--;
              }
            }

            if (missedPeaks <= maxNumMissedPeaks) {
              //cerr<<"["<<specIdx<<","<<peakIdx<<","<<charge<<"]: missedPeaks = "<<missedPeaks<<", tagsIdx = "<<tagsIdx<<", tagMatchScore = "<<tagMatchScore<<"\n"; cerr.flush();
              matchedTagsIdx[specIdx].push_back(tagsIdx);
              matchedTagsScores[specIdx].push_back(tagMatchScore);
              matchedTagsPos[specIdx].push_back(peakIdx);
            }
            
          } // for (peakIdx = 0;
          
        } // for (tagsIdx = 0;
        
      } // for (float charge = 1.0; 
      
    } // for (specIdx = 0;
    
  } // MatchTagsToSpectra()
  
  
} // namespace specnets

