#ifndef TAGS_H
#define TAGS_H

#include "spectrum.h"
#include "SpecSet.h"
#include <vector>
#include <list>
#include <string>

namespace specnets
{
	using namespace std;

	/**
	 * TODO: add description
	 */
	class Tag {
	public:

		/**
		 * Tag sequence (each entry indexes the corresponding amino
		 * acid in AAJumps::aaLetters).
		 */
		vector<char> sequence;

		/**
		 * Tag sequence as a string
		 */
                string strSequence;

		/**
		 * Prefix flanking masses.
		 */
		float flankingPrefix;

		/**
		 * Suffix flanking masses.
		 */
		float flankingSuffix;

		/**
		 * Score for this tag.
		 */
		float score;

		/**
		 * TODO: add description
		 *
		 *@param tagLen
		 */
		Tag(unsigned short tagLen = 0) {
			sequence.resize(tagLen);
		}

		/**
		 * TODO: add description
		 *
		 *@param putHere
		 *@return
		 */
		string &asString(string &putHere) {
			putHere.resize(sequence.size());
			for (unsigned int i = 0; i < sequence.size(); i++)
				putHere[i] = sequence[i];
			return putHere;
		}

		/**
		 * TODO: add description
		 *
		 *@param other
		 *@return
		 */
		const Tag &operator=(const Tag &other) {
			unsigned int seqLen = min(sequence.size(), other.sequence.size());
			for (unsigned int i = 0; i < seqLen; i++)
				sequence[i] = other.sequence[i];
			score = other.score;
			flankingPrefix = other.flankingPrefix;
			flankingSuffix = other.flankingSuffix;
			return *this;
		}
	};


        unsigned int ExtractTagsAllJumps(Spectrum &spec,
                           list<Tag> &tags,
                           float peakTol,
                           unsigned int tagLen,
                           unsigned int maxNumJumps,
                           unsigned int maxNumTags,
                           float peakSkipPenalty);
	/**
	 * TODO: add description
	 *
	 *@param spec
	 *@param tags
	 *@param peakTol
	 *@param tagLen
	 *@param maxNumJumps
	 *@param maxNumTags
	 *@return
	 */
	unsigned int ExtractTags(Spectrum &spec, list<Tag> &tags, float peakTol,
			unsigned int tagLen, unsigned int maxNumJumps, unsigned int maxNumTags = 0,
                        float peakSkipPenalty = 0.0);

	/**
	 * TODO: add description
	 *
	 *@param sequence
	 *@param tags
	 *@param tagLen
	 *@return
	 */
	unsigned int ExtractTags(char *sequence, vector<Tag> &tags, unsigned int tagLen);

	/**
	 * TODO: add description
	 *
	 *@param specs
	 *@param tags
	 *@param peakTol
	 *@param maxCharge
	 *@param missingPeakPenalty
	 *@param noisePeakPenalty
	 *@param maxNumMissedPeaks
	 *@param matchedTagsIdx
	 *@param matchedTagsScores
	 *@param matchedTagsPos
	 */
	void MatchTagsToSpectra(SpecSet &specs, list<Tag> &tags, float peakTol,
			float maxCharge, float missingPeakPenalty, float noisePeakPenalty,
			unsigned int maxNumMissedPeaks,
			vector<list<unsigned int> > &matchedTagsIdx,
			vector<list<float> > &matchedTagsScores,
			vector<list<unsigned int> > &matchedTagsPos);

}

#endif
