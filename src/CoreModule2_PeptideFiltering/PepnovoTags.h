#ifndef PEPNOVOTAGS_H_
#define PEPNOVOTAGS_H_

#include <string.h>
#include <stdlib.h>

/*
 * Mass-spectrometry Filtering Project
 * This program is property of University of California at San Diego
 *
 * Supervisor: Pavel Pevzner
 * @author Dumitru Brinza
 *
 * ppevzner@cs.ucsd.edu
 * dima@cs.ucsd.edu
 *
 */
struct TTag {

	/**
	 * TODO: add description
	 */
	int static NumberOfFlankingMassesToMatch;

	/**
	 * TODO: add description
	 */
	char sequence[3];

	/**
	 * TODO: add description
	 */
	unsigned short n_gap;

	/**
	 * TODO: add description
	 */
	unsigned short c_gap;

	/**
	 * TODO: add description
	 */
	unsigned char charge;

	/**
	 * TODO: add description
	 */
	float confedence;

	/**
	 * TODO: add description
	 *
	 *@param t
	 */
	void set(TTag & t) {
		confedence = t.confedence;
		n_gap = t.n_gap;
		c_gap = t.c_gap;
		charge = t.charge;
		memcpy(sequence, t.sequence, 3);
	}

	/**
	 * TODO: add description
	 *
	 *@param tag
	 *@return
	 */
	bool operator==(const TTag& tag) const {

		if (sequence[0] == tag.sequence[0] && sequence[1] == tag.sequence[1]
				&& sequence[2] == tag.sequence[2]
				&& (NumberOfFlankingMassesToMatch == 0
				|| (NumberOfFlankingMassesToMatch == 1
				&& (abs(n_gap - tag.n_gap) <= 15 || abs(c_gap - tag.c_gap) <= 30))
				|| (NumberOfFlankingMassesToMatch == 2
				&& (abs(n_gap - tag.n_gap) <= 15
				&& abs(c_gap - tag.c_gap) <= 30))))
			return true;

		if (sequence[0] == tag.sequence[2] && sequence[1] == tag.sequence[1]
				&& sequence[2] == tag.sequence[0]
				&& (NumberOfFlankingMassesToMatch == 0
				|| (NumberOfFlankingMassesToMatch == 1
				&& (abs(n_gap - (tag.c_gap - 10)) <= 15 || abs(c_gap - (tag.n_gap + 10)) <= 30))
				|| (NumberOfFlankingMassesToMatch == 2
				&& (abs(n_gap - (tag.c_gap - 10)) <= 15 || abs(c_gap - (tag.n_gap + 10)) <= 30))))
			return true;
		return false;
	}
};

/**
 * TODO: add description
 *
 *@param FileName
 *@param tags
 *@param numberOfFlankingMassesToMatch
 */
void LoadTags(std::string FileName, std::vector<std::vector<TTag> > &tags,
		int numberOfFlankingMassesToMatch);

/**
 * TODO: add description
 *
 *@param tags
 *@param spectrumId
 *@param candidates
 *@param numberOfTagsToMatch
 *@param numberOfFlankingMassesToMatch
 */
void IntersectTags(std::vector<std::vector<TTag> > &tags,
		unsigned int spectrumId, std::vector<unsigned int>& candidates,
		int numberOfTagsToMatch, int numberOfFlankingMassesToMatch);

#endif /*PEPNOVOTAGS_H_*/
