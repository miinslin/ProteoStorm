#ifndef LABEL_H
#define LABEL_H

#include "twovalues.h"
#include "spectrum.h"
#include "SpecSet.h"
#include <string>

namespace specnets
{
	using namespace std;

	/**
	 * TODO: add description
	 */
	class Label {
	public:

		/**
		 * TODO: add description
		 */
		string text;

		Label() {
			text.clear();
		}

		/**
		 * TODO: add description
		 *
		 *@param inText
		 */
		Label(char *inText) {
			text.assign(inText);
		}

		/**
		 * TODO: add description
		 *
		 *@param inText
		 */
		Label(string &inText) {
			text.assign(inText);
		}

		/**
		 * TODO: add description
		 *
		 *@return
		 */
		string &asString() {
			return text;
		}
	};


	/**
	 *
	 * PrmLabel - labels for PRM peaks' annotations of the form "b(nb) y(ny) m(nm) o(no) e(ne)"
	 *	   b/nb - number of b-ion peaks
	 *	   y/ny - number of y-ion peaks
	 *	   o/no - number of non-b and non-y peaks
	 *	   m/nm - number of b _and_ y ion peaks (mix)
	 *	   e/ne - number of endpoint peaks (entries above account for spectrum peaks only and exclude these)
	 *
	 */
	class PeakLabel: public Label {

		static char typeLabels[5];
	public:

		/**
		 * TODO: add description
		 * positions 0-4 are nb/ny/no/nm/ne
		 */
		vector<int> counts;

		/**
		 * TODO: add description
		 */
		PeakLabel(int cb = 0, int cy = 0, int co = 0, int cm = 0, int ce = 0) {
			counts.resize(5);
			set(cb, cy, co, cm, ce);
		}

		/**
		 * TODO: add description
		 *
		 *@return
		 */
		int &operator[](int idx) {
			return counts[idx];
		}

		/**
		 * TODO: add description
		 *
		 *@return
		 */
		PeakLabel &operator=(const PeakLabel &other) {
			return set(other.counts[0], other.counts[1], other.counts[2],
					other.counts[3], other.counts[4]);
		}

		/**
		 * b/nb - number of b-ion peaks.
		 *
		 *@return the number of b-ion peaks.
		 */
		int &nb() {
			return counts[0];
		}

		/**
		 * y/ny - number of y-ion peaks.
		 *
		 *@return the number of y-ion peaks.
		 */
		int &ny() {
			return counts[1];
		}

		/**
		 * o/no - number of non-b and non-y peaks.
		 *
		 *@return the number of non-b and non-y peaks.
		 */
		int &no() {
			return counts[2];
		}

		/**
		 * m/nm - number of b _and_ y ion peaks (mix)
		 *
		 *@return the number of b _and_ y ion peaks.
		 */
		int &nm() {
			return counts[3];
		}

		/**
		 * e/ne - number of endpoint peaks (entries above account for spectrum peaks only and exclude these).
		 *
		 *@return the number of endpoint peaks.
		 */
		int &ne() {
			return counts[4];
		}

		/**
		 * TODO: add description
		 *
		 *@return
		 */
		PeakLabel &merge(PeakLabel &other) {
			for (unsigned int i = 0; i < counts.size(); i++)
				counts[i] += other.counts[i];
			return *this;
		}

		/**
		 * TODO: add description
		 *
		 *@return true if counts is empty and false otherwise.
		 */
		bool isEmpty() {
			bool empty = true;
			for (unsigned int i = 0; i < counts.size(); i++)
				empty = empty and counts[i] == 0;
			return empty;
		}


		/**
		 * TODO: add description
		 *
		 *@return
		 */
		PeakLabel &set(int cb, int cy, int co, int cm, int ce) {
			counts[0] = cb;
			counts[1] = cy;
			counts[2] = co;
			counts[3] = cm;
			counts[4] = ce;
			return *this;
		}

		/**
		 * TODO: add description
		 *
		 *@return
		 */
		PeakLabel &reverse() {
			int tmp = counts[0];
			counts[0] = counts[1];
			counts[1] = tmp;
			return *this;
		}

		/**
		 * TODO: add description
		 *
		 *@return
		 */
		string &asString();
	};


	/**
	 * TODO: add description
	 */
	class SpectrumPeakLabels {
	public:

		/**
		 * TODO: add description
		 */
		vector<PeakLabel> peakLabels;

		/**
		 * TODO: add description
		 */
		SpectrumPeakLabels() {
			peakLabels.resize(0);
		}

		/**
		 * TODO: add description
		 */
		unsigned int size() {
			return peakLabels.size();
		}

		/**
		 * TODO: add description
		 */
		void resize(unsigned int newSize) {
			peakLabels.resize(newSize);
		}

		/**
		 * TODO: add description
		 */
		PeakLabel &operator[](unsigned int idx) {
			return peakLabels[idx];
		}

		/**
		 * TODO: add description
		 */
		SpectrumPeakLabels &operator=(const SpectrumPeakLabels &other) {
			peakLabels.resize(other.peakLabels.size());
			for (unsigned int i = 0; i < peakLabels.size(); i++)
				peakLabels[i] = other.peakLabels[i];
			return *this;
		}

		/**
		 * TODO: add description
		 */
		void reverse();
	};


	/**
	 * TODO: add description
	 */
	int LoadLabels(SpecSet &specSet, const char *filename,
			vector<SpectrumPeakLabels> &labels);
}

#endif
