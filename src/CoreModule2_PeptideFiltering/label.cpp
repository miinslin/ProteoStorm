#include "label.h"
//#include "batch.h"
#include <cstdio>

namespace specnets
{

	char PeakLabel::typeLabels[5] = {'b', 'y', 'o', 'm', 'e'};

	string &PeakLabel::asString() {
		text.clear();   char *buffer = (char*)malloc(128);
		bool firstEntry=true;
		for(unsigned int i=0; i<counts.size(); i++)
			if(counts[i]>0) {
				if(firstEntry) { sprintf(buffer,"%c(%d)",typeLabels[i],counts[i]); firstEntry=false; }
				else sprintf(buffer," %c(%d)",typeLabels[i],counts[i]);
				text.append(buffer);
			}
		free(buffer);
		return text;
	}

	void SpectrumPeakLabels::reverse() {
		vector<PeakLabel> oldLabels(peakLabels.size());
		for(unsigned int i=0; i<peakLabels.size(); i++) oldLabels[peakLabels.size()-1-i] = peakLabels[i];
		for(unsigned int i=0; i<peakLabels.size(); i++) peakLabels[i] = oldLabels[i].reverse();
	}

	int LoadLabels(SpecSet &specSet, const char *filename, vector<SpectrumPeakLabels> &labels) {
		vector<vector<short> > data;
		int res = Load_binArray(filename, data); if (res<=0) return res;
		if(data.size()<1 or data[0].size()!=1) { cerr<<"ERROR: missing values in "<<filename<<".\n"; return -1; }

		unsigned int dataIdx=0, specIdx, peakIdx;
		labels.resize(specSet.size());
		for(specIdx=0; specIdx<specSet.size(); specIdx++) {
			labels[specIdx].resize(specSet[specIdx].size());
			for(peakIdx=0; peakIdx<specSet[specIdx].size(); peakIdx++) {
				if(dataIdx>=data.size()) { cerr<<"ERROR: not enough labels in "<<filename<<" for all peaks in the spectrum set (dataIdx="<<dataIdx<<").\n"; return -1; }
				labels[specIdx][peakIdx][data[dataIdx][0]]++;
				dataIdx++;
			}
		}
		return 1;
	}
}
