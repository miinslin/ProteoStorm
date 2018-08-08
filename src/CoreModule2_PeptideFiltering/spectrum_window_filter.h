#ifndef SPECTRUMWINDOWFILTER_H
#define SPECTRUMWINDOWFILTER_H

#include "spectrum.h"
#include <iostream>
#include <algorithm>

using namespace std;
using namespace specnets;

bool mass_intensity_pair_comp(pair < float, float > mass_intensity_1, pair < float, float > mass_intensity_2);
bool mass_intensity_pair_mass_comp(pair < float, float > mass_intensity_1, pair < float, float > mass_intensity_2);
int filter_window(Spectrum * s, int window, int max_peaks);


#endif
