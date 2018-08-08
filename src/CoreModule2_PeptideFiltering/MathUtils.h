/*
 * combination.h
 *
 *  Created on: Jul 20, 2012
 *      Author: jsnedecor
 */

#ifndef COMBINATION_H_
#define COMBINATION_H_

#include <vector>
#include <iostream>

using namespace std;

namespace specnets
{
  namespace MathUtils
  {
    /** Generates all k length combinations of n.
        */
    void combinations(int n, int k, vector<vector<int> > &combinations);
  }
}

#endif /* COMBINATION_H_ */
