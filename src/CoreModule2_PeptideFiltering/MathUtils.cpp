/*
 * combination.cpp
 *
 *  Created on: Jul 20, 2012
 *      Author: jsnedecor
 */

/*
 * Algorithm T from Knuth volume 4, 7.2.1.3, for enumerating all combinations k choose N in lexicographic order.
 */
#include "MathUtils.h"

namespace specnets 
{
  void MathUtils::combinations(int n, int k, vector<vector<int> > &combinations)
  {
    // T1. Initialize.
    int j, x;
    int * c = NULL;
  
    c = new int[k + 3];
  
    for (j = 1; j <= k; j++)
    {
      c[j] = j - 1;
    }
    c[k + 1] = n;
    c[k + 2] = 0;
  
    j = k;
  
    do
    {
      //store result
      vector<int> temp(k);
      copy(c+1,c+k+1,temp.begin());

      combinations.push_back(temp);
  
      if (j > 0)
      {
        x = j;
      }
      else
      {
        //T3
        if (c[1] + 1 < c[2])
        {
          c[1]++;
          continue;
        }
        else
        {
          j = 2;
        }
      }
  
      //T4. Find j
      do
      {
        c[j - 1] = j - 2;
        x = c[j] + 1;
        if (x == c[j+1])
        {
          j++;
        }
      } while (x == c[j+1]);
  
      //T6. Increase c[j]
      c[j] = x;
      j--;
    } while (j < k); //T5
  
  }
}
