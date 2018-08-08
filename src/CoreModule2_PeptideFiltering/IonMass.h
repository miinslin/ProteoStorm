/*
 * IonMass.h
 *
 *  Created on: May 14, 2013
 *      Author: jsnedecor
 */

#ifndef IONMASS_H_
#define IONMASS_H_

#include "spectrum_scoring.h"

namespace specnets
{
  class IonMass
  {
  public:
    IonMass();

    inline bool operator<(const IonMass & rhs) const
    {
        return  m_mass < rhs.m_mass;
    }

    float m_mass;
    string m_name;
    unsigned int m_index;
  };
}

#endif /* IONMASS_H_ */
